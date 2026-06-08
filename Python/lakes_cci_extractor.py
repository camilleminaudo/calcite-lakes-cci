#!/usr/bin/env python
# coding: utf-8

"""
===============================================================================
Script Name   : lakes_cci_extractor.py
Author        : Camille Minaudo
Created       : 2026-03-17
Last Updated  : 2026-06-05
Description   : Extracts time series of water-leaving reflectance (Rw) from the
                lakes_cci database for a list of lakes from L2/L3 satellite data
                (OLCI / MERIS / MODIS), with parallel processing, mask caching,
                and statistics including bright pixel proportion.

                Masks and extraction are interleaved in batches to avoid
                spending days precomputing masks before any extraction starts.

Usage         : python lakes_cci_extractor.py \
                    --lake-file list_lakes.txt \
                    --output-folder ./results \
                    --sensor OLCI \
                    --level L3 \
                    --workers 6

                Optional arguments:
                    --lakes GLWD00000411
                    --level L2/L3
                    --sensor OLCI/MERIS/MODIS
                    --start-date YYYY-MM-DD
                    --end-date YYYY-MM-DD
                    --workers N
                    --batch-size N   (lakes per mask+extract cycle, default 5)

Notes         :
    - Mask is computed once per unique spatial grid (lake × grid extent hash)
      and cached on disk. This is safe: lat/lon grids are fixed by sensor
      geometry and do not vary with cloud cover or data quality.
    - Requires Python 3.11+, xarray, geopandas, rioxarray, tqdm
    - NetCDF engine auto-detected at startup: netCDF4 (preferred) → h5netcdf → scipy
===============================================================================
"""

# ----------------------------
# IMPORTS
# ----------------------------
import argparse
import hashlib
import xarray as xr
import geopandas as gpd
from shapely import wkt
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
import os
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

try:
    from rasterio.features import geometry_mask
    from rasterio.transform import from_bounds
    USE_RASTERIO = True
except ImportError:
    from shapely import contains_xy
    USE_RASTERIO = False
    logging.getLogger(__name__).warning(
        "rasterio not found — falling back to shapely contains_xy (slower). "
        "Install rasterio for significantly faster mask computation."
    )

# ----------------------------
# NETCDF ENGINE SELECTION
# ----------------------------
def _probe_nc_engine() -> str:
    """
    Return the best available xarray engine for reading NetCDF/HDF5 files.
    Preference order: netcdf4 → h5netcdf → scipy.
    Raises RuntimeError if none are available.
    """
    _candidates = ["netcdf4", "h5netcdf", "scipy"]
    _probe_file = None
    for engine in _candidates:
        try:
            import importlib
            # netcdf4 engine needs the 'netCDF4' package; h5netcdf needs 'h5netcdf'; scipy needs 'scipy'
            _pkg = {"netcdf4": "netCDF4", "h5netcdf": "h5netcdf", "scipy": "scipy"}[engine]
            importlib.import_module(_pkg)
            return engine
        except ImportError:
            continue
    raise RuntimeError(
        "No NetCDF reading engine found. Install at least one of: netCDF4, h5netcdf, scipy.\n"
        "  conda install -c conda-forge netcdf4      # recommended\n"
        "  pip install netCDF4"
    )

NC_ENGINE = _probe_nc_engine()
logging.getLogger(__name__).info(f"NetCDF engine: {NC_ENGINE}")


def open_nc(path, **kwargs) -> xr.Dataset:
    """Thin wrapper around xr.open_dataset using the detected engine."""
    return xr.open_dataset(path, engine=NC_ENGINE, **kwargs)


# ----------------------------
# CONSTANTS
# ----------------------------
BASE_DATA_FOLDER = Path("/data/datasets/Projects/lakes_cci/processing/data")

SENSOR_MAP = {
    "OLCI": "olci",
    "MERIS": "meris",
    "MODIS": "modis"
}

MASK_CACHE_DIR = Path("./mask_cache")
MASK_CACHE_DIR.mkdir(exist_ok=True)

# Realistic reflectance bounds (applied to Rw bands only, not chla/tsm/BGR)
RW_MIN, RW_MAX = 0.0, 1.0

# ----------------------------
# LOGGING SETUP
# ----------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(processName)s | %(message)s"
)
logger = logging.getLogger(__name__)

# ----------------------------
# ARGUMENT PARSER
# ----------------------------
parser = argparse.ArgumentParser(description="Extract Rw time series for lakes")

parser.add_argument("--lakes", type=str, help="Comma-separated lake IDs")
parser.add_argument("--lake-file", type=str, help="File with lake IDs")
parser.add_argument("--output-folder", type=str, required=True)
parser.add_argument("--sensor", type=str, required=True, choices=["MODIS", "MERIS", "OLCI"])
parser.add_argument("--level", type=str, default="L2", choices=["L2", "L3"])
parser.add_argument("--start-date", type=str, help="YYYY-MM-DD")
parser.add_argument("--end-date", type=str, help="YYYY-MM-DD")
parser.add_argument("--workers", type=int, default=4)
parser.add_argument("--batch-size", type=int, default=5,
                    help="Number of lakes to mask+extract per cycle (default: 5)")

args = parser.parse_args()

# ----------------------------
# INPUT HANDLING
# ----------------------------
sensor = args.sensor.upper()
sensor_str = SENSOR_MAP[sensor]
data_level = args.level.upper()
output_folder = Path(args.output_folder)

date_start = pd.to_datetime(args.start_date) if args.start_date else None
date_end = pd.to_datetime(args.end_date) if args.end_date else None

# lakes
if args.lakes:
    lake_ids = [l.strip() for l in args.lakes.split(",") if l.strip()]
elif args.lake_file:
    with open(args.lake_file) as f:
        lake_ids = [line.strip() for line in f if line.strip()]
else:
    raise ValueError("Provide --lakes or --lake-file")

# ----------------------------
# STATIC PATHS
# ----------------------------
lake_info_folder = Path("/users/rsg-new/cami/project_Rw_time_series/data/lake-polygons-PML/")
metadata_file = lake_info_folder / "cglops_4k_metadata_20210301.csv"
wkt_folder = lake_info_folder / "CLMS4k"

lake_meta = pd.read_csv(metadata_file)
meta_lookup = lake_meta.set_index("id_str")


# ----------------------------
# UTIL: DATA ROOT
# ----------------------------
def get_data_root(lake_id: str) -> Path | None:
    if sensor == "OLCI":
        if data_level == "L2":
            base = BASE_DATA_FOLDER / "lwlr_v3.0_olci_archive" / lake_id
        else:
            base = BASE_DATA_FOLDER / "lwlr_v3.0_olci_current" / lake_id
    else:
        if data_level == "L2":
            base = BASE_DATA_FOLDER / f"lwlr_v3.0_{sensor_str}" / "L2" / "v3.0" / lake_id
        else:
            base = BASE_DATA_FOLDER / f"lwlr_v3.0_{sensor_str}" / "L3_com" / "v3.0.0" / "1D" / "1km" / lake_id

    return base if base.exists() else None


# ----------------------------
# UTIL: DATE EXTRACTION
# ----------------------------
def extract_date_from_file(f: str, root: str) -> datetime | None:
    # OLCI (both L2 and L3): .../lake_id/YYYY/MM/DD/filename.nc
    # root is the *directory* containing f, so parse from the full path
    if sensor == "OLCI":
        parts = Path(root).parts
        for i in range(len(parts) - 2):
            y, m, d = parts[i], parts[i + 1], parts[i + 2]
            if len(y) == 4 and len(m) == 2 and len(d) == 2:
                try:
                    return datetime.strptime(f"{y}{m}{d}", "%Y%m%d")
                except ValueError:
                    continue

    # MERIS / MODIS: date embedded in filename
    try:
        if data_level == "L2":
            return datetime.strptime(f.split("____")[1][:8], "%Y%m%d")
        elif data_level == "L3":
            return datetime.strptime(f.split("-")[-2], "%Y%m%d")
    except Exception:
        pass

    return None


# ----------------------------
# UTIL: GRID COORDINATE EXTRACTION
# ----------------------------
def extract_grid_coords(ds: xr.Dataset) -> tuple[np.ndarray, np.ndarray]:
    """Return (lat2d, lon2d) from a dataset, handling both 1-D and 2-D coordinates."""
    lat = ds["lat"]
    lon = ds["lon"]
    if "time" in lat.dims:
        lat = lat.isel(time=0)
        lon = lon.isel(time=0)
    if lat.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon.values, lat.values)
    else:
        lat2d, lon2d = lat.values, lon.values
    return lat2d, lon2d


# ----------------------------
# UTIL: STABLE GRID HASH
# ----------------------------
def grid_hash(lat2d: np.ndarray, lon2d: np.ndarray) -> str:
    key_vals = np.array([
        lat2d.shape[0], lat2d.shape[1],
        lat2d[0, 0], lat2d[-1, -1],
        lon2d[0, 0], lon2d[-1, -1],
        lat2d[lat2d.shape[0] // 2, lat2d.shape[1] // 2],
        lon2d[lat2d.shape[0] // 2, lat2d.shape[1] // 2],
    ])
    return hashlib.md5(key_vals.tobytes()).hexdigest()[:12]


# ----------------------------
# UTIL: MASK COMPUTATION
# ----------------------------
def compute_mask(geom_shape, lat2d: np.ndarray, lon2d: np.ndarray) -> np.ndarray:
    rows, cols = lat2d.shape

    if USE_RASTERIO:
        transform = from_bounds(
            lon2d.min(), lat2d.min(), lon2d.max(), lat2d.max(), cols, rows
        )
        outside = geometry_mask(
            [geom_shape.__geo_interface__],
            transform=transform,
            invert=False,
            out_shape=(rows, cols),
        )
        return (~outside).flatten()
    else:
        mask_flat = contains_xy(
            geom_shape,
            lon2d.flatten(),
            lat2d.flatten()
        )
        return mask_flat


# ----------------------------
# UTIL: LOAD OR BUILD MASK
# ----------------------------
def get_or_build_mask(lake_id: str, geom_shape, lat2d: np.ndarray, lon2d: np.ndarray) -> np.ndarray:
    """Load the cached mask or compute and persist it."""
    ghash = grid_hash(lat2d, lon2d)
    mask_file = MASK_CACHE_DIR / f"{lake_id}_{ghash}.npy"

    if mask_file.exists():
        return np.load(mask_file)

    mask_flat = compute_mask(geom_shape, lat2d, lon2d)

    tmp_file = mask_file.with_suffix(".tmp.npy")
    np.save(tmp_file, mask_flat)
    tmp_file.rename(mask_file)

    return mask_flat


# ----------------------------
# UTIL: ALREADY-DONE CHECK
# ----------------------------
def is_already_done(lake_id: str) -> bool:
    """
    Return True if a non-empty output CSV already exists for this lake.

    The filename encodes sensor, level, lake name, lake id, and the date
    range. When --start-date / --end-date are supplied we can build the
    exact name. When they are not, we fall back to a glob on lake_id so we
    don't need to scan any NetCDF files just to decide whether to skip.

    Logs a clear reason so the user can see what was skipped and why.
    """
    if lake_id not in meta_lookup.index:
        # Metadata missing — let process_lake emit the proper error later.
        return False

    lake_name = meta_lookup.loc[lake_id, "name"]

    if date_start and date_end:
        # Exact filename is known without touching any data files.
        candidate = output_folder / (
            f"Rw_{sensor}_{data_level}_{lake_name}_{lake_id}_"
            f"{date_start.strftime('%Y%m%d')}_{date_end.strftime('%Y%m%d')}.csv"
        )
        if candidate.exists() and candidate.stat().st_size > 0:
            logger.info(f"Skipping {lake_id} — output already exists: {candidate.name}")
            return True
        if candidate.exists():
            logger.warning(f"{lake_id}: empty output file found, will reprocess")
        return False

    # No explicit date range → glob for any matching file for this lake/sensor/level.
    # Pattern: Rw_<SENSOR>_<LEVEL>_<name>_<id>_*.csv
    pattern = f"Rw_{sensor}_{data_level}_{lake_name}_{lake_id}_*.csv"
    matches = [p for p in output_folder.glob(pattern) if p.stat().st_size > 0]
    if matches:
        logger.info(
            f"Skipping {lake_id} — found existing output(s): "
            + ", ".join(m.name for m in matches)
        )
        return True
    return False


# ----------------------------
# UTIL: MASK PRE-COMPUTATION (single lake)
# ----------------------------
def precompute_masks_for_lake(lake_id: str) -> None:
    """
    Walk the data tree for one lake, find one representative file per unique
    grid, and pre-build masks so worker processes never need to compute them.
    """
    if lake_id not in meta_lookup.index:
        return

    wkt_file = wkt_folder / f"{lake_id}.txt"
    if not wkt_file.exists():
        return

    geom_shape = wkt.loads(wkt_file.read_text().strip())
    data_root = get_data_root(lake_id)

    if not data_root or not data_root.exists():
        return

    seen_grids: set[str] = set()

    for root, _, files in os.walk(data_root):
        for f in files:
            if not f.endswith(".nc"):
                continue
            nc_file = Path(root) / f
            try:
                with open_nc(nc_file) as ds:
                    lat2d, lon2d = extract_grid_coords(ds)

                ghash = grid_hash(lat2d, lon2d)
                if ghash in seen_grids:
                    continue
                seen_grids.add(ghash)

                get_or_build_mask(lake_id, geom_shape, lat2d, lon2d)

            except Exception as e:
                logger.warning(f"Mask precompute failed for {lake_id} / {f}: {e}")


# ----------------------------
# UTIL: VECTORISED BAND STATS
# ----------------------------
STAT_PERCENTILES = [5, 10, 25, 75, 90, 95]
STAT_NAMES = ["mean", "median", "std", "p5", "p10", "p25", "p75", "p90", "p95"]
_NAN_STATS = {k: np.nan for k in STAT_NAMES}


def compute_stats(pixels: np.ndarray) -> dict:
    if pixels.size == 0:
        return {**_NAN_STATS, "count": 0}
    return {
        "mean": np.mean(pixels),
        "median": np.median(pixels),
        "std": np.std(pixels),
        "p5": np.percentile(pixels, 5),
        "p10": np.percentile(pixels, 10),
        "p25": np.percentile(pixels, 25),
        "p75": np.percentile(pixels, 75),
        "p90": np.percentile(pixels, 90),
        "p95": np.percentile(pixels, 95),
        "count": pixels.size,
    }


def extract_pixels(ds: xr.Dataset, var_name: str, mask_flat: np.ndarray,
                   clip_min: float = -np.inf, clip_max: float = np.inf) -> np.ndarray:
    arr = ds[var_name]
    values = (arr.isel(time=0).values if "time" in arr.dims else arr.values).flatten()
    pix = values[mask_flat]
    pix = pix[np.isfinite(pix)]
    if np.isfinite(clip_min) or np.isfinite(clip_max):
        pix = pix[(pix >= clip_min) & (pix <= clip_max)]
    return pix


# ----------------------------
# CORE FUNCTION
# ----------------------------
def process_lake(lake_id: str) -> None:
    try:
        logger.info(f"Processing lake {lake_id}")

        if lake_id not in meta_lookup.index:
            logger.error(f"Lake {lake_id} not in metadata")
            return

        lake_name = meta_lookup.loc[lake_id, "name"]

        wkt_file = wkt_folder / f"{lake_id}.txt"
        if not wkt_file.exists():
            logger.error(f"WKT missing for {lake_id}")
            return

        geom_shape = wkt.loads(wkt_file.read_text().strip())

        data_root = get_data_root(lake_id)
        if not data_root or not data_root.exists():
            logger.warning(f"No data folder for {lake_id}")
            return

        # ----------------------------
        # COLLECT (file, date) PAIRS
        # ----------------------------
        file_date_pairs: list[tuple[Path, datetime]] = []
        for root, _, fs in os.walk(data_root):
            for f in fs:
                if not f.endswith(".nc"):
                    continue
                d = extract_date_from_file(f, root)
                if d:
                    file_date_pairs.append((Path(root) / f, d))

        if not file_date_pairs:
            logger.warning(f"No valid files for {lake_id}")
            return

        # ----------------------------
        # DATE RANGE FILTER
        # ----------------------------
        all_dates = [d for _, d in file_date_pairs]
        start_use = date_start or min(all_dates)
        end_use = date_end or max(all_dates)

        file_date_pairs = [
            (fp, d) for fp, d in file_date_pairs if start_use <= d <= end_use
        ]
        logger.info(f"{lake_id}: {len(file_date_pairs)} files in range")

        # ----------------------------
        # SKIP IF ALREADY DONE
        # ----------------------------
        output_folder.mkdir(parents=True, exist_ok=True)
        output_csv = output_folder / (
            f"Rw_{sensor}_{data_level}_{lake_name}_{lake_id}_"
            f"{start_use.strftime('%Y%m%d')}_{end_use.strftime('%Y%m%d')}.csv"
        )

        if output_csv.exists() and output_csv.stat().st_size > 0:
            logger.info(f"Skipping {lake_id} (already processed)")
            return
        if output_csv.exists():
            logger.warning(f"Empty file detected for {lake_id}, reprocessing")

        # ----------------------------
        # BAND NAMES
        # ----------------------------
        chla_var = "chla_top_3_weighted" if data_level == "L2" else "chla"
        tsm_var  = "tsm_top_3_weighted"  if data_level == "L2" else "tsm"

        # ----------------------------
        # PROCESS FILES
        # ----------------------------
        all_stats: list[dict] = []
        mask_cache: dict[str, np.ndarray] = {}

        for nc_file, file_date in sorted(file_date_pairs, key=lambda x: x[1]):

            try:
                with open_nc(nc_file) as ds:

                    date = pd.Timestamp(file_date)

                    lat2d, lon2d = extract_grid_coords(ds)
                    ghash = grid_hash(lat2d, lon2d)

                    if ghash not in mask_cache:
                        mask_cache[ghash] = get_or_build_mask(
                            lake_id, geom_shape, lat2d, lon2d
                        )
                    mask_flat = mask_cache[ghash]

                    if not mask_flat.any():
                        logger.debug(f"{lake_id}: empty mask for {nc_file.name}, skipping")
                        continue

                    rw_vars = [
                        v for v in ds.data_vars
                        if v.startswith("Rw") and "uncertainty" not in v
                    ]

                    def get_2d(name: str) -> np.ndarray:
                        arr = ds[name]
                        return (arr.isel(time=0).values if "time" in arr.dims else arr.values)

                    rw490_2d = get_2d("Rw490")
                    rw560_2d = get_2d("Rw560")
                    rw665_2d = get_2d("Rw665")

                    bgr_2d = 0.5 * np.abs(
                        490 * rw560_2d + 560 * rw665_2d + 665 * rw490_2d
                        - 560 * rw490_2d - 665 * rw560_2d - 490 * rw665_2d
                    )

                    if data_level == "L2":
                        try:
                            rw412_2d = get_2d("Rw412")
                            rw865_2d = get_2d("Rw865")
                            bright_2d = (rw412_2d > 0.4) | (rw560_2d > 0.4) | (rw865_2d > 0.4)
                            n_bright = int(np.sum(bright_2d.flatten()[mask_flat]))
                        except KeyError:
                            n_bright = np.nan
                    elif data_level == "L3" and "lwlr_quality_flags" in ds:
                        flags_2d = get_2d("lwlr_quality_flags")
                        n_bright = int(np.sum((flags_2d.flatten()[mask_flat]) == 8))
                    else:
                        n_bright = np.nan

                    bands_to_process = rw_vars + [chla_var, tsm_var, "BGR"]

                    for band in bands_to_process:
                        if band == "BGR":
                            raw = bgr_2d.flatten()[mask_flat]
                            pix = raw[np.isfinite(raw)]
                        elif band in (chla_var, tsm_var):
                            if band not in ds.data_vars:
                                continue
                            pix = extract_pixels(ds, band, mask_flat)
                        else:
                            pix = extract_pixels(ds, band, mask_flat, RW_MIN, RW_MAX)

                        stats = compute_stats(pix)
                        all_stats.append({
                            "date": date,
                            "band": band,
                            "n_bright_pixels": n_bright,
                            **stats,
                        })

            except Exception as e:
                logger.warning(f"{lake_id}: failed to process {nc_file.name}: {e}")
                continue

        if not all_stats:
            logger.warning(f"{lake_id}: no stats collected, writing empty file")

        pd.DataFrame(all_stats).to_csv(output_csv, index=False)
        logger.info(f"Saved {output_csv}")

    except Exception as e:
        logger.exception(f"Failed lake {lake_id}: {e}")


# ----------------------------
# BATCHED EXECUTION
# ----------------------------
def process_lakes_in_batches(lake_ids: list[str], batch_size: int, n_workers: int) -> None:
    """
    For each batch of `batch_size` lakes:
      1. Precompute masks (single process, sequential — mask I/O is the bottleneck)
      2. Run extraction in parallel with ProcessPoolExecutor

    This prevents the job from spending days on masks before any data is written,
    and keeps memory/disk pressure bounded at `batch_size` lakes at a time.
    """
    total = len(lake_ids)
    n_batches = (total + batch_size - 1) // batch_size

    for batch_idx in range(n_batches):
        batch = lake_ids[batch_idx * batch_size : (batch_idx + 1) * batch_size]
        logger.info(
            f"=== Batch {batch_idx + 1}/{n_batches}: "
            f"lakes {batch_idx * batch_size + 1}–{min((batch_idx + 1) * batch_size, total)} "
            f"of {total} ==="
        )

        # -- Step 1: skip lakes whose output already exists --
        todo = [lid for lid in batch if not is_already_done(lid)]
        n_skip = len(batch) - len(todo)
        if n_skip:
            logger.info(f"Batch {batch_idx + 1}: skipping {n_skip} already-done lake(s)")
        if not todo:
            logger.info(f"Batch {batch_idx + 1}: nothing to do, moving to next batch")
            continue

        # -- Step 2: masks only for lakes that still need processing --
        logger.info(f"Batch {batch_idx + 1}: precomputing masks for {len(todo)} lake(s)")
        for lake_id in tqdm(todo, desc=f"Masks (batch {batch_idx + 1}/{n_batches})"):
            precompute_masks_for_lake(lake_id)

        # -- Step 3: extraction for this batch in parallel --
        logger.info(f"Batch {batch_idx + 1}: starting extraction with {n_workers} workers")
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(process_lake, lid): lid for lid in todo}
            for future in tqdm(
                as_completed(futures),
                total=len(futures),
                desc=f"Extracting (batch {batch_idx + 1}/{n_batches}, step 3)",
            ):
                lid = futures[future]
                try:
                    future.result()
                except Exception as e:
                    logger.error(f"Unhandled exception for lake {lid}: {e}")

        logger.info(f"Batch {batch_idx + 1}/{n_batches} complete")


# ----------------------------
# ENTRY POINT
# ----------------------------
if __name__ == "__main__":
    logger.info(
        f"Starting: {len(lake_ids)} lakes, batch_size={args.batch_size}, "
        f"workers={args.workers}, sensor={sensor}, level={data_level}"
    )
    process_lakes_in_batches(lake_ids, batch_size=args.batch_size, n_workers=args.workers)
    logger.info("All lakes processed")
