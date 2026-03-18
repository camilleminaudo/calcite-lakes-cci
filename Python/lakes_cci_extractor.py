#!/usr/bin/env python
# coding: utf-8

"""
===============================================================================
Script Name   : lakes_cci_extractor.py
Author        : Camille Minaudo
Created       : 2026-03-17
Last Updated  : 2026-03-17
Description   : Extracts time series of water-leaving reflectance (Rw) from the lakes_cci
                database for a list of lakes from L2/L3 satellite data
                (OLCI / MERIS / MODIS), with parallel processing, mask caching,
                and statistics including bright pixel proportion.

Usage         : python lakes_cci_extractor.py   --lake-file list_lakes.txt   --output-folder ./results   --sensor OLCI   --level L3   --workers 6

                Optional arguments:
                    --lakes GLWD00000411
                    --level L2/L3
                    --sensor OLCI/MERIS/MODIS
                    --start-date YYYY-MM-DD
                    --end-date YYYY-MM-DD
                    --workers N

Notes         :
    - Automatically detects date range if start/end not provided
    - Caches masks per lake/grid for performance
    - Requires Python 3.11+, xarray, geopandas, rioxarray, tqdm
===============================================================================
"""

# ----------------------------
# IMPORTS
# ----------------------------
import argparse
import xarray as xr
import geopandas as gpd
from shapely import wkt
from shapely import contains_xy
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
import os
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm  # for progress bars

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
# UTIL: DATE EXTRACTION
# ----------------------------
def extract_date_from_file(f, root, level):
    try:
        if level == "L2":
            return datetime.strptime(f.split("____")[1][:8], "%Y%m%d")
        elif level == "L3":
            return datetime.strptime(f.split("-")[-2], "%Y%m%d")
    except Exception:
        pass

    for p in Path(root).parts:
        try:
            return datetime.strptime(p, "%Y%m%d")
        except Exception:
            continue
    return None


# ----------------------------
# UTIL: LOAD OR CREATE MASK
# ----------------------------
def get_mask(lake_id, geom, lat2d, lon2d):
    grid_sig = f"{lat2d.shape}_{lon2d.shape}"
    mask_file = MASK_CACHE_DIR / f"{lake_id}_{grid_sig}.npy"

    if mask_file.exists():
        return np.load(mask_file)

    mask_flat = contains_xy(
        geom.iloc[0].geometry,
        lon2d.flatten(),
        lat2d.flatten()
    )

    np.save(mask_file, mask_flat)
    return mask_flat


# ----------------------------
# CORE FUNCTION
# ----------------------------
def process_lake(lake_id):
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

        geom = gpd.GeoDataFrame(
            geometry=[wkt.loads(wkt_file.read_text().strip())],
            crs="EPSG:4326"
        )

        # ----------------------------
        # DATA ROOT
        # ----------------------------
        if data_level == "L2":
            data_root = BASE_DATA_FOLDER / f"lwlr_v3.0_{sensor_str}" / "L2" / "v3.0" / lake_id
        else:
            data_root = BASE_DATA_FOLDER / f"lwlr_v3.0_{sensor_str}" / "L3_com" / "v3.0.0" / "1D" / "1km" / lake_id

        if not data_root.exists():
            logger.warning(f"No data for {lake_id}")
            return

        # ----------------------------
        # FIND FILES
        # ----------------------------
        files = []
        dates = []

        for root, _, fs in os.walk(data_root):
            for f in fs:
                if not f.endswith(".nc"):
                    continue

                d = extract_date_from_file(f, root, data_level)
                if d:
                    files.append((Path(root) / f, d))
                    dates.append(d)

        if not files:
            logger.warning(f"No valid files for {lake_id}")
            return

        # ----------------------------
        # DATE RANGE
        # ----------------------------
        start_use = date_start or min(dates)
        end_use = date_end or max(dates)

        files = [f for f, d in files if start_use <= d <= end_use]

        logger.info(f"{lake_id}: {len(files)} files")

        # ----------------------------
        # OUTPUT
        # ----------------------------
        output_folder.mkdir(parents=True, exist_ok=True)

        output_csv = output_folder / (
            f"Rw_{sensor}_{data_level}_{lake_name}_{lake_id}_"
            f"{start_use.strftime('%Y%m%d')}_{end_use.strftime('%Y%m%d')}.csv"
        )

        # Proceed with extraction only if data was not extracted before or if the corresponding file is empty
        output_csv = output_folder / (
            f"Rw_{sensor}_{data_level}_{lake_name}_{lake_id}_"
            f"{start_use.strftime('%Y%m%d')}_{end_use.strftime('%Y%m%d')}.csv"
        )

        if output_csv.exists() and not args.force:
            if output_csv.stat().st_size > 0:
                logger.info(f"Skipping {lake_id} (already processed)")
                return
            else:
                logger.warning(f"Empty file detected for {lake_id}, reprocessing")

        # ----------------------------
        # PROCESS FILES
        # ----------------------------
        all_stats = []
        mask_cache = {}

        for nc_file in sorted(files):

            ds = xr.open_dataset(nc_file)

            if "time" in ds.coords:
                date = pd.to_datetime(str(ds.time.values[0]))
            else:
                date = pd.to_datetime(nc_file.name.split("-")[-2])

            lat = ds["lat"]
            lon = ds["lon"]

            if "time" in lat.dims:
                lat = lat.isel(time=0)
                lon = lon.isel(time=0)

            if lat.ndim == 1:
                lon2d, lat2d = np.meshgrid(lon.values, lat.values)
            else:
                lat2d, lon2d = lat.values, lon.values

            # Dealing with masks and grid: do not compute mask if already in cache
            grid_sig = f"{lat2d.shape}_{lon2d.shape}_{np.nanmin(lat2d):.4f}_{np.nanmax(lat2d):.4f}"

            if grid_sig not in mask_cache:
                mask_cache[grid_sig] = get_mask(lake_id, geom, lat2d, lon2d)

            mask_flat = mask_cache[grid_sig]

            rw_vars = [v for v in ds.data_vars if v.startswith("Rw") and "uncertainty" not in v]

            chla_var = "chla_top_3_weighted" if data_level == "L2" else "chla"
            tsm_var = "tsm_top_3_weighted" if data_level == "L2" else "tsm"

            # BGR
            def get_arr(name):
                arr = ds[name]
                return arr.isel(time=0).values if "time" in arr.dims else arr.values

            rw490, rw560, rw665 = get_arr("Rw490"), get_arr("Rw560"), get_arr("Rw665")

            bgr = 0.5 * np.abs(
                490 * rw560 + 560 * rw665 + 665 * rw490
                - 560 * rw490 - 665 * rw560 - 490 * rw665
            )

            # -------------------------
            # BRIGHT PIXELS COUNT
            # -------------------------
            if data_level == "L2":
                try:
                    rw412 = get_arr("Rw412")
                    rw560_b = get_arr("Rw560")  # avoid overwrite
                    rw865 = get_arr("Rw865")

                    bright_mask = (
                            (rw412 > 0.4) |
                            (rw560_b > 0.4) |
                            (rw865 > 0.4)
                    )

                    n_bright_pixels = np.sum(bright_mask.flatten()[mask_flat])

                except KeyError:
                    n_bright_pixels = np.nan

            elif data_level == "L3":
                if "lwlr_quality_flags" in ds:
                    flags = get_arr("lwlr_quality_flags")
                    bright_mask = (flags == 8)
                    n_bright_pixels = np.sum(bright_mask.flatten()[mask_flat])
                else:
                    n_bright_pixels = np.nan

            else:
                n_bright_pixels = np.nan

            for band in rw_vars + [chla_var, tsm_var, "BGR"]:

                arr = bgr if band == "BGR" else get_arr(band)

                pix = arr.flatten()[mask_flat]
                pix = pix[~np.isnan(pix)]
                pix = pix[(pix >= 0) & (pix <= 10)]  # constrain to realistic Rw values (which should be 0-1 anyway)

                if pix.size == 0:
                    stats = dict.fromkeys(
                        ["mean", "median", "std", "p5", "p10", "p25", "p75", "p90", "p95"], np.nan
                    )
                    count = 0
                else:
                    stats = {
                        "mean": np.mean(pix),
                        "median": np.median(pix),
                        "std": np.std(pix),
                        "p5": np.percentile(pix, 5),
                        "p10": np.percentile(pix, 10),
                        "p25": np.percentile(pix, 25),
                        "p75": np.percentile(pix, 75),
                        "p90": np.percentile(pix, 90),
                        "p95": np.percentile(pix, 95),
                    }
                    count = pix.size

                all_stats.append({
                    "date": date,
                    "band": band,
                    "count": count,
                    "n_bright_pixels": n_bright_pixels,
                    **stats
                })

            ds.close()

        pd.DataFrame(all_stats).to_csv(output_csv, index=False)

        logger.info(f"Saved {output_csv}")

    except Exception as e:
        logger.exception(f"Failed lake {lake_id}: {e}")


# ----------------------------
# PARALLEL EXECUTION
# ----------------------------
if __name__ == "__main__":

    logger.info(f"Starting processing with {args.workers} workers")

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(process_lake, lid) for lid in lake_ids]

        for _ in tqdm(as_completed(futures), total=len(futures), desc="Processing lakes"):
            try:
                _.result()
            except Exception:
                pass

    logger.info("All lakes processed")