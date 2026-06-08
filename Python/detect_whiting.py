#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Detect whiting events from lake time series (BGR / Rw + bright pixel counts)

Key improvements over v1:
- Baseline estimated from the LOWER half of the distribution (so whitings
  don't inflate the reference, making the z-score more sensitive)
- bright_pixel_fraction used as a corroborating signal alongside BGR/Rw
- Seasonal context: each observation is compared against its own month's
  baseline, so gradual seasonal whitings are handled separately from sharp
  episodic events
- Dual-criteria flagging: a point can be flagged by reflectance alone,
  bright-pixel fraction alone, or both (recorded in flag_source)
- Summary table written across all lakes at the end

How to run:
    python detect_whiting.py \
        --input-folder  ./results \
        --output-folder ./events \
        --band          BGR \
        --value-col     mean \
        --z-thresh      5 \
        --bright-thresh 0.10 \
        --min-obs       3 \
        --max-gap-days  10 \
        --min-duration  3

Arguments:
    --band            Band name to use for reflectance signal  [BGR]
    --value-col       Column to use within that band           [mean]
    --z-thresh        MAD z-score threshold for reflectance    [5]
    --bright-thresh   Fraction of bright pixels to flag (0–1)  [0.10]
    --min-obs         Minimum flagged observations per event   [3]
    --max-gap-days    Max gap (days) to merge nearby flags     [10]
    --min-duration    Minimum event duration in days           [3]
    --seasonal        Use month-level baseline (recommended for
                      lakes with regular summer whitings)      [False]
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=RuntimeWarning)


# ─────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────

def robust_zscore(values: np.ndarray, lower_only: bool = True) -> np.ndarray:
    """
    MAD-based z-score.

    When lower_only=True the median and MAD are computed on values ≤ the
    overall median, so that whiting peaks do not inflate the baseline.
    This makes the detector more sensitive to genuine anomalies.
    """
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.full_like(values, np.nan, dtype=float)

    if lower_only:
        reference = finite[finite <= np.nanmedian(finite)]
    else:
        reference = finite

    med = np.nanmedian(reference)
    mad = np.nanmedian(np.abs(reference - med))
    mad = max(mad, 1e-9)

    return (values - med) / mad


def bright_fraction(series_n_bright: pd.Series,
                    mask_flat: pd.Series) -> pd.Series:
    """
    Compute fraction of in-lake pixels that are bright.

    mask_flat is the total pixel count column ('count') for the same rows.
    Returns NaN where count is 0 or missing.
    """
    total = mask_flat.replace(0, np.nan)
    return series_n_bright / total


def merge_events(flagged_dates: pd.Series,
                 max_gap_days: int) -> pd.Series:
    """
    Assign an integer event_id to each flagged date.
    Consecutive flagged dates separated by ≤ max_gap_days are merged.
    Returns a Series aligned with flagged_dates.
    """
    sorted_dates = flagged_dates.sort_values()
    gaps = sorted_dates.diff().dt.days
    new_event = (gaps.isna()) | (gaps > max_gap_days)
    event_ids = new_event.cumsum()
    return event_ids


# ─────────────────────────────────────────────
# CORE DETECTOR
# ─────────────────────────────────────────────

def detect_whiting_events(
    df: pd.DataFrame,
    value_col: str = "mean",
    date_col: str = "date",
    count_col: str = "count",
    n_bright_col: str = "n_bright_pixels",
    z_thresh: float = 5.0,
    bright_thresh: float = 0.10,
    min_obs: int = 3,
    max_gap_days: int = 10,
    min_duration: int = 3,
    seasonal: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Detect whiting events in a single-lake, single-band time series.

    Parameters
    ----------
    df            : DataFrame with at least [date_col, value_col]
    value_col     : reflectance column (e.g. 'mean')
    date_col      : date column
    count_col     : total in-lake pixel count column (for bright fraction)
    n_bright_col  : number of bright pixels column
    z_thresh      : MAD z-score threshold for reflectance anomaly
    bright_thresh : fraction of bright pixels (0–1) to flag independently
    min_obs       : minimum flagged obs for a valid event
    max_gap_days  : days gap allowed within a single event
    min_duration  : minimum event duration in calendar days
    seasonal      : if True, z-score is computed per calendar month so that
                    gradual seasonal whitings are compared to their own
                    monthly baseline rather than the annual one

    Returns
    -------
    df_out  : full time series with z-score, bright_fraction, flags, event_id
    events  : one row per valid detected event with summary metrics
    """

    df = df.copy()
    df[date_col] = pd.to_datetime(df[date_col])
    df = df.sort_values(date_col).reset_index(drop=True)

    # ── Reflectance z-score ──────────────────────────────────────────────
    if seasonal:
        df["month"] = df[date_col].dt.month
        df["z"] = np.nan
        for month, grp in df.groupby("month"):
            df.loc[grp.index, "z"] = robust_zscore(grp[value_col].values)
        df.drop(columns="month", inplace=True)
    else:
        df["z"] = robust_zscore(df[value_col].values)

    # ── Bright pixel fraction ────────────────────────────────────────────
    has_bright = (n_bright_col in df.columns) and (count_col in df.columns)
    if has_bright:
        df["bright_frac"] = bright_fraction(
            df[n_bright_col], df[count_col]
        )
    else:
        df["bright_frac"] = np.nan

    # ── Dual-criteria flagging ───────────────────────────────────────────
    flag_refl  = df["z"] > z_thresh
    flag_bright = df["bright_frac"] > bright_thresh if has_bright else pd.Series(False, index=df.index)

    df["flag_refl"]   = flag_refl
    df["flag_bright"] = flag_bright
    df["event_flag"]  = flag_refl | flag_bright

    # Record which criterion triggered
    df["flag_source"] = "none"
    df.loc[flag_refl,                      "flag_source"] = "reflectance"
    df.loc[flag_bright,                    "flag_source"] = "bright_pixels"
    df.loc[flag_refl & flag_bright,        "flag_source"] = "both"

    # ── Early exit if nothing flagged ────────────────────────────────────
    df_flagged = df[df["event_flag"]].copy()
    if df_flagged.empty:
        df["event_id"]    = np.nan
        df["valid_event"] = False
        return df, pd.DataFrame()

    # ── Merge nearby flags into events ───────────────────────────────────
    event_ids = merge_events(df_flagged[date_col], max_gap_days)
    df_flagged = df_flagged.copy()
    df_flagged["event_id"] = event_ids.values

    # ── Summarise candidate events ───────────────────────────────────────
    agg_dict = {
        "start_date":    (date_col,    "min"),
        "end_date":      (date_col,    "max"),
        "n_obs":         (date_col,    "count"),
        "span_days":     (date_col,    lambda x: (x.max() - x.min()).days + 1),
        "peak_value":    (value_col,   "max"),
        "mean_value":    (value_col,   "mean"),
        "peak_z":        ("z",         "max"),
    }
    if has_bright:
        agg_dict["peak_bright_frac"] = ("bright_frac", "max")

    # flag_source majority vote
    def dominant_source(s):
        counts = s.value_counts()
        return counts.index[0] if not counts.empty else "unknown"

    events = (
        df_flagged
        .groupby("event_id")
        .agg(**{k: v for k, v in agg_dict.items()})
        .reset_index()
    )
    events["flag_source"] = (
        df_flagged.groupby("event_id")["flag_source"]
        .apply(dominant_source)
        .values
    )

    # ── Apply validity filters ────────────────────────────────────────────
    events = events[
        (events["n_obs"] >= min_obs) &
        (events["span_days"] >= min_duration)
    ].copy()

    # ── Tag back to full time series ─────────────────────────────────────
    valid_ids = set(events["event_id"])
    df_flagged["valid_event"] = df_flagged["event_id"].isin(valid_ids)

    df = df.merge(
        df_flagged[[date_col, "event_id", "valid_event"]].drop_duplicates(date_col),
        on=date_col,
        how="left",
    )
    df["valid_event"] = df["valid_event"].fillna(False)
    df["event_id"]    = df["event_id"].where(df["valid_event"])

    return df, events.reset_index(drop=True)


# ─────────────────────────────────────────────
# LAKE-LEVEL SUMMARY METRICS
# ─────────────────────────────────────────────

def compute_lake_summary(
    lake_id: str,
    lake_name: str,
    df_ts: pd.DataFrame,
    events: pd.DataFrame,
    value_col: str,
) -> dict:
    """
    Collapse all events for one lake into a single summary row.
    Includes time-series-level baseline stats as well as event metrics.
    """
    n_obs_total  = len(df_ts)
    n_obs_valid  = df_ts[value_col].notna().sum()
    date_min     = df_ts["date"].min()
    date_max     = df_ts["date"].max()
    record_years = (date_max - date_min).days / 365.25 if n_obs_valid > 1 else np.nan

    ts_median    = df_ts[value_col].median()
    ts_p95       = df_ts[value_col].quantile(0.95)
    ts_std       = df_ts[value_col].std()

    if events.empty:
        return {
            "lake_id":              lake_id,
            "lake_name":            lake_name,
            "record_start":         date_min,
            "record_end":           date_max,
            "record_years":         round(record_years, 2) if not np.isnan(record_years) else np.nan,
            "n_obs_total":          n_obs_total,
            "n_obs_valid":          n_obs_valid,
            "ts_median":            round(ts_median, 5),
            "ts_p95":               round(ts_p95, 5),
            "ts_std":               round(ts_std, 5),
            "n_events":             0,
            "total_event_days":     0,
            "mean_event_duration":  np.nan,
            "max_event_duration":   np.nan,
            "mean_peak_value":      np.nan,
            "max_peak_value":       np.nan,
            "mean_peak_z":          np.nan,
            "max_peak_z":           np.nan,
            "peak_bright_frac":     np.nan,
            "events_per_year":      0.0,
            "dominant_flag_source": "none",
            "first_event_date":     pd.NaT,
            "last_event_date":      pd.NaT,
            "peak_month":           np.nan,
        }

    n_events        = len(events)
    events_per_year = n_events / record_years if record_years and record_years > 0 else np.nan

    # Month with most event starts → seasonality indicator
    start_months    = pd.to_datetime(events["start_date"]).dt.month
    peak_month      = int(start_months.mode().iloc[0]) if not start_months.empty else np.nan

    row = {
        "lake_id":              lake_id,
        "lake_name":            lake_name,
        "record_start":         date_min,
        "record_end":           date_max,
        "record_years":         round(record_years, 2) if not np.isnan(record_years) else np.nan,
        "n_obs_total":          n_obs_total,
        "n_obs_valid":          n_obs_valid,
        "ts_median":            round(ts_median, 5),
        "ts_p95":               round(ts_p95, 5),
        "ts_std":               round(ts_std, 5),
        "n_events":             n_events,
        "total_event_days":     int(events["span_days"].sum()),
        "mean_event_duration":  round(events["span_days"].mean(), 1),
        "max_event_duration":   int(events["span_days"].max()),
        "mean_peak_value":      round(events["peak_value"].mean(), 5),
        "max_peak_value":       round(events["peak_value"].max(), 5),
        "mean_peak_z":          round(events["peak_z"].mean(), 2),
        "max_peak_z":           round(events["peak_z"].max(), 2),
        "peak_bright_frac":     round(events["peak_bright_frac"].max(), 4) if "peak_bright_frac" in events.columns else np.nan,
        "events_per_year":      round(events_per_year, 3) if not np.isnan(events_per_year) else np.nan,
        "dominant_flag_source": events["flag_source"].mode().iloc[0] if not events.empty else "none",
        "first_event_date":     events["start_date"].min(),
        "last_event_date":      events["end_date"].max(),
        "peak_month":           peak_month,
    }
    return row


# ─────────────────────────────────────────────
# FILENAME PARSING
# ─────────────────────────────────────────────

def parse_lake_info(stem: str) -> tuple[str, str]:
    """
    Extract (lake_id, lake_name) from the CSV filename produced by
    lakes_cci_extractor.py:
        Rw_<SENSOR>_<LEVEL>_<lake_name>_<lake_id>_<date1>_<date2>
    Falls back gracefully if the pattern doesn't match.
    """
    parts = stem.split("_")
    # Minimum expected tokens: Rw, SENSOR, LEVEL, name..., id, date, date
    if len(parts) >= 5 and parts[0] == "Rw":
        lake_id   = parts[-3]   # third from end (before two date tokens)
        lake_name = "_".join(parts[3:-3])
        return lake_id, lake_name
    return stem, stem


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():

    parser = argparse.ArgumentParser(
        description="Detect whiting events in lakes_cci time series"
    )
    parser.add_argument("--input-folder",  required=True,
                        help="Folder containing per-lake Rw CSV files")
    parser.add_argument("--output-folder", required=True,
                        help="Destination for per-lake and summary outputs")
    parser.add_argument("--band",          default="BGR",
                        help="Band to analyse  [BGR]")
    parser.add_argument("--value-col",     default="mean",
                        help="Statistic column to use  [mean]")
    parser.add_argument("--z-thresh",      type=float, default=5.0,
                        help="MAD z-score threshold  [5]")
    parser.add_argument("--bright-thresh", type=float, default=0.10,
                        help="Bright-pixel fraction threshold  [0.10]")
    parser.add_argument("--min-obs",       type=int,   default=3,
                        help="Min flagged observations per event  [3]")
    parser.add_argument("--max-gap-days",  type=int,   default=10,
                        help="Max gap (days) to merge flags  [10]")
    parser.add_argument("--min-duration",  type=int,   default=3,
                        help="Min event duration in days  [3]")
    parser.add_argument("--seasonal",      action="store_true",
                        help="Use per-month baseline (recommended for "
                             "lakes with regular seasonal whitings)")
    args = parser.parse_args()

    input_folder  = Path(args.input_folder)
    output_folder = Path(args.output_folder)
    output_folder.mkdir(exist_ok=True, parents=True)

    files = sorted(input_folder.glob("*.csv"))
    if not files:
        print(f"No CSV files found in {input_folder}")
        return

    summary_rows: list[dict] = []

    for f in files:
        print(f"→ {f.name}")
        lake_id, lake_name = parse_lake_info(f.stem)

        try:
            df = pd.read_csv(f)
        except Exception as e:
            print(f"  [WARN] Could not read file: {e}")
            continue

        if "date" not in df.columns:
            print("  [SKIP] No 'date' column")
            continue

        # ── Filter to requested band ──────────────────────────────────
        if "band" in df.columns:
            df_band = df[df["band"] == args.band].copy()
        else:
            df_band = df.copy()   # single-band file

        if df_band.empty or args.value_col not in df_band.columns:
            print(f"  [SKIP] Band '{args.band}' or column '{args.value_col}' not found")
            continue

        # ── Detect ────────────────────────────────────────────────────
        df_out, events = detect_whiting_events(
            df_band,
            value_col     = args.value_col,
            z_thresh      = args.z_thresh,
            bright_thresh = args.bright_thresh,
            min_obs       = args.min_obs,
            max_gap_days  = args.max_gap_days,
            min_duration  = args.min_duration,
            seasonal      = args.seasonal,
        )

        n_events = len(events)
        print(f"  {n_events} event(s) detected")

        # ── Per-lake outputs ──────────────────────────────────────────
        base    = f.stem
        ts_out  = output_folder / f"{base}_with_events.csv"
        ev_out  = output_folder / f"{base}_events.csv"

        df_out.to_csv(ts_out,  index=False)
        events.to_csv(ev_out, index=False)

        # ── Accumulate summary ────────────────────────────────────────
        summary_rows.append(
            compute_lake_summary(
                lake_id   = lake_id,
                lake_name = lake_name,
                df_ts     = df_out,
                events    = events,
                value_col = args.value_col,
            )
        )

    # ── Global summary table ──────────────────────────────────────────
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows).sort_values(
            ["n_events", "max_peak_value"], ascending=[False, False]
        )
        summary_path = output_folder / "whiting_summary.csv"
        summary_df.to_csv(summary_path, index=False)
        print(f"\n✓ Summary table saved → {summary_path}  ({len(summary_df)} lakes)")

        # Quick console overview
        print("\n── Top 10 lakes by event count ────────────────────────────────")
        cols_show = ["lake_name", "n_events", "events_per_year",
                     "max_peak_value", "max_peak_z", "peak_bright_frac",
                     "mean_event_duration", "peak_month", "dominant_flag_source"]
        cols_show = [c for c in cols_show if c in summary_df.columns]
        print(summary_df.head(10)[cols_show].to_string(index=False))
    else:
        print("No data processed.")

    print("\nDone.")


if __name__ == "__main__":
    main()
