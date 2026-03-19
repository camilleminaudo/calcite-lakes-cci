#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Detect whiting events from lake time series (BGR / Rw)

Features:
- Robust z-score (MAD-based)
- Threshold detection
- Merge events with temporal gaps
- Event duration filtering
- Outputs per-lake CSVs

How to run:
python detect_whiting.py \
  --input-folder ./results \
  --output-folder ./events \
  --band BGR \
  --value-col mean \
  --z-thresh 5 \
  --min-obs 3 \
  --max-gap-days 10

"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse


# ----------------------------
# CORE FUNCTION
# ----------------------------
def detect_whiting_events(df,
                          value_col="mean",
                          date_col="date",
                          z_thresh=5,
                          min_obs=3,
                          max_gap_days=10):

    df = df.sort_values(date_col).copy()
    df[date_col] = pd.to_datetime(df[date_col])

    values = df[value_col].values

    # ----------------------------
    # Robust stats
    # ----------------------------
    med = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - med))

    if mad == 0:
        mad = 1e-6

    df["z"] = (df[value_col] - med) / mad
    df["event_flag"] = df["z"] > z_thresh

    # ----------------------------
    # Keep only candidate points
    # ----------------------------
    df_events = df[df["event_flag"]].copy()

    if df_events.empty:
        df["event_id"] = np.nan
        df["valid_event"] = False
        return df, pd.DataFrame()

    df_events = df_events.sort_values(date_col)

    # ----------------------------
    # Compute gaps
    # ----------------------------
    df_events["gap"] = df_events[date_col].diff().dt.days
    df_events["new_event"] = (df_events["gap"].isna()) | (df_events["gap"] > max_gap_days)
    df_events["event_id"] = df_events["new_event"].cumsum()

    # ----------------------------
    # Summarise events
    # ----------------------------
    events = (
        df_events
        .groupby("event_id")
        .agg(
            start_date=(date_col, "min"),
            end_date=(date_col, "max"),
            n_obs=(date_col, "count"),
            span_days=(date_col, lambda x: (x.max() - x.min()).days + 1),
            peak_value=(value_col, "max"),
            mean_value=(value_col, "mean")
        )
        .reset_index()
    )

    # Filter valid events
    events = events[events["n_obs"] >= min_obs]

    # ----------------------------
    # Tag back to full dataset
    # ----------------------------
    df_events["valid_event"] = df_events["event_id"].isin(events["event_id"])

    df = df.merge(
        df_events[[date_col, "event_id", "valid_event"]],
        on=date_col,
        how="left"
    )

    df["valid_event"] = df["valid_event"].fillna(False)

    return df, events


# ----------------------------
# MAIN
# ----------------------------
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--input-folder", required=True)
    parser.add_argument("--output-folder", required=True)
    parser.add_argument("--band", default="BGR")
    parser.add_argument("--value-col", default="mean")
    parser.add_argument("--z-thresh", type=float, default=5)
    parser.add_argument("--min-obs", type=int, default=3)
    parser.add_argument("--max-gap-days", type=int, default=10)

    args = parser.parse_args()

    input_folder = Path(args.input_folder)
    output_folder = Path(args.output_folder)
    output_folder.mkdir(exist_ok=True, parents=True)

    files = list(input_folder.glob("*.csv"))

    for f in files:
        print(f"Processing {f.name}")

        df = pd.read_csv(f)

        if "date" not in df.columns:
            print("Skipping (no date column)")
            continue

        # Filter band
        df_band = df[df["band"] == args.band].copy()

        if df_band.empty:
            print("Skipping (no band data)")
            continue

        df_out, events = detect_whiting_events(
            df_band,
            value_col=args.value_col,
            z_thresh=args.z_thresh,
            min_obs=args.min_obs,
            max_gap_days=args.max_gap_days
        )

        # ----------------------------
        # OUTPUT FILES
        # ----------------------------
        base = f.stem

        ts_out = output_folder / f"{base}_with_events.csv"
        ev_out = output_folder / f"{base}_events.csv"

        df_out.to_csv(ts_out, index=False)
        events.to_csv(ev_out, index=False)

    print("Done.")


if __name__ == "__main__":
    main()