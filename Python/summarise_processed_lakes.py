#!/usr/bin/env python3
"""
Scan an output folder from lakes_cci_extractor.py and build a summary CSV
of all lakes already processed, based on filenames only.

Expected filename pattern:
    Rw_<SENSOR>_<LEVEL>_<lake_name>_<lake_id>_<YYYYMMDD>_<YYYYMMDD>.csv

Usage:
    python summarise_processed_lakes.py \
        --output-folder ./results \
        --summary-file  ./processed_lakes.csv
"""

import argparse
import re
import pandas as pd
from pathlib import Path

# Rw_OLCI_L2_Caspian Sea_GLWD00000001_20160401_20230101.csv
PATTERN = re.compile(
    r"^Rw_"
    r"(?P<sensor>[^_]+)_"
    r"(?P<level>L[23])_"
    r"(?P<lake_name>.+)_"
    r"(?P<lake_id>(?:GLWD|HYLA)\d+)_"
    r"(?P<date_start>\d{8})_"
    r"(?P<date_end>\d{8})"
    r"\.csv$"
)

parser = argparse.ArgumentParser()
parser.add_argument("--output-folder", required=True,
                    help="Folder containing extracted lake CSVs")
parser.add_argument("--summary-file",  default="processed_lakes.csv",
                    help="Path for the output summary CSV  [processed_lakes.csv]")
args = parser.parse_args()

rows = []
unmatched = []

for f in sorted(Path(args.output_folder).glob("*.csv")):
    m = PATTERN.match(f.name)
    if not m:
        unmatched.append(f.name)
        continue
    rows.append({
        "lake_name":  m["lake_name"],
        "lake_id":    m["lake_id"],
        "sensor":     m["sensor"],
        "level":      m["level"],
        "date_start": pd.to_datetime(m["date_start"], format="%Y%m%d").date(),
        "date_end":   pd.to_datetime(m["date_end"],   format="%Y%m%d").date(),
        "filename":   f.name,
    })

if not rows:
    print("No matching files found.")
else:
    df = pd.DataFrame(rows).sort_values(["sensor", "level", "lake_name"])
    df.to_csv(args.summary_file, index=False)
    print(f"✓ {len(df)} lake(s) summarised → {args.summary_file}")

if unmatched:
    print(f"  {len(unmatched)} file(s) skipped (pattern mismatch):")
    for name in unmatched:
        print(f"    {name}")
