# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Build the cached viewer data (cutouts of real SPHEREx exposures) for one challenge target.

  python scripts/build_viewer_cache.py Barnard [--half 150] [--workers 8] [--max-frames 60]
Needs scripts/scan_targets.py to have written data_cache/_work/hits_<target>.json first.
Reads only the needed pixels from the public AWS mirror of the IRSA SPHEREx archive.
"""
import argparse
import socket

socket.setdefaulttimeout(30)
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from astropy.time import Time

from threx import aws_mirror as aws, blink, spherex_io
from scan_targets import position_at

ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument("target")
p.add_argument("--half", type=float, default=150.0)
p.add_argument("--workers", type=int, default=8)
p.add_argument("--max-frames", type=int, default=70)
a = p.parse_args()

tg = json.loads((ROOT / "data_cache" / "targets_simbad.json").read_text())
name = next(n for n in tg["targets"] if n.lower().startswith(a.target.lower()))
t = tg["targets"][name]
ra, dec = position_at(t, 2026.0)
# every exposure that covers the target: IRSA SIA2 search (all epochs, all detectors); pixels from the AWS mirror
df = spherex_io.query_exposures(ra, dec, radius_arcsec=1.0, retries=1).sort_values("t_mid").reset_index(drop=True)
if df.empty:
    sys.exit(f"{name}: IRSA returned no exposure")
print(name, "SIA2 exposures:", len(df), "from", df["time_utc"].iloc[0][:10], "to", df["time_utc"].iloc[-1][:10])
total = len(df)
if len(df) > a.max_frames:
    df = df.iloc[np.unique(np.round(np.linspace(0, len(df) - 1, a.max_frames)).astype(int))]
rows = pd.DataFrame({"obs_id": df["access_url"].str.split("/").str[-1].str.replace(".fits", "", regex=False),
                     "obs_collection": df["obs_collection"], "t_mid": df["t_mid"], "time_utc": df["time_utc"],
                     "access_url": df["access_url"], "key": df["access_url"].str.replace(aws.IRSA_DATA, "", regex=False)})
print(name, len(rows), "frames; field centre", round(ra, 5), round(dec, 5))
bundle = blink.build_bundle(
    name, ra, dec, rows, half_size_arcsec=a.half, workers=a.workers,
    reader=lambda row: aws.read_cutout_s3(row["key"], ra, dec, a.half),
    progress=lambda k, n, r: print(f"  {k + 1}/{n}", end="\r"),
    extra={"source": "SPHEREx Level 2 spectral images from the public IRSA cloud mirror (AWS bucket "
                     "nasa-irsa-spherex), read with HTTP range requests; exposures found with the IRSA SIA2 search",
           "sia2_exposures_total": total, "simbad": {k: t[k] for k in ("simbad_id", "ra", "dec", "pmra", "pmdec", "plx_mas")},
           "field_centre_epoch": 2026.0})
path = blink.save_bundle(bundle)
fr = bundle["meta"]["frames"]
df = pd.DataFrame([{
    "obs_id": f["obs_id"], "obs_collection": f["collection"], "t_mid": f["t_mid"], "time_utc": f["time_utc"],
    "lam_min_um": f["lam_um"] - f["dlam_um"] / 2, "lam_max_um": f["lam_um"] + f["dlam_um"] / 2,
    "detector": f["detector"], "access_url": f["access_url"]} for f in fr])
spherex_io.save_exposures(name, df, ra, dec)
print(f"\nsaved {path.name} {path.stat().st_size / 1e6:.2f} MB; frames {len(fr)}; failed {len(bundle['meta']['failed'])}",
      bundle["meta"]["failed"][:2])
