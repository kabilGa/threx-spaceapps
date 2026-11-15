# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Find SPHEREx exposures that contain a KNOWN asteroid (Horizons positions) on the AWS mirror.

  python scripts/search_asteroid_aws.py 6 2025-11-25 2025-12-07     # id, first day, last day (UTC)
Writes data_cache/_work/asteroid_<id>_window0.csv, which scripts/build_asteroid_proof.py build reads.
"""
import datetime as dt
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from astropy.time import Time

from threx import aws_mirror as aws, jpl

ast, d0, d1 = sys.argv[1], dt.date.fromisoformat(sys.argv[2]), dt.date.fromisoformat(sys.argv[3])
jds = np.arange(Time(d0.isoformat()).jd, Time(d1.isoformat()).jd + 1, 1 / 24)
eph = []
for chunk in np.array_split(jds, int(np.ceil(len(jds) / 40))):
    eph += jpl.horizons_observer(f"{ast};", list(chunk))
mjd = np.array([e["jd"] for e in eph]) - 2400000.5
ra = np.degrees(np.unwrap(np.radians([e["ra"] for e in eph])))
dec = np.array([e["dec"] for e in eph])
pos = lambda m: (float(np.interp(m, mjd, ra)) % 360, float(np.interp(m, mjd, dec)))
periods = [p for p in aws.list_periods()
           if d0 - dt.timedelta(days=8) <= p[2] <= d1 + dt.timedelta(days=1)]
print(len(periods), "candidate period directories")
hits = []
for rel, name, w0, prefix in periods:
    h = aws.scan_period(prefix, pos, 0.0, stride=4)
    for x in h:
        x.update(release=rel, period=name)
    hits += h
    print(" ", name, len(h), "hits", flush=True)
hits = [h for h in hits if mjd[0] <= h["mjd"] <= mjd[-1]]
hits.sort(key=lambda h: (h["mjd"], h["det"]))
rows = []
for chunk in np.array_split(np.arange(len(hits)), max(1, int(np.ceil(len(hits) / 40)))):
    if len(chunk) == 0:
        continue
    e = jpl.horizons_observer(f"{ast};", [hits[i]["mjd"] + 2400000.5 for i in chunk])
    for i, ev in zip(chunk, e):
        h = hits[i]
        rows.append({"obs_id": h["key"].split("/")[-1][:-5], "obs_collection": h["release"], "t_mid": h["mjd"],
                     "time_utc": Time(h["mjd"], format="mjd", scale="utc").iso, "detector": h["det"],
                     "access_url": aws.BUCKET + h["key"], "key": h["key"], "pred_ra": ev["ra"], "pred_dec": ev["dec"],
                     "pred_dra_cosd": ev["dra_cosd"], "pred_ddec": ev["ddec"], "pred_vmag": ev["vmag"],
                     "pred_delta": ev["delta"], "pred_r": ev["r"], "pred_elong": ev["elong"]})
out = Path(__file__).resolve().parents[1] / "data_cache" / "_work" / f"asteroid_{ast}_window0.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print(len(rows), "exposure/detector hits ->", out.name)
if rows:
    df = pd.DataFrame(rows)
    print(df.groupby(df["time_utc"].str[:13]).size().to_string())
