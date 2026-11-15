# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Find the SPHEREx exposures that cover each challenge target, using the AWS mirror.

  python scripts/scan_targets.py            # all targets  -> data_cache/_work/hits_<slug>.json
  python scripts/scan_targets.py Barnard    # one target (name prefix)
Target positions: SIMBAD (CDS, via the CfA mirror), ICRS, epoch J2000, with proper motion;
moved to the epoch of 2026.0 before searching (the field moves by up to a few arcmin).
"""
import datetime as dt
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from astropy import units as u
from astropy.time import Time

from threx import aws_mirror as aws
from threx.observer import earth_heliocentric, unit_vector, separation_arcsec
from threx.spherex_io import _slug

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "data_cache" / "_work"


def position_at(t, year):
    """Proper-motion-corrected (ra, dec) in degrees at a decimal year (linear motion)."""
    dy = year - 2000.0
    dec = t["dec"] + t["pmdec"] * dy / 3.6e6
    ra = t["ra"] + t["pmra"] * dy / 3.6e6 / np.cos(np.radians(t["dec"]))
    return ra, dec


def candidate_periods(t, periods, tol=14.0):
    ra, dec = position_at(t, 2026.0)
    days = Time(dt.date(2025, 4, 1).isoformat()) + np.arange(0, 600) * u.day
    elong = separation_arcsec(unit_vector(ra, dec), -earth_heliocentric(days)) / 3600
    good = {d.date() for d, e in zip(days.to_datetime(), elong) if abs(e - 90) < tol}
    keep = []
    for rel, name, w0, prefix in periods:
        if any(w0 + dt.timedelta(days=k) in good for k in range(-1, 8)):
            keep.append((rel, name, w0, prefix))
    return keep


if __name__ == "__main__":
    targets = json.loads((ROOT / "data_cache" / "targets_simbad.json").read_text())["targets"]
    only = sys.argv[1].lower() if len(sys.argv) > 1 else ""
    WORK.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    periods = aws.list_periods()
    print(len(periods), "planning periods on the mirror", round(time.time() - t0), "s")
    for name, t in targets.items():
        if only and not name.lower().startswith(only):
            continue
        ra, dec = position_at(t, 2026.0)
        keep = candidate_periods(t, periods)
        print(f"{name}: RA {ra:.5f} Dec {dec:.5f} (2026.0); {len(keep)} candidate periods")
        hits = []
        for rel, pname, w0, prefix in keep:
            h = aws.scan_period(prefix, ra, dec)
            for x in h:
                x.update(release=rel, period=pname)
            hits += h
            print(f"  {pname}: {len(h)} hits ({round(time.time() - t0)} s)")
        (WORK / f"hits_{_slug(name)}.json").write_text(json.dumps(hits, indent=1))
        print("  total", len(hits), "exposures")
