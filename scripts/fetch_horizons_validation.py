"""Download the JPL Horizons data used to validate our orbit code (section 6.6 of the brief).

For Sedna (90377) and Eris (136199) we store:
 * their heliocentric osculating elements at a middle epoch (input of OUR code),
 * Horizons geocentric astrometric RA/Dec, rates and distances on 13 dates (the truth).
Output: tests/data/horizons_validation.json  (REAL data from JPL, with the download date).

Run:  python scripts/fetch_horizons_validation.py
"""
import json
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from threx import jpl

EPOCH_JD = 2461090.5                                   # 2026-Feb-19 (TDB), middle of the test window
JDS = [EPOCH_JD + 30.0 * k for k in range(-6, 7)]      # 13 dates, +-180 days
BODIES = {"Sedna (90377)": "90377;", "Eris (136199)": "136199;"}

out = {"source": "NASA/JPL Horizons API, https://ssd.jpl.nasa.gov/api/horizons.api",
       "downloaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
       "note": "RA/Dec are geocentric astrometric ICRF (quantity 1); rates are quantity 3.",
       "bodies": {}}
for name, cmd in BODIES.items():
    el = jpl.horizons_elements(cmd, EPOCH_JD, label=name)
    eph = jpl.horizons_observer(cmd, JDS)
    out["bodies"][name] = {"elements": asdict(el), "ephemeris": eph}
    print(name, "a =", round(el.a, 2), "AU, e =", round(el.e, 4), "rows:", len(eph))

path = Path(__file__).resolve().parents[1] / "tests" / "data" / "horizons_validation.json"
path.write_text(json.dumps(out, indent=1))
print("wrote", path)
