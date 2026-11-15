"""Exploration helper: check the sb_ident call on a known asteroid (not used by the app)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import requests
from astropy.time import Time

from threx import jpl

t = Time("2026-03-01T00:00:00", scale="utc")
eph = jpl.horizons_observer("4;", [t.jd])[0]
print("Vesta from Horizons:", eph)
for observer in ({"lat": "0", "lon": "0", "alt": "-6378.137"}, {"mpc-code": "500"}, {"mpc-code": "568"}):
    params = {"sb-kind": "a", "obs-time": t.strftime("%Y-%m-%d_%H:%M:%S"),
              "fov-ra-center": jpl.sexagesimal(eph["ra"], hours=True), "fov-dec-center": jpl.sexagesimal(eph["dec"]),
              "fov-ra-hwidth": "0.05", "fov-dec-hwidth": "0.05", "two-pass": "true", "mag-required": "true", "vmag-lim": "16", **observer}
    r = requests.get(jpl.SB_IDENT_URL, params=params, timeout=420)
    print(observer, r.status_code, r.text[:1200])
    print()
