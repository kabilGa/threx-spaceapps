"""Exploration helper: what does IRSA return for one sky position?  (not used by the app)

Run:  python scripts/explore_sia.py 60 20
"""
import sys
import time

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroquery.ipac.irsa import Irsa

ra, dec = (float(sys.argv[1]), float(sys.argv[2])) if len(sys.argv) > 2 else (60.0, 20.0)
print(Irsa.list_collections(servicetype="SIA", filter="spherex"))
t0 = time.time()
res = Irsa.query_sia(pos=(SkyCoord(ra, dec, unit="deg"), 1 * u.arcsec), collection=["spherex_dr1", "spherex_qr3"])
print("rows:", len(res), "in", round(time.time() - t0, 1), "s")
print(res.colnames)
res.sort("t_min")
t = Time(res["t_min"], format="mjd")
print("first", t[0].iso, "last", t[-1].iso)
for row in res[:: max(1, len(res) // 25)]:
    print(Time(row["t_min"], format="mjd").iso[:19], row["obs_collection"], row["energy_bandpassname"],
          f'{row["em_min"] * 1e6:.2f}-{row["em_max"] * 1e6:.2f}um', row["t_exptime"], row["access_estsize"],
          row["access_url"].split("/")[-1])
print(np.unique(res["obs_collection"], return_counts=True))
