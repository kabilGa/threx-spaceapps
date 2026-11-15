"""Compare our orbit code with JPL Horizons for real distant bodies (Sedna, Eris).

Input of our code : the heliocentric osculating elements Horizons gives at one middle epoch.
Truth             : Horizons geocentric astrometric RA/Dec, rates and distances at 13 dates.
The data file is REAL JPL output saved by scripts/fetch_horizons_validation.py.
"""
import json
from pathlib import Path

import numpy as np
from astropy.time import Time

from .observer import geocentric_vector, radec, unit_vector, separation_arcsec
from .orbits import OrbitalElements
from .signature import motion_rate

DATA_FILE = Path(__file__).resolve().parents[1] / "tests" / "data" / "horizons_validation.json"


def load(path=DATA_FILE):
    return json.loads(Path(path).read_text())


def compare(path=DATA_FILE):
    """Return {body: dict of residual arrays and summary numbers}."""
    data = load(path)
    results = {"downloaded_utc": data["downloaded_utc"], "bodies": {}}
    for name, body in data["bodies"].items():
        el = OrbitalElements(**body["elements"])
        eph = body["ephemeris"]
        t = Time([row["jd"] for row in eph], format="jd", scale="utc")
        ra_h = np.array([row["ra"] for row in eph])
        dec_h = np.array([row["dec"] for row in eph])
        vec = geocentric_vector(el, t)
        ra, dec = radec(vec)
        truth = unit_vector(ra_h, dec_h)
        pos_err = separation_arcsec(vec, truth)                         # arcsec
        # apparent rate: ours vs Horizons (dRA*cosDec, dDec in arcsec/hour)
        rate_h = np.hypot([row["dra_cosd"] for row in eph], [row["ddec"] for row in eph])
        rate = motion_rate(el, t)["arcsec_per_hour"]
        # shift between first and last date (360 days) and between dates 180 days apart
        shift_h = separation_arcsec(truth[:-6], truth[6:])
        shift = separation_arcsec(vec[:-6], vec[6:])
        delta_h = np.array([row["delta"] for row in eph])
        results["bodies"][name] = {
            "jd": t.jd, "ra": ra, "dec": dec, "ra_horizons": ra_h, "dec_horizons": dec_h,
            "pos_err_arcsec": pos_err, "rate": rate, "rate_horizons": rate_h,
            "shift_180d": shift, "shift_180d_horizons": shift_h,
            "delta_err_au": np.linalg.norm(vec, axis=-1) - delta_h,
            "max_pos_err_arcsec": float(pos_err.max()),
            "max_rate_err_arcsec_per_hour": float(np.abs(rate - rate_h).max()),
            "max_shift_err_arcsec": float(np.abs(shift - shift_h).max()),
            "max_shift_rel_err": float(np.abs(shift / shift_h - 1).max()),
            "distance_au": float(delta_h.mean()),
        }
    return results


if __name__ == "__main__":
    res = compare()
    print("Horizons data downloaded:", res["downloaded_utc"])
    for name, r in res["bodies"].items():
        print(f"{name}: distance ~{r['distance_au']:.1f} AU | max position error "
              f"{r['max_pos_err_arcsec']:.2f} arcsec | max rate error "
              f"{r['max_rate_err_arcsec_per_hour']:.4f} arcsec/h | max 180-day shift error "
              f"{r['max_shift_err_arcsec']:.2f} arcsec ({100 * r['max_shift_rel_err']:.3f} %)")
        print("   position errors (arcsec) by date:", np.round(r["pos_err_arcsec"], 2))
