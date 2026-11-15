# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Reproduce the back-of-envelope numbers of the project brief (section 6.5).

Run:  python scripts/sanity_checks.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from astropy import units as u
from astropy.time import Time

from threx import orbits, signature, photometry
from threx.observer import apparent_position

T0 = Time("2026-01-01T00:00:00", scale="utc")
el = orbits.circular_orbit(400.0, inc=0.0, epoch_jd=T0.tdb.jd)

# 1. brightness
m = photometry.reflected_mag(400, 400, 3.0, 0.4)
print(f"1. V-like AB magnitude at r = Delta = 400 AU, R = 3 R_earth, p = 0.4 : {m:.2f}  (expected ~19.7)")

# 2-4. motion: scan one year
yr = signature.rate_through_year(el, T0, n=1462)
i_opp = int(np.argmax(yr["elongation_deg"]))
t_opp = T0 + yr["days"][i_opp] * u.day
r = signature.motion_rate(el, t_opp)
print(f"3. opposition on {t_opp.iso[:10]}: {float(r['arcsec_per_day']):.2f} arcsec/day = "
      f"{float(r['arcsec_per_hour']):.3f} arcsec/hour = {float(r['pixels_per_day']):.2f} pixels/day, "
      f"retrograde = {bool(r['retrograde'])}  (expected ~9, ~0.36, ~1.4, True)")

s = signature.shift_between(el, t_opp - 91.3 * u.day, t_opp + 91.3 * u.day)
print(f"2. shift over 6 months centred on opposition: total {s['total_arcsec']:.0f} arcsec = "
      f"{s['total_pixels']:.0f} pixels; parallax part {s['parallax_arcsec']:.0f} arcsec "
      f"(expected ~1000 arcsec, ~170 pixels)")
print(f"4. orbital part over the same 6 months: {s['orbital_arcsec']:.0f} arcsec (expected ~80)")

pos = apparent_position(el, t_opp - 91.3 * u.day)
d = signature.distance_from_shift(s["total_arcsec"], t_opp - 91.3 * u.day, t_opp + 91.3 * u.day,
                                  float(pos["ra"]), float(pos["dec"]))
d2 = signature.distance_from_shift_refined(s["total_arcsec"], t_opp - 91.3 * u.day, t_opp + 91.3 * u.day,
                                           float(pos["ra"]), float(pos["dec"]))
print(f"5. distance recovered from that shift: first order {d['distance_au']:.0f} AU "
      f"(+-{d['d_err_meas_au']:.1f} AU for a 1-pixel error, systematic up to {100 * d['sys_frac']:.0f}%), "
      f"refined {d2:.0f} AU  (truth 400 AU)")

for lam in (0.55, 1.24, 1.65, 2.16, 3.4, 4.6):
    print(f"   Sun AB mag at 1 AU, {lam:4.2f} um: Willmer table {float(photometry.sun_ab_mag(lam)):.2f}, "
          f"5772 K blackbody {float(photometry.sun_ab_mag_blackbody(lam)):.2f}")
for T in (10, 30, 50):
    print(f"   thermal AB mag at 5 um, T = {T} K, R = 3 R_earth, 400 AU: "
          f"{float(photometry.thermal_mag(T, 3.0, 400.0, 5.0)):.0f}")
