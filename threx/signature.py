"""The "signature" of a distant body: how fast and how far it appears to move.

Plain words: a body hundreds of AU away barely moves on its own orbit. What we see is mostly
the reflection of OUR motion: Earth circles the Sun, so the body seems to swing back and forth
once a year (parallax). The closer the body, the bigger the swing. So the size of the swing
is a distance measurement.
"""
import numpy as np
from astropy import units as u

from .constants import ARCSEC_PER_RAD, GAUSS_K, SPHEREX_PIXEL_ARCSEC, OBLIQUITY_DEG
from .observer import (to_time, earth_heliocentric, body_heliocentric, geocentric_vector,
                       ecliptic_lonlat, separation_arcsec, unit_vector)


def motion_rate(el, t, dt_hours=1.0):
    """Apparent motion rate at time t (finite difference over dt_hours).

    Returns arcsec/hour, arcsec/day, SPHEREx pixels/day, and 'retrograde' (True when the
    ecliptic longitude is decreasing, i.e. the body moves westward against the stars).
    """
    t = to_time(t)
    half = 0.5 * dt_hours * u.hour
    v1 = geocentric_vector(el, t - half)
    v2 = geocentric_vector(el, t + half)
    sep = separation_arcsec(v1, v2)
    lon1, _ = ecliptic_lonlat(v1)
    lon2, _ = ecliptic_lonlat(v2)
    dlon = (lon2 - lon1 + 180.0) % 360.0 - 180.0
    per_hour = sep / dt_hours
    return {
        "arcsec_per_hour": per_hour,
        "arcsec_per_day": per_hour * 24.0,
        "pixels_per_day": per_hour * 24.0 / SPHEREX_PIXEL_ARCSEC,
        "retrograde": dlon < 0,
        "signed_arcsec_per_day": np.where(dlon < 0, -1.0, 1.0) * per_hour * 24.0,
    }


def solar_elongation_deg(el, t):
    """Angle Sun - Earth - body in degrees (180 = opposition, 90 = quadrature)."""
    t = to_time(t)
    to_body = geocentric_vector(el, t)
    to_sun = -earth_heliocentric(t)
    return separation_arcsec(to_body, to_sun) / 3600.0


def rate_through_year(el, t0, n=366):
    """Motion rate and solar elongation sampled over one year starting at t0."""
    t0 = to_time(t0)
    days = np.linspace(0.0, 365.25, n)
    t = t0 + days * u.day
    out = motion_rate(el, t)
    out["days"] = days
    out["elongation_deg"] = solar_elongation_deg(el, t)
    return out


def quadrature_times(el, t0, n=1462):
    """Dates within one year of t0 when the body is 90 deg from the Sun.

    Why it matters: SPHEREx scans along a great circle about 90 deg from the Sun
    (spherex.caltech.edu/page/survey), so it visits any sky position about every six
    months, exactly around these dates."""
    t0 = to_time(t0)
    days = np.linspace(0.0, 365.25, n)
    elong = solar_elongation_deg(el, t0 + days * u.day) - 90.0
    idx = np.where(np.sign(elong[:-1]) != np.sign(elong[1:]))[0]
    # linear interpolation of each zero crossing
    cross = [days[i] - elong[i] * (days[i + 1] - days[i]) / (elong[i + 1] - elong[i]) for i in idx]
    return [t0 + d * u.day for d in cross]


def opposition_time(el, t0, n=1462):
    """Date within one year of t0 when the body is opposite to the Sun (largest elongation)."""
    t0 = to_time(t0)
    days = np.linspace(0.0, 365.25, n)
    return t0 + days[int(np.argmax(solar_elongation_deg(el, t0 + days * u.day)))] * u.day


def sky_track(el, t0, n_days=365.25, n=400):
    """Apparent path on the sky (ecliptic lon/lat, RA/Dec in deg) over n_days from t0."""
    t0 = to_time(t0)
    days = np.linspace(0.0, n_days, n)
    vec = geocentric_vector(el, t0 + days * u.day)
    lon, lat = ecliptic_lonlat(vec)
    return {"days": days, "lon": lon, "lat": lat, "delta": np.linalg.norm(vec, axis=-1)}


def shift_between(el, t1, t2):
    """Apparent shift between two epochs, split into its two causes.

    total    : what an observer on Earth measures.
    parallax : only Earth moves (the body is frozen where it was at t1).
    orbital  : only the body moves (seen from the Sun).
    The two parts add as vectors on the sky, so total is not exactly parallax + orbital.
    """
    t1, t2 = to_time(t1), to_time(t2)
    total = separation_arcsec(geocentric_vector(el, t1), geocentric_vector(el, t2))
    p1, p2 = body_heliocentric(el, t1), body_heliocentric(el, t2)
    e1, e2 = earth_heliocentric(t1), earth_heliocentric(t2)
    parallax = separation_arcsec(p1 - e1, p1 - e2)
    orbital = separation_arcsec(p1, p2)
    px = SPHEREX_PIXEL_ARCSEC
    return {
        "total_arcsec": float(total), "parallax_arcsec": float(parallax),
        "orbital_arcsec": float(orbital),
        "total_pixels": float(total / px), "parallax_pixels": float(parallax / px),
        "orbital_pixels": float(orbital / px),
        "baseline_days": float((t2 - t1).to_value("day")),
    }


def earth_baseline_perp_au(t1, t2, ra_deg, dec_deg):
    """Length (AU) of the part of Earth's displacement between t1 and t2 that is
    perpendicular to the line of sight (ra, dec). Only this part produces parallax."""
    b = earth_heliocentric(t2) - earth_heliocentric(t1)
    los = unit_vector(ra_deg, dec_deg)
    b_perp = b - np.sum(b * los, axis=-1)[..., None] * los
    return float(np.linalg.norm(b_perp, axis=-1))


def distance_from_shift(shift_arcsec, t1, t2, ra_deg, dec_deg, sigma_arcsec=SPHEREX_PIXEL_ARCSEC):
    """First-order distance from a measured shift:  d = B_perp / theta.

    Assumptions: (1) the whole shift is parallax (the body's own motion is ignored);
    (2) small angles. Errors returned:
      d_err_meas : from the position measurement error sigma_arcsec (default one pixel),
                   d_err/d = sigma/theta.
      sys_frac   : how wrong assumption (1) can be: orbital shift of a circular orbit at
                   distance d divided by the measured shift.
    """
    b_perp = earth_baseline_perp_au(t1, t2, ra_deg, dec_deg)
    theta = shift_arcsec / ARCSEC_PER_RAD
    d = b_perp / theta
    dt_days = abs(float((to_time(t2) - to_time(t1)).to_value("day")))
    orbital_arcsec = GAUSS_K / d ** 1.5 * dt_days * ARCSEC_PER_RAD
    return {
        "distance_au": d, "b_perp_au": b_perp,
        "d_err_meas_au": d * sigma_arcsec / shift_arcsec,
        "sys_frac": orbital_arcsec / shift_arcsec,
    }


def _model_shift_arcsec(d, t1, t2, ra_deg, dec_deg):
    """Forward model used by the refined estimate: body at distance d from Earth at t1 in
    direction (ra, dec), on a circular prograde orbit parallel to the ecliptic plane."""
    e1, e2 = earth_heliocentric(t1), earth_heliocentric(t2)
    p1 = e1 + d * unit_vector(ra_deg, dec_deg)
    r = np.linalg.norm(p1)
    eps = np.radians(OBLIQUITY_DEG)
    ecl_pole = np.array([0.0, -np.sin(eps), np.cos(eps)])     # ecliptic north pole, equatorial axes
    vdir = np.cross(ecl_pole, p1)
    vdir /= np.linalg.norm(vdir)
    dt_days = float((to_time(t2) - to_time(t1)).to_value("day"))
    p2 = p1 + vdir * (GAUSS_K / np.sqrt(r)) * dt_days         # circular speed k/sqrt(r), AU/day
    return float(separation_arcsec(p1 - e1, p2 - e2))


def track_from_position(ra_deg, dec_deg, distance_au, t_ref, times):
    """'What if a planet sat exactly HERE?'  Put a body at distance_au from Earth in the
    direction (ra, dec) at time t_ref, on a circular prograde orbit parallel to the ecliptic,
    and return its apparent offset from that position at each of `times`.

    Returns dict: offset_arcsec (total), offset_pixels, and the direction vectors.
    Used with the REAL SPHEREx exposure dates of a field."""
    t_ref, times = to_time(t_ref), to_time(times)
    e_ref = earth_heliocentric(t_ref)
    p_ref = e_ref + distance_au * unit_vector(ra_deg, dec_deg)
    r = np.linalg.norm(p_ref)
    eps = np.radians(OBLIQUITY_DEG)
    ecl_pole = np.array([0.0, -np.sin(eps), np.cos(eps)])
    vdir = np.cross(ecl_pole, p_ref)
    vdir /= np.linalg.norm(vdir)
    ang = GAUSS_K / r ** 1.5 * (times - t_ref).to_value("day")          # orbital angle [rad]
    ang = np.atleast_1d(ang)[:, None]
    p = r * (np.cos(ang) * (p_ref / r) + np.sin(ang) * vdir)             # exact circular motion
    vec = p - np.atleast_2d(earth_heliocentric(times))
    off = separation_arcsec(vec, p_ref - e_ref)
    return {"offset_arcsec": off, "offset_pixels": off / SPHEREX_PIXEL_ARCSEC, "vectors": vec}


def distance_from_shift_refined(shift_arcsec, t1, t2, ra_deg, dec_deg):
    """Distance that reproduces the measured shift when the body's own motion is included
    (assumed: circular, prograde, parallel to the ecliptic). Falls back to the first-order
    value if no solution is bracketed."""
    t1, t2 = to_time(t1), to_time(t2)
    f = lambda d: _model_shift_arcsec(d, t1, t2, ra_deg, dec_deg) - shift_arcsec
    lo, hi = 20.0, 20000.0
    f_lo, f_hi = f(lo), f(hi)
    if f_lo * f_hi > 0:                              # no solution bracketed
        return distance_from_shift(shift_arcsec, t1, t2, ra_deg, dec_deg)["distance_au"]
    for _ in range(60):                              # plain bisection: simple and safe
        mid = 0.5 * (lo + hi)
        f_mid = f(mid)
        if f_lo * f_mid <= 0:
            hi = mid
        else:
            lo, f_lo = mid, f_mid
    return 0.5 * (lo + hi)
