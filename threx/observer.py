# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Where the body appears in the sky as seen from Earth.

APPROXIMATION (stated on purpose): SPHEREx is in low Earth orbit (roughly 700 km up). We put
the observer at the centre of the Earth. At 400 AU the largest error this can cause is
about (6400 + 700) km / 400 AU = 0.02 arcsec, far below one 6.2 arcsec SPHEREx pixel.
"""
import numpy as np
from astropy.coordinates import get_body_barycentric, solar_system_ephemeris
from astropy.time import Time

from .constants import ARCSEC_PER_RAD, C_AU_PER_DAY
from .orbits import heliocentric_ecliptic, ecliptic_to_equatorial, equatorial_to_ecliptic


def to_time(t):
    """Accept an astropy Time, an ISO string (UTC), or a Julian Date number (UTC)."""
    if isinstance(t, Time):
        return t
    arr = np.asarray(t)
    if arr.dtype.kind in "US":
        return Time(t, scale="utc")
    return Time(t, format="jd", scale="utc")


def earth_heliocentric(t):
    """Centre of the Earth relative to the Sun, equatorial axes, AU. Uses the ephemeris
    built into astropy (no download); it is good to a few km, far better than we need."""
    t = to_time(t)
    with solar_system_ephemeris.set("builtin"):
        earth = get_body_barycentric("earth", t).xyz.to_value("au")
        sun = get_body_barycentric("sun", t).xyz.to_value("au")
    return np.moveaxis(earth - sun, 0, -1)          # shape (..., 3)


def body_heliocentric(el, t):
    """Body relative to the Sun, equatorial axes, AU."""
    return ecliptic_to_equatorial(heliocentric_ecliptic(el, to_time(t).tdb.jd))


def geocentric_vector(el, t, light_time=True):
    """Vector Earth -> body (equatorial axes, AU). With light_time=True we use the position
    the body had when the light left it (what JPL Horizons calls 'astrometric')."""
    t = to_time(t)
    jd = t.tdb.jd
    earth = earth_heliocentric(t)
    vec = ecliptic_to_equatorial(heliocentric_ecliptic(el, jd)) - earth
    if light_time:
        for _ in range(2):                           # two iterations are plenty
            tau = np.linalg.norm(vec, axis=-1) / C_AU_PER_DAY
            vec = ecliptic_to_equatorial(heliocentric_ecliptic(el, jd - tau)) - earth
    return vec


def radec(vec):
    """RA, Dec in degrees from an equatorial vector."""
    x, y, z = np.moveaxis(np.asarray(vec, dtype=float), -1, 0)
    return np.degrees(np.arctan2(y, x)) % 360.0, np.degrees(np.arctan2(z, np.hypot(x, y)))


def ecliptic_lonlat(vec):
    """Ecliptic longitude and latitude in degrees from an equatorial vector."""
    return radec(equatorial_to_ecliptic(vec))


def unit_vector(ra_deg, dec_deg):
    """Equatorial unit vector for a sky position given in degrees."""
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)


def separation_arcsec(v1, v2):
    """Angle between two direction vectors, in arcsec (numerically stable for small angles)."""
    v1 = np.asarray(v1, dtype=float)
    v2 = np.asarray(v2, dtype=float)
    cross = np.linalg.norm(np.cross(v1, v2), axis=-1)
    dot = np.sum(v1 * v2, axis=-1)
    return np.arctan2(cross, dot) * ARCSEC_PER_RAD


def apparent_position(el, t):
    """Dict with RA/Dec, ecliptic lon/lat (deg) and distance from Earth (AU)."""
    vec = geocentric_vector(el, t)
    ra, dec = radec(vec)
    lon, lat = ecliptic_lonlat(vec)
    return {"ra": ra, "dec": dec, "lon": lon, "lat": lat, "delta": np.linalg.norm(vec, axis=-1)}
