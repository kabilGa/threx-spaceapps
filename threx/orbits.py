"""Two-body (Kepler) orbit of a body around the Sun.

Plain words: an orbit is an ellipse fixed in space. Six numbers describe it:
size (a), shape (e), tilt (i), orientation (node, argp) and where the body is at a
reference date (M0). From them we compute the position at any date.
"""
from dataclasses import dataclass, replace

import numpy as np

from .constants import (GAUSS_K, OBLIQUITY_DEG, P9_A_AU, P9_Q_AU, P9_INC_DEG,
                        DEMO_NODE_DEG, DEMO_ARGP_DEG)

J2000_JD = 2451545.0


@dataclass(frozen=True)
class OrbitalElements:
    """Heliocentric elements, ecliptic J2000 frame. Angles in degrees, a in AU."""
    a: float            # semi-major axis [AU]
    e: float            # eccentricity
    inc: float          # inclination [deg]
    node: float         # longitude of ascending node [deg]
    argp: float         # argument of perihelion [deg]
    M0: float           # mean anomaly at epoch [deg]
    epoch_jd: float = J2000_JD   # epoch of M0 [JD, TDB]
    label: str = "user-defined"

    @property
    def period_years(self):
        return self.a ** 1.5          # Kepler's third law (Sun's mass only)

    @property
    def mean_motion_rad_per_day(self):
        return GAUSS_K / self.a ** 1.5


def solve_kepler(M, e, tol=1e-13, max_iter=50):
    """Solve Kepler's equation M = E - e sin(E) for E (radians) with Newton's method."""
    M = np.asarray(M, dtype=float)
    Mw = (M + np.pi) % (2 * np.pi) - np.pi          # wrap to [-pi, pi)
    E = Mw + e * np.sin(Mw)                         # good first guess for e < 1
    for _ in range(max_iter):
        dE = (E - e * np.sin(E) - Mw) / (1.0 - e * np.cos(E))
        E = E - dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def heliocentric_ecliptic(el: OrbitalElements, jd):
    """Position of the body relative to the Sun, ecliptic J2000 axes. Returns (..., 3) in AU."""
    jd = np.asarray(jd, dtype=float)
    M = np.radians(el.M0) + el.mean_motion_rad_per_day * (jd - el.epoch_jd)
    E = solve_kepler(M, el.e)
    # position inside the orbital plane (x axis points to perihelion)
    xp = el.a * (np.cos(E) - el.e)
    yp = el.a * np.sqrt(1.0 - el.e ** 2) * np.sin(E)
    # rotate by argp (inside the plane), inc (tilt), node (around the ecliptic pole)
    w, i, O = np.radians([el.argp, el.inc, el.node])
    cw, sw, ci, si, cO, sO = np.cos(w), np.sin(w), np.cos(i), np.sin(i), np.cos(O), np.sin(O)
    x = (cO * cw - sO * sw * ci) * xp + (-cO * sw - sO * cw * ci) * yp
    y = (sO * cw + cO * sw * ci) * xp + (-sO * sw + cO * cw * ci) * yp
    z = (sw * si) * xp + (cw * si) * yp
    return np.stack([x, y, z], axis=-1)


def ecliptic_to_equatorial(xyz):
    """Rotate ecliptic J2000 axes to equatorial (ICRF-like) axes."""
    eps = np.radians(OBLIQUITY_DEG)
    x, y, z = np.moveaxis(np.asarray(xyz, dtype=float), -1, 0)
    return np.stack([x, y * np.cos(eps) - z * np.sin(eps), y * np.sin(eps) + z * np.cos(eps)], axis=-1)


def equatorial_to_ecliptic(xyz):
    """Inverse of ecliptic_to_equatorial."""
    eps = np.radians(OBLIQUITY_DEG)
    x, y, z = np.moveaxis(np.asarray(xyz, dtype=float), -1, 0)
    return np.stack([x, y * np.cos(eps) + z * np.sin(eps), -y * np.sin(eps) + z * np.cos(eps)], axis=-1)


def mean_anomaly_for_distance(a, e, r, outbound=True):
    """Mean anomaly [deg] at which the heliocentric distance equals r (needs q <= r <= Q)."""
    if e < 1e-9:
        return 0.0
    cosE = np.clip((1.0 - r / a) / e, -1.0, 1.0)
    E = np.arccos(cosE)
    if not outbound:
        E = -E
    return float(np.degrees(E - e * np.sin(E)))


def brown_batygin_2021(r_au=None, node=DEMO_NODE_DEG, argp=DEMO_ARGP_DEG, epoch_jd=J2000_JD):
    """HYPOTHETICAL Planet Nine, Brown & Batygin 2021 best fit: a=380 AU, q=300 AU, i=16 deg.

    e is derived from a and q. node/argp are NOT in the abstract: demonstration defaults.
    r_au sets where on the orbit the planet is at epoch (default: at distance a = 380 AU);
    it is clipped to the range the orbit allows (300 to 460 AU).
    """
    e = 1.0 - P9_Q_AU / P9_A_AU
    r = P9_A_AU if r_au is None else float(np.clip(r_au, P9_Q_AU, P9_A_AU * (1 + e)))
    M0 = mean_anomaly_for_distance(P9_A_AU, e, r)
    return OrbitalElements(P9_A_AU, e, P9_INC_DEG, node, argp, M0, epoch_jd,
                           "Brown & Batygin 2021 best fit (hypothetical)")


def circular_orbit(r_au, inc=P9_INC_DEG, node=DEMO_NODE_DEG, lon_in_orbit=150.0, epoch_jd=J2000_JD):
    """User-defined circular orbit at distance r_au (lets us scan any distance)."""
    return OrbitalElements(float(r_au), 0.0, inc, node, 0.0, lon_in_orbit, epoch_jd,
                           f"circular orbit at {r_au:g} AU (user-defined)")


def with_epoch(el: OrbitalElements, epoch_jd):
    """Same orbit, with M0 re-expressed at another epoch."""
    M = el.M0 + np.degrees(el.mean_motion_rad_per_day) * (epoch_jd - el.epoch_jd)
    return replace(el, M0=float(M % 360.0), epoch_jd=float(epoch_jd))
