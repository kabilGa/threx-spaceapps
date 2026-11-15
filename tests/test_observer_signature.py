# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Geocentric transform, parallax formula, motion rate, distance estimate."""
import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord, BarycentricMeanEcliptic
from astropy.time import Time

from threx import orbits, signature
from threx.constants import ARCSEC_PER_RAD, GAUSS_K
from threx.observer import (earth_heliocentric, geocentric_vector, radec, ecliptic_lonlat,
                            apparent_position, separation_arcsec, unit_vector)

T0 = Time("2026-01-01T00:00:00", scale="utc")


def test_earth_orbit_is_sane():
    t = T0 + np.linspace(0, 365.25, 200) * u.day
    r = np.linalg.norm(earth_heliocentric(t), axis=-1)
    assert 0.983 < r.min() < 0.984          # perihelion distance in AU
    assert 1.016 < r.max() < 1.017          # aphelion distance in AU
    # in early January the Earth is near perihelion and the Sun is seen near ecliptic longitude 280 deg
    lon_sun, lat_sun = ecliptic_lonlat(-earth_heliocentric(Time("2026-01-03T00:00:00")))
    assert lon_sun == pytest.approx(282.5, abs=1.0)
    assert abs(lat_sun) < 0.01


def test_radec_matches_astropy_frames():
    # a far-away direction given in ecliptic coordinates must convert like astropy does
    lon, lat = 123.4, 16.0
    c = SkyCoord(lon * u.deg, lat * u.deg, frame=BarycentricMeanEcliptic(equinox="J2000")).icrs
    ecl = unit_vector(lon, lat)                      # same formula works for (lon, lat)
    ra, dec = radec(orbits.ecliptic_to_equatorial(ecl))
    assert separation_arcsec(unit_vector(ra, dec), unit_vector(c.ra.deg, c.dec.deg)) < 0.1


def test_geocentric_minus_heliocentric_is_earth():
    el = orbits.circular_orbit(400.0, epoch_jd=T0.tdb.jd)
    geo = geocentric_vector(el, T0, light_time=False)
    helio = orbits.ecliptic_to_equatorial(orbits.heliocentric_ecliptic(el, T0.tdb.jd))
    assert np.allclose(helio - geo, earth_heliocentric(T0), atol=1e-12)
    assert 399.0 < np.linalg.norm(geo) < 401.1


def _opposition(el):
    yr = signature.rate_through_year(el, T0, n=1462)
    return T0 + yr["days"][int(np.argmax(yr["elongation_deg"]))] * u.day


def test_parallax_shift_matches_small_angle_formula():
    """Shift over 6 months centred on opposition ~ 2 AU / distance (brief: ~1000 arcsec at 400 AU)."""
    el = orbits.circular_orbit(400.0, inc=0.0, epoch_jd=T0.tdb.jd)
    t_opp = _opposition(el)
    s = signature.shift_between(el, t_opp - 91.3 * u.day, t_opp + 91.3 * u.day)
    expected = 2.0 / 400.0 * ARCSEC_PER_RAD                       # 1031 arcsec
    assert s["parallax_arcsec"] == pytest.approx(expected, rel=0.03)
    assert s["parallax_arcsec"] == pytest.approx(1000.0, rel=0.10)
    assert s["parallax_pixels"] == pytest.approx(170.0, rel=0.10)
    # the body's own motion: mean motion x time (brief: ~80 arcsec)
    orbital_expected = GAUSS_K / 400.0 ** 1.5 * 182.6 * ARCSEC_PER_RAD
    assert s["orbital_arcsec"] == pytest.approx(orbital_expected, rel=0.01)
    assert s["orbital_arcsec"] == pytest.approx(80.0, rel=0.10)
    # around opposition the two parts point in opposite directions
    assert s["total_arcsec"] == pytest.approx(s["parallax_arcsec"] - s["orbital_arcsec"], rel=0.01)
    assert s["parallax_arcsec"] > 10 * s["orbital_arcsec"]         # parallax dominates


def test_parallax_scales_as_one_over_distance():
    shifts = []
    for d in (200.0, 400.0, 800.0):
        el = orbits.circular_orbit(d, inc=0.0, epoch_jd=T0.tdb.jd)
        t_opp = _opposition(el)
        shifts.append(signature.shift_between(el, t_opp - 91.3 * u.day, t_opp + 91.3 * u.day)["parallax_arcsec"])
    assert shifts[0] / shifts[1] == pytest.approx(2.0, rel=0.02)
    assert shifts[1] / shifts[2] == pytest.approx(2.0, rel=0.02)


def test_opposition_rate_is_retrograde_and_right_size():
    """Brief: ~9 arcsec/day, ~0.36 arcsec/hour, ~1.4 pixels/day, retrograde at 400 AU."""
    el = orbits.circular_orbit(400.0, inc=0.0, epoch_jd=T0.tdb.jd)
    r = signature.motion_rate(el, _opposition(el))
    assert bool(r["retrograde"])
    assert float(r["arcsec_per_day"]) == pytest.approx(9.0, rel=0.10)
    assert float(r["arcsec_per_hour"]) == pytest.approx(0.36, rel=0.10)
    assert float(r["pixels_per_day"]) == pytest.approx(1.4, rel=0.10)


def test_rate_is_small_at_quadrature():
    """90 deg from the Sun (where SPHEREx looks) Earth moves along the line of sight,
    so the apparent motion almost stops."""
    el = orbits.circular_orbit(400.0, inc=0.0, epoch_jd=T0.tdb.jd)
    yr = signature.rate_through_year(el, T0, n=1462)
    i_quad = int(np.argmin(np.abs(yr["elongation_deg"] - 90.0)))
    assert yr["arcsec_per_day"][i_quad] < 0.1 * yr["arcsec_per_day"].max()


def test_distance_from_shift_recovers_truth():
    el = orbits.circular_orbit(400.0, inc=0.0, epoch_jd=T0.tdb.jd)
    t_opp = _opposition(el)
    t1, t2 = t_opp - 91.3 * u.day, t_opp + 91.3 * u.day
    s = signature.shift_between(el, t1, t2)
    pos = apparent_position(el, t1)
    simple = signature.distance_from_shift(s["total_arcsec"], t1, t2, float(pos["ra"]), float(pos["dec"]))
    assert simple["distance_au"] == pytest.approx(400.0, rel=0.10)
    assert abs(simple["distance_au"] / 400.0 - 1) <= simple["sys_frac"] * 1.3     # error budget holds
    assert simple["b_perp_au"] == pytest.approx(2.0, rel=0.03)
    refined = signature.distance_from_shift_refined(s["total_arcsec"], t1, t2, float(pos["ra"]), float(pos["dec"]))
    assert refined == pytest.approx(400.0, rel=0.02)
