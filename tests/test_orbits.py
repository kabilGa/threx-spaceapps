"""Kepler solver and orbit geometry."""
import numpy as np
import pytest

from threx import orbits


@pytest.mark.parametrize("e", [0.0, 0.1, 0.2105, 0.6, 0.86, 0.97])
def test_kepler_equation_is_satisfied(e):
    M = np.linspace(-10, 10, 2001)
    E = orbits.solve_kepler(M, e)
    Mw = (M + np.pi) % (2 * np.pi) - np.pi
    assert np.max(np.abs(E - e * np.sin(E) - Mw)) < 1e-11


def test_kepler_known_values():
    assert orbits.solve_kepler(1.234, 0.0) == pytest.approx(1.234)          # circle: E = M
    assert orbits.solve_kepler(0.0, 0.5) == pytest.approx(0.0, abs=1e-12)   # perihelion
    # aphelion: E = pi (the solver may return -pi, which is the same point on the orbit)
    assert np.cos(orbits.solve_kepler(np.pi, 0.5)) == pytest.approx(-1.0)
    # textbook example: M = 27 deg, e = 0.5  ->  E = 48.43418 deg
    assert np.degrees(orbits.solve_kepler(np.radians(27.0), 0.5)) == pytest.approx(48.43418, abs=1e-4)


def test_distance_at_perihelion_and_aphelion():
    el = orbits.OrbitalElements(a=380.0, e=1 - 300 / 380, inc=16, node=100, argp=150, M0=0.0)
    r_peri = np.linalg.norm(orbits.heliocentric_ecliptic(el, el.epoch_jd))
    half_period_days = 0.5 * el.period_years * 365.25
    r_apo = np.linalg.norm(orbits.heliocentric_ecliptic(el, el.epoch_jd + half_period_days))
    assert r_peri == pytest.approx(300.0, rel=1e-9)
    assert r_apo == pytest.approx(460.0, rel=1e-6)


def test_period_and_inclination():
    el = orbits.circular_orbit(400.0, inc=16.0)
    assert el.period_years == pytest.approx(8000.0)                # 400^1.5
    jd = el.epoch_jd + np.linspace(0, el.period_years * 365.25, 500)
    xyz = orbits.heliocentric_ecliptic(el, jd)
    lat = np.degrees(np.arcsin(xyz[:, 2] / np.linalg.norm(xyz, axis=1)))
    assert lat.max() == pytest.approx(16.0, abs=0.01)              # highest latitude = inclination
    assert np.allclose(np.linalg.norm(xyz, axis=1), 400.0)


def test_brown_batygin_preset():
    el = orbits.brown_batygin_2021(r_au=400.0)
    assert el.e == pytest.approx(1 - 300 / 380)
    assert np.linalg.norm(orbits.heliocentric_ecliptic(el, el.epoch_jd)) == pytest.approx(400.0, rel=1e-9)
    # distances outside what the orbit allows are clipped to perihelion/aphelion
    el_far = orbits.brown_batygin_2021(r_au=900.0)
    assert np.linalg.norm(orbits.heliocentric_ecliptic(el_far, el_far.epoch_jd)) == pytest.approx(460.0, rel=1e-6)


def test_frame_rotation_roundtrip():
    v = np.random.default_rng(1).normal(size=(50, 3))
    back = orbits.equatorial_to_ecliptic(orbits.ecliptic_to_equatorial(v))
    assert np.allclose(back, v)
    # the ecliptic north pole is at RA 270 deg, Dec 90 - obliquity
    pole = orbits.ecliptic_to_equatorial(np.array([0.0, 0.0, 1.0]))
    assert np.degrees(np.arctan2(pole[1], pole[0])) % 360 == pytest.approx(270.0)
    assert np.degrees(np.arcsin(pole[2])) == pytest.approx(90 - 23.4392911, abs=1e-6)
