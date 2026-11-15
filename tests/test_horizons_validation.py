"""Our Kepler + parallax code against JPL Horizons for Sedna and Eris.

Offline test: uses tests/data/horizons_validation.json (REAL Horizons output, saved by
scripts/fetch_horizons_validation.py). Set THREX_LIVE=1 to also query Horizons live.

Stated tolerance: 1 arcsec in position (0.16 SPHEREx pixel) over +-180 days around the
element epoch. The remaining error (about 0.3 arcsec) comes from using a fixed two-body
orbit around the Sun: we ignore the pull of the giant planets.
"""
import os

import pytest

from threx import validation, jpl

POS_TOL_ARCSEC = 1.0
RATE_TOL_ARCSEC_PER_HOUR = 0.002
SHIFT_REL_TOL = 0.001


@pytest.fixture(scope="module")
def results():
    return validation.compare()


@pytest.mark.parametrize("body", ["Sedna (90377)", "Eris (136199)"])
def test_positions_rates_and_shifts(results, body):
    r = results["bodies"][body]
    assert len(r["jd"]) == 13
    assert r["max_pos_err_arcsec"] < POS_TOL_ARCSEC
    assert r["max_rate_err_arcsec_per_hour"] < RATE_TOL_ARCSEC_PER_HOUR
    assert r["max_shift_rel_err"] < SHIFT_REL_TOL
    assert abs(r["delta_err_au"]).max() < 1e-3          # Earth-body distance within 0.001 AU


@pytest.mark.skipif(os.environ.get("THREX_LIVE") != "1", reason="live network test (set THREX_LIVE=1)")
def test_live_horizons_still_agrees():
    from astropy.time import Time
    from threx.observer import geocentric_vector, unit_vector, separation_arcsec
    el = jpl.horizons_elements("90377;", 2461090.5)
    eph = jpl.horizons_observer("90377;", [2461090.5])[0]
    vec = geocentric_vector(el, Time(eph["jd"], format="jd", scale="utc"))
    assert float(separation_arcsec(vec, unit_vector(eph["ra"], eph["dec"]))) < POS_TOL_ARCSEC
