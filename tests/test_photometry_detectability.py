# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Brightness formula and detectability logic."""
import numpy as np
import pytest

from threx import photometry, detectability
from threx.constants import SPHEREX_BANDS


def test_brief_sanity_value():
    """Brief: Delta = r = 400 AU, R = 3 R_earth, p = 0.4 gives m ~ 19.7 (V-like, AB)."""
    assert float(photometry.reflected_mag(400, 400, 3.0, 0.4)) == pytest.approx(19.7, abs=0.1)


def test_scaling_laws():
    m0 = float(photometry.reflected_mag(400, 400, 3.0, 0.4))
    assert float(photometry.reflected_mag(800, 800, 3.0, 0.4)) - m0 == pytest.approx(10 * np.log10(2))   # r^-4
    assert float(photometry.reflected_mag(400, 400, 6.0, 0.4)) - m0 == pytest.approx(-5 * np.log10(2))   # R^2
    assert float(photometry.reflected_mag(400, 400, 3.0, 0.2)) - m0 == pytest.approx(2.5 * np.log10(2))  # p


def test_sun_magnitudes():
    assert float(photometry.sun_ab_mag(0.5511)) == pytest.approx(-26.77)
    assert float(photometry.sun_ab_mag(1.2393)) == pytest.approx(-27.03)
    # a 5772 K blackbody must agree with the tabulated Sun to within ~0.25 mag over 0.55-5 micron
    lam = np.linspace(0.55, 5.0, 40)
    assert np.max(np.abs(photometry.sun_ab_mag(lam) - photometry.sun_ab_mag_blackbody(lam))) < 0.25
    # the Sun gets fainter (in AB) towards 5 micron
    assert float(photometry.sun_ab_mag(5.0)) > float(photometry.sun_ab_mag(3.0)) > float(photometry.sun_ab_mag(1.2))


def test_thermal_emission_is_negligible():
    """Blackbody glow at 10-50 K is more than 30 magnitudes fainter than reflected light at 3-5 micron."""
    for T in (10, 30, 50):
        for lam in (3.0, 4.0, 5.0):
            thermal = float(photometry.thermal_mag(T, 3.0, 400.0, lam))
            reflected = float(photometry.reflected_mag(400, 400, 3.0, 0.4, lam_um=lam))
            assert thermal - reflected > 30
    assert photometry.equilibrium_temperature(400.0) == pytest.approx(12.7, abs=0.5)


def test_thresholds_are_ordered():
    for b in SPHEREX_BANDS:
        th = detectability.thresholds(b["band"])
        assert th["one_channel_one_visit"] < th["one_channel_four_visits"] < th["band_stack_four_visits"]
        assert th["one_channel_one_visit"] < th["band_stack_one_visit"]
        assert th["one_channel_four_visits"] == b["depth"][0]
    assert detectability.SINGLE_VISIT_PENALTY == pytest.approx(0.753, abs=0.001)
    assert detectability.stacking_gain(100) == pytest.approx(2.5)


def test_classification_gets_harder_with_distance():
    d = np.linspace(200, 1000, 50)
    grid = detectability.magnitude_grid(d, [3.0], 0.4, band=2)
    cls = detectability.classify(grid[0], 2)
    assert np.all(np.diff(grid[0]) > 0)          # fainter when farther
    assert np.all(np.diff(cls) >= 0)             # class never improves with distance
    assert cls[0] == 0 and cls[-1] == 3          # bright at 200 AU, out of reach at 1000 AU


def test_verdict_text_mentions_inputs():
    txt = detectability.verdict_text(400.0, 3.0, 0.4)
    assert "400" in txt and "Verdict" in txt
