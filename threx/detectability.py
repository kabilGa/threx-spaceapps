# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Could SPHEREx see it?  Compare the predicted magnitude with the SPHEREx depth.

Published number (MEASURED by the mission, Bock et al. 2025 abstract): 5-sigma depth per
spectral channel of 19.5-19.9 AB (0.75-3.8 micron) and 17.8-18.8 AB (3.8-5.0 micron) AFTER the
four all-sky passes of the two-year mission.

Our ASSUMPTIONS on top of that:
 * one pass (one "visit") is shallower by 2.5 log10(sqrt(4)) = 0.75 mag (noise averages
   down as the square root of the number of passes);
 * adding K channels with equal signal-to-noise gains 2.5 log10(sqrt(K)) mag. This is an
   upper bound: a real planet does not have the same brightness in every channel.
"""
import numpy as np

from .constants import SPHEREX_BANDS, SPHEREX_N_CHANNELS, SPHEREX_N_SURVEYS
from .photometry import reflected_mag

CHANNELS_PER_BAND = SPHEREX_N_CHANNELS // len(SPHEREX_BANDS)       # 17
SINGLE_VISIT_PENALTY = 2.5 * np.log10(np.sqrt(SPHEREX_N_SURVEYS))   # 0.75 mag (assumption)


def stacking_gain(n):
    """Magnitudes gained by adding n equally good measurements: 2.5 log10(sqrt(n))."""
    return 2.5 * np.log10(np.sqrt(n))


def band_info(band):
    return next(b for b in SPHEREX_BANDS if b["band"] == band)


def magnitude_grid(distances_au, radii_rearth, albedo, band):
    """Predicted AB magnitude on a (radius, distance) grid for one SPHEREx detector band.
    Geometry: zero phase, Earth-planet distance = Sun-planet distance (error < 0.01 mag
    beyond 200 AU)."""
    b = band_info(band)
    d, R = np.meshgrid(np.asarray(distances_au, dtype=float), np.asarray(radii_rearth, dtype=float))
    return reflected_mag(d, d, R, albedo, lam_um=b["lam_c"])


def thresholds(band, pessimistic=True):
    """Limiting magnitudes for the four ways of looking, for one band.

    pessimistic=True uses the bright end of the published depth range (19.5 / 17.8)."""
    depth = band_info(band)["depth"][0 if pessimistic else 1]
    return {
        "one_channel_one_visit": depth - SINGLE_VISIT_PENALTY,
        "band_stack_one_visit": depth - SINGLE_VISIT_PENALTY + stacking_gain(CHANNELS_PER_BAND),
        "one_channel_four_visits": depth,
        "band_stack_four_visits": depth + stacking_gain(CHANNELS_PER_BAND),
    }


CLASS_LABELS = [
    "Seen in ONE channel in ONE visit",
    "Needs the ~17 channels of the band added (one visit)",
    "Needs channels added AND the 4 visits shift-and-stacked",
    "Out of reach",
]


def classify(mag, band, pessimistic=True):
    """0..3 index into CLASS_LABELS for each magnitude."""
    th = thresholds(band, pessimistic)
    mag = np.asarray(mag, dtype=float)
    out = np.full(mag.shape, 3, dtype=int)
    out[mag <= th["band_stack_four_visits"]] = 2
    out[mag <= th["band_stack_one_visit"]] = 1
    out[mag <= th["one_channel_one_visit"]] = 0
    return out


def summary_table(distance_au, radius_rearth, albedo):
    """One row per band: predicted magnitude, thresholds and verdicts (computed)."""
    rows = []
    for b in SPHEREX_BANDS:
        m = float(reflected_mag(distance_au, distance_au, radius_rearth, albedo, lam_um=b["lam_c"]))
        rows.append({
            "band": b["band"], "range_um": f'{b["lam_min"]:.2f}-{b["lam_max"]:.2f}',
            "mag_ab": m, "depth_published": f'{b["depth"][0]}-{b["depth"][1]}',
            "verdict_pessimistic": CLASS_LABELS[int(classify(m, b["band"], True))],
            "verdict_optimistic": CLASS_LABELS[int(classify(m, b["band"], False))],
        })
    return rows


def verdict_text(distance_au, radius_rearth, albedo):
    """Two plain sentences for the app (computed from the model, not a measurement)."""
    rows = summary_table(distance_au, radius_rearth, albedo)
    best = min(rows, key=lambda r: r["mag_ab"] - thresholds(r["band"], True)["one_channel_one_visit"])
    cls_pess = int(classify(best["mag_ab"], best["band"], True))
    cls_opt = int(classify(best["mag_ab"], best["band"], False))
    s1 = (f"A hypothetical planet of {radius_rearth:g} Earth radii with albedo {albedo:g} at "
          f"{distance_au:g} AU would have AB magnitude {best['mag_ab']:.1f} in SPHEREx band "
          f"{best['band']} ({best['range_um']} micron), its most favourable band.")
    if cls_pess == cls_opt:
        s2 = f"Verdict: {CLASS_LABELS[cls_pess].lower()}."
    else:
        s2 = (f"Verdict: between '{CLASS_LABELS[cls_opt].lower()}' (optimistic depth) and "
              f"'{CLASS_LABELS[cls_pess].lower()}' (pessimistic depth).")
    return s1 + " " + s2
