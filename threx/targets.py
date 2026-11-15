# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Preset sky positions for the real-data page. Each one says where its coordinates come from."""

import json
from pathlib import Path

_SIMBAD = Path(__file__).resolve().parents[1] / "data_cache" / "targets_simbad.json"


def _load():
    """The five targets named by the challenge, from SIMBAD (see data_cache/targets_simbad.json).
    Coordinates are moved to 2026.0 with the catalogue proper motion (the SPHEREx epoch)."""
    out = {}
    try:
        data = json.loads(_SIMBAD.read_text())
    except OSError:
        return out
    for name, t in data["targets"].items():
        dy = 26.0
        out[name] = {
            "ra": t["ra"] + t["pmra"] * dy / 3.6e6 / __import__("math").cos(__import__("math").radians(t["dec"])),
            "dec": t["dec"] + t["pmdec"] * dy / 3.6e6,
            "note": f"{name} ({t['simbad_id']}): SIMBAD position, J2000, moved to 2026.0 with its proper motion "
                    f"({t['pmra']:.0f}, {t['pmdec']:.0f} mas/yr).",
        }
    return out


TARGETS = _load()


def cached_proof_targets():
    """Presets created by scripts/build_asteroid_proof.py: any cached bundle with a 'proof'."""
    from .spherex_io import CACHE_DIR
    out = {}
    for js in sorted(CACHE_DIR.glob("bundle_*.json")):
        meta = json.loads(js.read_text())
        if meta.get("proof"):
            out[meta["name"]] = {"ra": meta["ra_deg"], "dec": meta["dec_deg"],
                                 "note": meta["proof"].get("target_note", "Known asteroid, used as a proof.")}
    return out


def all_targets():
    """The five challenge targets first, then the asteroid proof."""
    return {**TARGETS, **cached_proof_targets()}


# What we honestly see in the cached frames (written after looking at the pixels and at
# data_cache/target_checks.json; positions: SIMBAD, proper motion from SIMBAD).
COMMENTS = {
    "Barnard's Star": "A very bright star (AB about 7-9 mag in every detector, so its core is saturated and flagged). The cached "
                      "frames span 13 months (Aug 2025 to Sep 2026, three visits). Its proper motion (10.4 arcsec/year) "
                      "moves it about 11 arcsec = 1.8 pixels over that time, and 0.4 arcsec inside one visit: **within a "
                      "visit it does not move visibly; across visits the shift is only about two pixels** (the orange "
                      "circle shows where the catalogue says it is). We did not measure the shift ourselves.",
    "Epsilon Indi B": "The brown-dwarf pair (not resolved at 6.2 arcsec per pixel) is clearly detected. The frames span "
                      "17 months; its proper motion (4.7 arcsec/year) moves it about 6.8 arcsec = 1.1 pixel in that time: "
                      "**too small to see by eye**. The orange circle follows the source.",
    "WISE 0855-0714": "Extremely cold brown dwarf moving 8.2 arcsec/year (1.3 pixel over the 12 months of these frames). "
                      "WISE 0855-0714 is extremely faint below about 1.6 micron (Hubble photometry in Vega magnitudes: 26.31 +/- 0.10 at 1.153 micron and 23.22 +/- 0.35 at 1.575 micron; Luhman and Esplin 2016, arXiv:1605.08620), far below SPHEREx's ~19.5 AB single-band depth, so a source measured near its position at those wavelengths is most likely a background blend and we claim no detection.",
    "RW Aurigae": "A bright variable star. Frames at different wavelengths cannot show its variability by themselves: "
                  "to see real changes you would compare frames of the same detector and wavelength at different dates. "
                  "We did not do this analysis.",
    "R Draconis": "A bright long-period variable (AB about 8 at 1-2 µm: the core is saturated or flagged in most "
                  "short-wavelength frames, so the simple photometry often fails there). Same remark as RW Aurigae: the "
                  "wavelength changes from frame to frame, so brightness differences are not variability. Not analysed.",
}
