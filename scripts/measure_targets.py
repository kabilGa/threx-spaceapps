"""Measure what is (and is not) visible at the catalogue position of each target, and make the ready-made GIFs.

  python scripts/measure_targets.py
Writes data_cache/target_checks.json and data_cache/gif_<target>.gif.  Aperture photometry (2-pixel radius,
sky ring 4-7 pixels, no aperture correction) at the SIMBAD position moved to each exposure date.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from threx import blink, spherex_io, targets

ROOT = Path(__file__).resolve().parents[1]
out = {}
for name in targets.TARGETS:
    try:
        b = blink.load_bundle(name)
    except spherex_io.DataUnavailable:
        print(name, "no bundle yet")
        continue
    cat = blink.catalogue_offsets(b)
    n, g = b["frames"].shape[-1], b["meta"]["grid_pix_arcsec"]
    rows = []
    for k, (e, nn) in enumerate(cat):
        col, row = n // 2 - e / g, n // 2 + nn / g
        inside = 8 < col < n - 8 and 8 < row < n - 8
        ph = blink.photometry(b, k, col, row) if inside else None
        f = b["meta"]["frames"][k]
        rows.append({"detector": f["detector"], "lam_um": f["lam_um"], "inside": bool(inside),
                     "snr": None if not ph or ph["snr"] is None else float(ph["snr"]),
                     "mag_ab": None if not ph or ph["mag_ab"] is None else float(ph["mag_ab"]),
                     "frac_flagged": None if not ph else ph["n_flagged"]})
    res = {"n_frames": len(rows), "per_detector": {}}
    for d in sorted({r["detector"] for r in rows}):
        sn = [r["snr"] for r in rows if r["detector"] == d and r["snr"] is not None]
        mg = [r["mag_ab"] for r in rows if r["detector"] == d and r["snr"] is not None and r["snr"] > 5 and r["mag_ab"] is not None]
        res["per_detector"][str(d)] = {"n": sum(r["detector"] == d for r in rows),
                                       "median_snr": float(np.median(sn)) if sn else None,
                                       "n_snr_gt5": int(sum(s > 5 for s in sn)),
                                       "median_ab_where_snr_gt5": float(np.median(mg)) if mg else None}
    out[name] = res
    gif = blink.make_gif(b, markers=cat, duration_ms=450, scale=3)
    (ROOT / "data_cache" / f"gif_{spherex_io._slug(name)}.gif").write_bytes(gif)
    print(name, json.dumps(res["per_detector"]), f"gif {len(gif) / 1e6:.2f} MB")
(ROOT / "data_cache" / "target_checks.json").write_text(json.dumps(out, indent=1))
