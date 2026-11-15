# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Save the Tier 1 figures as PNG files in screenshots/ (same code as the app).

Run:  python scripts/make_figures.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from astropy.time import Time

from threx import orbits, plots, signature, validation

OUT = Path(__file__).resolve().parents[1] / "screenshots"
OUT.mkdir(exist_ok=True)
T0 = Time("2026-01-01T00:00:00", scale="utc")
el = orbits.brown_batygin_2021(r_au=400.0, epoch_jd=T0.tdb.jd)

plots.fig_sky_track(el, T0).savefig(OUT / "figure_signature_sky_track.png")
plots.fig_rate(el, T0).savefig(OUT / "figure_motion_rate.png")
q = signature.quadrature_times(el, T0)
plots.fig_shift_split(signature.shift_between(el, q[0], q[1])).savefig(OUT / "figure_shift_split.png")
plots.fig_heatmaps(0.4, mark=(400, 3.0)).savefig(OUT / "figure_detectability_heatmap.png")
plots.fig_heatmaps(0.03, mark=(400, 3.0)).savefig(OUT / "figure_detectability_methane_dark.png")
plots.fig_validation(validation.compare()).savefig(OUT / "figure_horizons_validation.png")
print("saved figures in", OUT)

# ---- real-data views (same functions as the app's viewer) ----
from threx import blink, spherex_io, targets  # noqa: E402

b = blink.load_bundle("Barnard's Star")
cat = blink.catalogue_offsets(b)
n = len(b["frames"])
plots.fig_frames(b, [0, 1, n // 3, n // 3 + 1, 2 * n // 3, 2 * n // 3 + 1, n - 2, n - 1], markers=cat, ncols=4).savefig(
    OUT / "figure_viewer_frames.png")
w = blink.load_bundle("WISE 0855-0714")
idx = [q for q, f in enumerate(w["meta"]["frames"]) if f["detector"] == 5][:4] + \
      [q for q, f in enumerate(w["meta"]["frames"]) if f["detector"] == 1][:4]
plots.fig_frames(w, idx, markers=blink.catalogue_offsets(w), ncols=4).savefig(OUT / "figure_viewer_wise0855.png")
df, _ = spherex_io.get_exposures("Barnard's Star", 0, 0, live=False)
vis = spherex_io.find_visits(df)
plots.fig_timeline(df, vis, (df["t_mid"].min(), df["t_max"].max() if "t_max" in df else df["t_mid"].max()),
                   "Barnard's Star: real SPHEREx exposures (cached subset)").savefig(OUT / "figure_schedule_timeline.png")
a = blink.load_bundle("Asteroid proof: 6 Hebe (A847 NA)")
pr = a["meta"]["proof"]
pred = {q: (f["pred_east_arcsec"], f["pred_north_arcsec"]) for q, f in enumerate(pr["frames"])}
meas = {q: (f["meas_east_arcsec"], f["meas_north_arcsec"]) for q, f in enumerate(pr["frames"]) if f.get("meas_east_arcsec") is not None}
plots.fig_pair(a, 0, 9, pred, meas).savefig(OUT / "figure_asteroid_proof_hebe.png")
print("saved viewer figures")
