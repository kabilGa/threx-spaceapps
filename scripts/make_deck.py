"""Make slides/threx_deck.pdf (7 slides, English) from the project figures. Run: python scripts/make_deck.py"""
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from threx import blink, plots

R = Path(__file__).resolve().parents[1]
W, H = 13.33, 7.5
INK, MUTED, BLUE = "#0b0b0b", "#52514e", "#2a78d6"


def slide(pdf, title, bullets, images=(), wrap=None, foot="THREX · NASA Space Apps 2026 · Planet X is hypothetical: no detection is claimed"):
    fig = plt.figure(figsize=(W, H), facecolor="#fcfcfb")
    fig.text(0.05, 0.9, title, fontsize=28, fontweight="bold", color=INK)
    wide = bool(images)
    y = 0.8
    for b in bullets:
        lines = textwrap.wrap(b, wrap or (34 if wide else 95))
        fig.text(0.05, y, "• " + lines[0], fontsize=16, color=INK, va="top")
        for ln in lines[1:]:
            y -= 0.045
            fig.text(0.062, y, ln, fontsize=16, color=INK, va="top")
        y -= 0.085
    for (path, box) in images:
        ax = fig.add_axes(box)
        ax.imshow(mpimg.imread(path) if isinstance(path, (str, Path)) else path)
        ax.axis("off")
    fig.text(0.05, 0.03, foot, fontsize=10, color=MUTED)
    pdf.savefig(fig)
    plt.close(fig)


tmp = R / "data_cache" / "_work"
b = blink.load_bundle("Barnard's Star")
plots.fig_frames(b, [0, 1, 20, 21, 40, 41, 56, 57], ncols=4).savefig(tmp / "deck_frames.png", dpi=110)
a = blink.load_bundle("Asteroid proof: 6 Hebe (A847 NA)")
pr = a["meta"]["proof"]
pred = {q: (f["pred_east_arcsec"], f["pred_north_arcsec"]) for q, f in enumerate(pr["frames"])}
meas = {q: (f["meas_east_arcsec"], f["meas_north_arcsec"]) for q, f in enumerate(pr["frames"]) if f.get("meas_east_arcsec") is not None}
plots.fig_pair(a, 0, 9, pred, meas).savefig(tmp / "deck_asteroid.png", dpi=110)
d = pr["displacement"]
S = R / "screenshots"
with PdfPages(R / "slides" / "threx_deck.pdf") as pdf:
    slide(pdf, "THREX: see how the sky changes in SPHEREx images",
          ["A web viewer for real SPHEREx images, plus a prediction of what a hypothetical Planet X would look like to SPHEREx.",
           "NASA Space Apps Challenge 2026: Planet X and SPHEREx", "Team: Kabil Gannouni (Tunisia)",
           "Everything is labelled measured, computed or assumed."])
    slide(pdf, "The problem", ["SPHEREx maps the whole sky in 102 near-infrared colours, about every six months.",
          "Asteroids, comets, nearby stars, brown dwarfs and maybe an unseen planet reveal themselves by moving between images.",
          "No person can blink through that much data: we need easy tools."])
    slide(pdf, "1 · The viewer (real SPHEREx pixels)",
          ["Type RA/Dec or pick a challenge target; the app lists the real exposures and plays aligned cutouts as a time-lapse.",
           "Each frame shows date, wavelength, exposure ID; flagged pixels can be shown.",
           "Frames differ in wavelength: colour is not variability."], [(tmp / "deck_frames.png", [0.42, 0.12, 0.56, 0.72])])
    slide(pdf, "2 · Proof on a known asteroid",
          [f"6 Hebe, predicted by JPL Horizons. Two groups of frames {d['dt_hours']:.1f} h apart.",
           f"Measured displacement {d['measured_arcsec']:.0f} arcsec; predicted {d['predicted_arcsec']:.0f} arcsec.",
           "JPL sb_ident lists 6 Hebe in the field. Hebe saturates: only its positions are meaningful. Already known: nothing discovered."], [(tmp / "deck_asteroid.png", [0.03, 0.04, 0.94, 0.42])], wrap=100)
    slide(pdf, "3 · What would Planet X look like?",
          ["A planet at 400 AU shifts by about 1000 arcsec (160 pixels) in six months, mostly Earth's parallax.",
           "Shift = distance measurement: d ≈ B⊥ / θ.", "Orbit code agrees with JPL Horizons to 0.4 arcsec (Sedna, Eris)."],
          [(S / "figure_signature_sky_track.png", [0.45, 0.15, 0.53, 0.65])])
    slide(pdf, "4 · Could SPHEREx see it?",
          ["Reflected sunlight at 400 AU, 3 Earth radii, albedo 0.4: AB ≈ 19.4-21.5 depending on band.",
           "Only within reach by adding channels in bands 1-3; out of reach beyond 3.8 µm or with a methane-dark albedo.",
           "Assumptions (albedo, radius, one-visit depth) are labelled."], [(S / "figure_detectability_heatmap.png", [0.42, 0.1, 0.57, 0.75])])
    slide(pdf, "Honesty, sources and AI use",
          ["Planet X is hypothetical. No detection, no discovery.",
           "Data: SPHEREx at IRSA, JPL Horizons / sb_ident / SBDB, SIMBAD. Papers: Bock 2025, Brown & Batygin 2021, Fortney 2016.",
           "Limits: saturated bright stars (Hebe, R Dra), only a subset of the exposures is cached (e.g. 68 of 358 for Barnard), simple photometry, no WISE 0855 detection claimed.",
           "AI: the concept was developed with Claude (Anthropic) in a chat; Claude Code wrote the code, tests, app and first drafts. The team member (Kabil Gannouni) directed the work and is responsible for all content.",
           "Code, data and tests: https://github.com/kabilGa/threx-spaceapps"])
print("slides/threx_deck.pdf written")
