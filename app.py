"""THREX - Streamlit app.   Run:  streamlit run app.py

Pages: 1 SPHEREx sky viewer | 2 Planet X signature | 3 Can SPHEREx see it? | 4 About and honesty
Planet X / Planet Nine is HYPOTHETICAL. Nothing here is a detection.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from astropy.time import Time

from threx import orbits, signature, photometry, detectability, plots, validation
from threx.constants import (SPHEREX_PIXEL_ARCSEC, P9_A_AU, P9_Q_AU, P9_INC_DEG, DEMO_NODE_DEG,
                             DEMO_ARGP_DEG)
from threx.observer import apparent_position

ROOT = Path(__file__).parent
st.set_page_config(page_title="THREX - Planet X and SPHEREx", page_icon="🔭", layout="wide")

PAGES = ["1 · SPHEREx sky viewer", "2 · Extra: Planet X signature", "3 · Extra: Can SPHEREx see it?", "4 · About and honesty"]
KIND = {  # the three kinds of statements, shown as small labels everywhere
    "measured": "🟢 **MEASURED** (real data)",
    "computed": "🔵 **COMPUTED** (our model)",
    "assumed": "🟠 **ASSUMED**",
}


def show(fig, caption):
    """Draw a matplotlib figure with its caption, then free the memory."""
    st.pyplot(fig, width="stretch")
    st.caption(caption)
    plt.close(fig)


def orbit_controls(key):
    """Sidebar-like controls shared by pages 1 and 2. Returns (elements, distance, t0)."""
    c1, c2, c3 = st.columns([1.3, 1.2, 1])
    preset = c1.selectbox(
        "Orbit preset", ["Brown & Batygin 2021 best fit (hypothetical)", "User-defined circular orbit"],
        key=f"{key}_preset",
        help="Brown & Batygin 2021 (arXiv:2108.09868): a = 380 AU, perihelion 300 AU, inclination 16°. "
             "That orbit only reaches 300-460 AU from the Sun; use the circular option for other distances.")
    bb = preset.startswith("Brown")
    lo, hi = (300, 460) if bb else (200, 1000)
    dist = c2.slider("Distance from the Sun today [AU]", lo, hi, 400, step=10, key=f"{key}_dist_{bb}")
    date = c3.date_input("Reference date (UTC)", value=pd.Timestamp("2026-01-01"), key=f"{key}_date")
    t0 = Time(f"{date}T00:00:00", scale="utc")
    with st.expander("Orbit orientation (free inputs: demonstration defaults, not published values)"):
        e1, e2, e3 = st.columns(3)
        inc = e1.slider("Inclination i [deg]", 0.0, 40.0, P9_INC_DEG, key=f"{key}_inc",
                        help="16 ± 5 deg in Brown & Batygin 2021.")
        node = e2.slider("Longitude of ascending node Ω [deg]", 0.0, 360.0, DEMO_NODE_DEG, key=f"{key}_node")
        argp = e3.slider("Argument of perihelion ω [deg]", 0.0, 360.0, DEMO_ARGP_DEG, key=f"{key}_argp")
        st.caption(f"{KIND['assumed']}: Ω and ω are not given in the abstract of Brown & Batygin (2021). "
                   "They change where on the sky the planet is, not how fast it moves.")
    if bb:
        el = orbits.brown_batygin_2021(r_au=dist, node=node, argp=argp, epoch_jd=t0.tdb.jd)
        el = orbits.OrbitalElements(el.a, el.e, inc, node, argp, el.M0, el.epoch_jd, el.label)
    else:
        el = orbits.circular_orbit(dist, inc=inc, node=node, epoch_jd=t0.tdb.jd)
    return el, dist, t0


# ------------------------------------------------------------------------------ page 1
def page_signature():
    st.title("What would Planet X look like to SPHEREx?")
    st.markdown(
        "A planet hundreds of AU away hardly moves on its own orbit. What we would see is mostly the "
        "reflection of **our own motion**: as Earth circles the Sun, the planet seems to swing back and "
        "forth once a year (parallax). The size of the swing tells us its distance. "
        "*Planet X is hypothetical; this page is a prediction, not an observation.*")
    el, dist, t0 = orbit_controls("sig")

    pos = apparent_position(el, t0)
    t_opp = signature.opposition_time(el, t0, n=732)
    quads = signature.quadrature_times(el, t0)
    r_opp = signature.motion_rate(el, t_opp)
    r_quad = signature.motion_rate(el, quads[0]) if quads else r_opp

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Distance from Earth", f"{float(pos['delta']):.0f} AU")
    m2.metric("Fastest motion (opposition)", f"{float(r_opp['arcsec_per_day']):.1f} arcsec/day",
              f"{float(r_opp['pixels_per_day']):.2f} px/day, retrograde", delta_color="off")
    m3.metric("Motion at a SPHEREx visit", f"{float(r_quad['arcsec_per_day']):.1f} arcsec/day",
              f"{float(r_quad['pixels_per_day']):.2f} px/day, 90° from the Sun", delta_color="off")
    m4.metric("Orbital period", f"{el.period_years:,.0f} years")
    st.caption(f"{KIND['computed']} from a Kepler orbit and the position of Earth (astropy built-in ephemeris). "
               f"Sky position on the reference date: RA {float(pos['ra']):.2f}°, Dec {float(pos['dec']):.2f}° "
               f"(ecliptic longitude {float(pos['lon']):.1f}°, latitude {float(pos['lat']):.1f}°): this depends on "
               f"the assumed Ω and ω. {KIND['assumed']}: the observer is at the centre of the Earth "
               "(SPHEREx is ~700 km up; the error is ~0.02 arcsec at 400 AU). 1 pixel = 6.2 arcsec.")

    st.subheader("Sky track")
    show(plots.fig_sky_track(el, t0),
         "Figure 1. Apparent path over two years, as offsets from the first position (ecliptic frame, arcsec). "
         "The loops are parallax; the dashed line is the planet's own slow orbital drift. Orange dots: dates "
         "when the planet is 90° from the Sun, which is where SPHEREx scans. Computed.")

    st.subheader("Motion rate")
    show(plots.fig_rate(el, t0),
         "Figure 2. Left: apparent motion in arcsec/day through one year for the chosen orbit. It is fastest "
         "at opposition (retrograde) and behind the Sun, and slowest 90° from the Sun, where it drops to the "
         "fastest rate times the sine of the ecliptic latitude. Right: fastest rate versus distance for "
         "circular orbits (right-hand scale: SPHEREx pixels/day). Computed.")

    st.subheader("Shift between two observations")
    mode = st.radio("Which two dates?", ["Two SPHEREx-like visits (90° from the Sun, about 6 months apart)",
                                         "Custom dates"], horizontal=True)
    if mode.startswith("Two") and len(quads) >= 2:
        t1, t2 = quads[0], quads[1]
    else:
        d1, d2 = st.columns(2)
        t1 = Time(f"{d1.date_input('First date', value=pd.Timestamp(t0.iso[:10]))}T00:00:00", scale="utc")
        t2 = Time(f"{d2.date_input('Second date', value=pd.Timestamp(t0.iso[:10]) + pd.Timedelta(days=183))}T00:00:00",
                  scale="utc")
    if abs((t2 - t1).to_value("day")) < 1:
        st.warning("Pick two different dates.")
        return
    s = signature.shift_between(el, t1, t2)
    show(plots.fig_shift_split(s),
         f"Figure 3. Shift between {t1.iso[:10]} and {t2.iso[:10]}, split into the part caused by Earth's motion "
         "(parallax) and the part caused by the planet's own orbital motion. They add as vectors, so the total "
         "can be smaller or larger than the parallax part. Computed.")

    st.subheader("Distance from the shift")
    p1 = apparent_position(el, t1)
    d = signature.distance_from_shift(s["total_arcsec"], t1, t2, float(p1["ra"]), float(p1["dec"]))
    d_ref = signature.distance_from_shift_refined(s["total_arcsec"], t1, t2, float(p1["ra"]), float(p1["dec"]))
    st.latex(r"d \;\approx\; \frac{B_\perp}{\theta}\qquad\text{(}B_\perp\text{: Earth's displacement across the "
             r"line of sight, }\theta\text{: measured shift in radians)}")
    c1, c2, c3 = st.columns(3)
    c1.metric("Earth baseline B⊥", f"{d['b_perp_au']:.2f} AU")
    c2.metric("First-order distance", f"{d['distance_au']:.0f} AU",
              f"± {d['d_err_meas_au']:.1f} AU for a 1-pixel error", delta_color="off")
    c3.metric("Refined distance", f"{d_ref:.0f} AU", f"true model value: {float(p1['delta']):.0f} AU",
              delta_color="off")
    st.caption(
        f"{KIND['computed']}. If someone measured this shift ({s['total_arcsec']:.0f} arcsec) and applied the formula, "
        f"they would get {d['distance_au']:.0f} AU. {KIND['assumed']} in the first-order value: the whole shift is "
        f"parallax. Ignoring the planet's own motion can bias it by up to {100 * d['sys_frac']:.0f} %. The refined "
        "value removes most of that bias by assuming a circular, prograde orbit near the ecliptic plane. "
        "With three visits (one year) the parallax and the orbital motion can be separated without that assumption.")


# ------------------------------------------------------------------------------ page 2
def page_detectability():
    st.title("Can SPHEREx see it?")
    st.markdown(
        "The planet shines only by **reflected sunlight** at these wavelengths. Its brightness falls roughly "
        "as distance to the 4th power: sunlight fades on the way out (1/r²) and again on the way back (1/Δ²).")
    st.latex(r"m = m_\odot(1\,\mathrm{AU}) - 2.5\log_{10}\!\left[\frac{p\,(R/1\,\mathrm{AU})^2}{r^2\,\Delta^2}\right]")
    c1, c2, c3, c4 = st.columns([1.5, 1, 1, 1.1])
    scen = c1.selectbox("Albedo scenario (assumed)", list(photometry.ALBEDO_SCENARIOS), index=1)
    albedo = photometry.ALBEDO_SCENARIOS[scen]
    radius = c2.slider("Radius [Earth radii]", 1.5, 6.0, 3.0, step=0.1,
                       help="Assumed from the mass (~6 Earth masses): 2-3.7 Earth radii in our brief; "
                            "Fortney et al. 2016 quote 3-6.")
    dist = c3.slider("Distance [AU]", 200, 1000, 400, step=10)
    pess = c4.radio("Published SPHEREx depth", ["pessimistic end", "optimistic end"], horizontal=False,
                    help="Bock et al. 2025 quote a range: 19.5-19.9 AB (0.75-3.8 µm), 17.8-18.8 AB (3.8-5.0 µm).")
    pessimistic = pess.startswith("pess")

    st.info("**Conclusion.** " + detectability.verdict_text(dist, radius, albedo) +
            "\n\n**Assumptions behind it:** reflected sunlight at zero phase angle; albedo the same at all "
            "wavelengths of a band; published 5σ depth is for four all-sky passes, so one visit is taken as "
            "0.75 mag shallower; adding the ~17 channels of a band gains at most 1.5 mag.")

    show(plots.fig_heatmaps(albedo, pessimistic, mark=(dist, radius)),
         "Figure 4. Predicted AB magnitude (colour) versus distance and radius, one panel per SPHEREx detector "
         "band. Lines are detection limits; everything to the LEFT of a line is detectable that way. The white dot "
         "is the planet chosen with the sliders. Magnitudes: computed. Albedo and radius: assumed. "
         "SPHEREx depth: published by the mission (Bock et al. 2025), scaled by us to one visit (assumed).")

    st.subheader("Band by band")
    rows = detectability.summary_table(dist, radius, albedo)
    table = pd.DataFrame([{
        "Band": r["band"], "Wavelength [µm]": r["range_um"], "Predicted AB mag (computed)": round(r["mag_ab"], 2),
        "Published 5σ depth, 4 passes [AB]": r["depth_published"],
        "One visit, one channel [AB] (assumed scaling)": "{:.2f}".format(
            detectability.thresholds(r["band"], pessimistic)["one_channel_one_visit"]),
        "Verdict": r["verdict_pessimistic"] if pessimistic else r["verdict_optimistic"],
    } for r in rows])
    st.dataframe(table, hide_index=True, width="stretch")

    st.subheader("What this means in plain words")
    st.markdown(f"""
- **One channel, one visit** is the easy case: the planet would simply be a faint dot in single SPHEREx images.
- **Adding channels within one visit** helps (up to ~1.5 mag per band, ~2.5 mag for all 102 channels if the planet
  were equally bright in all of them, which a real planet is not). During a visit the planet drifts slowly
  (see page 1: SPHEREx looks 90° from the Sun, where the motion is slowest), but it can still cross a few pixels
  in a week, so the images must be shifted along a trial motion before adding.
- **Adding the four passes** is the hard part for a moving source: between two passes six months apart the planet
  jumps by hundreds of pixels (page 1), in a direction and by an amount that depend on its unknown distance and
  orbit. One has to try many possible tracks ("shift-and-stack"), which raises the false-alarm rate.
- **Out of reach** is a valid result. With a methane-dark near-infrared albedo the planet is about
  {2.5 * np.log10(0.4 / 0.03):.1f} mag fainter than with p = 0.4 and SPHEREx cannot see it at any plausible distance.
""")
    with st.expander("Does the planet's own heat matter at 3-5 µm? (thermal emission check)"):
        temps = [10, 30, 50]
        therm = pd.DataFrame({
            "Temperature [K]": temps,
            "Blackbody AB mag at 5 µm": [round(float(photometry.thermal_mag(T, radius, dist, 5.0)), 0) for T in temps],
            "Reflected AB mag at 5 µm": [round(float(photometry.reflected_mag(dist, dist, radius, albedo, lam_um=5.0)), 1)] * 3,
        })
        st.dataframe(therm, hide_index=True)
        st.caption(
            f"{KIND['computed']}: a blackbody at 10-50 K is tens to hundreds of magnitudes fainter than the reflected "
            f"light at 5 µm, so it is negligible. Equilibrium temperature at {dist} AU: "
            f"{photometry.equilibrium_temperature(dist):.0f} K. Caveat: Fortney et al. 2016 (arXiv:1604.07424) find that "
            "real cold-planet atmospheres emit many orders of magnitude more than a blackbody at 3-5 µm (and predict "
            "internal-heat temperatures up to ~35-50 K). We do not model that.")


# ------------------------------------------------------------------------------ page 4
def page_about():
    st.title("About and honesty")
    st.markdown("""
**THREX** asks: *what would a hypothetical Planet X look like to SPHEREx, and could SPHEREx see it?*
It combines a small physics model (pages 1-2) with a viewer of real SPHEREx images that makes moving
objects visible (page 3). Built for the NASA Space Apps Challenge 2026, challenge "Planet X and SPHEREx".

### Three kinds of statements
- 🟢 **MEASURED**: real data (SPHEREx images and exposure times from IRSA; positions from JPL).
- 🔵 **COMPUTED**: output of our model (orbits, parallax, magnitudes).
- 🟠 **ASSUMED**: inputs nobody knows (albedo, radius, orbit orientation, noise scaling).

**Planet X / Planet Nine is hypothetical. This tool reports no detection and no discovery.**

### Sources
- SPHEREx mission and sensitivity: Bock et al. 2025, [arXiv:2511.02985](https://arxiv.org/abs/2511.02985);
  [SPHEREx survey page](https://spherex.caltech.edu/page/survey).
- Planet Nine orbit: Brown & Batygin 2021, [arXiv:2108.09868](https://arxiv.org/abs/2108.09868).
- Planet Nine atmosphere, radius, albedo: Fortney et al. 2016, [arXiv:1604.07424](https://arxiv.org/abs/1604.07424).
- Sun's AB magnitudes: Willmer 2018, ApJS 236, 47, [online table](http://mips.as.arizona.edu/~cnaw/sun.html).
- SPHEREx data: NASA/IPAC Infrared Science Archive (IRSA),
  [data access guide](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access),
  [tutorial](https://caltech-ipac.github.io/irsa-tutorials/spherex-intro),
  [cutout service](https://irsa.ipac.caltech.edu/ibe/cutouts.html).
- Solar-system positions: [JPL Horizons](https://ssd-api.jpl.nasa.gov/doc/horizons.html),
  [Small-Body Identification](https://ssd-api.jpl.nasa.gov/doc/sb_ident.html),
  [Small-Body Database](https://ssd-api.jpl.nasa.gov/doc/sbdb.html).
- Challenge page: [Planet X and SPHEREx](https://www.spaceappschallenge.org/2026/challenges/planet-x-and-spherex/).

### What the data do and do not show
- Frames in the viewer are an evenly spaced subset of the exposures IRSA lists: Barnard's Star 68 of 358, Epsilon Indi B 69 of 458, WISE 0855-0714 69 of 258, RW Aurigae 66 of 335, R Draconis 61 of 539.
- R Draconis is saturated in most short-wavelength frames and its variability was not analysed.
- 6 Hebe saturates the detector: only its positions are meaningful (766 arcsec measured vs 765 arcsec predicted by Horizons).
- Proper motions are catalogue (SIMBAD) values, not measured by us. Photometry is rough and not independently calibrated.
- Custom RA/Dec queries work only when IRSA's search answers; presets use cached data.
- Assumptions: one visit is 0.75 mag shallower than the four-pass depth; the methane-dark albedo 0.03; the demonstration values of Omega and omega.
- WISE 0855-0714 is extremely faint below about 1.6 micron (Hubble, Vega magnitudes: 26.31 +/- 0.10 at 1.153 micron and 23.22 +/- 0.35 at 1.575 micron; Luhman and Esplin 2016, arXiv:1605.08620), far below SPHEREx's ~19.5 AB depth, so a source near its position at those wavelengths is most likely a background blend: we claim no detection.

### Assumptions and limitations
- The observer is the centre of the Earth, not the SPHEREx spacecraft (error ≤ 0.02 arcsec at 400 AU,
  but up to a few arcsec for a main-belt asteroid).
- Two-body orbit around the Sun (no planetary perturbations). Checked against JPL Horizons below.
- Brightness: reflected sunlight only, zero phase angle, one albedo per band. Real ice giants have deep
  methane absorption bands in the near infrared; the "methane-dark" scenario is only illustrative.
- The published SPHEREx depth is for a fixed source after four passes. Our scaling to one visit
  (0.75 mag) and the gain from adding channels are simple square-root-of-N estimates.
- Confusion with stars and galaxies (many per SPHEREx pixel near the Galactic plane) is not modelled.
- The image viewer does simple alignment and subtraction. Different frames have different wavelengths, so
  a difference image also shows every star whose brightness changes with wavelength.
""")
    st.subheader("Validation against JPL Horizons")
    try:
        res = validation.compare()
        show(plots.fig_validation(res),
             f"Figure 5. Difference between our computed positions and JPL Horizons for two real distant bodies, "
             f"using their orbital elements at one epoch. Horizons data downloaded {res['downloaded_utc'][:10]}. "
             "Measured: Horizons positions. Computed: ours.")
        st.dataframe(pd.DataFrame([{
            "Body": name, "Distance [AU]": round(r["distance_au"], 1),
            "Max position error [arcsec]": round(r["max_pos_err_arcsec"], 2),
            "Max rate error [arcsec/hour]": round(r["max_rate_err_arcsec_per_hour"], 4),
            "Max 180-day shift error [%]": round(100 * r["max_shift_rel_err"], 3),
        } for name, r in res["bodies"].items()]), hide_index=True)
    except Exception as exc:
        st.warning(f"Validation data not available: {exc}")
    st.subheader("Use of Artificial Intelligence")
    ai = ROOT / "AI_USE.md"
    st.markdown(ai.read_text(encoding="utf-8") if ai.exists() else "See AI_USE.md in the repository.")


# ------------------------------------------------------------------------------- router
def main():
    params = st.query_params
    default = 0                                                      # the viewer is the front page; ?ra=..&dec=.. is read there
    st.sidebar.title("THREX")
    st.sidebar.caption("Planet X and SPHEREx · NASA Space Apps 2026")
    page = st.sidebar.radio("Page", PAGES, index=default)
    st.sidebar.divider()
    st.sidebar.markdown("Planet X is **hypothetical**. Nothing here is a detection.")
    st.sidebar.markdown("\n\n".join(KIND.values()))
    if page == PAGES[0]:
        try:
            import app_realdata
            app_realdata.page_real_data(show, KIND)
        except Exception as exc:                       # never crash the whole app
            st.error(f"The viewer could not be shown: {type(exc).__name__}: {exc}")
    elif page == PAGES[1]:
        page_signature()
    elif page == PAGES[2]:
        page_detectability()
    else:
        page_about()


main()
