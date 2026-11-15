# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Page 3 of the THREX app: real SPHEREx data.

A. Schedule : when SPHEREx really looked at a sky position (IRSA SIA2), and what a
              hypothetical Planet X at that position would do on those exact dates.
B. Blink    : aligned cutouts of that position shown as a short sequence (the short animation
              asked for by the challenge), with flags and a difference image.
C. Proof    : the same idea catching a KNOWN asteroid, checked against JPL.
Live data when IRSA answers, cached data (with its date) otherwise. Never invented data.
"""
import os

import numpy as np
import pandas as pd
import streamlit as st
from astropy.time import Time

from threx import blink, plots, signature, spherex_io, targets
from threx.constants import SPHEREX_PIXEL_ARCSEC

OFFLINE = os.environ.get("THREX_OFFLINE") == "1"


@st.cache_data(ttl=3600, show_spinner=False)
def _exposures(name, ra, dec, live):
    return spherex_io.get_exposures(name, ra, dec, live=live)


@st.cache_data(show_spinner=False)
def _checks():
    import json
    f = spherex_io.CACHE_DIR / "target_checks.json"
    return json.loads(f.read_text()) if f.exists() else {}


def _status(meta, what):
    if meta.get("mode") == "live":
        st.success(f"🟢 **LIVE** {what} from IRSA, fetched {meta.get('queried_utc', meta.get('built_utc', ''))[:19]} UTC.")
    else:
        when = meta.get("queried_utc", meta.get("built_utc", "unknown date"))[:10]
        why = (" Live service unavailable, showing cached targets (" + meta["live_error"][:90] + ").") if meta.get("live_error") else ""
        st.warning(f"🟡 **CACHED on {when}**: {what} saved earlier from the IRSA archive (SPHEREx files on the IRSA cloud mirror).{why}")


def _pick_target():
    params = st.query_params
    presets = targets.all_targets()
    names = list(presets) + ["Custom RA/Dec"]
    custom_from_url = "ra" in params and "dec" in params
    c1, c2, c3 = st.columns([1.6, 1, 1])
    name = c1.selectbox("Sky position", names, index=len(names) - 1 if custom_from_url else 0,
                        help="Presets have cached data. 'Custom' queries IRSA live. "
                             "You can also open the app with ?ra=..&dec=.. in the URL.")
    if name == "Custom RA/Dec":
        try:
            ra0, dec0 = float(params.get("ra", 60.0)), float(params.get("dec", 20.0))
        except ValueError:
            ra0, dec0 = 60.0, 20.0
        ra = c2.number_input("RA [deg, ICRS]", 0.0, 360.0, ra0, format="%.5f")
        dec = c3.number_input("Dec [deg, ICRS]", -90.0, 90.0, dec0, format="%.5f")
        return f"custom_{ra:.4f}_{dec:.4f}", ra, dec, {"note": "User-chosen position."}
    t = presets[name]
    c2.metric("RA [deg]", f"{t['ra']:.5f}")
    c3.metric("Dec [deg]", f"{t['dec']:.5f}")
    return name, t["ra"], t["dec"], t


def _schedule(name, ra, dec, df, show, KIND):
    st.header("A · When did SPHEREx look here?")
    visits = spherex_io.find_visits(df)
    t = df["t_mid"].to_numpy(dtype=float)
    pair = (float(t.min()), float(t.max()))
    show(plots.fig_timeline(df, visits, pair),
         f"Figure 6. Each dot is one real SPHEREx exposure covering this position ({len(df)} exposures in the cached subset; the archive lists more). "
         "x: exposure date (UTC, converted from MJD). y: wavelength range centre of that exposure's detector "
         "band at this position. Grey bands: visits. The arrow marks the pair of exposures with the longest "
         "time baseline. Measured (IRSA exposure metadata).")
    vt = pd.DataFrame([{
        "Visit": k + 1, "First exposure (UTC)": Time(v["start_mjd"], format="mjd").iso[:16],
        "Last exposure (UTC)": Time(v["end_mjd"], format="mjd").iso[:16],
        "Exposures": v["n"], "Duration [days]": round(v["duration_days"], 1),
    } for k, v in enumerate(visits)])
    st.dataframe(vt, hide_index=True, width="stretch")
    # solar elongation of the field at the real exposure times (measured dates, computed angle)
    from threx.observer import earth_heliocentric, unit_vector, separation_arcsec
    tt = Time(t, format="mjd", scale="utc")
    elong = separation_arcsec(unit_vector(ra, dec), -earth_heliocentric(tt)) / 3600
    st.caption(f"{KIND['measured']}: exposure dates. {KIND['computed']}: at those dates this field was "
               f"{elong.min():.0f}° to {elong.max():.0f}° from the Sun (median {np.median(elong):.0f}°), "
               "as expected from the SPHEREx scan strategy (about 90° from the Sun).")

    st.subheader("A hypothetical Planet X at this position, on these real dates")
    dist = st.slider("Assumed distance of the hypothetical planet [AU]", 200, 1000, 400, step=10, key="real_dist")
    show(plots.fig_real_shift(t, ra, dec, distances=(max(200, dist - 100), dist, dist + 200)),
         "Figure 7. Offset (in SPHEREx pixels) that a hypothetical planet sitting at this position at the first "
         "exposure would show at each later real exposure, for three assumed distances. Dates: measured. "
         "Offsets: computed (circular prograde orbit near the ecliptic plane: assumed).")
    tr = signature.track_from_position(ra, dec, dist, tt[np.argmin(t)], tt[np.argmax(t)])
    drifts = []
    for v in visits:
        a = Time(v["start_mjd"], format="mjd", scale="utc")
        b = Time(v["end_mjd"], format="mjd", scale="utc")
        drifts.append(float(signature.track_from_position(ra, dec, dist, a, b)["offset_pixels"][0]))
    c1, c2, c3 = st.columns(3)
    c1.metric("Longest baseline", f"{pair[1] - pair[0]:.0f} days")
    c2.metric(f"Shift over that baseline at {dist} AU", f"{float(tr['offset_pixels'][0]):.0f} px",
              f"{float(tr['offset_arcsec'][0]):.0f} arcsec", delta_color="off")
    c3.metric("Drift inside one visit", f"{max(drifts):.1f} px" if drifts else "n/a",
              "largest of the visits", delta_color="off")
    st.caption(f"{KIND['computed']} for {KIND['assumed'].split(' ', 1)[1]} distance {dist} AU. A drift of more than "
               "about one pixel inside a visit means the channels of that visit cannot be added without "
               "shifting them first. No planet is known at this position: this is a what-if.")


def _blink(name, ra, dec, df, live, show, KIND):
    st.header("B · Blink the real images")
    st.markdown(
        "SPHEREx images of the same patch, aligned north-up and shown one after the other. "
        "A fixed star stays put; something that moves jumps. Each exposure has a **different wavelength** "
        "at a given sky position, so stars also change brightness from frame to frame: that is colour, not "
        "variability.")
    key = f"bundle_{name}"
    bundle = st.session_state.get(key)
    if bundle is None:
        try:
            bundle = blink.load_bundle(name)
        except spherex_io.DataUnavailable:
            bundle = None
    if live and not OFFLINE:
        with st.expander("Fetch images live from IRSA" + ("" if bundle is None else " (replace the cached set)"),
                         expanded=bundle is None):
            c1, c2, c3 = st.columns(3)
            dets = c1.multiselect("Detectors (1-6, short to long wavelength)", [1, 2, 3, 4, 5, 6], default=[2],
                                  help="Blinking works best inside one detector band.")
            n = c2.slider("Number of frames", 2, 24, 8)
            half = c3.slider("Half field size [arcsec]", 60, 300, 150, step=30)
            if st.button("Fetch cutouts from IRSA", type="primary"):
                rows = blink.select_frames(df, n, dets)
                if len(rows) < 2:
                    st.warning("Fewer than two exposures match. Choose other detectors.")
                else:
                    bar = st.progress(0.0, text="Downloading cutouts from IRSA ...")
                    try:
                        bundle = blink.build_bundle(
                            name, ra, dec, rows, half_size_arcsec=float(half),
                            progress=lambda k, tot, row: bar.progress(k / tot, text=f"Cutout {k + 1} of {tot} ..."))
                        st.session_state[key] = bundle
                    except Exception as exc:
                        st.error(f"IRSA cutouts could not be downloaded: {type(exc).__name__}: {str(exc)[:200]}")
                    bar.empty()
    if bundle is None:
        st.info("No images for this position yet. Use 'Fetch images live from IRSA' above "
                "(needs a network connection and IRSA being up), or pick a preset with cached images.")
        return None
    _status(bundle["meta"], f"{len(bundle['frames'])} image cutouts")
    frames = bundle["meta"]["frames"]
    proof = bundle["meta"].get("proof")
    markers_pix = markers_arc = None
    if proof:
        markers_arc = [(f["pred_east_arcsec"], f["pred_north_arcsec"]) for f in proof["frames"]]
        markers_pix = [blink.offset_to_pixel(bundle, f["pred_ra"], f["pred_dec"]) for f in proof["frames"]]
    cat = blink.catalogue_offsets(bundle)
    c0, c1_, c2_ = st.columns(3)
    dets = sorted({f["detector"] for f in frames})
    pick = c0.multiselect("Detectors shown", dets, default=dets,
                          help="Detectors 1-3 are short wavelength (0.75-5 µm in three bands); 4-6 the longer ones. "
                               "Detectors 1&4, 2&5, 3&6 look at the same sky at the same time but at different wavelengths.")
    stretch = c1_.radio("Brightness stretch", ["per frame", "fixed for all frames"], horizontal=True,
                        help="Per frame: every frame scaled to its own range (best to see faint stars). "
                             "Fixed: one scale for all frames (best to see real brightness changes).")
    mark_cat = c2_.checkbox("Mark the catalogue position", value=cat is not None, disabled=cat is None,
                            help="Orange circle: where SIMBAD (position + proper motion) puts the target at that date.")
    idx = [q for q, f in enumerate(frames) if f["detector"] in pick] or list(range(len(frames)))
    if cat and mark_cat and not proof:
        markers_arc = cat
        n_, g_ = bundle["frames"].shape[-1], bundle["meta"]["grid_pix_arcsec"]
        markers_pix = [(n_ // 2 - e / g_, n_ // 2 + n / g_) for e, n in cat]       # east is to the left
    lims = None
    if stretch.startswith("fixed"):
        good = np.concatenate([bundle["frames"][q][~bundle["bad"][q]].ravel() for q in idx])
        lims = tuple(np.nanpercentile(good, [30, 99.7]))
    c1, c2 = st.columns([1, 1.25])
    with c1:
        gif_path = spherex_io.CACHE_DIR / f"gif_{spherex_io._slug(name)}.gif"
        if gif_path.exists() and mark_cat and not proof and stretch.startswith("per") and len(pick) == len(dets):
            st.image(str(gif_path), width="stretch")                    # ready-made animation
        else:
            st.image(blink.make_gif(bundle, indices=idx, markers=markers_pix if (proof or mark_cat) else None),
                     width="stretch")
        st.caption("Animation of the aligned frames (north up, east left), in time order. Label: frame number, "
                   "exposure time, wavelength at the field centre, detector. "
                   + ("Orange circle: expected position. " if (proof or mark_cat) else "") + "Measured.")
    with c2:
        pos = st.slider("Frame (press the arrows or drag)", 1, len(idx), 1) - 1
        k = idx[pos]
        show_bad = st.checkbox("Show flagged pixels (FLAGS extension) in colour", value=False)
        show(plots.fig_frames(bundle, [k], markers=markers_arc if (proof or mark_cat) else None,
                              show_bad=show_bad, ncols=1, lims=lims),
             f"Frame {pos + 1} of {len(idx)}: {blink.frame_label(frames[k])}; exposure ID {frames[k]['obs_id']}. "
             f"{100 * frames[k]['frac_flagged']:.1f} % of its pixels are flagged. Black areas inside a bright star "
             "are saturated/flagged pixels. Measured.")
    t_all = [frames[q]["t_mid"] for q in idx]
    st.caption(f"Time span covered by the frames shown: {Time(min(t_all), format='mjd').iso[:16]} to "
               f"{Time(max(t_all), format='mjd').iso[:16]} UTC ({max(t_all) - min(t_all):.1f} days); wavelengths "
               f"{min(frames[q]['lam_um'] for q in idx):.2f}-{max(frames[q]['lam_um'] for q in idx):.2f} µm.")
    sim = bundle["meta"].get("simbad")
    if sim:
        em = blink.expected_motion(bundle, idx[0], idx[-1])
        c1, c2, c3 = st.columns(3)
        c1.metric("Catalogue proper motion (SIMBAD)", f"{np.hypot(sim['pmra'], sim['pmdec']) / 1000:.2f} arcsec/yr")
        c2.metric("Expected motion across these frames", f"{em[0]:.2f} arcsec", f"{em[1]:.2f} pixel", delta_color="off")
        c3.metric("Frames span", f"{em[2]:.0f} days")
        verdict = ("**Not visible as motion**: the expected shift is below one SPHEREx pixel (6.2 arcsec), so within "
                   "these frames the target stays in the same place. Seeing it move needs frames about "
                   f"{6.2 / (np.hypot(sim['pmra'], sim['pmdec']) / 1000):.1f} year(s) apart or more."
                   if em[1] < 1 else
                   "**Visible as motion** if the target is bright enough: the expected shift exceeds one pixel.")
        st.markdown(verdict)
    if name in targets.COMMENTS:
        st.info("**What we see (and do not see):** " + targets.COMMENTS[name])
    chk = _checks().get(name)
    if chk:
        st.caption("Rough aperture photometry at the catalogue position (2-pixel radius, no aperture correction; "
                   "flux calibration not independently verified): " + "; ".join(
            f"detector {d}: {v['n_snr_gt5']}/{v['n']} frames with S/N > 5" for d, v in chk["per_detector"].items()) + ".")
    st.warning("Colour is not variability: every exposure samples a different wavelength at a given position, so a "
               "star can look brighter or fainter from frame to frame only because of its colour. "
               "Genuine variability (for example RW Aur, R Dra) needs frames of the same wavelength.")
    st.subheader("Subtract two frames")
    a, b = st.columns(2)
    i = a.selectbox("First frame", range(len(frames)), index=0, format_func=lambda q: f"{q + 1}: {blink.frame_label(frames[q])}")
    j = b.selectbox("Second frame", range(len(frames)), index=len(frames) - 1,
                    format_func=lambda q: f"{q + 1}: {blink.frame_label(frames[q])}")
    if i != j:
        pred = {q: markers_arc[q] for q in (i, j)} if markers_arc else None
        meas = None
        if proof:
            meas = {q: (proof["frames"][q]["meas_east_arcsec"], proof["frames"][q]["meas_north_arcsec"])
                    for q in (i, j) if proof["frames"][q].get("meas_east_arcsec") is not None}
        show(plots.fig_pair(bundle, i, j, pred, meas),
             "Figure 8. Two real SPHEREx frames of the same patch and their difference (each frame scaled to its "
             "own brightness range first). Fixed stars cancel only approximately because the two frames have "
             "different wavelengths. Measured (pixels); simple processing.")
    with st.expander("Frame list and archive links"):
        st.dataframe(pd.DataFrame([{
            "Frame": q + 1, "Exposure mid-time (UTC)": str(f["time_utc"])[:19], "Wavelength [µm]": round(f["lam_um"], 3),
            "Detector": f["detector"], "Flagged pixels [%]": round(100 * f["frac_flagged"], 1),
            "Release": f["collection"], "IRSA file": f["access_url"],
        } for q, f in enumerate(frames)]), hide_index=True, width="stretch")
    return bundle


def _proof(bundle, show, KIND):
    proof = bundle["meta"].get("proof") if bundle else None
    if not proof:
        return
    st.header("C · Proof on real pixels: a known asteroid")
    st.markdown(
        f"**{proof['name']}** is a known asteroid. JPL Horizons predicts where it was at the time of each SPHEREx "
        "exposure (orange circles above). If the moving-object idea works, a point of light must sit inside each "
        "circle and move from frame to frame, while the stars stay put.")
    rows = []
    for q, f in enumerate(proof["frames"]):
        rows.append({
            "Frame": q + 1, "Exposure (UTC)": str(bundle["meta"]["frames"][q]["time_utc"])[:19],
            "Wavelength [µm]": round(bundle["meta"]["frames"][q]["lam_um"], 3),
            "Predicted RA [deg]": round(f["pred_ra"], 5), "Predicted Dec [deg]": round(f["pred_dec"], 5),
            "Measured − predicted [arcsec]": None if f.get("miss_arcsec") is None else round(f["miss_arcsec"], 1),
            "Flux [µJy]": None if f.get("flux_ujy") is None else round(f["flux_ujy"]),
            "S/N": None if f.get("snr") is None else round(f["snr"], 1),
            "AB mag": None if f.get("mag_ab") is None else round(f["mag_ab"], 2),
        })
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    d = proof.get("displacement")
    if d:
        c1, c2, c3 = st.columns(3)
        c1.metric("Time between first and last frame", f"{d['dt_hours']:.2f} h")
        c2.metric("Measured displacement", f"{d['measured_arcsec']:.1f} arcsec",
                  f"{d['measured_arcsec'] / SPHEREX_PIXEL_ARCSEC:.1f} px", delta_color="off")
        c3.metric("Predicted by JPL Horizons", f"{d['predicted_arcsec']:.1f} arcsec",
                  f"{d['predicted_arcsec'] / SPHEREX_PIXEL_ARCSEC:.1f} px", delta_color="off")
    st.caption(
        f"{KIND['measured']}: SPHEREx pixels (IRSA), centres and fluxes measured by us with a simple centroid and "
        "a 2-pixel-radius aperture (no aperture correction). Predicted positions: JPL Horizons, geocentric "
        "(SPHEREx is in low Earth orbit, which can shift a main-belt asteroid by a few arcsec). "
        "6 Hebe is bright enough to saturate the detector (tens of aperture pixels are flagged per frame), so only its "
        "positions are meaningful: the flux and AB columns are unreliable.")
    ident = proof.get("sb_ident")
    if ident:
        st.subheader("Independent check: JPL Small-Body Identification")
        st.markdown(ident["summary"])
        if ident.get("rows"):
            st.dataframe(pd.DataFrame(ident["rows"]), hide_index=True, width="stretch")
    if proof.get("sbdb"):
        s = proof["sbdb"]
        st.caption("JPL Small-Body Database: " + "; ".join(f"{k}: {v}" for k, v in s.items()))
    st.info("This asteroid was already known. Nothing was discovered here: the point is that the method "
            "(predict, align, blink, subtract) recovers a real moving object in real SPHEREx pixels.")


def page_real_data(show, KIND):
    st.title("SPHEREx sky viewer")
    st.markdown("Type a sky position, see when SPHEREx observed it, and blink the real images to spot what moves. "
                "Data: SPHEREx spectral images from the NASA/IPAC Infrared Science Archive (IRSA).")
    name, ra, dec, info = _pick_target()
    if info.get("note"):
        st.caption(info["note"])
    live = st.toggle("Ask IRSA live (falls back to cached data if IRSA does not answer)",
                     value=not OFFLINE and name.startswith("custom"), disabled=OFFLINE)
    with st.spinner("Asking IRSA which exposures cover this position ..." if live else "Loading cached data ..."):
        try:
            df, meta = _exposures(name, ra, dec, live)
        except spherex_io.DataUnavailable as exc:
            try:                                    # no exposure list: rebuild it from the cached frames
                bd = blink.load_bundle(name)
                df = pd.DataFrame([{**{k: f[k] for k in ("obs_id", "t_mid", "time_utc", "detector", "access_url")},
                                    "obs_collection": f["collection"], "lam_min_um": f["lam_um"] - f["dlam_um"] / 2,
                                    "lam_max_um": f["lam_um"] + f["dlam_um"] / 2} for f in bd["meta"]["frames"]])
                meta = {"mode": "cached", "queried_utc": bd["meta"]["built_utc"]}
            except Exception:
                st.error(f"No data for this position: {exc}. IRSA's live search is not answering right now; pick a "
                         "preset target (cached data) or try again later.")
                return
    _status(meta, f"exposure list ({len(df)} exposures)")
    if len(df) == 0:
        st.info("IRSA returned no SPHEREx exposure for this position (not yet observed or not yet released).")
        return
    _schedule(name, ra, dec, df, show, KIND)
    bundle = _blink(name, ra, dec, df, live, show, KIND)
    _proof(bundle, show, KIND)
