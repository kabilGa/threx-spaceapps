"""All figures of the app and of the README (matplotlib, so the app and the saved PNGs match).

Every figure has axis labels with units. Captions are returned by the app next to each figure.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from matplotlib.colors import LinearSegmentedColormap

from . import signature, detectability
from .constants import SPHEREX_BANDS, SPHEREX_PIXEL_ARCSEC
from .observer import to_time

# One small palette, used the same way everywhere.
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"
# sequential map, one hue: dark blue = bright planet, pale = faint planet
SEQ = LinearSegmentedColormap.from_list("threx_blue", ["#0d366b", "#256abf", "#6da7ec", "#cde2fb", "#f4f8fd"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "lines.linewidth": 2.0, "legend.frameon": False, "figure.dpi": 110, "savefig.dpi": 160,
})


def fig_sky_track(el, t0, years=2.0):
    """Apparent path on the sky over `years`, as offsets from the first position."""
    tr = signature.sky_track(el, t0, n_days=365.25 * years, n=600)
    lat0 = np.radians(tr["lat"][0])
    dx = ((tr["lon"] - tr["lon"][0] + 180) % 360 - 180) * 3600 * np.cos(lat0)   # arcsec, along ecliptic
    dy = (tr["lat"] - tr["lat"][0]) * 3600                                       # arcsec, across ecliptic
    fig, ax = plt.subplots(figsize=(9.0, 3.9))
    ax.plot(dx, dy, color=BLUE, label="seen from Earth (parallax loops + slow orbital drift)")
    # after a whole number of years Earth is back at the same place, so first -> last
    # position shows the planet's own orbital drift alone
    ax.plot([dx[0], dx[-1]], [dy[0], dy[-1]], color=MUTED, ls="--", lw=1.2,
            label="planet's own orbital drift alone")
    # mark the SPHEREx-like visits (body 90 deg from the Sun)
    t0 = to_time(t0)
    first = True
    for yr in range(int(np.ceil(years))):
        for tq in signature.quadrature_times(el, t0 + yr * 365.25 * u.day):
            day = (tq - t0).to_value("day")
            if day <= 365.25 * years:
                i = int(np.argmin(np.abs(tr["days"] - day)))
                ax.plot(dx[i], dy[i], "o", ms=8, color=ORANGE, mec=SURFACE, mew=1.5,
                        label="90 deg from the Sun (when SPHEREx visits)" if first else None)
                first = False
    ax.plot(dx[0], dy[0], "s", ms=7, color=INK, mec=SURFACE, mew=1.5, label="start")
    ax.set_xlabel("offset along the ecliptic [arcsec]  (east to the left)")
    ax.set_ylabel("offset across the ecliptic [arcsec]")
    ax.invert_xaxis()
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=8.5)
    ax.set_title(f"Apparent path over {years:g} years")
    fig.tight_layout()
    return fig


def fig_rate(el, t0):
    """Left: motion rate through one year. Right: rate at opposition versus distance."""
    from . import orbits
    yr = signature.rate_through_year(el, t0, n=366)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.4), gridspec_kw={"width_ratios": [1.5, 1]})
    a1.plot(yr["days"], yr["arcsec_per_day"], color=BLUE)
    t0 = to_time(t0)
    for k, tq in enumerate(signature.quadrature_times(el, t0)):
        d = (tq - t0).to_value("day")
        a1.axvline(d, color=ORANGE, lw=1.4, ls="--", label="90 deg from the Sun (SPHEREx visit)" if k == 0 else None)
    i_opp = int(np.argmax(yr["elongation_deg"]))
    i_con = int(np.argmin(yr["elongation_deg"]))
    for i, txt in ((i_opp, "opposition (retrograde)"), (i_con, "behind the Sun")):
        a1.annotate(txt, (yr["days"][i], yr["arcsec_per_day"][i]), xytext=(0, 6),
                    textcoords="offset points", ha="center", fontsize=8.5, color=INK2)
    a1.set_xlabel("days after the reference date")
    a1.set_ylabel("apparent motion [arcsec / day]")
    a1.set_ylim(0, yr["arcsec_per_day"].max() * 1.18)
    a1.legend(loc="lower center", fontsize=8.5)
    a1.set_title("Motion rate through the year")
    # rate at opposition versus distance (circular orbits)
    dist = np.array([200, 250, 300, 400, 500, 600, 800, 1000.0])
    rates = []
    for d in dist:
        e2 = orbits.circular_orbit(d, inc=el.inc, node=el.node, epoch_jd=t0.tdb.jd)
        rates.append(float(signature.motion_rate(e2, signature.opposition_time(e2, t0, n=366))["arcsec_per_day"]))
    a2.plot(dist, rates, "o-", color=BLUE, ms=5)
    a2.set_xlabel("distance [AU]")
    a2.set_ylabel("rate at opposition [arcsec / day]")
    a2.set_ylim(bottom=0)
    sec = a2.secondary_yaxis("right", functions=(lambda x: x / SPHEREX_PIXEL_ARCSEC, lambda x: x * SPHEREX_PIXEL_ARCSEC))
    sec.set_ylabel("SPHEREx pixels / day")          # same quantity in other units (not a second scale)
    a2.set_title("Fastest rate versus distance")
    fig.tight_layout()
    return fig


def fig_shift_split(shift):
    """Horizontal bars: parallax part, orbital part, total, in SPHEREx pixels."""
    labels = ["Parallax (Earth moved)", "Orbital (planet moved)", "Total measured shift"]
    vals = [shift["parallax_pixels"], shift["orbital_pixels"], shift["total_pixels"]]
    arcs = [shift["parallax_arcsec"], shift["orbital_arcsec"], shift["total_arcsec"]]
    colors = [BLUE, ORANGE, INK2]
    fig, ax = plt.subplots(figsize=(7.4, 2.2))
    ax.set_axisbelow(True)
    ax.barh(labels[::-1], vals[::-1], color=colors[::-1], height=0.55)
    for i, (v, a) in enumerate(zip(vals[::-1], arcs[::-1])):
        ax.text(v, i, f"  {v:.0f} px  ({a:.0f} arcsec)", va="center", fontsize=9, color=INK)
    ax.set_xlim(0, max(vals) * 1.35)
    ax.set_xlabel("shift between the two epochs [SPHEREx pixels, 6.2 arcsec each]")
    ax.grid(axis="y", visible=False)
    ax.set_title(f"Shift over {shift['baseline_days']:.0f} days")
    fig.tight_layout()
    return fig


def fig_heatmaps(albedo, pessimistic=True, mark=None):
    """Six panels (one per SPHEREx detector band): predicted AB magnitude versus distance
    and radius, with the SPHEREx limits drawn as contours. mark = (distance, radius)."""
    dist = np.linspace(200, 1000, 161)
    radii = np.linspace(1.5, 6.0, 91)
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.0), sharex=True, sharey=True)
    vmin, vmax = 15.0, 27.0
    for ax, b in zip(axes.ravel(), SPHEREX_BANDS):
        grid = detectability.magnitude_grid(dist, radii, albedo, b["band"])
        im = ax.pcolormesh(dist, radii, grid, cmap=SEQ, vmin=vmin, vmax=vmax, shading="auto", rasterized=True)
        th = detectability.thresholds(b["band"], pessimistic)
        for key, ls, col in (("one_channel_one_visit", "-", ORANGE),
                             ("band_stack_one_visit", "--", ORANGE),
                             ("band_stack_four_visits", ":", INK)):
            if grid.min() < th[key] < grid.max():
                ax.contour(dist, radii, grid, levels=[th[key]], colors=[col], linestyles=[ls], linewidths=1.8)
        if mark:
            ax.plot(mark[0], mark[1], "o", ms=7, color="white", mec=INK, mew=1.5)
        ax.set_title(f"Band {b['band']}: {b['lam_min']:.2f}-{b['lam_max']:.2f} µm", fontsize=10)
        ax.set_xticks([200, 400, 600, 800])
        ax.grid(False)
    for ax in axes[-1]:
        ax.set_xlabel("distance [AU]")
    for ax in axes[:, 0]:
        ax.set_ylabel("planet radius [Earth radii]")
    fig.subplots_adjust(left=0.07, right=0.88, top=0.90, bottom=0.17, wspace=0.08, hspace=0.22)
    cax = fig.add_axes([0.90, 0.17, 0.018, 0.73])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("predicted AB magnitude (darker = brighter)")
    cb.ax.invert_yaxis()
    handles = [plt.Line2D([], [], color=ORANGE, ls="-", lw=1.8, label="limit: one channel, one visit"),
               plt.Line2D([], [], color=ORANGE, ls="--", lw=1.8, label="limit: ~17 channels of the band added, one visit"),
               plt.Line2D([], [], color=INK, ls=":", lw=1.8, label="limit: channels added + 4 visits shift-and-stacked")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=8.5, bbox_to_anchor=(0.48, 0.02))
    fig.suptitle(f"Predicted AB magnitude for albedo p = {albedo:g}.  Left of a line = detectable that way.",
                 x=0.07, ha="left", fontsize=10.5, fontweight="bold")
    return fig


def fig_validation(results):
    """Position error of our code against JPL Horizons for each validation body."""
    fig, ax = plt.subplots(figsize=(7.4, 3.0))
    for (name, r), col in zip(results["bodies"].items(), (BLUE, ORANGE)):
        days = r["jd"] - r["jd"][len(r["jd"]) // 2]
        ax.plot(days, r["pos_err_arcsec"], "o-", ms=5, color=col, label=f"{name}, {r['distance_au']:.0f} AU")
    ax.axhline(SPHEREX_PIXEL_ARCSEC / 10, color=MUTED, ls="--", lw=1.2)
    ax.text(ax.get_xlim()[0], SPHEREX_PIXEL_ARCSEC / 10, " 1/10 of a SPHEREx pixel", va="bottom", fontsize=8.5, color=INK2)
    ax.set_xlabel("days from the epoch of the orbital elements (2026-02-19)")
    ax.set_ylabel("position error [arcsec]")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8.5, loc="upper center")
    ax.set_title("Our code minus JPL Horizons")
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------- real SPHEREx data
DETECTOR_COLORS = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a", 4: "#eda100", 5: "#e87ba4", 6: "#6250d6"}


def _dates(mjd):
    from astropy.time import Time
    return Time(np.atleast_1d(mjd), format="mjd", scale="utc").to_datetime()


def fig_timeline(df, visits=None, pair=None, title="When did SPHEREx look here?"):
    """Real exposure times (x) against the wavelength of each exposure (y), by detector."""
    fig, ax = plt.subplots(figsize=(9.6, 3.6))
    lam = 0.5 * (df["lam_min_um"] + df["lam_max_um"])
    for det, col in DETECTOR_COLORS.items():
        m = df["detector"] == det
        if m.any():
            ax.scatter(_dates(df.loc[m, "t_mid"]), lam[m], s=14, color=col, label=f"detector {det}", zorder=3)
    for k, v in enumerate(visits or []):
        a, b = _dates([v["start_mjd"] - 2, v["end_mjd"] + 2])
        ax.axvspan(a, b, color=GRID, alpha=0.6, zorder=1)
        ax.text(_dates(v["mid_mjd"])[0], 5.25, f"visit {k + 1}\n{v['n']} exp.", ha="center", va="bottom",
                fontsize=8, color=INK2)
    if pair is not None:
        a, b = _dates(pair)
        ax.annotate("", xy=(b, 0.45), xytext=(a, 0.45), arrowprops=dict(arrowstyle="<->", color=INK, lw=1.3))
        ax.text(a + (b - a) / 2, 0.52, f"longest baseline: {pair[1] - pair[0]:.0f} days", ha="center",
                va="bottom", fontsize=8.5, color=INK)
    ax.set_ylim(0.3, 6.0)
    ax.set_ylabel("wavelength of the exposure\nat this position [µm]")
    ax.set_xlabel("exposure date (UTC)")
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8.5, markerscale=1.6)
    ax.set_title(title)
    fig.autofmt_xdate(rotation=0, ha="center")
    fig.tight_layout()
    return fig


def fig_real_shift(t_mid, ra_deg, dec_deg, distances=(300.0, 400.0, 600.0)):
    """Offset a hypothetical planet at this position would show at the REAL exposure dates."""
    from astropy.time import Time
    t = Time(np.sort(np.asarray(t_mid, dtype=float)), format="mjd", scale="utc")
    fine = Time(np.linspace(t.mjd.min(), t.mjd.max(), 400), format="mjd", scale="utc")
    fig, ax = plt.subplots(figsize=(9.6, 3.4))
    for d, col in zip(distances, (BLUE, ORANGE, AQUA)):
        curve = signature.track_from_position(ra_deg, dec_deg, d, t[0], fine)["offset_pixels"]
        pts = signature.track_from_position(ra_deg, dec_deg, d, t[0], t)["offset_pixels"]
        ax.plot(fine.to_datetime(), curve, color=col, lw=1.5, label=f"if at {d:g} AU")
        ax.plot(t.to_datetime(), pts, "o", ms=4, color=col, mec=SURFACE, mew=0.5)
    ax.set_ylabel("offset from the position at the\nfirst exposure [SPHEREx pixels]")
    ax.set_xlabel("date (UTC); dots = real SPHEREx exposures of this field")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8.5, loc="upper left")
    ax.set_title("A hypothetical Planet X at this position, on the real SPHEREx dates")
    fig.autofmt_xdate(rotation=0, ha="center")
    fig.tight_layout()
    return fig


def _show_frame(ax, bundle, k, show_bad=False, lims=None):
    from . import blink
    m = bundle["meta"]
    h = m["half_size_arcsec"]
    if lims is None:
        img = blink.normalise(bundle["frames"][k], bundle["bad"][k]) ** 0.5
    else:
        img = np.clip((np.nan_to_num(bundle["frames"][k], nan=lims[0]) - lims[0]) / max(lims[1] - lims[0], 1e-9), 0, 1) ** 0.5
    ax.imshow(img, origin="lower", cmap="gray", extent=[h, -h, -h, h], interpolation="nearest")
    if show_bad:
        bad = np.ma.masked_where(~bundle["bad"][k], np.ones_like(img))
        ax.imshow(bad, origin="lower", cmap="autumn", alpha=0.55, extent=[h, -h, -h, h], interpolation="nearest")
    ax.grid(False)
    ax.set_title(blink.frame_label(m["frames"][k]), fontsize=8, fontweight="normal")


def fig_frames(bundle, indices=None, markers=None, show_bad=False, ncols=4, lims=None):
    """Grid of aligned, north-up frames. markers: per-frame (east, north) offsets in arcsec."""
    idx = list(range(len(bundle["frames"]))) if indices is None else list(indices)
    ncols = min(ncols, len(idx))
    nrows = int(np.ceil(len(idx) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.9 * ncols, 3.05 * nrows), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, k in zip(axes.ravel(), idx):
        ax.axis("on")
        _show_frame(ax, bundle, k, show_bad, lims)
        if markers is not None and markers[k] is not None:
            ax.plot(*markers[k], "o", ms=16, mfc="none", mec=ORANGE, mew=1.6)
        ax.tick_params(labelsize=7)
    for ax in axes[-1]:
        ax.set_xlabel("east offset [arcsec]", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("north offset [arcsec]", fontsize=8)
    fig.tight_layout()
    return fig


def fig_pair(bundle, i, j, predicted=None, measured=None):
    """Frame i, frame j and their difference. predicted / measured: dicts {frame index:
    (east, north) offset in arcsec} drawn as an orange circle / a blue cross."""
    from . import blink
    h = bundle["meta"]["half_size_arcsec"]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.9))
    for ax, k in zip(axes[:2], (i, j)):
        _show_frame(ax, bundle, k)
    diff = blink.difference(bundle, i, j)
    axes[2].imshow(diff, origin="lower", cmap="RdBu_r", vmin=-0.6, vmax=0.6, extent=[h, -h, -h, h],
                   interpolation="nearest")
    axes[2].grid(False)
    axes[2].set_title("second minus first (blue = only in first, red = only in second)", fontsize=8,
                      fontweight="normal")
    for ax, ks in zip(axes, ((i,), (j,), (i, j))):
        for k in ks:
            if predicted and k in predicted:
                ax.plot(*predicted[k], "o", ms=18, mfc="none", mec=ORANGE, mew=1.6,
                        label="JPL Horizons prediction" if k == ks[0] else None)
            if measured and k in measured:
                ax.plot(*measured[k], "+", ms=11, color=BLUE, mew=1.8,
                        label="measured centre" if k == ks[0] else None)
        ax.set_xlabel("east offset [arcsec]", fontsize=8.5)
        ax.tick_params(labelsize=7.5)
    axes[0].set_ylabel("north offset [arcsec]", fontsize=8.5)
    if predicted or measured:
        axes[0].legend(loc="lower left", fontsize=7.5, frameon=True, facecolor="white", framealpha=0.85)
    fig.tight_layout()
    return fig
