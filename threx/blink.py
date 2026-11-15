"""Blink viewer data: a short, aligned sequence of real SPHEREx cutouts of one sky patch.

A "bundle" is everything the viewer needs for one target:
  frames[i]  : north-up image in micro-Jy/pixel (zodiacal light and median removed)
  sigma[i]   : noise image from the VARIANCE extension, same grid
  bad[i]     : flagged pixels (FLAGS extension), same grid
  meta       : target name/position, when and how it was made, and per-frame facts
               (exposure time, wavelength at the centre, detector, archive URL, ...)
Bundles are saved as one .npz + one .json in data_cache/ so the app works without IRSA.
"""
import io
import json
from pathlib import Path

import numpy as np

from . import spherex_io
from .constants import SPHEREX_PIXEL_ARCSEC

GRID_PIX_ARCSEC = SPHEREX_PIXEL_ARCSEC / 2.0       # resample at half a native pixel


def select_frames(df, n_frames=12, detectors=None, t_min=None, t_max=None):
    """Pick up to n_frames exposures, evenly spread in time, optionally from given detectors
    and a time window (MJD). Returns a DataFrame sorted by time."""
    sel = df
    if detectors:
        sel = sel[sel["detector"].isin(list(detectors))]
    if t_min is not None:
        sel = sel[sel["t_mid"] >= t_min]
    if t_max is not None:
        sel = sel[sel["t_mid"] <= t_max]
    sel = sel.sort_values("t_mid").reset_index(drop=True)
    if len(sel) > n_frames:
        idx = np.unique(np.round(np.linspace(0, len(sel) - 1, n_frames)).astype(int))
        sel = sel.iloc[idx].reset_index(drop=True)
    return sel


def build_bundle(name, ra_deg, dec_deg, rows, half_size_arcsec=150.0, progress=None, extra=None, reader=None,
                 workers=1):
    """Download one cutout per row of `rows` (LIVE from IRSA) and align them north-up on
    (ra, dec). Frames that fail to download are skipped and listed in meta['failed']."""
    frames, sigmas, bads, info, failed = [], [], [], [], []
    size_arcmin = 2 * half_size_arcsec * 1.5 / 60.0        # margin: cutouts are rotated
    from concurrent.futures import ThreadPoolExecutor
    records = rows.to_dict("records")

    def get(row):
        try:
            if reader is not None:
                return reader(row), None
            hdul = spherex_io.fetch_cutout(row["access_url"], ra_deg, dec_deg, size_arcmin)
            return spherex_io.read_cutout(hdul), None
        except Exception as exc:
            return None, f"{type(exc).__name__}: {str(exc)[:120]}"
    with ThreadPoolExecutor(workers) as ex:
        results = list(ex.map(get, records))
    for k, (row, (d, err)) in enumerate(zip(records, results)):
        if progress:
            progress(k, len(records), row)
        if d is None:
            failed.append({"obs_id": row.get("obs_id"), "error": err})
            continue
        args = (d["wcs"], ra_deg, dec_deg, half_size_arcsec, GRID_PIX_ARCSEC)
        frames.append(spherex_io.north_up(d["image"], *args).astype(np.float32))
        sigmas.append(spherex_io.north_up(d["sigma"], *args).astype(np.float32))
        bads.append(np.nan_to_num(spherex_io.north_up(d["bad"].astype(float), *args, order=0), nan=1.0) > 0.5)
        info.append({
            "obs_id": row.get("obs_id"), "collection": row.get("obs_collection"),
            "time_utc": row.get("time_utc"), "t_mid": float(row["t_mid"]),
            "detector": int(d["detector"]) if d["detector"] is not None else None,
            "lam_um": d["lam_um"], "dlam_um": d["dlam_um"],
            "access_url": row["access_url"],
            "cutout_url": spherex_io.cutout_url(row["access_url"], ra_deg, dec_deg, size_arcmin),
            "frac_flagged": float(bads[-1].mean()),
        })
    if not frames:
        raise spherex_io.DataUnavailable(f"no cutout could be downloaded for '{name}': {failed[:2]}")
    meta = {"name": name, "ra_deg": ra_deg, "dec_deg": dec_deg, "built_utc": spherex_io._now_utc(),
            "half_size_arcsec": half_size_arcsec, "grid_pix_arcsec": GRID_PIX_ARCSEC,
            "unit": "micro-Jy per native 6.2 arcsec pixel; ZODI model and median subtracted",
            "source": "SPHEREx spectral image cutouts, NASA/IPAC Infrared Science Archive (IRSA)",
            "frames": info, "failed": failed, "mode": "live"}
    if extra:
        meta.update(extra)
    return {"frames": np.stack(frames), "sigma": np.stack(sigmas), "bad": np.stack(bads), "meta": meta}


def _paths(name):
    slug = spherex_io._slug(name)
    return spherex_io.CACHE_DIR / f"bundle_{slug}.npz", spherex_io.CACHE_DIR / f"bundle_{slug}.json"


def save_bundle(bundle):
    npz, js = _paths(bundle["meta"]["name"])
    npz.parent.mkdir(exist_ok=True)
    np.savez_compressed(npz, frames=bundle["frames"], bad=bundle["bad"])      # sigma is not stored: keeps the cache small
    js.write_text(json.dumps(bundle["meta"], indent=1))
    return npz


def load_bundle(name):
    npz, js = _paths(name)
    if not npz.exists():
        raise spherex_io.DataUnavailable(f"no cached image bundle for '{name}'")
    with np.load(npz) as z:
        out = {"frames": z["frames"], "bad": z["bad"]}
        out["sigma"] = z["sigma"] if "sigma" in z.files else np.zeros(z["frames"].shape, dtype=np.float32)
    out["meta"] = json.loads(js.read_text())
    out["meta"]["mode"] = "cached"
    return out


def offset_to_pixel(bundle, ra_deg, dec_deg):
    """(column, row) on the north-up grid for a sky position (east is to the left)."""
    m = bundle["meta"]
    n = bundle["frames"].shape[-1]
    dra = ((ra_deg - m["ra_deg"] + 180) % 360 - 180) * 3600 * np.cos(np.radians(m["dec_deg"]))
    ddec = (dec_deg - m["dec_deg"]) * 3600
    return n // 2 - dra / m["grid_pix_arcsec"], n // 2 + ddec / m["grid_pix_arcsec"]


def pixel_to_offset_arcsec(bundle, col, row):
    """(east, north) offset in arcsec from the bundle centre for a grid position."""
    n = bundle["frames"].shape[-1]
    g = bundle["meta"]["grid_pix_arcsec"]
    return -(col - n // 2) * g, (row - n // 2) * g


def normalise(frame, bad=None, lo=30.0, hi=99.7):
    """Scale one frame to 0..1 between its lo and hi percentiles (robust, per frame), so that
    frames taken at different wavelengths can be blinked."""
    f = np.array(frame, dtype=float)
    good = np.isfinite(f) if bad is None else (np.isfinite(f) & ~bad)
    if good.sum() < 10:
        return np.zeros_like(f)
    a, b = np.percentile(f[good], [lo, hi])
    return np.clip((np.nan_to_num(f, nan=a) - a) / max(b - a, 1e-9), 0, 1)


def frame_label(info):
    lam = f"{info['lam_um']:.2f} µm" if info.get("lam_um") == info.get("lam_um") else "λ n/a"
    return f"{str(info['time_utc'])[:16]} UTC · {lam} · det {info.get('detector')}"


def make_gif(bundle, indices=None, scale=4, duration_ms=600, markers=None):
    """Animated GIF (bytes) of the frames: the 'little video' asked for by the challenge.
    markers: optional list (one per frame) of (col, row) grid positions to circle."""
    from PIL import Image, ImageDraw
    idx = list(range(len(bundle["frames"]))) if indices is None else list(indices)
    images = []
    for k in idx:
        g = (normalise(bundle["frames"][k], bundle["bad"][k]) ** 0.5 * 255).astype(np.uint8)[::-1]   # north up
        im = Image.fromarray(g, mode="L").resize((g.shape[1] * scale, g.shape[0] * scale), Image.NEAREST).convert("RGB")
        draw = ImageDraw.Draw(im)
        draw.rectangle([0, 0, im.width, 16], fill=(0, 0, 0))
        draw.text((4, 2), f"{k + 1}/{len(bundle['frames'])}  " + frame_label(bundle["meta"]["frames"][k]).replace("µ", "u"),
                  fill=(255, 255, 255))
        if markers is not None and markers[k] is not None:
            c, r = markers[k]
            x, y = c * scale, (g.shape[0] - 1 - r) * scale
            draw.ellipse([x - 14, y - 14, x + 14, y + 14], outline=(235, 104, 52), width=2)
        images.append(im)
    buf = io.BytesIO()
    images[0].save(buf, format="GIF", save_all=True, append_images=images[1:], duration=duration_ms, loop=0)
    return buf.getvalue()


def difference(bundle, i, j):
    """Normalised frame j minus normalised frame i. Fixed stars cancel (only roughly when the
    two wavelengths differ); something that moved shows as a dark/bright pair."""
    a = normalise(bundle["frames"][i], bundle["bad"][i])
    b = normalise(bundle["frames"][j], bundle["bad"][j])
    d = b - a
    d[bundle["bad"][i] | bundle["bad"][j]] = np.nan
    return d


def find_peaks(image, sigma_image=None, threshold=5.0, max_peaks=20, border=3):
    """Very simple detector: local maxima brighter than threshold x robust noise.
    Returns a list of (col, row, signal-to-noise), brightest first."""
    img = np.nan_to_num(np.array(image, dtype=float))
    med = np.median(img)
    noise = 1.4826 * np.median(np.abs(img - med)) or 1e-9
    snr = (img - med) / noise
    ny, nx = img.shape
    core = snr[1:-1, 1:-1]
    is_max = np.ones_like(core, dtype=bool)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dx or dy:
                is_max &= core >= snr[1 + dy:ny - 1 + dy, 1 + dx:nx - 1 + dx]
    rows, cols = np.where(is_max & (core > threshold))
    rows, cols = rows + 1, cols + 1
    keep = (rows >= border) & (rows < ny - border) & (cols >= border) & (cols < nx - border)
    rows, cols = rows[keep], cols[keep]
    order = np.argsort(-snr[rows, cols])[:max_peaks]
    return [(int(cols[k]), int(rows[k]), float(snr[rows[k], cols[k]])) for k in order]


def photometry(bundle, k, col, row, radius_native_pix=2.0, sky_in=4.0, sky_out=7.0):
    """Aperture photometry on frame k of a bundle at grid position (col, row).

    Sum inside radius_native_pix (in native 6.2 arcsec pixels) minus the median of a sky
    ring. Frame values are micro-Jy per NATIVE pixel sampled on a finer grid, so the sum is
    multiplied by (grid pixel / native pixel)^2. No aperture correction is applied."""
    g = bundle["meta"]["grid_pix_arcsec"] / SPHEREX_PIXEL_ARCSEC          # 0.5
    img, sig, bad = bundle["frames"][k], bundle["sigma"][k], bundle["bad"][k]
    yy, xx = np.indices(img.shape)
    r = np.hypot(xx - col, yy - row) * g                                  # radius in native pixels
    ap = (r <= radius_native_pix) & np.isfinite(img)
    ring = (r >= sky_in) & (r <= sky_out) & np.isfinite(img) & ~bad
    if ap.sum() < 3 or ring.sum() < 10:
        return None
    sky = float(np.median(img[ring]))
    flux = float(np.sum(img[ap] - sky) * g ** 2)
    if np.any(sig > 0):
        err = float(np.sqrt(np.sum(sig[ap] ** 2) * g ** 2))             # neighbouring grid pixels are correlated: rough
    else:                                                               # cached bundles: noise measured in the sky ring
        err = float(1.4826 * np.median(np.abs(img[ring] - sky)) * np.sqrt(ap.sum()) * g ** 2)
    mag = float(-2.5 * np.log10(flux * 1e-6 / 3631.0)) if flux > 0 else None
    return {"flux_ujy": flux, "err_ujy": err, "snr": flux / err if err > 0 else None, "mag_ab": mag,
            "n_flagged": int(bad[ap].sum())}


def centroid(bundle, k, col, row, box=4):
    """Brightness-weighted centre near (col, row) on frame k (search box of +-box grid pixels)."""
    return spherex_io.centroid(bundle["frames"][k], col, row, box=box)


def catalogue_offsets(bundle):
    """Where the catalogue says the target is in each frame: (east, north) offsets in arcsec from the
    field centre, from the SIMBAD position and proper motion moved to the exposure date.
    None when the bundle has no catalogue entry (custom position)."""
    sim = bundle["meta"].get("simbad")
    if not sim:
        return None
    m = bundle["meta"]
    out = []
    for f in bundle["meta"]["frames"]:
        dy = (f["t_mid"] - 51544.5) / 365.25                          # years since J2000
        ra = sim["ra"] + sim["pmra"] * dy / 3.6e6 / np.cos(np.radians(sim["dec"]))
        dec = sim["dec"] + sim["pmdec"] * dy / 3.6e6
        out.append((float(((ra - m["ra_deg"] + 180) % 360 - 180) * 3600 * np.cos(np.radians(m["dec_deg"]))),
                    float((dec - m["dec_deg"]) * 3600)))
    return out


def expected_motion(bundle, i, j):
    """Catalogue proper motion between frames i and j: (arcsec, pixels, days). None without catalogue."""
    sim = bundle["meta"].get("simbad")
    if not sim:
        return None
    days = abs(bundle["meta"]["frames"][j]["t_mid"] - bundle["meta"]["frames"][i]["t_mid"])
    arc = float(np.hypot(sim["pmra"], sim["pmdec"]) / 1000.0 * days / 365.25)
    return arc, arc / SPHEREX_PIXEL_ARCSEC, float(days)
