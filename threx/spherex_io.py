"""Real SPHEREx data from the NASA/IPAC Infrared Science Archive (IRSA).

Docs followed:
  https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access
  https://caltech-ipac.github.io/irsa-tutorials/spherex-intro      (FITS handling)
  https://irsa.ipac.caltech.edu/ibe/cutouts.html                    (cutout service)

Two things are fetched:
  1. the list of exposures covering a sky position (SIA2 query; one row per exposure),
  2. small cutouts of those exposures around the position (cutout service, on-premises).
Everything fetched can be saved in data_cache/ so the app still works when IRSA is down.
Nothing here ever invents data: if IRSA fails and there is no cache, an error is raised.
"""
import http.client
import io
import json
import time
import urllib.error
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.utils.data import conf
from astropy.wcs import WCS

from .constants import SPHEREX_PIXEL_ARCSEC

conf.remote_timeout = 120          # SPHEREx files can be slow to read (IRSA tutorial)

CACHE_DIR = Path(__file__).resolve().parents[1] / "data_cache"
COLLECTIONS = ["spherex_dr1", "spherex_qr3"]     # wide survey, current releases (QR2 is not searchable)
KEEP_COLUMNS = ["obs_id", "obs_collection", "energy_bandpassname", "em_min", "em_max", "t_min", "t_max",
                "t_exptime", "s_ra", "s_dec", "access_url", "access_estsize"]
# 1 MJy/sr on one 6.2 arcsec pixel, in micro-Jansky (pixel solid angle = (6.2/206265)^2 sr)
UJY_PER_MJYSR_PIXEL = 1e12 * (SPHEREX_PIXEL_ARCSEC / 206264.806) ** 2
# flags that make a pixel unusable (names from the FLAGS header; SOURCE just means "known source")
BAD_FLAG_NAMES = ["TRANSIENT", "OVERFLOW", "SUR_ERROR", "NONFUNC", "DICHROIC", "MISSING_DATA", "HOT",
                  "COLD", "NONLINEAR", "PERSIST", "OUTLIER", "STREAK", "BLOOM", "SNOWBALL", "HALO"]
# bit numbers printed in the IRSA tutorial, used only if a header lacks the MSKN keywords
DEFAULT_FLAG_BITS = {"TRANSIENT": 0, "OVERFLOW": 1, "SUR_ERROR": 2, "NONFUNC": 6, "DICHROIC": 7,
                     "MISSING_DATA": 9, "HOT": 10, "COLD": 11, "FULLSAMPLE": 12, "PHANMISS": 14,
                     "NONLINEAR": 15, "PERSIST": 17, "OUTLIER": 19, "CROSSTALK": 20, "SOURCE": 21,
                     "GHOST": 22, "GHOST_FPA": 23, "GHOST_EXT": 24, "STREAK": 25, "BLOOM": 26,
                     "SNOWBALL": 27, "HALO": 28}


class DataUnavailable(RuntimeError):
    """IRSA could not be reached and no cached copy exists."""


def _now_utc():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# ----------------------------------------------------------------------------- exposures
def query_exposures(ra_deg, dec_deg, radius_arcsec=1.0, collections=COLLECTIONS, retries=1):
    """LIVE SIA2 query: every SPHEREx spectral image whose footprint covers the position.
    Returns a pandas DataFrame sorted by time, with wavelengths in micron and UTC strings."""
    from astroquery.ipac.irsa import Irsa
    coord = SkyCoord(ra_deg, dec_deg, unit="deg")
    frames = []
    for coll in collections:                       # one collection per call: works with any astroquery
        for attempt in range(3):                   # the archive answers 502 now and then
            try:
                tab = Irsa.query_sia(pos=(coord, radius_arcsec * u.arcsec), collection=coll)
                break
            except Exception:
                if attempt == 2 or retries == 0:
                    raise
                time.sleep(8 * (attempt + 1))
        if len(tab):
            frames.append(tab[[c for c in KEEP_COLUMNS if c in tab.colnames]].to_pandas())
    if not frames:
        return pd.DataFrame(columns=KEEP_COLUMNS + ["time_utc", "lam_min_um", "lam_max_um"])
    df = pd.concat(frames, ignore_index=True)
    for col in df.columns:                          # bytes -> str
        if df[col].dtype == object:
            df[col] = df[col].apply(lambda v: v.decode() if isinstance(v, bytes) else v)
    df = df.sort_values("t_min").drop_duplicates("access_url", keep="last").reset_index(drop=True)
    return add_derived_columns(df)


def add_derived_columns(df):
    """UTC date strings, wavelengths in micron, detector number."""
    df = df.copy()
    t_mid = 0.5 * (df["t_min"].astype(float) + df["t_max"].astype(float))
    df["t_mid"] = t_mid
    df["time_utc"] = Time(t_mid.to_numpy(), format="mjd", scale="utc").iso if len(df) else []
    df["lam_min_um"] = df["em_min"].astype(float) * 1e6        # SIA gives metres
    df["lam_max_um"] = df["em_max"].astype(float) * 1e6
    # file names end in ..._<obs>D<detector>_spx_...  -> detector 1..6
    df["detector"] = df["access_url"].str.extract(r"D(\d)_spx")[0].astype(float)
    return df


def _slug(name):
    return "".join(ch if ch.isalnum() else "_" for ch in name).strip("_").lower()


def save_exposures(name, df, ra_deg, dec_deg):
    CACHE_DIR.mkdir(exist_ok=True)
    df.to_csv(CACHE_DIR / f"exposures_{_slug(name)}.csv", index=False)
    meta = {"name": name, "ra_deg": ra_deg, "dec_deg": dec_deg, "queried_utc": _now_utc(),
            "rows": int(len(df)), "collections": COLLECTIONS,
            "source": "IRSA SIA2, https://irsa.ipac.caltech.edu/SIA"}
    (CACHE_DIR / f"exposures_{_slug(name)}.json").write_text(json.dumps(meta, indent=1))
    return meta


def load_exposures(name):
    csv = CACHE_DIR / f"exposures_{_slug(name)}.csv"
    if not csv.exists():
        raise DataUnavailable(f"no cached exposure list for '{name}'")
    meta = json.loads((CACHE_DIR / f"exposures_{_slug(name)}.json").read_text())
    return pd.read_csv(csv), meta


def get_exposures(name, ra_deg, dec_deg, live=True):
    """Exposure list for a position: live from IRSA if possible, otherwise the cached copy.
    Returns (DataFrame, meta) where meta['mode'] is 'live' or 'cached'."""
    error = None
    if live:
        try:
            df = query_exposures(ra_deg, dec_deg)
            meta = {"name": name, "ra_deg": ra_deg, "dec_deg": dec_deg, "queried_utc": _now_utc(),
                    "rows": int(len(df)), "mode": "live"}
            return df, meta
        except Exception as exc:                     # network down, IRSA down, ...
            error = f"{type(exc).__name__}: {str(exc)[:200]}"
    try:
        df, meta = load_exposures(name)
    except DataUnavailable:
        raise DataUnavailable(f"IRSA not reachable ({error}) and no cached list for '{name}'")
    meta["mode"] = "cached"
    meta["live_error"] = error
    return df, meta


def find_visits(df, gap_days=30.0):
    """Group exposures into 'visits' (runs of exposures separated by less than gap_days).
    Returns a list of dicts: start/end (MJD), n exposures, duration in days."""
    if len(df) == 0:
        return []
    t = np.sort(df["t_mid"].to_numpy(dtype=float))
    breaks = np.where(np.diff(t) > gap_days)[0]
    starts = np.concatenate([[0], breaks + 1])
    ends = np.concatenate([breaks, [len(t) - 1]])
    return [{"start_mjd": float(t[a]), "end_mjd": float(t[b]), "n": int(b - a + 1),
             "duration_days": float(t[b] - t[a]), "mid_mjd": float(0.5 * (t[a] + t[b]))}
            for a, b in zip(starts, ends)]


# ------------------------------------------------------------------------------- cutouts
def cutout_url(access_url, ra_deg, dec_deg, size_arcmin):
    """URL of a cutout: the image URL plus ?center=RA,Dec&size=... (IRSA cutout service)."""
    return f"{access_url}?center={ra_deg:.6f},{dec_deg:.6f}deg&size={size_arcmin:.3f}arcmin"


def fetch_cutout(access_url, ra_deg, dec_deg, size_arcmin=6.0, max_retries=3, save_to=None):
    """Download one cutout (all extensions) and return it as an in-memory HDUList.
    Retries on transient errors, as in the IRSA tutorial. Optionally saves the FITS file."""
    url = cutout_url(access_url, ra_deg, dec_deg, size_arcmin)
    for attempt in range(max_retries):
        try:
            with fits.open(url, memmap=False, cache=False) as hdul:
                for hdu in hdul:
                    hdu.data                          # force the read while the file is open
                copy = fits.HDUList([h.copy() for h in hdul])
            break
        except (TimeoutError, urllib.error.URLError, http.client.IncompleteRead, OSError):
            if attempt == max_retries - 1:
                raise
            time.sleep(5 * (attempt + 1))
    if save_to is not None:
        Path(save_to).parent.mkdir(parents=True, exist_ok=True)
        copy.writeto(save_to, overwrite=True)
    return copy


def open_cutout(path_or_bytes):
    """Open a cached cutout file fully in memory."""
    if isinstance(path_or_bytes, (bytes, bytearray)):
        path_or_bytes = io.BytesIO(path_or_bytes)
    with fits.open(path_or_bytes, memmap=False) as hdul:
        for hdu in hdul:
            hdu.data
        return fits.HDUList([h.copy() for h in hdul])


def bad_pixel_mask(flags, flags_header):
    """True where a pixel carries one of the BAD_FLAG_NAMES flags."""
    bits = dict(DEFAULT_FLAG_BITS)
    for key in flags_header.get("MSKN*", {}):        # flag names are stored in the header
        bits[str(flags_header[key]).strip()] = int(key[4:])
    mask_value = sum(1 << bits[name] for name in BAD_FLAG_NAMES if name in bits)
    return (np.asarray(flags).astype(np.int64) & mask_value) != 0


def read_cutout(hdul):
    """Turn a SPHEREx cutout into the arrays we use.

    Returns a dict:
      image   : IMAGE minus the ZODI model minus the remaining median, in micro-Jy per pixel
                (the zodiacal light is NOT subtracted in the archive files),
      sigma   : sqrt(VARIANCE) in the same unit,      bad : flagged pixels (bool),
      wcs     : celestial WCS (with distortion),      lam_um / dlam_um : wavelength and
                bandwidth at the cutout centre,       detector, header.
    """
    hdr = hdul["IMAGE"].header
    img = np.asarray(hdul["IMAGE"].data, dtype=float)
    zodi = np.asarray(hdul["ZODI"].data, dtype=float) if "ZODI" in hdul else 0.0
    var = np.asarray(hdul["VARIANCE"].data, dtype=float)
    bad = bad_pixel_mask(hdul["FLAGS"].data, hdul["FLAGS"].header) | ~np.isfinite(img)
    sub = img - zodi
    good = ~bad & np.isfinite(sub)
    sub = sub - (np.median(sub[good]) if good.any() else 0.0)
    wcs = WCS(hdr)
    ny, nx = img.shape
    lam = dlam = np.nan
    try:                                              # wavelength lookup table (WCS "W")
        wwcs = WCS(header=hdr, fobj=hdul, key="W")
        wwcs.sip = None                               # distortion must be off for this WCS
        lam_q, dlam_q = wwcs.pixel_to_world((nx - 1) / 2, (ny - 1) / 2)
        lam, dlam = float(lam_q.to_value("um")), float(dlam_q.to_value("um"))
    except Exception:
        pass
    return {"image": sub * UJY_PER_MJYSR_PIXEL, "sigma": np.sqrt(np.clip(var, 0, None)) * UJY_PER_MJYSR_PIXEL,
            "bad": bad, "wcs": wcs, "lam_um": lam, "dlam_um": dlam,
            "detector": hdr.get("DETECTOR"), "header": hdr, "shape": img.shape}


def north_up(array, wcs, ra_deg, dec_deg, half_size_arcsec=120.0, pix_arcsec=3.1, order=1):
    """Resample onto a north-up, east-left grid centred on (ra, dec).

    SPHEREx images are rotated differently in every exposure, so frames must be put on a
    common grid before they can be compared or blinked. Bilinear interpolation (order=1),
    or nearest neighbour (order=0, used for masks). Pixels outside the cutout become NaN.
    """
    n = int(2 * half_size_arcsec / pix_arcsec) + 1
    off = (np.arange(n) - n // 2) * pix_arcsec                      # arcsec
    xi, eta = np.meshgrid(-off, off)                                # east to the left, north up
    xi_r, eta_r = np.radians(xi / 3600), np.radians(eta / 3600)
    ra0, dec0 = np.radians(ra_deg), np.radians(dec_deg)
    # inverse gnomonic (tangent plane) projection
    den = np.cos(dec0) - eta_r * np.sin(dec0)
    ra = ra0 + np.arctan2(xi_r, den)
    dec = np.arctan2(np.sin(dec0) + eta_r * np.cos(dec0), np.hypot(xi_r, den))
    x, y = wcs.all_world2pix(np.degrees(ra), np.degrees(dec), 0)
    arr = np.asarray(array, dtype=float)
    ny, nx = arr.shape
    if order == 0:
        xi_, yi_ = np.rint(x).astype(int), np.rint(y).astype(int)
        inside = (xi_ >= 0) & (xi_ < nx) & (yi_ >= 0) & (yi_ < ny)
        out = np.full(x.shape, np.nan)
        out[inside] = arr[yi_[inside], xi_[inside]]
        return out
    x0, y0 = np.floor(x).astype(int), np.floor(y).astype(int)
    fx, fy = x - x0, y - y0
    inside = (x0 >= 0) & (x0 < nx - 1) & (y0 >= 0) & (y0 < ny - 1)
    x0c, y0c = np.clip(x0, 0, nx - 2), np.clip(y0, 0, ny - 2)
    out = (arr[y0c, x0c] * (1 - fx) * (1 - fy) + arr[y0c, x0c + 1] * fx * (1 - fy)
           + arr[y0c + 1, x0c] * (1 - fx) * fy + arr[y0c + 1, x0c + 1] * fx * fy)
    out[~inside] = np.nan
    return out


def aperture_photometry(data, ra_deg, dec_deg, radius_pix=2.0, sky_in=4.0, sky_out=7.0):
    """Simple aperture photometry at a sky position on a read_cutout() dict.

    Sum inside radius_pix (native pixels) minus the median of a sky ring. No aperture
    correction is applied, so fluxes are slightly underestimated. Returns flux and error in
    micro-Jy, AB magnitude, signal-to-noise, pixel position, and how many aperture pixels
    are flagged."""
    x, y = data["wcs"].all_world2pix(ra_deg, dec_deg, 0)
    yy, xx = np.indices(data["image"].shape)
    r = np.hypot(xx - x, yy - y)
    ap, ring = r <= radius_pix, (r >= sky_in) & (r <= sky_out) & ~data["bad"]
    if not ap.any() or not ring.any():
        return None
    sky = np.nanmedian(data["image"][ring])
    flux = float(np.nansum(data["image"][ap] - sky))
    err = float(np.sqrt(np.nansum(data["sigma"][ap] ** 2)))
    mag = -2.5 * np.log10(flux * 1e-6 / 3631.0) if flux > 0 else np.nan
    return {"flux_ujy": flux, "err_ujy": err, "snr": flux / err if err > 0 else np.nan,
            "mag_ab": float(mag), "x": float(x), "y": float(y), "n_flagged": int(data["bad"][ap].sum())}


def centroid(image, x0, y0, box=3):
    """Brightness-weighted centre in a small box around (x0, y0) on a north-up image."""
    ny, nx = image.shape
    xa, xb = max(int(x0) - box, 0), min(int(x0) + box + 1, nx)
    ya, yb = max(int(y0) - box, 0), min(int(y0) + box + 1, ny)
    sub = np.nan_to_num(image[ya:yb, xa:xb].copy())
    sub = np.clip(sub - np.median(sub), 0, None)
    if sub.sum() <= 0:
        return float(x0), float(y0)
    yy, xx = np.indices(sub.shape)
    return float(xa + (sub * xx).sum() / sub.sum()), float(ya + (sub * yy).sum() / sub.sum())
