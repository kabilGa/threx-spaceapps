# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Find and read SPHEREx exposures on the public IRSA cloud mirror (AWS bucket nasa-irsa-spherex).

Used when IRSA's SIA2/TAP search is down. The bucket has NO search by sky position, so we
build the search ourselves from the file headers:
  1. keep only the planning periods (weeks) in which a target was about 90 deg from the Sun
     (the SPHEREx scan strategy),
  2. read the first ~40 kB of every Nth exposure of detector 2 (HTTP Range request) to get its
     pointing and time,
  3. around the pointings that come close, read the headers of all 6 detectors and test with
     the file's own WCS whether the target falls on the detector.
Pixels are then read with astropy + fsspec (only the needed rows).  Nothing is invented.
Docs: https://irsa.ipac.caltech.edu/cloud_access/#spherex ,
      https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access
"""
import datetime as dt
import re
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import requests
from astropy.io import fits
from astropy.wcs import WCS

BUCKET = "https://nasa-irsa-spherex.s3.us-east-1.amazonaws.com/"
IRSA_DATA = "https://irsa.ipac.caltech.edu/ibe/data/spherex/"       # same paths as the bucket
SESSION = requests.Session()
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
SESSION.mount("https://", HTTPAdapter(max_retries=Retry(total=6, backoff_factor=1.0, status_forcelist=(500, 502, 503, 504)),
                                      pool_connections=32, pool_maxsize=32))


def list_prefix(prefix, delimiter="/"):
    """All sub-prefixes and (key, size) of an S3 prefix (follows continuation tokens)."""
    prefixes, keys, token = [], [], None
    while True:
        params = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if delimiter:
            params["delimiter"] = delimiter
        if token:
            params["continuation-token"] = token
        xml = SESSION.get(BUCKET, params=params, timeout=60).text
        prefixes += re.findall(r"<CommonPrefixes><Prefix>([^<]+)</Prefix>", xml)
        keys += list(zip(re.findall(r"<Key>([^<]+)</Key>", xml), map(int, re.findall(r"<Size>(\d+)</Size>", xml))))
        nxt = re.findall(r"<NextContinuationToken>([^<]+)<", xml)
        if not nxt:
            return prefixes, keys
        token = nxt[0]


def period_week_start(name):
    """'2026W21_2A' -> date of the Monday of ISO week 21 of 2026."""
    m = re.match(r"(\d{4})W(\d{2})_", name)
    return dt.date.fromisocalendar(int(m.group(1)), int(m.group(2)), 1)


def list_periods(releases=("dr1", "qr3")):
    """[(release, period_name, week_start_date, prefix_with_version)] for every planning period."""
    out = []
    for rel in releases:
        for p in list_prefix(f"{rel}/level2/")[0]:
            name = p.rstrip("/").split("/")[-1]
            if not re.match(r"\d{4}W\d{2}_", name):
                continue
            for v in list_prefix(p)[0]:                      # version directory, e.g. l2b-v27-2026-268/
                out.append((rel, name, period_week_start(name), v))
    return out


def read_header(key, nbytes=40000):
    """IMAGE-extension header of an exposure file from its first bytes (HTTP Range request)."""
    r = SESSION.get(BUCKET + key, headers={"Range": f"bytes=0-{nbytes - 1}"}, timeout=60)
    r.raise_for_status()
    raw = r.content[2880:]                                    # skip the short primary header
    for n in range(0, len(raw), 2880):
        blk = raw[n:n + 2880]
        if any(blk[j:j + 80].startswith(b"END" + b" " * 77) for j in range(0, 2880, 80)):
            return fits.Header.fromstring(raw[:n + 2880].decode("ascii", "replace"))
    raise ValueError("header longer than the bytes read: " + key)


def key_for_detector(key, det):
    """Same exposure on another detector: swap the directory number and the 'D<n>' tag."""
    parts = key.split("/")
    parts[-2] = str(det)
    parts[-1] = re.sub(r"D\d_spx", f"D{det}_spx", parts[-1])
    return "/".join(parts)


def sep_deg(ra1, dec1, ra2, dec2):
    a1, d1, a2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(a1 - a2)
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


def covers(header, ra, dec, margin_pix=25):
    """Pixel (x, y) of (ra, dec) in this exposure and whether it is inside the detector."""
    w = WCS(header)
    x, y = w.all_world2pix(ra, dec, 0)
    n = header["NAXIS1"]
    return bool(margin_pix <= x < n - margin_pix and margin_pix <= y < n - margin_pix), float(x), float(y)


def scan_period(prefix, ra, dec, stride=6, near_deg=8.0, workers=16):
    """Hits [{key, det, mjd, x, y}] for one (release, period, version) prefix.
    ra may also be a function mjd -> (ra, dec) for a MOVING target (dec is then ignored)."""
    pos = ra if callable(ra) else (lambda m: (ra, dec))
    keys = sorted(k for k, _ in list_prefix(prefix + "2/")[1] if k.endswith(".fits"))
    if not keys:
        return []
    idx = list(range(0, len(keys), stride))

    def probe(i):
        try:
            h = read_header(keys[i])
            r_, d_ = pos(float(h["MJD-AVG"]))
            return i, sep_deg(h["CRVAL1"], h["CRVAL2"], r_, d_)
        except Exception:
            return i, 999.0
    with ThreadPoolExecutor(workers) as ex:
        probes = list(ex.map(probe, idx))
    near = sorted({j for i, s in probes if s < near_deg for j in range(max(0, i - stride), min(len(keys), i + stride + 1))})
    cand = [(keys[j], d) for j in near for d in range(1, 7)]

    def test(item):
        key, det = item
        k2 = key_for_detector(key, det)
        try:
            h = read_header(k2)
            r_, d_ = pos(float(h["MJD-AVG"]))
            ok, x, y = covers(h, r_, d_)
            return {"key": k2, "det": det, "mjd": float(h["MJD-AVG"]), "x": x, "y": y} if ok else None
        except Exception:
            return None
    with ThreadPoolExecutor(workers) as ex:
        return [h for h in ex.map(test, cand) if h]


def read_cutout_s3(key, ra, dec, half_size_arcsec=150.0):
    """Read a small piece of one exposure (only the needed pixel rows are transferred).

    Returns the dict of spherex_io.read_cutout() plus 'lam_um'/'dlam_um' at the target position
    and 'pix' (x, y of the target on the full detector)."""
    from . import spherex_io
    url = BUCKET + key
    with fits.open(url, use_fsspec=True, lazy_load_hdus=True,
                   fsspec_kwargs={"block_size": 1 << 18}) as hdul:
        hdr = hdul["IMAGE"].header.copy()
        w_full = WCS(hdr)
        x, y = w_full.all_world2pix(ra, dec, 0)
        lam = dlam = float("nan")
        try:
            ww = WCS(header=hdr, fobj=hdul, key="W")
            ww.sip = None
            l, b = ww.pixel_to_world(float(x), float(y))
            lam, dlam = float(l.to_value("um")), float(b.to_value("um"))
        except Exception:
            pass
        n = hdr["NAXIS1"]
        if not (30 <= x < n - 30 and 30 <= y < n - 30):
            raise ValueError("target too close to the detector edge")
        half = int(half_size_arcsec * 1.5 / spherex_io.SPHEREX_PIXEL_ARCSEC) + 4
        x0, x1 = max(int(x) - half, 0), min(int(x) + half + 1, n)
        y0, y1 = max(int(y) - half, 0), min(int(y) + half + 1, n)
        img = np.asarray(hdul["IMAGE"].section[y0:y1, x0:x1], dtype=float)
        var = np.asarray(hdul["VARIANCE"].section[y0:y1, x0:x1], dtype=float)
        zodi = np.asarray(hdul["ZODI"].section[y0:y1, x0:x1], dtype=float)
        flags = np.asarray(hdul["FLAGS"].section[y0:y1, x0:x1])
        fhdr = hdul["FLAGS"].header.copy()
    chdr = hdr.copy()
    chdr["CRPIX1"] = hdr["CRPIX1"] - x0
    chdr["CRPIX2"] = hdr["CRPIX2"] - y0
    bad = spherex_io.bad_pixel_mask(flags, fhdr) | ~np.isfinite(img)
    sub = img - zodi
    good = ~bad & np.isfinite(sub)
    sub = sub - (np.median(sub[good]) if good.any() else 0.0)
    u = spherex_io.UJY_PER_MJYSR_PIXEL
    return {"image": sub * u, "sigma": np.sqrt(np.clip(var, 0, None)) * u, "bad": bad, "wcs": WCS(chdr),
            "lam_um": lam, "dlam_um": dlam, "detector": hdr.get("DETECTOR"), "header": chdr,
            "shape": img.shape, "pix": (float(x), float(y))}
