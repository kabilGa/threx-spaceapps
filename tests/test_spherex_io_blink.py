"""Image handling (alignment, masks, photometry, peak finding) on a tiny SYNTHETIC cutout.

The synthetic image exists only inside this test (never saved, never shown in the app): it
has the same extensions as a SPHEREx cutout and a rotated WCS, so we can check that a star
put at a known sky position ends up where it should.
"""
import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

from threx import spherex_io, blink

RA0, DEC0 = 150.0, 2.0
STAR = (150.004, 2.003)          # sky position of the synthetic star
FLUX_MJYSR_PEAK = 5.0


def synthetic_cutout(n=81, rot_deg=37.0):
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [RA0, DEC0]
    w.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    s, c = np.sin(np.radians(rot_deg)), np.cos(np.radians(rot_deg))
    scale = 6.2 / 3600
    w.wcs.cd = np.array([[-scale * c, scale * s], [scale * s, scale * c]])
    x, y = w.all_world2pix(*STAR, 0)
    yy, xx = np.indices((n, n))
    star = FLUX_MJYSR_PEAK * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * 1.3 ** 2))
    zodi = np.full((n, n), 0.4)
    image = star + zodi + 0.02                     # zodiacal light is NOT subtracted in real files
    flags = np.zeros((n, n), dtype=np.int32)
    flags[5, 7] = 1 << 10                          # HOT
    flags[20, 20] = 1 << 21                        # SOURCE only: not a bad pixel
    hdr = w.to_header()
    hdr["DETECTOR"] = 2
    fhdr = fits.Header()
    fhdr["MSKN0010"] = "HOT"
    fhdr["MSKN0021"] = "SOURCE"
    hdul = fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(image.astype(np.float32), header=hdr, name="IMAGE"),
        fits.ImageHDU(flags, header=fhdr, name="FLAGS"),
        fits.ImageHDU(np.full((n, n), 1e-4, dtype=np.float32), name="VARIANCE"),
        fits.ImageHDU(zodi.astype(np.float32), name="ZODI"),
    ])
    return hdul, (x, y)


def test_read_cutout_subtracts_background_and_masks():
    hdul, _ = synthetic_cutout()
    d = spherex_io.read_cutout(hdul)
    assert d["bad"][5, 7] and not d["bad"][20, 20]
    assert d["bad"].sum() == 1
    assert abs(np.median(d["image"])) < 1e-6 * spherex_io.UJY_PER_MJYSR_PIXEL + 1e-3   # background removed
    assert d["detector"] == 2
    assert spherex_io.UJY_PER_MJYSR_PIXEL == pytest.approx(903.5, rel=1e-3)


def test_north_up_puts_star_at_the_right_offset():
    hdul, _ = synthetic_cutout()
    d = spherex_io.read_cutout(hdul)
    img = spherex_io.north_up(d["image"], d["wcs"], RA0, DEC0, half_size_arcsec=60, pix_arcsec=3.1)
    n = img.shape[0]
    row, col = np.unravel_index(np.nanargmax(img), img.shape)
    east = -(col - n // 2) * 3.1
    north = (row - n // 2) * 3.1
    assert east == pytest.approx((STAR[0] - RA0) * 3600 * np.cos(np.radians(DEC0)), abs=3.1)
    assert north == pytest.approx((STAR[1] - DEC0) * 3600, abs=3.1)


def test_aperture_photometry_recovers_flux():
    hdul, _ = synthetic_cutout()
    d = spherex_io.read_cutout(hdul)
    phot = spherex_io.aperture_photometry(d, *STAR, radius_pix=4.0, sky_in=7.0, sky_out=10.0)
    total = FLUX_MJYSR_PEAK * 2 * np.pi * 1.3 ** 2 * spherex_io.UJY_PER_MJYSR_PIXEL     # Gaussian integral
    assert phot["flux_ujy"] == pytest.approx(total, rel=0.03)
    assert phot["snr"] > 50 and phot["n_flagged"] == 0
    assert phot["mag_ab"] == pytest.approx(-2.5 * np.log10(total * 1e-6 / 3631), abs=0.05)


def test_cutout_url_and_visits():
    url = spherex_io.cutout_url("https://irsa.ipac.caltech.edu/ibe/data/spherex/x.fits", 10.5, -3.25, 6.0)
    assert url.endswith("x.fits?center=10.500000,-3.250000deg&size=6.000arcmin")
    import pandas as pd
    df = pd.DataFrame({"t_mid": [60800.0, 60801.5, 60803.0, 60985.0, 60986.0]})
    visits = spherex_io.find_visits(df)
    assert [v["n"] for v in visits] == [3, 2]
    assert visits[0]["duration_days"] == pytest.approx(3.0)


def _bundle():
    hdul, _ = synthetic_cutout()
    d = spherex_io.read_cutout(hdul)
    args = (d["wcs"], RA0, DEC0, 60.0, blink.GRID_PIX_ARCSEC)
    f = spherex_io.north_up(d["image"], *args).astype(np.float32)
    bad = np.zeros(f.shape, dtype=bool)
    info = {"time_utc": "2026-01-01 00:00:00", "t_mid": 61041.0, "detector": 2, "lam_um": 1.2, "dlam_um": 0.03}
    meta = {"name": "SYNTHETIC test", "ra_deg": RA0, "dec_deg": DEC0, "half_size_arcsec": 60.0,
            "grid_pix_arcsec": blink.GRID_PIX_ARCSEC, "frames": [info, dict(info)]}
    return {"frames": np.stack([f, np.roll(f, 4, axis=1)]), "sigma": np.ones((2,) + f.shape, dtype=np.float32),
            "bad": np.stack([bad, bad]), "meta": meta}


def test_blink_helpers():
    b = _bundle()
    col, row = blink.offset_to_pixel(b, *STAR)
    peaks = blink.find_peaks(b["frames"][0])
    assert peaks and abs(peaks[0][0] - col) <= 1.5 and abs(peaks[0][1] - row) <= 1.5
    east, north = blink.pixel_to_offset_arcsec(b, col, row)
    assert east == pytest.approx((STAR[0] - RA0) * 3600 * np.cos(np.radians(DEC0)), abs=0.01)
    assert north == pytest.approx((STAR[1] - DEC0) * 3600, abs=0.01)
    # the star was shifted by 4 grid pixels between the frames: the difference image shows a pair
    d = blink.difference(b, 0, 1)
    assert np.nanmin(d) < -0.5 and np.nanmax(d) > 0.5
    gif = blink.make_gif(b)
    assert gif[:6] in (b"GIF87a", b"GIF89a")


def test_select_frames():
    import pandas as pd
    df = pd.DataFrame({"t_mid": np.arange(100.0), "detector": np.tile([1.0, 2.0], 50)})
    sel = blink.select_frames(df, n_frames=10, detectors=[2])
    assert len(sel) == 10 and set(sel["detector"]) == {2.0}
    assert sel["t_mid"].is_monotonic_increasing
