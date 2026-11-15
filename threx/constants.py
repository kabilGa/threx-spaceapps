"""Physical constants and SPHEREx numbers used everywhere (each with its source)."""

AU_KM = 149_597_870.7             # km, IAU 2012 definition of the astronomical unit
R_EARTH_KM = 6378.1               # km, IAU 2015 nominal equatorial Earth radius
C_AU_PER_DAY = 173.144632674      # speed of light in AU/day
GAUSS_K = 0.01720209895           # Gaussian gravitational constant, AU^1.5/day (Sun only)
OBLIQUITY_DEG = 84381.448 / 3600  # J2000 obliquity (IAU 1976), the one JPL Horizons uses
ARCSEC_PER_RAD = 206264.806247

# --- SPHEREx (Bock et al. 2025, arXiv:2511.02985, abstract) -------------------
SPHEREX_PIXEL_ARCSEC = 6.2        # arcsec per pixel
SPHEREX_N_CHANNELS = 102          # spectral channels, 0.75-5.0 micron
SPHEREX_N_SURVEYS = 4             # all-sky passes in the 2-year mission

# Detector bands (IRSA tutorial "Introduction to SPHEREx Spectral Images").
# depth = published 5-sigma point-source depth range (bright end, faint end) in AB mag,
# valid per spectral channel AFTER the four all-sky passes (Bock et al. abstract).
SPHEREX_BANDS = [
    {"band": 1, "lam_min": 0.75, "lam_max": 1.09, "depth": (19.5, 19.9)},
    {"band": 2, "lam_min": 1.10, "lam_max": 1.62, "depth": (19.5, 19.9)},
    {"band": 3, "lam_min": 1.63, "lam_max": 2.41, "depth": (19.5, 19.9)},
    {"band": 4, "lam_min": 2.42, "lam_max": 3.82, "depth": (19.5, 19.9)},
    {"band": 5, "lam_min": 3.83, "lam_max": 4.41, "depth": (17.8, 18.8)},
    {"band": 6, "lam_min": 4.42, "lam_max": 5.00, "depth": (17.8, 18.8)},
]
for _b in SPHEREX_BANDS:
    _b["lam_c"] = 0.5 * (_b["lam_min"] + _b["lam_max"])   # band centre, micron

# --- Planet Nine (Brown & Batygin 2021, arXiv:2108.09868, abstract). HYPOTHETICAL. ---
P9_A_AU = 380.0         # semi-major axis (+140/-80)
P9_Q_AU = 300.0         # perihelion distance (+85/-60)
P9_INC_DEG = 16.0       # inclination (+-5)
P9_MASS_EARTH = 6.2     # mass in Earth masses (+2.2/-1.3)
# Not given in the abstract -> free inputs. OUR demonstration defaults, not published values.
DEMO_NODE_DEG = 100.0
DEMO_ARGP_DEG = 150.0
