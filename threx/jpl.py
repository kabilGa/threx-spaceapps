"""Small clients for NASA/JPL Solar System Dynamics services (no API key needed).

Horizons    : https://ssd-api.jpl.nasa.gov/doc/horizons.html   (positions, orbital elements)
sb_ident    : https://ssd-api.jpl.nasa.gov/doc/sb_ident.html   (which known bodies are in a field)
SBDB        : https://ssd-api.jpl.nasa.gov/doc/sbdb.html       (facts about one body)
Every function raises on network errors; callers decide whether to fall back to a cache.
"""
import requests

from .orbits import OrbitalElements

HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
SB_IDENT_URL = "https://ssd-api.jpl.nasa.gov/sb_ident.api"
SBDB_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api"
TIMEOUT = 90


def _horizons_table(params):
    """Call Horizons (text format) and return (header names, rows) of the CSV block."""
    base = {"format": "text", "OBJ_DATA": "NO", "MAKE_EPHEM": "YES", "CSV_FORMAT": "YES"}
    r = requests.get(HORIZONS_URL, params={**base, **params}, timeout=TIMEOUT)
    r.raise_for_status()
    lines = r.text.splitlines()
    if "$$SOE" not in lines:
        raise RuntimeError("Horizons returned no ephemeris:\n" + r.text[-800:])
    i0, i1 = lines.index("$$SOE"), lines.index("$$EOE")
    header = [h.strip() for h in lines[i0 - 2].split(",")]
    rows = [[c.strip() for c in ln.split(",")] for ln in lines[i0 + 1:i1]]
    return header, rows


def horizons_observer(command, jd_list, center="500@399", quantities="1,3,9,19,20,23"):
    """Geocentric (default) astrometric ephemeris at the given Julian Dates (UT).

    Returns a list of dicts with: jd, ra, dec [deg, ICRF, astrometric], dra_cosd, ddec
    [arcsec/hour], vmag (when available), r, delta [AU], elong [deg].
    command examples: "90377;" (asteroid number), "DES=2003 VB12;".
    """
    header, rows = _horizons_table({
        "COMMAND": f"'{command}'", "EPHEM_TYPE": "OBSERVER", "CENTER": f"'{center}'",
        "TLIST": ",".join(f"'{jd:.6f}'" for jd in jd_list), "QUANTITIES": f"'{quantities}'",
        "ANG_FORMAT": "DEG", "CAL_FORMAT": "JD", "EXTRA_PREC": "YES",
    })
    names = {"dRA*cosD": "dra_cosd", "d(DEC)/dt": "ddec", "APmag": "vmag", "r": "r",
             "delta": "delta", "S-O-T": "elong", "T-mag": "vmag"}

    def key_for(h):
        # RA/Dec column titles change with the precision option, so match their start
        if h.startswith("Date"):
            return "jd"
        if h.startswith("R.A."):
            return "ra"
        if h.startswith("DEC"):
            return "dec"
        return names.get(h)

    out = []
    for row in rows:
        d = {}
        for h, c in zip(header, row):
            k = key_for(h)
            if k:
                try:
                    d[k] = float(c)
                except ValueError:
                    d[k] = None
        out.append(d)
    return out


def horizons_elements(command, epoch_jd, label=None):
    """Heliocentric osculating elements (ecliptic J2000, AU, deg) at epoch_jd (TDB)."""
    header, rows = _horizons_table({
        "COMMAND": f"'{command}'", "EPHEM_TYPE": "ELEMENTS", "CENTER": "'500@10'",
        "REF_PLANE": "ECLIPTIC", "OUT_UNITS": "AU-D", "TLIST": f"'{epoch_jd:.6f}'",
    })
    v = dict(zip(header, rows[0]))
    return OrbitalElements(a=float(v["A"]), e=float(v["EC"]), inc=float(v["IN"]),
                           node=float(v["OM"]), argp=float(v["W"]), M0=float(v["MA"]),
                           epoch_jd=float(v["JDTDB"]), label=label or command)


def sexagesimal(value_deg, hours=False):
    """Degrees -> 'dd-mm-ss.ss' (or hours 'hh-mm-ss.ss'); minus sign written 'M' (sb_ident format)."""
    sign = "M" if value_deg < 0 else ""
    total = round(abs(value_deg) / (15.0 if hours else 1.0) * 360000)      # hundredths of a second
    a, rest = divmod(total, 360000)
    b, c = divmod(rest, 6000)
    return f"{sign}{a:02d}-{b:02d}-{c / 100:05.2f}"


def sb_ident(obs_time_utc, ra_deg, dec_deg, half_width_deg=0.1, kind="a"):
    """Known small bodies inside a box around (ra, dec) at a UTC time ('YYYY-MM-DD_hh:mm:ss').

    Observer: the centre of the Earth, given as latitude 0, longitude 0, altitude -6378.137 km
    (the service does not accept MPC code 500). Returns the decoded JSON."""
    params = {
        "sb-kind": kind, "lat": "0", "lon": "0", "alt": "-6378.137", "obs-time": obs_time_utc,
        "fov-ra-center": sexagesimal(ra_deg, hours=True), "fov-dec-center": sexagesimal(dec_deg),
        "fov-ra-hwidth": f"{half_width_deg:.4f}", "fov-dec-hwidth": f"{half_width_deg:.4f}",
        "two-pass": "true", "mag-required": "true", "req-elem": "false",
    }
    r = requests.get(SB_IDENT_URL, params=params, timeout=TIMEOUT)
    r.raise_for_status()
    return r.json()


def sbdb(sstr, phys=True):
    """Small-Body Database record for one object (name, orbit class, H magnitude, diameter...)."""
    r = requests.get(SBDB_URL, params={"sstr": sstr, "phys-par": "1" if phys else "0"}, timeout=TIMEOUT)
    r.raise_for_status()
    return r.json()
