"""Check that every compiled library THREX needs can be loaded on this machine.

Some Windows machines run an "Application Control" policy that blocks individual DLLs of
Python packages. This script tells you which ones (if any) are blocked.

Run:  python scripts/check_environment.py
"""
import importlib
import sys

MODULES = [
    "numpy", "pandas", "matplotlib.pyplot", "PIL.Image", "streamlit", "requests",
    "astropy.time", "astropy.coordinates", "astropy.io.fits",
    "astropy.io.fits.hdu.compressed._codecs", "astropy.io.fits._tiled_compression._compression",
    "astropy.wcs", "astropy.io.votable", "astropy.utils.xml._iterparser", "astropy.table",
    "pyvo", "astroquery.ipac.irsa",
]
print("Python", sys.version.split()[0])
blocked = 0
for name in MODULES:
    try:
        importlib.import_module(name)
        print(f"  ok       {name}")
    except Exception as exc:
        msg = str(exc).splitlines()[0][:110]
        optional = "No module named" in msg
        blocked += 0 if optional else 1
        print(f"  {'absent ' if optional else 'BLOCKED'}  {name}: {msg}")
print("All required libraries load." if blocked == 0 else f"{blocked} module(s) could not be loaded.")
