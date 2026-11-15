# data_cache

Small copies of REAL results, so the app still works when IRSA or JPL cannot be reached.
The app always says whether it shows live data or data "cached on <date>".
Nothing in this folder is synthetic.

| File | What it is | Made by |
|---|---|---|
| `exposures_<target>.csv` + `.json` | List of SPHEREx exposures covering one sky position (one row per exposure: time in MJD and UTC, wavelength range, detector, archive URL). Source: IRSA SIA2, collections `spherex_dr1` and `spherex_qr3`. The `.json` holds the position and the query date. | `python scripts/build_cache.py schedule` |
| `bundle_<target>.npz` + `.json` | A short sequence of SPHEREx image cutouts of one sky patch, resampled north-up on a 3.1 arcsec grid: `frames` (micro-Jy per native pixel, ZODI model and median subtracted), `sigma` (from the VARIANCE extension), `bad` (flagged pixels from the FLAGS extension). The `.json` lists, for every frame, the exposure time, wavelength, detector and the exact IRSA URL it came from. | `python scripts/build_cache.py blink ...` |
| `bundle_asteroid_proof_*.npz` + `.json` | Same, centred on the track of a known asteroid, plus the JPL Horizons predictions, our measured centres and fluxes, and the JPL sb_ident answer. | `python scripts/build_asteroid_proof.py search` then `build` |
| `_work/` | Intermediate tables of the asteroid search (candidate exposures with Horizons positions). | `build_asteroid_proof.py search` |

The Horizons data used by the tests is in `tests/data/horizons_validation.json`
(`python scripts/fetch_horizons_validation.py`).

To rebuild everything from scratch, delete the files and run the commands above
(needs a network connection; IRSA must be up).
