# Text for the Space Apps submission form

Copy the fields in this order. The character count is under each block. The images to upload are listed at the end.

## 1. Project name

```text
THREX: SPHEREx time-lapse viewer and Planet X toolkit
```

53 characters.

## 2. High-level summary

```text
THREX is a web tool where you type a sky position and watch real SPHEREx images of it as a time-lapse to spot what moves, checked on a known asteroid (6 Hebe, 766 arcsec measured vs 765 predicted by JPL). An extra toolkit predicts how a hypothetical Planet X would move and whether SPHEREx could see it.
```

303 characters.

## 3. Demo link

```text
https://github.com/kabilGa/threx-spaceapps/blob/main/slides/threx_deck.pdf
```

74 characters.

## 4. Link to final project

```text
https://github.com/kabilGa/threx-spaceapps
```

42 characters.

## 5. Detailed description

```text
THREX answers the challenge in two parts.

1) The viewer (front page of the Streamlit app). The user types RA/Dec or picks one of the five challenge targets: Barnard's Star, Epsilon Indi B, WISE 0855-0714, RW Aurigae and R Draconis. The app lists the real SPHEREx exposures that covered the position (IRSA SIA2 search), aligns small cutouts of 5 arcmin on a north-up grid and plays them as a time-lapse with a frame slider. Every frame shows its UTC date, wavelength and exposure ID; flagged pixels (FLAGS extension) can be shown; the zodiacal-light model is removed; two frames can be subtracted; the SIMBAD catalogue position is marked; the brightness stretch can be per frame or fixed. The app warns that every exposure samples a different wavelength, so brightness changes between frames are colour, not variability.

Cached data (so the app works offline): evenly spaced subsets of the exposures IRSA lists, namely Barnard's Star 68 of 358, Epsilon Indi B 69 of 458, WISE 0855-0714 69 of 258, RW Aurigae 66 of 335 and R Draconis 61 of 539 (2025-04 to 2026-10). Custom RA/Dec queries work only when IRSA's search answers.

Proof on a moving object: the known asteroid 6 Hebe (2025-12-01). JPL Horizons predicts its position in each frame; two groups of frames 13.2 hours apart show it at the predicted places: 766 arcsec measured vs 765 arcsec predicted. JPL sb_ident lists 6 Hebe in the field. Hebe saturates the detector, so only its positions are meaningful. It was already known: nothing was discovered.

What we do not see: the catalogue proper motions of the targets are about 0 to 1.8 pixels over the cached frames (values from SIMBAD, not measured by us), too small to see by eye. R Draconis is saturated in most short-wavelength frames and its variability was not analysed. WISE 0855-0714 is extremely faint below about 1.6 micron (Hubble photometry in Vega magnitudes: 26.31 +/- 0.10 at 1.153 micron and 23.22 +/- 0.35 at 1.575 micron; Luhman and Esplin 2016, arXiv:1605.08620), far below SPHEREx's ~19.5 AB single-band depth, so a source measured near its position at those wavelengths is most likely a background blend and we claim no detection. Photometry is rough and not independently calibrated.

2) Extra pages: what would a hypothetical Planet X look like? From a Kepler orbit (Brown and Batygin 2021 best fit, labelled hypothetical) and the motion of Earth we compute that a planet at 400 AU appears to move about 8 arcsec/day at opposition and shifts by about 1000 arcsec (about 160 SPHEREx pixels) in six months, mostly Earth's parallax, which also gives its distance. From reflected sunlight versus the published SPHEREx depth we find it could be seen only marginally (AB 19.4 to 21.5 depending on band for 3 Earth radii and albedo 0.4) and not at all with a methane-dark near-infrared albedo. Assumptions are labelled as such: one visit is 0.75 mag shallower than the four-pass depth, the albedo 0.03 scenario, and the demonstration values of the orbit angles Omega and omega. Our orbit code agrees with JPL Horizons to 0.4 arcsec for Sedna and Eris. Planet X is hypothetical and no detection is claimed.

Everything is labelled measured, computed or assumed. The code has automated tests. To run it: pip install -r requirements.txt, then streamlit run app.py.
```

3287 characters.

## 6. NASA data used

```text
- SPHEREx Level 2 spectral images, cutouts and exposure metadata (DR1 and QR3), NASA/IPAC Infrared Science Archive (IRSA): SIA2 search and the public IRSA cloud copy of the SPHEREx files.
- JPL Horizons API: ephemerides of 6 Hebe, orbital elements and positions of Sedna and Eris.
- JPL Small-Body Identification (sb_ident) and Small-Body Database (SBDB) APIs: independent check of the asteroid and its properties.
- Published SPHEREx mission numbers (depth, bands, pixel scale): Bock et al. 2025.
```

497 characters.

## 7. Other data and tools

```text
SIMBAD astronomical database (CDS): positions and proper motions of the five targets. Sun's AB magnitudes: Willmer 2018 table. Planet Nine orbit: Brown and Batygin 2021. Atmosphere and albedo context: Fortney et al. 2016. WISE 0855 near-infrared photometry: Luhman and Esplin 2016. Software: Python, NumPy, pandas, Matplotlib, Astropy (including its built-in ephemeris for Earth and the Sun), astroquery, Pillow, Streamlit, pytest, Playwright (screenshots only). Claude and Claude Code (Anthropic) as described in the AI field.
```

527 characters.

## 8. References

```text
Bock, J. J. et al. 2025, The SPHEREx Satellite Mission, arXiv:2511.02985.
Brown, M. E. and Batygin, K. 2021, The orbit of Planet Nine, AJ 162, 219, arXiv:2108.09868.
Fortney, J. J. et al. 2016, The Hunt for Planet Nine: Atmosphere, Spectra, Evolution, and Detectability, ApJL 824, L25, arXiv:1604.07424.
Luhman, K. L. and Esplin, T. L. 2016, The Spectral Energy Distribution of the Coldest Known Brown Dwarf, arXiv:1605.08620.
Willmer, C. N. A. 2018, The Absolute Magnitude of the Sun in Several Filters, ApJS 236, 47.
SPHEREx Archive at IRSA, user guide and tutorials: https://caltech-ipac.github.io/spherex-archive-documentation/ and https://caltech-ipac.github.io/irsa-tutorials/spherex-intro
JPL Solar System Dynamics APIs: https://ssd-api.jpl.nasa.gov/
SIMBAD: https://simbad.cds.unistra.fr/
```

796 characters.

## 9. Use of Artificial Intelligence

```text
Yes, AI tools were used. The team member (Kabil Gannouni) chose the challenge and approved the project concept, which was developed with Claude (Anthropic) in a chat. Claude Code (Anthropic's coding agent) wrote the Python code, the automated tests, the Streamlit app and the first drafts of the documentation and slides, and it queried and processed the public NASA/IPAC SPHEREx and JPL data described in this submission. The team member directed the work and is responsible for all content. All data shown are real SPHEREx and JPL data; synthetic images exist only inside automated tests and are labelled as such.
```

615 characters.

## Images to upload (in this order)

All four are real browser screenshots of the running app (Google Chrome, 1400 px wide), each under 10 MB.

1. `screenshots/01_viewer_barnards_star.png`: the viewer on Barnard's Star
2. `screenshots/02_viewer_asteroid_6_hebe.png`: the viewer on the 6 Hebe asteroid proof
3. `screenshots/03_planet_x_signature.png`: extra page, Planet X signature
4. `screenshots/04_planet_x_detectability.png`: extra page, can SPHEREx see it?
