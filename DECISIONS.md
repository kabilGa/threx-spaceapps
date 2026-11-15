# DECISIONS

Decisions taken by Claude Code while Kabil was away on Sunday 15 Nov 2026 (start 11:13 Tunis
time). Each entry says what was decided and why. Kabil owns the idea and can overrule any of them.

## D0. The official challenge text asks for something different from our brief (READ THIS FIRST)

I read the official page (https://www.spaceappschallenge.org/2026/challenges/planet-x-and-spherex/,
Details and Resources tabs) at 11:25 through the in-app browser. `curl` could not open it (the
certificate chain was not trusted on this machine; I did not bypass certificate checks).

What the page asks for, in my words: a **public web tool that shows SPHEREx images of a patch
of sky chosen by RA/Dec as a short sequence ("like a little ____ (the word for a short moving picture, as on the page)")** so that moving or
changing objects are easy to spot. It names test targets: Barnard's Star, Epsilon Indi B,
WISE 0855 (RA 08h55m10.831s, Dec -07d14m42.53s), and the variable stars RW Aurigae and
R Draconis. Optional ideas: adding/subtracting images, stacking per detector, spectra on the
fly, coordinates in the URL, showing flags, Gaia/WISE overlays. Listed resources: IRSA
SPHEREx pages, the IRSA SPHEREx tutorial notebooks, WiseView, Gaia DR3.

Our brief (THREX) is a physics model of a hypothetical Planet X plus a real-data proof. It
does not, by itself, deliver the blink viewer the page asks for.

**Decision:** keep the THREX concept and the tier order exactly as briefed (Kabil owns the
idea), and build the "Real SPHEREx data" page as a blink viewer for any RA/Dec, so the tool
also answers the official objective:
- the user types RA/Dec (or picks a preset, including the challenge test targets),
- the app lists the real SPHEREx exposures (schedule, Tier 2 C),
- it shows north-up, aligned cutouts as a frame-by-frame sequence with a difference image
  (Tier 2 D generalised), with wavelength, date and flagged pixels shown for each frame,
- the asteroid proof is one worked example of that viewer,
- RA/Dec can be passed in the URL (`?ra=..&dec=..`).

**Kabil must check:** that the framing in `form_text.md` ("a viewer plus a physics model that
tells you what a Planet X would look like in it") is how he wants to present the project.

## D1. Where the project lives

The session started inside `Sentinel_Starter_Kit` (another project, its own git repository).
The rule was "do not touch files outside this folder", so THREX lives in the sub-folder
`threx-spaceapps/` with **its own git repository** (`git init` inside it). Nothing in the
Sentinel repository was modified; it only sees one new untracked folder. To publish, Kabil
can push `threx-spaceapps/` as it is, or move the folder anywhere first.

## D2. Environment

- Python 3.14 virtual environment in `.venv/` (3.15 is also installed on this machine but is
  too new for some wheels).
- **scipy removed from the dependencies.** On this machine Windows "Application Control"
  blocks one scipy DLL (`_superlu`), so `scipy.optimize` cannot be imported. I did not try to
  bypass the policy. The two places that needed scipy now use 20 lines of NumPy (a bisection
  root finder and a bilinear image resampler). Fewer dependencies is also better for
  Streamlit Cloud.
- Horizons is called through its HTTP API with `requests` (simple and transparent) instead
  of `astroquery.jplhorizons`. IRSA is called through `astroquery.ipac.irsa`, as in the
  official tutorial.

## D3. Science choices (Tier 1)

- **Observer = centre of the Earth** (SPHEREx is ~700 km up). Error at 400 AU is about
  0.02 arcsec, 300 times smaller than a pixel. Stated in the code, app and README.
- **Earth and Sun positions:** astropy built-in ephemeris (no download).
- **Light-time:** included (two iterations), so our positions are "astrometric" like
  Horizons quantity 1. Aberration and light bending are not modelled (they move all sources
  in a field together and do not matter for a shift measured against the stars).
- **Omega, omega, M0 of Planet Nine** are not in the Brown & Batygin (2021) abstract. They
  are free inputs. Defaults: node = 100 deg, argument of perihelion = 150 deg, labelled in
  the app as "demonstration default, not a published value". M0 is set from the distance
  slider.
- **Distance slider beyond the orbit:** the Brown & Batygin orbit only spans 300-460 AU. For
  other distances (200-1000 AU) the app uses a circular orbit at that distance ("user-defined").
- **Sun's brightness:** apparent AB magnitudes at 1 AU from Willmer (2018, ApJS 236, 47),
  online table http://mips.as.arizona.edu/~cnaw/sun.html, read on 2026-11-15, interpolated in
  log wavelength. A 5772 K blackbody is kept as a labelled cross-check (agrees within 0.25 mag).
- **Albedo scenarios:** 0.2, 0.4, 0.6 (brief) plus "methane-dark near-IR" p = 0.03. The 0.03
  is an illustrative order of magnitude (Uranus and Neptune reflect a few percent in methane
  bands); it is labelled "illustrative" everywhere. Not a fitted number.
- **Radius:** 2-3.7 Earth radii is the brief's assumption. Fortney et al. (2016,
  arXiv:1604.07424 abstract) quote 3-6 Earth radii and a visible albedo up to ~0.75; both are
  cited as context.
- **Thermal emission:** blackbody at 10-50 K is more than 30 magnitudes below the reflected
  light at 3-5 micron, so negligible. Caveat written everywhere: Fortney et al. (2016) say real
  atmospheres can be many orders of magnitude above a blackbody at 3-5 micron. We do not
  model that.
- **Depth for a moving source (my assumption, important):** the published SPHEREx depth
  (19.5-19.9 AB and 17.8-18.8 AB) is "after mapping the full sky four times". A planet at
  400 AU moves about 150 pixels between two passes, so the four passes cannot simply be
  averaged. I assume one pass is 0.75 mag shallower (2.5 log10 sqrt(4)) and that adding K
  channels gains 2.5 log10 sqrt(K) at best. Both are labelled "assumed".
- **Heat maps:** six panels (one per SPHEREx detector band) rather than two, because the
  Sun is 2 magnitudes fainter (AB) at 4.7 micron than at 1.2 micron.

## D4. Sanity checks of the brief (section 6.5): all reproduced

| Quantity | Brief | Our code | Comment |
|---|---|---|---|
| V-like AB mag, 400 AU, 3 R_earth, p = 0.4 | ~19.7 | 19.71 | |
| Parallax shift over 6 months at 400 AU | ~1000 arcsec, ~170 px | 1031 arcsec (166 px) | |
| Total measured shift (parallax minus own motion) | - | 950 arcsec (153 px) | the two parts are opposite around opposition |
| Rate at opposition | ~9 arcsec/day, 0.36 arcsec/h, 1.4 px/day, retrograde | 8.3 arcsec/day, 0.347 arcsec/h, 1.34 px/day, retrograde | 8% lower: the planet's own motion (0.44 arcsec/day) is subtracted and Earth is slower than average in late May |
| Planet's own motion over 6 months | ~80 arcsec | 81 arcsec | |

## D5. Validation against JPL Horizons (section 6.6)

Sedna (83 AU) and Eris (95.5 AU), 13 dates over +-180 days around 2026-02-19, using their
heliocentric osculating elements at that date. Maximum position error 0.33 arcsec (Sedna)
and 0.36 arcsec (Eris), i.e. 0.06 pixel; 180-day shift error 0.01 %; rate error
0.0002 arcsec/hour. Stated test tolerance: 1 arcsec. The residual grows with time from the
epoch because we ignore the giant planets (fixed two-body orbit).

## D6. STATUS WHEN WORK STOPPED (about 12:05, usage limit reached)

- Tier 1 (A + B) is complete: package, 38 passing tests, Streamlit pages 1, 2 and 4, figures.
- Tier 2 (C + D) is NOT done. IRSA was down from about 11:36 until work stopped: SIA2 and TAP
  returned "502 Proxy Error", the data directories returned 404. No SPHEREx exposure list
  and no SPHEREx pixel was downloaded. `data_cache/` holds no data. Page 3 of the app shows
  an honest "no data" message. The code for it (threx/spherex_io.py, threx/blink.py,
  app_realdata.py, scripts/build_cache.py, scripts/build_asteroid_proof.py) is written and
  tested only on a synthetic in-test image, never against real IRSA data.
- Works now: JPL Horizons, sb_ident (slow, up to several minutes; latitude 0 / longitude 0 /
  altitude -6378.137 km works as geocentre), SBDB. Asteroid windows near 90 deg elongation
  were found (for example 9 Metis: 2025-08-03..13 and 2026-05-22..06-01).
- The IRSA cloud mirror (https://nasa-irsa-spherex.s3.us-east-1.amazonaws.com/) was up; real
  headers there carry MJD-AVG and the spacecraft position (X_SC, Y_SC, Z_SC). It has no
  position search, so it was only explored (scripts/explore_aws_mirror.py), not used.
- Not written: demo_script.md, form_text.md, SUBMIT_CHECKLIST.md. README section 6.4 is empty.
- Not done: clean-venv test run, Tier 3. Nothing was committed, pushed or published.
- Correction to D2: the Application Control block is transient (a second run loads scipy and
  astropy.wcs). scipy stays out of the dependencies anyway.
- One file was added outside this folder: `Sentinel_Starter_Kit/.claude/launch.json`
  (preview launcher for the app). Safe to delete.
- Gaia archive was unreachable, so no star preset besides WISE 0855 (challenge page coordinates).

Next steps when IRSA answers:
  python scripts/build_cache.py schedule
  python scripts/build_asteroid_proof.py search --ids 9 15 7 3
  python scripts/build_asteroid_proof.py build --id 9 --window 1

## D7. Official objective (verbatim, read 2026-11-15 on the challenge page, Details tab) and new priority (16:35)

> "Your challenge is to create a public-facing web tool to display images of the sky from the SPHEREx mission, making it easy for anyone to quickly see how they change over time."

> "Users will enter the coordinates of a patch of sky (Right Ascension and declination). Your tool will display a series of images of that patch like a little ____ (the word for a short moving picture, as on the page), making it easy to spot moving objects."

Test targets named: Barnard's Star, Epsilon Indi B, WISE 0855 (RA 08h55m10.831s, Dec -7d14m42.53s), RW Aurigae, R Draconis.
Decision: the app's front page is now the viewer ("1 · SPHEREx sky viewer"); the Planet X pages are "Extra" pages 2 and 3. IRSA SIA2/TAP still returned 502 at 16:33, so data access goes through the public AWS mirror (bucket nasa-irsa-spherex) if the viewer is to have real data.

## D8. Phase 1 done (16:45)
Docs written for the current state (form_text, SUBMIT_CHECKLIST). Viewer is now front page.

## D9. Phase 2 progress (17:00)
Real SPHEREx pixels obtained via the AWS mirror (IRSA SIA2/TAP still 502): Barnard bundle built (59 frames, 2026-03-11..25, 4 MB). Others scanning. Asteroid proof search (6 Hebe) started.

## D10. State at 17:45 (user ordered: no waiting on IRSA, 30 s network timeouts)
- Phase 1 files exist (demo_script, form_text, SUBMIT_CHECKLIST, README 6.4, README typo fixed, slides/threx_deck.pdf via scripts/make_deck.py).
- Real pixels cached: Barnard (59, scan-based, Mar 2026 only), Epsilon Indi B (SIA2-based), WISE 0855 (scan-based), RW Aur (SIA2-based), asteroid proof 6 Hebe. R Draconis NOT cached (scan killed; SIA2 returns 502 intermittently under load).
- SIA2 is intermittent (502). Background jobs were killed at 17:42 on request. IRSA SIA2 cross-check showed my scan misses visits not near 90 deg elongation.

## D11. Block A results (18:00)
- All five bundles rebuilt/checked from IRSA SIA2 exposure lists (all epochs); pixels from the AWS mirror. Frames are an even-in-time subset (max 70 requested; some dropped because the target is within 30 px of the detector edge).
- R Draconis: succeeded at the first attempt (539 SIA2 exposures). So 5 of 5 targets are cached.
- sigma (VARIANCE) is no longer stored in the cache (size); S/N figures were computed before and stored in target_checks.json.
- data_cache (without _work) is about 25 MB, no file over 10 MB.

## D12. Block B results (18:25)
- Full suite: 38 passed + 1 skipped (live test) plus the app smoke tests that now also fail on st.error and test every preset; same result in a clean venv (.venv-clean, pinned requirements-dev.txt) and in a copy of the repo without .venv and data_cache/_work.
- Streamlit started from that copy: health check 200.
- README links all exist (scripts/check_links.py). No data placeholders left. Deck re-rendered and checked (7 slides, overlap fixed).
- Screenshots (generated by the same plotting code as the app, not browser captures): screenshots/10..13 plus 01, 02, 05.
- Lesson: a long-running Streamlit server keeps old imports of threx modules; restart it after code edits.

## D13. Final pass (18:50)
Form text rewritten per form field with character counts; the slide deck is the demo; real Chrome screenshots 01-04; figures renamed figure_*; WISE 0855 sentence now cites Luhman and Esplin 2016.
