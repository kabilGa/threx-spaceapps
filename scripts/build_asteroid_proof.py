"""Tier 2 D: catch a KNOWN asteroid in real SPHEREx pixels and check it against JPL.

Steps (all with real services, nothing invented):
 1. JPL Horizons: daily positions of a bright asteroid; keep the days when it is about 90 deg
    from the Sun, because that is where SPHEREx scans.
 2. IRSA SIA2: which SPHEREx exposures cover the asteroid's position on those days.
 3. JPL Horizons again: exact predicted position at each exposure mid-time.
 4. IRSA cutouts around the predicted position; align north-up; measure the centre and the
    flux of the source nearest to the prediction.
 5. JPL sb_ident: independent list of known bodies in that small field at that time.
 6. Save everything as a cached bundle (data_cache/bundle_<name>.npz/.json).

  python scripts/build_asteroid_proof.py search  --ids 9 15 7 3        # steps 1-3, prints candidates
  python scripts/build_asteroid_proof.py build   --id 9 --window 0 --detectors 1 2 3 --frames 6
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from astropy.time import Time

from threx import blink, jpl, spherex_io
from threx.observer import unit_vector, separation_arcsec

START, STOP = "2025-05-01", "2026-11-01"       # span of public SPHEREx data (approx.)
SCRATCH = Path(__file__).resolve().parents[1] / "data_cache" / "_work"


def quadrature_windows(ast_id, tol_deg=4.0):
    """Runs of days when the asteroid is within tol_deg of 90 deg solar elongation."""
    jds = np.arange(Time(START).jd, Time(STOP).jd, 1.0)
    eph = []
    for chunk in np.array_split(jds, int(np.ceil(len(jds) / 40))):
        eph += jpl.horizons_observer(f"{ast_id};", list(chunk))
    df = pd.DataFrame(eph)
    near = df[np.abs(df["elong"] - 90.0) < tol_deg].reset_index(drop=True)
    if near.empty:
        return []
    groups = (near["jd"].diff() > 1.5).cumsum()
    return [g.reset_index(drop=True) for _, g in near.groupby(groups)]


def matching_exposures(ast_id, window):
    """SPHEREx exposures that should contain the asteroid during one window."""
    found = []
    for row in window.to_dict("records"):
        mjd = row["jd"] - 2400000.5
        df = spherex_io.query_exposures(row["ra"], row["dec"])
        df = df[np.abs(df["t_mid"] - mjd) < 0.5]            # asteroid moves << field size in half a day
        found.append(df)
    if not found:
        return pd.DataFrame()
    df = pd.concat(found, ignore_index=True)
    if df.empty:
        return df
    df = df.drop_duplicates("obs_id").sort_values("t_mid").reset_index(drop=True)
    # exact predicted position at each exposure mid-time
    eph = []
    jd = (df["t_mid"] + 2400000.5).to_numpy()
    for chunk in np.array_split(jd, int(np.ceil(len(jd) / 40))):
        eph += jpl.horizons_observer(f"{ast_id};", list(chunk))
    for key in ("ra", "dec", "dra_cosd", "ddec", "vmag", "delta", "r", "elong"):
        df["pred_" + key] = [e.get(key) for e in eph]
    df["sep_from_image_centre_deg"] = separation_arcsec(
        unit_vector(df["pred_ra"].to_numpy(), df["pred_dec"].to_numpy()),
        unit_vector(df["s_ra"].to_numpy(dtype=float), df["s_dec"].to_numpy(dtype=float))) / 3600
    return df


def do_search(args):
    SCRATCH.mkdir(parents=True, exist_ok=True)
    for ast_id in args.ids:
        windows = quadrature_windows(ast_id)
        print(f"asteroid {ast_id}: {len(windows)} window(s) near 90 deg elongation")
        for k, w in enumerate(windows):
            print(f"  window {k}: {Time(w['jd'].iloc[0], format='jd').iso[:10]} .. "
                  f"{Time(w['jd'].iloc[-1], format='jd').iso[:10]}  V = {w['vmag'].mean():.1f}")
            try:
                df = matching_exposures(ast_id, w)
            except Exception as exc:
                print("    IRSA query failed:", type(exc).__name__, str(exc)[:120])
                continue
            df.to_csv(SCRATCH / f"asteroid_{ast_id}_window{k}.csv", index=False)
            print(f"    {len(df)} matching exposures")
            if len(df):
                print(df[["time_utc", "detector", "lam_min_um", "lam_max_um", "pred_ra", "pred_dec",
                          "sep_from_image_centre_deg"]].to_string(max_rows=40))


def do_build(args):
    df = pd.read_csv(SCRATCH / f"asteroid_{args.id}_window{args.window}.csv")
    df = df[df["detector"].isin(args.detectors)].sort_values("t_mid").reset_index(drop=True)
    if args.t0 is not None:                                   # keep a run of frames starting at index t0
        df = df.iloc[args.t0:].reset_index(drop=True)
    if args.max_hours:
        df = df[(df["t_mid"] - df["t_mid"].iloc[0]) * 24 <= args.max_hours].reset_index(drop=True)
    rows = blink.select_frames(df, args.frames)
    info = jpl.sbdb(str(args.id))
    full_name = info["object"]["fullname"].strip()
    # field centre: middle of the predicted positions, so the whole track fits
    ra_c = float(0.5 * (rows["pred_ra"].min() + rows["pred_ra"].max()))
    dec_c = float(0.5 * (rows["pred_dec"].min() + rows["pred_dec"].max()))
    name = f"Asteroid proof: {full_name}"
    from threx import aws_mirror as aws
    reader = (lambda row: aws.read_cutout_s3(row["key"], ra_c, dec_c, args.half)) if "key" in rows else None
    bundle = blink.build_bundle(name, ra_c, dec_c, rows, half_size_arcsec=args.half, reader=reader, workers=8,
                                progress=lambda k, n, row: print(f"  cutout {k + 1}/{n} {row['time_utc']}"))
    kept = {f["obs_id"] for f in bundle["meta"]["frames"]}
    rows = rows[rows["obs_id"].isin(kept)].reset_index(drop=True)
    frames = []
    for k, row in enumerate(rows.to_dict("records")):
        col, rw = blink.offset_to_pixel(bundle, row["pred_ra"], row["pred_dec"])
        east, north = blink.pixel_to_offset_arcsec(bundle, col, rw)
        f = {"pred_ra": row["pred_ra"], "pred_dec": row["pred_dec"], "pred_east_arcsec": float(east),
             "pred_north_arcsec": float(north), "vmag_horizons": row["pred_vmag"],
             "rate_arcsec_per_hour": float(np.hypot(row["pred_dra_cosd"], row["pred_ddec"])),
             "delta_au": row["pred_delta"], "r_au": row["pred_r"], "elong_deg": row["pred_elong"]}
        phot = blink.photometry(bundle, k, col, rw)
        if phot:                                   # nothing to measure when the target is outside this frame
            c, r = blink.centroid(bundle, k, col, rw, box=5)
            me, mn = blink.pixel_to_offset_arcsec(bundle, c, r)
            f.update({"meas_east_arcsec": float(me), "meas_north_arcsec": float(mn),
                      "miss_arcsec": float(np.hypot(me - east, mn - north))})
            f.update(phot)
        frames.append(f)
        print(f"  frame {k + 1}: miss {f['miss_arcsec']:.1f} arcsec, S/N {f.get('snr')}, AB {f.get('mag_ab')}, "
              f"flagged in aperture {f.get('n_flagged')}")
    got = [f for f in frames if f.get("miss_arcsec") is not None]
    first, last = got[0], got[-1]
    dt_h = (rows["t_mid"].iloc[-1] - rows["t_mid"].iloc[0]) * 24
    proof = {
        "name": full_name, "frames": frames,
        "displacement": {
            "dt_hours": float(dt_h),
            "predicted_arcsec": float(np.hypot(last["pred_east_arcsec"] - first["pred_east_arcsec"],
                                               last["pred_north_arcsec"] - first["pred_north_arcsec"])),
            "measured_arcsec": float(np.hypot(last["meas_east_arcsec"] - first["meas_east_arcsec"],
                                              last["meas_north_arcsec"] - first["meas_north_arcsec"])),
        },
        "sbdb": {"full name": full_name, "orbit class": info["object"]["orbit_class"]["name"],
                 **{p["name"]: f"{p['value']} {p.get('units') or ''}".strip()
                    for p in info.get("phys_par", []) if p["name"] in ("H", "diameter", "albedo", "rot_per")}},
        "target_note": f"Field centred on the track of {full_name} during the SPHEREx exposures. "
                       "A known asteroid: nothing new was found.",
        "note": "",
    }
    # independent check with sb_ident at the time of the first frame
    t_first = Time(rows["t_mid"].iloc[0], format="mjd", scale="utc")
    try:
        ident = jpl.sb_ident(t_first.strftime("%Y-%m-%d_%H:%M:%S"), first["pred_ra"], first["pred_dec"],
                             half_width_deg=args.ident_hw)
        fields = ident.get("fields_second") or ident.get("fields_first") or []
        data = ident.get("data_second_pass") or ident.get("data_first_pass") or []
        table = [dict(zip(fields, r)) for r in data]
        proof["sb_ident"] = {
            "summary": f"JPL sb_ident, {t_first.iso[:19]} UTC, box of +-{args.ident_hw} deg around the predicted "
                       f"position, observer at the centre of the Earth: {len(table)} known small "
                       f"bod{'y' if len(table) == 1 else 'ies'} with a magnitude estimate.",
            "rows": table[:15], "n": len(table),
        }
        print("sb_ident:", len(table), "objects;", [r.get(fields[0]) for r in table[:5]])
    except Exception as exc:
        proof["sb_ident"] = {"summary": f"sb_ident could not be queried ({type(exc).__name__}: {str(exc)[:100]}).",
                             "rows": []}
        print("sb_ident failed:", exc)
    bundle["meta"]["proof"] = proof
    path = blink.save_bundle(bundle)
    try:
        spherex_io.save_exposures(name, spherex_io.query_exposures(ra_c, dec_c), ra_c, dec_c)
    except Exception as exc:                       # IRSA search down: not needed for the proof
        print("exposure list not saved:", type(exc).__name__)
    print("saved", path, f"({path.stat().st_size / 1e6:.2f} MB)")
    print(json.dumps(proof["displacement"], indent=1))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("search")
    s.add_argument("--ids", nargs="+", default=["9", "15", "7", "3"])
    b = sub.add_parser("build")
    b.add_argument("--id", required=True)
    b.add_argument("--window", type=int, default=0)
    b.add_argument("--detectors", type=int, nargs="*", default=[1, 2, 3])
    b.add_argument("--frames", type=int, default=6)
    b.add_argument("--half", type=float, default=150.0)
    b.add_argument("--t0", type=int, default=None)
    b.add_argument("--max-hours", type=float, default=None)
    b.add_argument("--ident-hw", type=float, default=0.05)
    a = p.parse_args()
    {"search": do_search, "build": do_build}[a.cmd](a)
