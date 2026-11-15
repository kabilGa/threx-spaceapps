# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Fill data_cache/ with REAL data from IRSA so the app also works offline.

  python scripts/build_cache.py schedule                 # exposure lists for all preset positions
  python scripts/build_cache.py blink "<preset name>" --detectors 2 --frames 8 --half 150
  python scripts/build_cache.py blink custom --name "My field" --ra 60 --dec 20

Needs a network connection and IRSA being up. Nothing is invented: on failure it stops.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from threx import blink, spherex_io, targets


def do_schedule(args):
    for name, t in targets.TARGETS.items():
        df = spherex_io.query_exposures(t["ra"], t["dec"])
        meta = spherex_io.save_exposures(name, df, t["ra"], t["dec"])
        visits = spherex_io.find_visits(df)
        print(f"{name}: {meta['rows']} exposures, {len(visits)} visits "
              f"{[(v['n'], round(v['duration_days'], 1)) for v in visits]}")


def do_blink(args):
    if args.target == "custom":
        name, ra, dec = args.name, args.ra, args.dec
    else:
        name = next(n for n in targets.TARGETS if n.lower().startswith(args.target.lower()))
        ra, dec = targets.TARGETS[name]["ra"], targets.TARGETS[name]["dec"]
    df = spherex_io.query_exposures(ra, dec)
    spherex_io.save_exposures(name, df, ra, dec)
    rows = blink.select_frames(df, args.frames, args.detectors)
    print(f"{name}: {len(df)} exposures, downloading {len(rows)} cutouts ...")
    bundle = blink.build_bundle(name, ra, dec, rows, half_size_arcsec=args.half,
                                progress=lambda k, n, row: print(f"  {k + 1}/{n} {row['time_utc']} det {row['detector']:.0f}"))
    path = blink.save_bundle(bundle)
    print("saved", path, f"({path.stat().st_size / 1e6:.2f} MB); failed:", bundle["meta"]["failed"])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("schedule")
    b = sub.add_parser("blink")
    b.add_argument("target")
    b.add_argument("--name", default="Custom field")
    b.add_argument("--ra", type=float)
    b.add_argument("--dec", type=float)
    b.add_argument("--detectors", type=int, nargs="*", default=[2])
    b.add_argument("--frames", type=int, default=8)
    b.add_argument("--half", type=float, default=150.0)
    a = p.parse_args()
    {"schedule": do_schedule, "blink": do_blink}[a.cmd](a)
