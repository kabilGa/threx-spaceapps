"""Exploration helper: look at the IRSA cloud mirror of SPHEREx (public AWS bucket) without
downloading whole files. Lists one planning period and reads only the FITS headers of one
file with HTTP range requests.  (not used by the app)

Run:  python scripts/explore_aws_mirror.py qr3/level2/2026W21_2A/
"""
import re
import sys

import requests
from astropy.io import fits

BUCKET = "https://nasa-irsa-spherex.s3.us-east-1.amazonaws.com/"


def list_prefix(prefix, delimiter="/", max_keys=1000, token=None):
    params = {"list-type": "2", "prefix": prefix, "max-keys": str(max_keys)}
    if delimiter:
        params["delimiter"] = delimiter
    if token:
        params["continuation-token"] = token
    xml = requests.get(BUCKET, params=params, timeout=60).text
    prefixes = re.findall(r"<CommonPrefixes><Prefix>([^<]+)</Prefix>", xml)
    keys = list(zip(re.findall(r"<Key>([^<]+)</Key>", xml), map(int, re.findall(r"<Size>(\d+)</Size>", xml))))
    nxt = re.findall(r"<NextContinuationToken>([^<]+)<", xml)
    return prefixes, keys, (nxt[0] if nxt else None)


def read_range(key, start, length):
    r = requests.get(BUCKET + key, headers={"Range": f"bytes={start}-{start + length - 1}"}, timeout=60)
    r.raise_for_status()
    return r.content


def read_header_at(key, offset):
    """Read one FITS header starting at byte `offset`. Returns (Header, bytes used)."""
    raw = b""
    while True:
        raw += read_range(key, offset + len(raw), 2880 * 8)
        blocks = [raw[i:i + 2880] for i in range(0, len(raw), 2880)]
        for n, blk in enumerate(blocks):
            cards = [blk[j:j + 80] for j in range(0, 2880, 80)]
            if any(c.startswith(b"END" + b" " * 77) for c in cards):
                used = (n + 1) * 2880
                return fits.Header.fromstring(raw[:used].decode("ascii", "replace")), used


if __name__ == "__main__":
    prefix = sys.argv[1] if len(sys.argv) > 1 else "qr3/level2/2026W21_2A/"
    p, k, _ = list_prefix(prefix)
    print("sub-prefixes:", p)
    p2, k2, _ = list_prefix(p[0])
    print("detector dirs:", p2)
    det2 = [x for x in p2 if x.rstrip("/").endswith("/2")][0]
    keys, token = [], None
    while True:
        _, kk, token = list_prefix(det2, token=token)
        keys += kk
        if not token:
            break
    print(len(keys), "files in", det2, "sizes MB:", round(keys[0][1] / 1e6, 1), "..", round(max(s for _, s in keys) / 1e6, 1))
    print(keys[0][0], keys[-1][0])
    key = keys[len(keys) // 2][0]
    h0, n0 = read_header_at(key, 0)
    h1, n1 = read_header_at(key, n0)
    print("PRIMARY cards:", len(h0), "IMAGE header bytes:", n1)
    for card in h1.cards:
        if not card.keyword.startswith(("A_", "B_", "AP_", "BP_", "HISTORY", "COMMENT")) and "L2 N_" not in card.keyword:
            print(str(card).rstrip())
