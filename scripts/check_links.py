# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Check that every local link/image in the Markdown files points to an existing file. Run: python scripts/check_links.py"""
import re
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
bad = 0
for md in root.glob("*.md"):
    for target in re.findall(r"\]\(([^)\s]+)\)", md.read_text(encoding="utf8")):
        if target.startswith(("http", "#", "mailto")):
            continue
        if not (root / target.split("#")[0]).exists():
            print(f"{md.name}: missing {target}")
            bad += 1
print("all local links exist" if not bad else f"{bad} broken link(s)")
sys.exit(1 if bad else 0)
