# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
import sys
from pathlib import Path

# make "import threx" work when pytest is started from the repo root
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
