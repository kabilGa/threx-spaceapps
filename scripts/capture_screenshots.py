# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Take real browser screenshots of the running app (needs: pip install playwright, and Google Chrome installed).

  python scripts/capture_screenshots.py
Starts the app on port 8791, captures four pages into screenshots/01_ .. 04_, stops the app.
"""
import subprocess
import sys
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "screenshots"
PORT = 8791
CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
proc = subprocess.Popen([sys.executable, "-m", "streamlit", "run", "app.py", "--server.port", str(PORT),
                         "--server.headless", "true", "--browser.gatherUsageStats", "false"], cwd=ROOT,
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
try:
    time.sleep(8)
    with sync_playwright() as p:
        b = p.chromium.launch(executable_path=CHROME, headless=True)
        pg = b.new_page(viewport={"width": 1400, "height": 3300})
        pg.set_default_timeout(30000)
        pg.goto(f"http://localhost:{PORT}")
        pg.wait_for_selector("text=SPHEREx sky viewer")
        pg.wait_for_selector("text=Frame (press the arrows")
        time.sleep(4)
        pg.screenshot(path=OUT / "01_viewer_barnards_star.png")
        pg.get_by_role("combobox").first.click()
        pg.get_by_role("option", name="Asteroid proof: 6 Hebe (A847 NA)").click()
        pg.wait_for_selector("text=Proof on real pixels")
        pg.set_viewport_size({"width": 1400, "height": 5600})      # tall: the proof section is at the bottom
        time.sleep(6)
        pg.screenshot(path=OUT / "02_viewer_asteroid_6_hebe.png")
        pg.set_viewport_size({"width": 1400, "height": 3300})
        pg.get_by_text("Extra: Planet X signature").click()
        pg.wait_for_selector("text=Sky track")
        time.sleep(6)
        pg.screenshot(path=OUT / "03_planet_x_signature.png")
        pg.get_by_text("Extra: Can SPHEREx see it?").click()
        pg.wait_for_selector("text=What this means in plain words")
        time.sleep(5)
        pg.screenshot(path=OUT / "04_planet_x_detectability.png")
        b.close()
finally:
    proc.terminate()
print("done")
