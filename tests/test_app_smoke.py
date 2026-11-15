# Written with SI assistance (Claude Code, Anthropic); see AI_USE.md.
"""Smoke test: every page of the Streamlit app renders without an exception (offline:
page 3 uses the cached data, never the network)."""
import os

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest

APP = os.path.join(os.path.dirname(__file__), "..", "app.py")


@pytest.mark.parametrize("page_index", [0, 1, 2, 3])
def test_page_renders(page_index, monkeypatch):
    monkeypatch.setenv("THREX_OFFLINE", "1")          # page 3: cached data only
    at = AppTest.from_file(APP, default_timeout=120).run()
    assert not at.exception, at.exception
    at.sidebar.radio[0].set_value(at.sidebar.radio[0].options[page_index]).run()
    assert not at.exception, at.exception
    assert not at.error, [e.value for e in at.error]
    assert len(at.title) >= 1


@pytest.mark.parametrize("target", ["Barnard's Star", "Epsilon Indi B", "WISE 0855-0714", "RW Aurigae",
                                    "R Draconis", "Asteroid proof: 6 Hebe (A847 NA)"])
def test_viewer_each_preset(target, monkeypatch):
    monkeypatch.setenv("THREX_OFFLINE", "1")
    at = AppTest.from_file(APP, default_timeout=180).run()
    at.selectbox[0].select(target).run()
    assert not at.exception, at.exception
    assert not at.error, [e.value for e in at.error]
