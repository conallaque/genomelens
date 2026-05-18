"""
Unit tests for `bloodwork.py`.

Critical invariants:
  1. `_classify` thresholds are 0.75 / 1.5 SD (loosened from 0.5 / 1.0 to match
     PRS variance-explained reality). Test the boundaries.
  2. Lab-key synonyms (ldl_c, hgb, b12, e2, …) all map correctly.
  3. Null / non-numeric values in the JSON are dropped, not zero-coerced.
  4. The interpretation text reflects which side dominates.
"""

from __future__ import annotations

import json

import pytest

import bloodwork as bw


# ── Threshold boundaries ─────────────────────────────────────────────────────

@pytest.mark.parametrize("delta_sd,expected", [
    (0.0,    "Confirmed"),
    (0.5,    "Confirmed"),
    (0.74,   "Confirmed"),
    (0.75,   "Confirmed"),     # inclusive upper bound
    (0.76,   "Partial"),
    (1.0,    "Partial"),
    (1.49,   "Partial"),
    (1.5,    "Partial"),       # inclusive upper bound
    (1.51,   "Diverged"),
    (-1.51,  "Diverged"),
    (3.0,    "Diverged"),
])
def test_classify_thresholds(delta_sd: float, expected: str) -> None:
    label, _ = bw._classify(delta_sd)
    assert label == expected


# ── Synonym mapping ──────────────────────────────────────────────────────────

@pytest.mark.parametrize("key,trait", [
    ("ldl",                 "LDL cholesterol"),
    ("LDL_C",               "LDL cholesterol"),   # case-insensitive
    ("ldl_cholesterol",     "LDL cholesterol"),
    ("hgb",                 "Hemoglobin"),
    ("hb",                  "Hemoglobin"),
    ("b12",                 "Serum vitamin B12"),
    ("vitamin_b12",         "Serum vitamin B12"),
    ("e2",                  "Estradiol (women)"),
    ("hs_crp",              "C-Reactive Protein"),
    ("a1c",                 "HbA1c (predicted)"),
])
def test_lab_key_synonyms(key: str, trait: str) -> None:
    normalised = bw._normalize_key(key)
    assert bw._LAB_TO_PHEWAS.get(normalised) == trait


def test_normalize_key_handles_dashes_and_spaces() -> None:
    assert bw._normalize_key(" Ferritin ") == "ferritin"
    assert bw._normalize_key("HDL-C") == "hdl_c"


# ── load_bloodwork ───────────────────────────────────────────────────────────

def test_load_bloodwork_drops_nulls_and_strings(tmp_path) -> None:
    p = tmp_path / "labs.json"
    p.write_text(json.dumps({
        "ldl": 142, "hdl": None, "crp": "not measured",
        "vitamin_d": 22.5, "ferritin": "180",   # numeric string is OK
    }))
    cleaned = bw.load_bloodwork(str(p))
    assert "ldl" in cleaned and cleaned["ldl"] == 142.0
    assert "vitamin_d" in cleaned
    assert "ferritin" in cleaned and cleaned["ferritin"] == 180.0
    assert "hdl" not in cleaned
    assert "crp" not in cleaned


def test_load_bloodwork_missing_file_raises(tmp_path) -> None:
    with pytest.raises(FileNotFoundError):
        bw.load_bloodwork(str(tmp_path / "nope.json"))


def test_load_bloodwork_non_object_raises(tmp_path) -> None:
    p = tmp_path / "wrong.json"
    p.write_text(json.dumps([1, 2, 3]))
    with pytest.raises(ValueError, match="object at the top level"):
        bw.load_bloodwork(str(p))


# ── End-to-end comparison ────────────────────────────────────────────────────

def _phewas_stub() -> dict:
    """Minimal PheWAS-shaped result for comparison tests."""
    return {"traits": {
        "LDL cholesterol": {
            "mean": 110, "sd": 32, "unit": "mg/dL", "category": "Lipids",
            "result": {"status": "ok", "predicted_value": 125,
                       "tier": "Above average", "callability_pct": 70,
                       "n_used": 6, "n_total": 8},
        },
        "C-Reactive Protein": {
            "mean": 1.5, "sd": 1.5, "unit": "mg/L", "category": "Inflammation",
            "result": {"status": "ok", "predicted_value": 1.8,
                       "tier": "Average", "callability_pct": 75,
                       "n_used": 3, "n_total": 4},
        },
        "25-OH Vitamin D": {
            "mean": 25, "sd": 10, "unit": "ng/mL", "category": "Vitamins",
            "result": {"status": "ok", "predicted_value": 19,
                       "tier": "Below average", "callability_pct": 75,
                       "n_used": 3, "n_total": 4},
        },
    }}


def test_compare_bloodwork_confirms_close_match(tmp_path) -> None:
    p = tmp_path / "labs.json"
    p.write_text(json.dumps({"ldl": 130, "crp": 1.9}))   # both within 0.75 SD
    result = bw.compare_bloodwork(str(p), _phewas_stub())
    assert result["status"] == "ok"
    assert result["n_confirmed"] == 2
    assert result["n_diverged"] == 0


def test_compare_bloodwork_flags_diverged(tmp_path) -> None:
    p = tmp_path / "labs.json"
    p.write_text(json.dumps({"crp": 6.5}))   # 1.5 SD = 4.05; 6.5 is >>1.5 SD above 1.8
    result = bw.compare_bloodwork(str(p), _phewas_stub())
    assert result["n_diverged"] == 1
    row = result["rows"][0]
    assert row["verdict"] == "Diverged"
    assert row["delta_sd"] > 1.5


def test_unmatched_labs_listed(tmp_path) -> None:
    p = tmp_path / "labs.json"
    p.write_text(json.dumps({"ldl": 110, "homocysteine": 9}))  # latter not in PheWAS stub
    result = bw.compare_bloodwork(str(p), _phewas_stub())
    assert "homocysteine" in result["unmatched"]


def test_no_phewas_short_circuits(tmp_path) -> None:
    p = tmp_path / "labs.json"
    p.write_text(json.dumps({"ldl": 120}))
    result = bw.compare_bloodwork(str(p), phewas_result=None)
    assert result["status"] == "no_phewas"


# ── HTML rendering ───────────────────────────────────────────────────────────

def test_html_contains_explainer(tmp_path) -> None:
    """The 'how to read this table' explainer is the key UX upgrade — make sure
    it renders so users don't misinterpret Diverged as 'genetics wrong'."""
    p = tmp_path / "labs.json"
    p.write_text(json.dumps({"crp": 1.9}))
    result = bw.compare_bloodwork(str(p), _phewas_stub())
    html = bw.render_bloodwork_html(result)
    assert "How to read this table" in html
    assert "0.75" in html and "1.5" in html
    assert "5–20%" in html or "5-20%" in html
