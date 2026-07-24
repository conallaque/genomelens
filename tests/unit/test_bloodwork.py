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


# ── V6.1 comprehensive clinical engine ───────────────────────────────────────

def _bm_by_id(cid):
    return bw._BM_BY_ID[cid]


def test_classify_clinical_status_tiers() -> None:
    ldl = _bm_by_id("ldl")           # high-bad, clinical <100, optimal 50-99
    assert bw.classify_clinical(75, ldl) == "optimal"
    assert bw.classify_clinical(120, ldl) == "high"
    assert bw.classify_clinical(200, ldl) == "critical_high"   # critical_high=190
    vitd = _bm_by_id("vitamin_d")    # low-bad, clinical 30-100, optimal 40-60
    assert bw.classify_clinical(50, vitd) == "optimal"
    assert bw.classify_clinical(34, vitd) == "borderline"      # in clinical, below optimal
    assert bw.classify_clinical(18, vitd) == "low"


def test_sex_specific_ranges() -> None:
    hdl = _bm_by_id("hdl")
    # 45 is optimal for a man (>=40 clinical, >=55 optimal → borderline) but LOW for a woman (<50)
    assert bw.classify_clinical(45, hdl, sex="M") == "borderline"
    assert bw.classify_clinical(45, hdl, sex="F") == "low"


def test_derived_markers_computed() -> None:
    vals = {"total_cholesterol": 200, "hdl": 40, "ldl": 130, "triglycerides": 160,
            "fasting_glucose": 100, "fasting_insulin": 10, "iron": 150, "tibc": 300,
            "creatinine": 1.0, "systolic_bp": 130, "diastolic_bp": 80}
    d = bw.compute_derived_markers(vals, {"sex": "M", "age": 40})
    assert d["non_hdl"] == 160
    assert d["tg_hdl_ratio"] == 4.0
    assert d["homa_ir"] == round(100 * 10 / 405, 2)
    assert d["transferrin_sat"] == 50.0
    assert "egfr" in d and d["egfr"] > 0
    assert d["map"] == 97


def test_analyze_clinical_scores_and_flags() -> None:
    labs = {"ldl": 180, "hdl": 65, "fasting_glucose": 85, "vitamin_d": 55}
    res = bw.analyze_clinical_bloodwork(labs, snps_df=None, meta={"sex": "M"})
    assert res["available"] and res["n_markers"] >= 4
    # LDL 180 must be a flag; systems must include Lipids
    names = {f["name"] for f in res["flags"]}
    assert "LDL Cholesterol" in names
    assert any(s["system"] == "Lipids & Cardiovascular" for s in res["systems"])
    assert 0 <= res["overall_score"] <= 100


def test_genotype_note_hfe_iron_overload() -> None:
    import pandas as pd
    # C282Y homozygous (rs1800562 AA) + high ferritin → hemochromatosis note
    df = pd.DataFrame({"genotype": {"rs1800562": "AA"}})
    note = bw._genotype_note("ferritin", "high", df)
    assert "C282Y" in note and "hemochromatosis" in note.lower()
    # no genotype data → no note
    assert bw._genotype_note("ferritin", "high", None) == ""


def test_clinical_attached_to_compare_output(tmp_path) -> None:
    p = tmp_path / "labs.json"
    p.write_text(json.dumps({"ldl": 160, "hdl": 40, "triglycerides": 200}))
    result = bw.compare_bloodwork(str(p), phewas_result=None)
    assert "clinical" in result and result["clinical"]["available"]
    assert result["clinical"]["n_derived"] >= 1   # tg_hdl_ratio at least
