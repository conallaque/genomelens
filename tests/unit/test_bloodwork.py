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

from wellness import bloodwork as bw

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
    normalized = bw._normalize_key(key)
    assert bw._LAB_TO_PHEWAS.get(normalized) == trait


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


# ── V6.2 advanced composite indices & biological age ──────────────────────────

def _adv(labs, meta=None):
    return bw.compute_advanced_indices(labs, {}, meta or {"sex": "M", "age": 41})


def _idx(adv, iid):
    return next((i for i in adv["indices"] if i["id"] == iid), None)


def test_phenoage_healthy_is_younger() -> None:
    # Healthy 40yo → biological age materially younger (validated hand-calc ~32.7)
    pa = bw._phenoage(4.5, 1.0, 90, 0.5, 30, 90, 13.0, 65, 5.5, 40)
    assert 31.5 < pa < 34.0
    # Unhealthy 40yo → older
    pa2 = bw._phenoage(4.0, 1.1, 130, 5.0, 20, 92, 15.5, 120, 9.0, 40)
    assert pa2 > 48


def test_phenoage_requires_si_conversion() -> None:
    # Sanity: feeding glucose raw (mg/dL) without conversion would explode xb;
    # our implementation converts, so a normal panel stays near chronological age.
    pa = bw._phenoage(4.4, 1.0, 100, 1.0, 28, 90, 13.5, 75, 6.0, 45)
    assert 38 < pa < 52   # plausible, not hundreds


def test_tyg_index() -> None:
    adv = _adv({"triglycerides": 180, "fasting_glucose": 104})
    i = _idx(adv, "tyg")
    assert i is not None and abs(i["value"] - 9.14) < 0.02
    assert i["status"] == "cl-high"   # >=8.75


def test_sampson_ldl() -> None:
    adv = _adv({"total_cholesterol": 210, "hdl": 42, "triglycerides": 180})
    i = _idx(adv, "ldl_sampson")
    assert i is not None and abs(i["value"] - 136) <= 2


def test_fib4_and_sii() -> None:
    adv = _adv({"ast": 30, "alt": 42, "platelets": 240, "neutrophils": 4.0,
                "lymphocytes": 1.6, "monocytes": 0.5}, {"sex": "M", "age": 41})
    assert abs(_idx(adv, "fib4_adv")["value"] - 0.79) < 0.03
    assert _idx(adv, "sii")["value"] == 600
    assert abs(_idx(adv, "siri")["value"] - 1.25) < 0.02


def test_metabolic_syndrome_flags_present() -> None:
    adv = _adv({"triglycerides": 180, "hdl": 38, "fasting_glucose": 104,
                "systolic_bp": 135, "diastolic_bp": 88}, {"sex": "M", "age": 41})
    i = _idx(adv, "metsyn")
    assert i["status"] == "cl-high"   # >=3 criteria
    assert i["status_label"] == "Present"


def test_corrected_calcium_and_anion_gap() -> None:
    adv = _adv({"calcium": 9.4, "albumin": 4.4, "sodium": 140,
                "chloride": 102, "co2": 25})
    assert abs(_idx(adv, "corr_ca")["value"] - 9.1) < 0.05
    assert _idx(adv, "anion_gap")["value"] == 13


def test_advanced_attached_to_clinical() -> None:
    res = bw.analyze_clinical_bloodwork(
        {"triglycerides": 150, "fasting_glucose": 95, "hdl": 45, "total_cholesterol": 190},
        snps_df=None, meta={"sex": "M", "age": 40})
    assert "advanced" in res and res["advanced"]["available"]
    assert res["advanced"]["n_indices"] >= 3


def test_mets_ir_and_fli_and_pni_aisi() -> None:
    adv = _adv({"fasting_glucose": 100, "triglycerides": 150, "hdl": 45, "bmi": 28,
                "ggt": 40, "waist": 95,
                "albumin": 4.4, "lymphocytes": 1.8, "neutrophils": 4.0,
                "monocytes": 0.5, "platelets": 250})
    assert _idx(adv, "mets_ir") is not None
    assert _idx(adv, "fli") is not None
    pni = _idx(adv, "pni")
    assert pni is not None and abs(pni["value"] - (10*4.4 + 0.005*1800)) < 0.1
    assert _idx(adv, "aisi") is not None


def test_phenoage_levers_and_mortality() -> None:
    pa, mort = bw._phenoage_core(4.4, 1.0, 104, 3.4, 30, 90, 13.8, 75, 6.1, 41)
    assert 0.0 < mort < 1.0
    lev = bw.phenoage_levers(4.4, 1.0, 104, 3.4, 30, 90, 13.8, 75, 6.1, 41)
    assert lev["recoverable_years"] > 0
    # RDW (largest CBC coefficient) should be a top lever here
    assert lev["levers"] and lev["levers"][0]["marker"] in ("RDW", "Fasting glucose", "hs-CRP")
    # each lever quantifies years recovered by optimizing that marker
    assert all(v["years_cost"] > 0 for v in lev["levers"])


def test_genetic_longevity_reads_variants() -> None:
    import pandas as pd
    df = pd.DataFrame({"genotype": {
        "rs2802292": "GG",   # FOXO3 longevity (favorable)
        "rs7412": "CT",      # APOE ε2 present (favorable)
        "rs429358": "TT",    # no ε4
    }})
    gl = bw._genetic_longevity(df)
    assert gl is not None and gl["n_favorable"] >= 2
    genes = {v["gene"] for v in gl["variants"]}
    assert "FOXO3" in genes and "APOE" in genes
    assert bw._genetic_longevity(None) is None


def test_bioage_simulator_renders_sliders() -> None:
    inputs = {"albumin": 4.4, "creatinine": 1.0, "glucose": 100, "crp": 1.0,
              "lymph_pct": 30, "mcv": 90, "rdw": 13.5, "alp": 70, "wbc": 6.0, "age": 41}
    html = bw._render_bioage_simulator(inputs)
    assert 'type="range"' in html and "phenoAge" in html and "baUpdate" in html


def test_prevent_ascvd_matches_reference() -> None:
    # Validated against the cross-checked open-source PREVENT implementation.
    assert abs(bw.prevent_ascvd_10yr("F", 50, 200, 50, 120, 90, 0, 0, 0, 0) - 1.3) < 0.3
    assert abs(bw.prevent_ascvd_10yr("M", 55, 210, 42, 132, 90, 0, 0, 0, 0) - 3.7) < 0.3
    # diabetic smoker with high BP → materially higher
    hi = bw.prevent_ascvd_10yr("F", 65, 240, 40, 150, 70, 1, 1, 0, 0)
    assert hi > 10
    assert bw.prevent_ascvd_10yr("X", 50, 200, 50, 120, 90) is None


def test_prevent_index_via_clinical() -> None:
    res = bw.analyze_clinical_bloodwork(
        {"total_cholesterol": 210, "hdl": 42, "systolic_bp": 132, "creatinine": 1.0},
        snps_df=None, meta={"sex": "M", "age": 55})
    ids = {i["id"] for i in res["advanced"]["indices"]}
    assert "prevent_ascvd" in ids


def test_longitudinal_trajectory(tmp_path) -> None:
    p = tmp_path / "labs.json"
    p.write_text(json.dumps({"sex": "M", "age": 41, "history": [
        {"date": "2024-01-01", "ldl": 160, "hdl": 40, "crp": 4.0,
         "albumin": 4.2, "creatinine": 1.0, "fasting_glucose": 110, "lymphocytes": 1.5,
         "wbc": 6.5, "mcv": 91, "rdw": 14.2, "alp": 85},
        {"date": "2025-01-01", "ldl": 110, "hdl": 52, "crp": 1.0,
         "albumin": 4.6, "creatinine": 0.95, "fasting_glucose": 88, "lymphocytes": 2.0,
         "wbc": 5.4, "mcv": 89, "rdw": 13.0, "alp": 68},
    ]}))
    res = bw.compare_bloodwork(str(p), phewas_result=None)
    tr = res["clinical"].get("trajectory")
    assert tr is not None and tr["n_timepoints"] == 2
    ldl = next(m for m in tr["metrics"] if m["key"] == "ldl")
    assert ldl["first"] == 160 and ldl["last"] == 110 and ldl["improving"]
    bio = next((m for m in tr["metrics"] if m["key"] == "bioage"), None)
    assert bio is not None and bio["improving"]
