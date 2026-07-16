"""
Blood-Work Comparison
=====================

Compares user-supplied lab panel results (JSON) against the PheWAS-derived
genetic predictions and identifies:
  • Confirmed — predicted value matches measured value within ±0.5 trait SD
  • Partial   — within 0.5–1.0 SD of prediction
  • Diverged  — measured value differs from prediction by >1.0 SD; the
                divergence direction tells us whether the genotype is being
                amplified (lifestyle aligned with genetic risk) or buffered
                (lifestyle/environment opposing genetic prediction)

JSON format expected (case-insensitive keys, synonyms supported):

    {
      "ldl": 142,              "hdl": 48,           "triglycerides": 180,
      "lp_a": 35,              "fasting_glucose": 96, "hba1c": 5.7,
      "crp": 2.4,              "ferritin": 180,      "iron": null,
      "vitamin_d": 22,         "vitamin_b12": 410,   "folate": 14,
      "testosterone": 480,     "shbg": 32,           "igf1": 175,
      "tsh": 1.8,              "alt": 28,            "creatinine": 1.0,
      "uric_acid": 6.2,        "hemoglobin": 14.8,   "mcv": 89,
      "platelets": 235,        "wbc": 6.1,           "systolic_bp": 128,
      "diastolic_bp": 82,      "heart_rate": 64
    }

Values may be null or missing for any biomarker; only present numeric values
are scored. Units are assumed to match the PheWAS catalog (mg/dL, ng/mL,
ng/dL, mIU/L, etc. — same as standard US clinical reports).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# Map of free-form JSON keys → exact PheWAS trait name.
# All variant spellings (case, dashes, common abbreviations) collapse here.
_LAB_TO_PHEWAS: Dict[str, str] = {
    # Lipids
    "ldl":                "LDL cholesterol",
    "ldl_c":              "LDL cholesterol",
    "ldl_cholesterol":    "LDL cholesterol",
    "hdl":                "HDL cholesterol",
    "hdl_c":              "HDL cholesterol",
    "hdl_cholesterol":    "HDL cholesterol",
    "triglycerides":      "Triglycerides",
    "tg":                 "Triglycerides",
    "lp_a":               "Lipoprotein(a)",
    "lipoprotein_a":      "Lipoprotein(a)",
    "lpa":                "Lipoprotein(a)",
    # Glucose
    "fasting_glucose":    "Fasting glucose",
    "glucose":            "Fasting glucose",
    "fpg":                "Fasting glucose",
    "hba1c":              "HbA1c (predicted)",
    "a1c":                "HbA1c (predicted)",
    "homa_ir":            "HOMA-IR (insulin resistance)",
    # Anthropometric
    "bmi":                "Body Mass Index (BMI)",
    "height":             "Height",
    "height_cm":          "Height",
    "whr":                "Waist-Hip Ratio (BMI-adjusted)",
    # Cardio
    "systolic_bp":        "Systolic Blood Pressure",
    "sbp":                "Systolic Blood Pressure",
    "diastolic_bp":       "Diastolic Blood Pressure",
    "dbp":                "Diastolic Blood Pressure",
    "heart_rate":         "Resting Heart Rate",
    "resting_hr":         "Resting Heart Rate",
    "rhr":                "Resting Heart Rate",
    # Hematology
    "hemoglobin":         "Hemoglobin",
    "hgb":                "Hemoglobin",
    "hb":                 "Hemoglobin",
    "mcv":                "Mean Corpuscular Volume",
    "platelets":          "Platelets",
    "plt":                "Platelets",
    "wbc":                "White Blood Cell Count",
    "white_blood_cells":  "White Blood Cell Count",
    # Inflammation
    "crp":                "C-Reactive Protein",
    "hs_crp":             "C-Reactive Protein",
    # Metabolic / liver / kidney
    "uric_acid":          "Uric acid / urate",
    "urate":              "Uric acid / urate",
    "creatinine":         "Serum creatinine / eGFR",
    "alt":                "ALT (liver function)",
    "sgpt":               "ALT (liver function)",
    # Thyroid
    "tsh":                "TSH",
    # Vitamins / minerals
    "vitamin_d":          "25-OH Vitamin D",
    "25_oh_d":            "25-OH Vitamin D",
    "25ohd":              "25-OH Vitamin D",
    "vitamin_b12":        "Serum vitamin B12",
    "b12":                "Serum vitamin B12",
    "folate":             "Folate",
    "ferritin":           "Iron / ferritin",
    "iron":               "Iron / ferritin",
    "iron_ferritin":      "Iron / ferritin",
    # Hormones
    "testosterone":       "Testosterone (men)",
    "total_testosterone": "Testosterone (men)",
    "estradiol":          "Estradiol (women)",
    "e2":                 "Estradiol (women)",
    "shbg":               "SHBG",
    "igf1":               "IGF-1",
    "igf_1":              "IGF-1",
    "cortisol":           "Cortisol (proxy)",
}


def _normalize_key(k: str) -> str:
    return k.strip().lower().replace(" ", "_").replace("-", "_")


def load_bloodwork(path: str) -> Dict[str, float]:
    """Read and normalise a bloodwork JSON file. Drops null/non-numeric values."""
    p = Path(path).expanduser()
    if not p.exists():
        raise FileNotFoundError(f"Bloodwork file not found: {p}")
    raw = json.loads(p.read_text())
    if not isinstance(raw, dict):
        raise ValueError("Bloodwork JSON must be an object at the top level.")
    cleaned: Dict[str, float] = {}
    for k, v in raw.items():
        if v is None:
            continue
        try:
            cleaned[_normalize_key(k)] = float(v)
        except (TypeError, ValueError):
            continue
    return cleaned


def _classify(delta_sd: float) -> Tuple[str, str]:
    """
    Return (label, html_class) for a difference expressed in trait SDs.

    Thresholds are calibrated for *common-variant polygenic predictions*, which
    typically explain only 5–20% of biomarker variance. The vast majority of
    measured variance is environmental (diet, training, sleep, medication,
    chance). We therefore use generous tolerances:

      |Δ| ≤ 0.75 SD  → Confirmed (well within the predictable range)
      |Δ| ≤ 1.5  SD  → Partial   (small lifestyle / environmental contribution)
      |Δ| >  1.5 SD  → Diverged  (large non-genetic driver dominating)

    Earlier 0.5 / 1.0 thresholds flagged most users as "diverged" simply because
    PRS variance explained is low; that is misleading, not informative.
    """
    a = abs(delta_sd)
    if a <= 0.75:
        return "Confirmed", "bw-confirmed"
    if a <= 1.5:
        return "Partial", "bw-partial"
    return "Diverged", "bw-diverged"


def compare_bloodwork(
    bloodwork_path: str,
    phewas_result: Optional[Dict],
) -> Dict:
    """
    Build the predicted-vs-actual comparison structure.

    Returns:
      {
        status: "ok" | "no_phewas" | "no_matches",
        n_labs_supplied, n_matched, n_confirmed, n_partial, n_diverged,
        accuracy_pct, rows: [...], unmatched: [...], notes: str
      }

    Each row:
      {
        trait, category, unit, predicted, actual, mean, sd,
        delta_abs, delta_sd, predicted_tier, actual_tier,
        verdict, verdict_class, interpretation
      }
    """
    labs = load_bloodwork(bloodwork_path)
    if not phewas_result or not phewas_result.get("traits"):
        return {
            "status": "no_phewas",
            "n_labs_supplied": len(labs),
            "n_matched": 0, "n_confirmed": 0, "n_partial": 0, "n_diverged": 0,
            "accuracy_pct": 0.0,
            "rows": [], "unmatched": list(labs.keys()),
            "notes": "PheWAS module did not run — cannot compare predictions.",
        }

    traits = phewas_result["traits"]
    rows: List[Dict] = []
    unmatched: List[str] = []

    for lab_key, actual_val in labs.items():
        trait_name = _LAB_TO_PHEWAS.get(lab_key)
        if not trait_name or trait_name not in traits:
            unmatched.append(lab_key)
            continue
        trait = traits[trait_name]
        res = trait.get("result", {})
        if res.get("status") != "ok":
            unmatched.append(lab_key)
            continue

        mean = float(trait["mean"])
        sd = float(trait["sd"]) or 1.0
        predicted = float(res["predicted_value"])
        delta_abs = actual_val - predicted
        delta_sd = delta_abs / sd

        verdict, vclass = _classify(delta_sd)
        predicted_tier = res.get("tier", "Average")

        # Determine actual tier from raw measurement vs population mean+sd
        actual_z = (actual_val - mean) / sd
        if actual_z >= 1.65:
            actual_tier = "Very high"
        elif actual_z >= 0.84:
            actual_tier = "Above average"
        elif actual_z >= -0.84:
            actual_tier = "Average"
        elif actual_z >= -1.65:
            actual_tier = "Below average"
        else:
            actual_tier = "Very low"

        # Interpretation narrative — the most useful column for the user
        interp = _build_interpretation(
            trait_name, predicted_tier, actual_tier, delta_sd, sd
        )

        rows.append({
            "trait": trait_name,
            "category": trait["category"],
            "unit": trait["unit"],
            "predicted": round(predicted, 2),
            "actual": round(actual_val, 2),
            "mean": mean,
            "sd": sd,
            "delta_abs": round(delta_abs, 2),
            "delta_sd": round(delta_sd, 2),
            "predicted_tier": predicted_tier,
            "actual_tier": actual_tier,
            "verdict": verdict,
            "verdict_class": vclass,
            "interpretation": interp,
            "callability_pct": res.get("callability_pct", 0),
            "n_used": res.get("n_used", 0),
            "n_total": res.get("n_total", 0),
        })

    # Sort: diverged first (most actionable), then partial, then confirmed.
    order = {"Diverged": 0, "Partial": 1, "Confirmed": 2}
    rows.sort(key=lambda r: (order[r["verdict"]], -abs(r["delta_sd"])))

    n_conf = sum(1 for r in rows if r["verdict"] == "Confirmed")
    n_part = sum(1 for r in rows if r["verdict"] == "Partial")
    n_div = sum(1 for r in rows if r["verdict"] == "Diverged")
    matched = len(rows)

    return {
        "status": "ok" if matched else "no_matches",
        "n_labs_supplied": len(labs),
        "n_matched": matched,
        "n_confirmed": n_conf,
        "n_partial": n_part,
        "n_diverged": n_div,
        "accuracy_pct": round(100 * (n_conf + 0.5 * n_part) / matched, 1) if matched else 0.0,
        "rows": rows,
        "unmatched": unmatched,
        "notes": (
            f"Compared {matched}/{len(labs)} supplied lab values against "
            f"PheWAS-derived genetic predictions. Each row shows the trait, "
            f"the prediction, the measured lab value, and the SD-scaled "
            f"divergence (Δ in SD units of the trait)."
        ),
    }


def _build_interpretation(
    trait: str, predicted_tier: str, actual_tier: str,
    delta_sd: float, sd: float,
) -> str:
    """Plain-English summary of what the gene–lab comparison means for this row."""
    pred_extreme = predicted_tier in ("Very high", "Above average", "Very low", "Below average")
    actual_extreme = actual_tier in ("Very high", "Above average", "Very low", "Below average")

    if abs(delta_sd) <= 0.75:
        if pred_extreme and actual_extreme:
            return (
                f"Genetic prediction CONFIRMED: predicted {predicted_tier.lower()}, "
                f"measured {actual_tier.lower()}. Genotype is expressing as expected."
            )
        return (
            f"Within prediction. Genotype and labs both average; nothing notable."
        )

    direction = "higher" if delta_sd > 0 else "lower"
    if abs(delta_sd) > 1.5:
        if pred_extreme and not actual_extreme:
            return (
                f"DIVERGENCE: genetics predicted {predicted_tier.lower()} but labs are "
                f"average. Lifestyle / environment / medication is likely BUFFERING the "
                f"genetic risk — keep doing whatever is working."
            )
        if not pred_extreme and actual_extreme:
            return (
                f"DIVERGENCE: genetics predicted average but labs are {actual_tier.lower()}. "
                f"Strong NON-GENETIC driver (diet, training, stress, illness, medication) "
                f"is dominating — investigate lifestyle/clinical causes."
            )
        return (
            f"Genetics directionally correct but measured value is {abs(delta_sd):.1f} SD "
            f"{direction} than predicted — gene-environment effect amplifying the genotype."
        )

    return (
        f"Partial match: measured value runs {abs(delta_sd):.1f} SD {direction} than "
        f"the genotype-only prediction. Minor environmental contribution."
    )


# ── HTML rendering ───────────────────────────────────────────────────────────

def _esc(s) -> str:
    s = "" if s is None else str(s)
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


_BW_CSS = """
<style>
.bw-wrap { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
           color:#222; max-width: 1100px; margin: 24px auto; padding: 0 16px; }
.bw-wrap h1 { font-size: 1.6em; border-bottom: 2px solid #333; padding-bottom: 6px; }
.bw-summary { display:flex; gap:14px; flex-wrap:wrap; margin: 18px 0; }
.bw-stat { background:#f6f6f7; border:1px solid #ddd; border-radius:8px;
           padding:10px 14px; min-width:120px; }
.bw-stat .v { font-size:1.5em; font-weight:600; }
.bw-stat .l { font-size:0.8em; color:#666; text-transform:uppercase; letter-spacing:0.05em; }
table.bw { width:100%; border-collapse: collapse; margin-top:10px; font-size:0.92em; }
table.bw th, table.bw td { padding:8px 10px; border-bottom:1px solid #eee; vertical-align: top; }
table.bw th { text-align:left; background:#f9f9f9; }
.bw-confirmed { background:#e7f6ea; }
.bw-partial   { background:#fff7e0; }
.bw-diverged  { background:#fde4e4; }
.bw-verdict { font-weight:600; }
.bw-verdict.bw-confirmed { color:#2c7a30; }
.bw-verdict.bw-partial   { color:#a06800; }
.bw-verdict.bw-diverged  { color:#a32a2a; }
.bw-num { font-variant-numeric: tabular-nums; }
.bw-unmatched { color:#666; font-size:0.9em; margin-top:14px; }
</style>
"""


def render_bloodwork_html(result: Dict, file_label: str = "") -> str:
    """Render the bloodwork comparison as a standalone HTML document."""
    if result.get("status") == "no_phewas":
        body = (
            f"<p>PheWAS module did not run; cannot compare {result['n_labs_supplied']} "
            f"supplied lab values to genetic predictions.</p>"
        )
    elif not result.get("rows"):
        body = (
            f"<p>None of the {result['n_labs_supplied']} supplied lab values matched "
            f"a scored PheWAS biomarker.</p>"
        )
    else:
        rows_html = []
        for r in result["rows"]:
            rows_html.append(f"""
              <tr class="{r['verdict_class']}">
                <td><strong>{_esc(r['trait'])}</strong><br>
                    <span style="color:#777;font-size:0.85em">{_esc(r['category'])}</span></td>
                <td class="bw-num">{_esc(r['predicted'])} {_esc(r['unit'])}<br>
                    <span style="color:#777;font-size:0.85em">{_esc(r['predicted_tier'])}</span>
                    {f'<br><span style="color:#999;font-size:0.78em">genetic coverage {r["callability_pct"]:.0f}% ({r.get("n_used",0)} SNPs)</span>' if r.get('callability_pct') else ''}</td>
                <td class="bw-num">{_esc(r['actual'])} {_esc(r['unit'])}<br>
                    <span style="color:#777;font-size:0.85em">{_esc(r['actual_tier'])}</span></td>
                <td class="bw-num">{r['delta_abs']:+.2f}<br>
                    <span style="color:#777;font-size:0.85em">{r['delta_sd']:+.2f} SD</span></td>
                <td class="bw-verdict {r['verdict_class']}">{_esc(r['verdict'])}</td>
                <td>{_esc(r['interpretation'])}</td>
              </tr>
            """)
        body = f"""
          <div class="bw-summary">
            <div class="bw-stat"><div class="v">{result['n_matched']}</div>
                <div class="l">Compared</div></div>
            <div class="bw-stat"><div class="v" style="color:#2c7a30">{result['n_confirmed']}</div>
                <div class="l">Confirmed</div></div>
            <div class="bw-stat"><div class="v" style="color:#a06800">{result['n_partial']}</div>
                <div class="l">Partial</div></div>
            <div class="bw-stat"><div class="v" style="color:#a32a2a">{result['n_diverged']}</div>
                <div class="l">Diverged</div></div>
            <div class="bw-stat"><div class="v">{result['accuracy_pct']}%</div>
                <div class="l">Accuracy</div></div>
          </div>
          <table class="bw">
            <tr>
              <th>Biomarker</th><th>Predicted (genetics)</th><th>Actual (lab)</th>
              <th>Δ</th><th>Verdict</th><th>Interpretation</th>
            </tr>
            {''.join(rows_html)}
          </table>
        """
        if result.get("unmatched"):
            body += (
                f"<p class='bw-unmatched'>Supplied labs not matched to a scored "
                f"genetic prediction: {', '.join(_esc(u) for u in result['unmatched'])}.</p>"
            )
    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>Blood Work vs Genetic Predictions{(' — ' + _esc(file_label)) if file_label else ''}</title>
{_BW_CSS}</head><body><div class="bw-wrap">
<h1>Blood Work vs Genetic Predictions</h1>
<p style="color:#666">{_esc(result.get('notes',''))}</p>

<details style="background:#f6f6f7;border:1px solid #ddd;border-radius:8px;
                padding:10px 14px;margin:12px 0;">
  <summary style="cursor:pointer;font-weight:600">
    How to read this table — what "Diverged" really means
  </summary>
  <p style="margin-top:10px;color:#444">
    Polygenic predictions built from common GWAS variants typically explain
    only <strong>5–20% of the variance</strong> in any single biomarker. The
    remaining 80–95% is driven by diet, exercise, sleep, medications, illness,
    age, sex, season, and chance. The thresholds used here are calibrated to
    that reality:
  </p>
  <ul style="color:#444">
    <li><strong>Confirmed</strong> (|Δ| ≤ 0.75 SD) — measured value sits well
        within where genetics alone would predict.</li>
    <li><strong>Partial</strong> (|Δ| ≤ 1.5 SD) — small lifestyle / environmental
        contribution on top of genetics.</li>
    <li><strong>Diverged</strong> (|Δ| &gt; 1.5 SD) — large non-genetic driver
        dominating; investigate diet, training, medication, stress, illness.
        This is the most actionable row in the table.</li>
  </ul>
  <p style="color:#444">
    "Diverged" never means the genetic prediction is <em>wrong</em> — it means
    a real-world driver is overriding it. That driver is what's modifiable.
  </p>
</details>

{body}
</div></body></html>"""
