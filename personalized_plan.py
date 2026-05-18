"""
Personalised Plan — Master Dashboard
====================================

Single-page synthesis of the V6 personalisation outputs:

  supplements + exercise + nutrition + bloodwork + headline PheWAS predictions

The page is meant to be the "executive briefing" version of the report:
one screen that tells the user what to do tomorrow morning, with the
underlying modules linked for the detail.

This module is deliberately a pure renderer — it does not re-derive any
genotype calls. Every input is the structured dict produced by the
respective analysis module.

Output keys in the returned dict (also rendered as HTML):
  headline:  most actionable single insight per pillar
  pillars:   {supplements, exercise, nutrition, biomarkers, bloodwork}
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional


# ── Source-of-truth: headline picker per pillar ─────────────────────────────

def _supp_headline(supp: Optional[Dict]) -> str:
    if not supp or supp.get("status") != "ok":
        return "—"
    essential = supp["tiers"].get("essential", [])
    if essential:
        names = [r["name"] for r in essential[:3]]
        return f"{len(essential)} essential supplement{'s' if len(essential) != 1 else ''}: " + ", ".join(names)
    recommended = supp["tiers"].get("recommended", [])
    if recommended:
        return f"{len(recommended)} recommended supplement{'s' if len(recommended) != 1 else ''}"
    return "No genotype-driven essentials"


def _ex_headline(ex: Optional[Dict]) -> str:
    if not ex or ex.get("status") != "ok":
        return "—"
    pe = ex["power_endurance"]
    recov = ex["recovery"]["speed"].lower()
    chrono = ex["chronotype"]["chronotype"].lower()
    window = ex["chronotype"]["optimal_window"]
    return (
        f"{pe['bias']} ({pe['ratio_pct_power']}/{pe['ratio_pct_endurance']}) · "
        f"{recov} recovery · train {window}"
    )


def _nu_headline(nu: Optional[Dict]) -> str:
    if not nu or nu.get("status") != "ok":
        return "—"
    m = nu["macros"]
    caf = nu["caffeine"]["metabolism"]
    alc = nu["alcohol"]["risk"]
    return f"{m['pct_carbs']}C/{m['pct_fat']}F/{m['pct_protein']}P · caffeine {caf.lower()} · alcohol {alc.lower()}"


def _bw_headline(bw: Optional[Dict]) -> str:
    if not bw or bw.get("status") != "ok":
        return "Bloodwork not supplied"
    return (
        f"{bw['n_matched']} labs compared · "
        f"{bw['n_confirmed']} confirmed · {bw['n_partial']} partial · "
        f"{bw['n_diverged']} diverged · {bw['accuracy_pct']}% accuracy"
    )


def _phewas_headlines(phewas: Optional[Dict]) -> List[Dict]:
    """Top 5 most extreme biomarker predictions (very-high / very-low)."""
    if not phewas:
        return []
    hl = phewas.get("headline", [])
    out: List[Dict] = []
    for h in hl[:5]:
        out.append({
            "trait": h.get("trait", ""),
            "tier": h.get("tier", ""),
            "percentile": h.get("percentile"),
            "predicted_value": h.get("predicted_value"),
            "unit": h.get("unit", ""),
        })
    return out


# ── Cross-module synthesis: "morning routine" extraction ────────────────────

def _morning_actions(supp: Optional[Dict], nu: Optional[Dict],
                     ex: Optional[Dict]) -> List[str]:
    """Pull the morning-timed actions from all three modules into one list."""
    actions: List[str] = []
    if supp and supp.get("status") == "ok":
        for tier_key in ("essential", "recommended"):
            for r in supp["tiers"].get(tier_key, []):
                t = (r.get("timing", "") or "").lower()
                if "morning" in t or "breakfast" in t or "am" in t:
                    actions.append(
                        f"{r['name']} — {r['dose']} ({tier_key})"
                    )
    if nu and nu.get("status") == "ok":
        cf = nu["caffeine"]
        if cf.get("limit_mg"):
            actions.append(
                f"Coffee cap {cf['limit_mg']} mg/day; last cup ≤ {cf['cutoff_time']}"
            )
        for meal in nu.get("daily_template", []):
            if meal["meal"] == "Breakfast":
                actions.append(f"Breakfast: {meal['example']}")
                break
    if ex and ex.get("status") == "ok":
        chrono = ex["chronotype"]
        if chrono["chronotype"] == "Morning":
            actions.append(f"Train in morning window {chrono['optimal_window']}")
    return actions


def _evening_actions(supp: Optional[Dict], nu: Optional[Dict],
                     ex: Optional[Dict]) -> List[str]:
    actions: List[str] = []
    if supp and supp.get("status") == "ok":
        for tier_key in ("essential", "recommended"):
            for r in supp["tiers"].get(tier_key, []):
                t = (r.get("timing", "") or "").lower()
                if "evening" in t or "night" in t or "dinner" in t or "pm" in t:
                    actions.append(f"{r['name']} — {r['dose']}")
    if ex and ex.get("status") == "ok":
        chrono = ex["chronotype"]
        if chrono["chronotype"].startswith("Evening") or chrono["chronotype"] == "Slight evening":
            actions.append(f"Train in evening window {chrono['optimal_window']}")
    return actions


# ── Bloodwork ↔ supplement reconciliation ────────────────────────────────────

def _reconciliation(supp: Optional[Dict], bw: Optional[Dict]) -> List[Dict]:
    """
    For each Diverged bloodwork row, find supplement recommendations that
    address the same biomarker and surface them as a "highest-leverage" list.
    """
    if not bw or bw.get("status") != "ok":
        return []
    out: List[Dict] = []
    diverged = [r for r in bw.get("rows", []) if r["verdict"] == "Diverged"]
    if not diverged or not supp:
        return out

    # Map biomarker → relevant supplement substring(s)
    targets = {
        "C-Reactive Protein": ["Curcumin", "EPA/DHA"],
        "25-OH Vitamin D":    ["Vitamin D3"],
        "LDL cholesterol":    ["EPA/DHA"],
        "Triglycerides":      ["EPA/DHA"],
        "HDL cholesterol":    ["EPA/DHA"],
        "Iron / ferritin":    ["AVOID iron", "Gentle Iron"],
        "Serum vitamin B12":  ["Methylcobalamin"],
        "Folate":             ["L-Methylfolate"],
    }

    flat: List[Dict] = []
    for tier_key in ("essential", "recommended", "optional", "avoid"):
        for r in (supp.get("tiers", {}).get(tier_key) or []):
            flat.append(r)

    for row in diverged:
        keys = targets.get(row["trait"], [])
        matching = [r for r in flat if any(k in r["name"] for k in keys)]
        out.append({
            "biomarker": row["trait"],
            "predicted": row["predicted"],
            "actual":    row["actual"],
            "unit":      row["unit"],
            "delta_sd":  row["delta_sd"],
            "interpretation": row["interpretation"],
            "supplements_in_play": [
                {"name": r["name"], "tier": r["tier"], "dose": r["dose"]}
                for r in matching
            ],
        })
    return out


# ── Public API ──────────────────────────────────────────────────────────────

def build_personalized_plan(
    supplement_result: Optional[Dict] = None,
    exercise_result: Optional[Dict] = None,
    nutrition_result: Optional[Dict] = None,
    bloodwork_result: Optional[Dict] = None,
    phewas_result: Optional[Dict] = None,
) -> Dict:
    headlines = {
        "supplements": _supp_headline(supplement_result),
        "exercise":    _ex_headline(exercise_result),
        "nutrition":   _nu_headline(nutrition_result),
        "bloodwork":   _bw_headline(bloodwork_result),
    }
    pillars = {
        "supplements": supplement_result,
        "exercise":    exercise_result,
        "nutrition":   nutrition_result,
        "bloodwork":   bloodwork_result,
        "biomarkers":  _phewas_headlines(phewas_result),
    }
    morning = _morning_actions(supplement_result, nutrition_result, exercise_result)
    evening = _evening_actions(supplement_result, nutrition_result, exercise_result)
    reconciliation = _reconciliation(supplement_result, bloodwork_result)

    return {
        "status": "ok",
        "headlines": headlines,
        "pillars": pillars,
        "morning_actions": morning,
        "evening_actions": evening,
        "reconciliation": reconciliation,
    }


# ── HTML rendering ──────────────────────────────────────────────────────────

def _esc(s) -> str:
    s = "" if s is None else str(s)
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


_CSS = """
<style>
  :root {
    --bg: #f7f7f9; --card: #ffffff; --border: #e2e2e6; --ink: #1c1c1f;
    --ink-soft: #5a5a62; --accent: #1f4f93; --good: #2c7a30; --warn: #b48a00;
    --bad: #a32a2a;
  }
  body { background: var(--bg); margin: 0; }
  .pp-wrap { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
             color: var(--ink); max-width: 1180px; margin: 0 auto;
             padding: 28px 20px; }
  .pp-wrap h1 { font-size: 1.9em; margin: 0 0 8px 0; letter-spacing: -0.01em; }
  .pp-wrap h2 { font-size: 1.15em; margin: 28px 0 10px 0;
                color: var(--accent); border-bottom: 1px solid var(--border);
                padding-bottom: 4px; }
  .pp-sub { color: var(--ink-soft); margin: 0 0 24px 0; }

  /* Pillar overview grid */
  .pp-grid { display: grid; grid-template-columns: repeat(4, 1fr);
             gap: 14px; margin-bottom: 28px; }
  @media (max-width: 900px) { .pp-grid { grid-template-columns: 1fr 1fr; } }
  @media (max-width: 540px) { .pp-grid { grid-template-columns: 1fr; } }
  .pp-pillar { background: var(--card); border: 1px solid var(--border);
               border-radius: 12px; padding: 16px;
               box-shadow: 0 1px 2px rgba(0,0,0,0.04); }
  .pp-pillar .icon { font-size: 1.4em; margin-bottom: 6px; }
  .pp-pillar .name { font-weight: 600; color: var(--ink); margin-bottom: 4px; }
  .pp-pillar .body { font-size: 0.92em; color: var(--ink-soft);
                     line-height: 1.4; }
  .pp-pillar a { color: var(--accent); text-decoration: none; font-weight: 500; }
  .pp-pillar a:hover { text-decoration: underline; }

  /* Routine columns */
  .pp-routine { display: grid; grid-template-columns: 1fr 1fr; gap: 18px; }
  @media (max-width: 700px) { .pp-routine { grid-template-columns: 1fr; } }
  .pp-routine .col { background: var(--card); border: 1px solid var(--border);
                     border-radius: 12px; padding: 16px; }
  .pp-routine ul { margin: 0; padding-left: 18px; }
  .pp-routine li { padding: 3px 0; color: var(--ink); }
  .pp-routine h3 { margin: 0 0 8px 0; font-size: 1.05em; }

  /* Reconciliation table */
  .pp-recon { background: var(--card); border: 1px solid var(--border);
              border-radius: 12px; padding: 14px 16px; }
  .pp-recon table { width: 100%; border-collapse: collapse; }
  .pp-recon th, .pp-recon td { padding: 8px 10px; border-bottom: 1px solid #f0f0f3;
                               text-align: left; vertical-align: top; }
  .pp-recon th { background: #fafafc; font-weight: 600; font-size: 0.9em; }
  .pp-recon .delta { color: var(--bad); font-weight: 600; }
  .pp-recon .supp-pill { display: inline-block; background: #eef2f8;
                         color: var(--accent); padding: 2px 8px; border-radius: 8px;
                         font-size: 0.82em; margin: 2px 4px 0 0; }

  /* Biomarker headlines */
  .pp-bio { display: grid; grid-template-columns: repeat(auto-fill, minmax(220px, 1fr));
            gap: 10px; }
  .pp-bio .item { background: var(--card); border: 1px solid var(--border);
                  border-radius: 10px; padding: 10px 12px; font-size: 0.9em; }
  .pp-bio .trait { font-weight: 600; color: var(--ink); margin-bottom: 2px; }
  .pp-bio .tier { font-size: 0.8em; }
  .pp-bio .tier.high { color: var(--bad); }
  .pp-bio .tier.low  { color: var(--accent); }
  .pp-bio .val { color: var(--ink-soft); font-variant-numeric: tabular-nums; }

  .pp-footer { color: var(--ink-soft); font-size: 0.82em; margin-top: 34px;
               border-top: 1px solid var(--border); padding-top: 16px; }
</style>
"""


_PILLAR_META = [
    ("supplements", "💊", "Supplement stack",  "supplements.html"),
    ("exercise",    "🏋",  "Exercise programming", "exercise.html"),
    ("nutrition",   "🥗", "Nutrition plan",     "nutrition.html"),
    ("bloodwork",   "🩸", "Blood-work check",   "bloodwork.html"),
]


def render_plan_html(plan: Dict, file_label: str = "",
                     report_link: str = "report.html") -> str:
    if not plan or plan.get("status") != "ok":
        return f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Personalised Plan</title>{_CSS}</head><body><div class="pp-wrap">
<h1>Personalised Plan</h1><p>Insufficient data to build the dashboard.</p>
</div></body></html>"""

    # Pillar cards
    headlines = plan["headlines"]
    pillar_cards: List[str] = []
    pillars = plan["pillars"]
    for key, icon, name, link in _PILLAR_META:
        has_data = pillars.get(key) is not None and (
            pillars[key].get("status") == "ok"
            if isinstance(pillars[key], dict) else False
        )
        body_text = headlines.get(key, "—")
        link_html = (f'<a href="{link}">Open detailed page →</a>'
                     if has_data else
                     '<span style="color:#888">(not generated)</span>')
        pillar_cards.append(f"""
<div class="pp-pillar">
  <div class="icon">{icon}</div>
  <div class="name">{_esc(name)}</div>
  <div class="body">{_esc(body_text)}</div>
  <div style="margin-top:10px;font-size:0.85em">{link_html}</div>
</div>""")

    # Morning / evening routine
    morning_items = (
        "".join(f"<li>{_esc(a)}</li>" for a in plan["morning_actions"])
        or "<li><em>No morning-timed actions identified.</em></li>"
    )
    evening_items = (
        "".join(f"<li>{_esc(a)}</li>" for a in plan["evening_actions"])
        or "<li><em>No evening-timed actions identified.</em></li>"
    )

    # Reconciliation
    if plan["reconciliation"]:
        recon_rows = []
        for row in plan["reconciliation"]:
            pills = "".join(
                f'<span class="supp-pill">{_esc(s["name"])} · {_esc(s["dose"])}</span>'
                for s in row["supplements_in_play"]
            ) or "<span style='color:#888'>No supplement directly addresses this</span>"
            recon_rows.append(f"""
<tr>
  <td><strong>{_esc(row['biomarker'])}</strong></td>
  <td class="delta">{row['delta_sd']:+.2f} SD<br>
      <span style="color:#5a5a62;font-weight:400">
        predicted {row['predicted']} {_esc(row['unit'])}, actual {row['actual']} {_esc(row['unit'])}
      </span></td>
  <td>{_esc(row['interpretation'])}<br>{pills}</td>
</tr>""")
        recon_block = f"""
<h2>Where genetics and labs disagree — highest-leverage fixes</h2>
<div class="pp-recon">
<table>
  <tr><th>Biomarker</th><th>Δ from prediction</th><th>What's happening + relevant supplements</th></tr>
  {''.join(recon_rows)}
</table>
</div>"""
    else:
        recon_block = ""

    # Biomarker headlines
    bio_items = plan["pillars"].get("biomarkers", [])
    if bio_items:
        bio_html = []
        for b in bio_items:
            tier_class = "high" if "high" in b["tier"].lower() else (
                          "low" if "low" in b["tier"].lower() else "")
            val_str = (f"{b['predicted_value']} {b['unit']}"
                       if b.get("predicted_value") is not None else b["unit"])
            pct_str = (f"{b['percentile']:.0f}th pct"
                       if b.get("percentile") is not None else "")
            bio_html.append(f"""
<div class="item">
  <div class="trait">{_esc(b['trait'])}</div>
  <div class="tier {tier_class}">{_esc(b['tier'])} · {_esc(pct_str)}</div>
  <div class="val">{_esc(val_str)}</div>
</div>""")
        bio_block = f"""
<h2>Headline biomarker predictions (PheWAS)</h2>
<div class="pp-bio">{''.join(bio_html)}</div>"""
    else:
        bio_block = ""

    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>Personalised Plan{(' — ' + _esc(file_label)) if file_label else ''}</title>
{_CSS}</head><body><div class="pp-wrap">

<h1>Your Personalised Plan</h1>
<p class="pp-sub">
  A one-screen synthesis of your supplement stack, exercise programming,
  nutrition plan, and (when supplied) blood-work comparison. Each pillar
  links to a full detail page.
</p>

<h2>Pillars at a glance</h2>
<div class="pp-grid">
  {''.join(pillar_cards)}
</div>

<h2>Daily routine extracted from your genotype</h2>
<div class="pp-routine">
  <div class="col">
    <h3>🌅 Morning</h3>
    <ul>{morning_items}</ul>
  </div>
  <div class="col">
    <h3>🌙 Evening</h3>
    <ul>{evening_items}</ul>
  </div>
</div>

{recon_block}

{bio_block}

<p class="pp-footer">
  Not medical advice. This dashboard summarises the V6 personalisation
  modules of the local DNA analysis pipeline; underlying SNPs, doses, and
  caveats are documented on the linked detail pages and in
  <a href="{_esc(report_link)}">the full report</a>. Refine with measured
  labs and clinical judgement.
</p>
</div></body></html>"""
