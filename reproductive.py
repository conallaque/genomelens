"""
Reproductive Genetics Simulator
===============================

For each carrier finding, estimate:
  * Probability that a partner of a given ancestry also carries the variant
  * Probability of having an affected child given two carriers (25%
    for autosomal recessive, 50% for autosomal dominant, etc.)
  * Specific variants the partner should test for before conception

Also flags:
  * Dominant pathogenic variants (50% inheritance per child)
  * X-linked variants (sons of carrier mothers — 50% affected)
  * Elevated ROH context (consanguinity risk for rare recessives)

Carrier-frequency strings come from carrier.py entries. We parse them
loosely into per-population probabilities; where no data is given we fall
back to a generic "1 in 100" assumption.
"""

from __future__ import annotations

import re


def _parse_frequency_string(freq_str: str) -> dict[str, float]:
    """Parse a free-form 'European ~1:25, Asian ~1:90' string into a dict."""
    if not freq_str:
        return {}
    out: dict[str, float] = {}
    # Find each population: freq pair
    parts = re.split(r"[,;]", freq_str)
    for p in parts:
        p = p.strip()
        # Try '1:X' or 'X%'
        m_ratio = re.search(r"([A-Za-z][A-Za-z\s]*?)\s*[~≈]?\s*1\s*[:/]\s*(\d+)", p)
        if m_ratio:
            pop = m_ratio.group(1).strip()
            denom = int(m_ratio.group(2))
            if denom > 0:
                out[pop] = 1.0 / denom
            continue
        m_pct = re.search(r"([A-Za-z][A-Za-z\s]*?)\s*[~≈]?\s*(\d+(?:\.\d+)?)\s*%", p)
        if m_pct:
            pop = m_pct.group(1).strip()
            pct = float(m_pct.group(2))
            out[pop] = pct / 100
    return out


def _inheritance_class(entry: dict) -> str:
    inh = entry.get("inheritance", "").lower()
    if "x-linked" in inh:
        return "x-linked"
    if "dominant" in inh:
        return "dominant"
    if "recessive" in inh:
        return "recessive"
    return "unknown"


def analyze_reproductive(carrier_result: dict | None,
                          roh_result: dict | None = None) -> dict:
    """Generate reproductive scenarios for each carrier/affected finding."""
    if not carrier_result:
        return {"available": False, "scenarios": []}

    scenarios: list[dict] = []
    all_findings = (carrier_result.get("carriers", []) +
                    carrier_result.get("affected", []))

    for f in all_findings:
        klass = _inheritance_class(f)
        carrier_freq_text = f.get("carrier_frequency", "")
        freqs = _parse_frequency_string(carrier_freq_text)
        # If we couldn't parse but text has data, present text
        if not freqs and carrier_freq_text:
            freqs = {"General population": 0.01}  # generic 1 in 100 fallback

        dosage = f.get("dosage")

        # Build scenario
        scenario = {
            "gene": f["gene"],
            "variant": f["variant"],
            "disease": f["disease"],
            "inheritance": f["inheritance"],
            "inheritance_class": klass,
            "user_dosage": dosage,
            "carrier_freq_text": carrier_freq_text,
            "by_population": [],
            "partner_testing_advice": f.get("partner_testing")
                                       or "Discuss with a genetic counsellor.",
            "treatment_outlook": f.get("treatment_outlook", ""),
            "chip_caveat": f.get("chip_caveat", ""),
        }

        if klass == "recessive":
            # User as carrier (dose 1) — partner-carrier × ¼ chance of affected child
            for pop, pcf in freqs.items():
                # Probability partner is a carrier
                p_partner_carrier = pcf
                # Probability affected child given both carriers = 1/4
                p_affected_child = p_partner_carrier * 0.25
                p_carrier_child = p_partner_carrier * 0.50 + (1 - p_partner_carrier) * 0.50
                scenario["by_population"].append({
                    "population": pop,
                    "partner_carrier_freq": p_partner_carrier,
                    "p_affected_child": p_affected_child,
                    "p_carrier_child": p_carrier_child,
                    "p_affected_str": _odds_str(p_affected_child),
                    "p_carrier_str": _odds_str(p_carrier_child),
                })
            # If user is HOMOZYGOUS (affected), every child gets ≥1 copy
            if dosage == 2:
                for sc in scenario["by_population"]:
                    sc["p_affected_child"] = sc["partner_carrier_freq"] * 0.5
                    sc["p_carrier_child"] = 1.0  # every child carries one
                    sc["p_affected_str"] = _odds_str(sc["p_affected_child"])
                    sc["p_carrier_str"] = "Every child"

        elif klass == "dominant":
            scenario["dominant_note"] = (
                "Each child has a 50% probability of inheriting this dominant "
                "variant. Penetrance — actual likelihood of developing disease — "
                "varies by condition and is often less than 100%."
            )

        elif klass == "x-linked":
            scenario["xlinked_note"] = (
                "X-linked: each son of a carrier mother has 50% chance of being "
                "affected (no functional copy); each daughter has 50% chance of "
                "being a carrier. If the user is male and carries, all daughters "
                "will be obligate carriers; sons unaffected (only one X comes from mum)."
            )

        scenarios.append(scenario)

    # ROH-based consanguinity context
    roh_context = ""
    if roh_result and roh_result.get("f_roh") is not None:
        f_roh = roh_result["f_roh"]
        if f_roh > 0.04:
            roh_context = (
                f"Elevated F_ROH ({f_roh:.4f}) suggests recent parental relatedness. "
                "Offspring have elevated risk of rare-recessive conditions. "
                "Expanded carrier screening (sequencing-based panel covering 250+ "
                "recessive conditions) is particularly valuable before conception."
            )
        elif f_roh > 0.015:
            roh_context = (
                f"F_ROH ({f_roh:.4f}) suggests founder-population background or "
                "distant consanguinity. Expanded carrier screening is reasonable "
                "if planning children with a partner from the same background."
            )

    return {
        "available": True,
        "scenarios": scenarios,
        "n_scenarios": len(scenarios),
        "roh_context": roh_context,
    }


def _odds_str(p: float) -> str:
    """Format probability as '1 in N' string."""
    if p <= 0:
        return "0"
    if p >= 1:
        return "Every child"
    n = round(1.0 / p)
    return f"1 in {n:,}"
