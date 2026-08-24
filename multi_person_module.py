"""
Multi-Person Comparison — two genomes side by side
==================================================

Compares two genotype files: genotype concordance, a **properly-computed
relationship estimate** (KING-robust kinship + IBS0), and a couple carrier-risk
summary for family planning.

Relationship inference — done properly, not naively
---------------------------------------------------
A naive "percent of identical genotype calls" will call two unrelated
Northern Europeans near-siblings, because they share a high fraction of common-
allele genotypes by chance. We therefore use the **KING-robust kinship
estimator** (Manichaikul et al., Bioinformatics 2010), which is robust to
population structure and needs no allele frequencies:

        φ = ( N_HetHet − 2·N_IBS0 ) / ( N_Het_i + N_Het_j )

where over the shared, biallelic, autosomal SNPs:
  * N_HetHet = loci where *both* individuals are heterozygous,
  * N_IBS0   = loci where they are *opposite homozygotes* (AA vs GG),
  * N_Het_i  = loci where individual i is heterozygous.

Standard KING kinship thresholds:
  * φ > 0.354           → duplicate sample / monozygotic twin
  * 0.177 < φ ≤ 0.354   → 1st-degree (parent-child or full siblings)
  * 0.0884 < φ ≤ 0.177  → 2nd-degree (grandparent, aunt/uncle, half-sib)
  * 0.0442 < φ ≤ 0.0884 → 3rd-degree (first cousin)
  * φ ≤ 0.0442          → unrelated / more distant

Parent-child vs full-sibling within the 1st-degree band is refined using
**IBS0**: parent-child pairs share ≥1 allele at every locus, so their IBS0
(opposite-homozygote) rate is ≈0; full siblings have a small but non-zero IBS0.

This is a research-grade estimate, not a legal relationship test — consumer-chip
overlap, genotyping error, and low shared-SNP counts across different vendors
all reduce confidence, which is reported.

Privacy
-------
This module never writes a second person's genotypes to disk. Comparison output
is returned to the caller and rendered to a STANDALONE page — it is deliberately
not threaded into the single-person report, because a second person's results do
not belong in the first person's document. An earlier version of this paragraph
claimed the opposite ("rendered into the main report"), contradicting the note
above ``build_comparison_html`` fifty lines below and describing wiring that has
never existed. Raw genome files stay wherever the user put them; the project
.gitignore excludes *.csv/*.txt/*.vcf/*.zip so no genome is ever committed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

_VALID = set("ACGT")


def _norm_gt(g) -> str | None:
    """Return a normalised 2-allele biallelic genotype as a sorted string, or
    None for no-calls / indels / non-ACGT."""
    if g is None:
        return None
    s = str(g).upper().replace(" ", "")
    if len(s) != 2:
        return None
    if s[0] not in _VALID or s[1] not in _VALID:
        return None
    return "".join(sorted(s))


def _is_het(gt: str) -> bool:
    return gt[0] != gt[1]


def _prep(df: pd.DataFrame) -> pd.Series:
    """Return a genotype Series indexed by rsid: deduplicated (first wins),
    autosomal-only when a chrom column is present, upper-cased/whitespace-
    stripped. Vectorised — safe for 600k+ markers."""
    d = df
    if not d.index.is_unique:
        d = d[~d.index.duplicated(keep="first")]
    if "chrom" in d.columns:
        chrom = (d["chrom"].astype(str).str.upper()
                 .str.replace("CHR", "", regex=False))
        d = d[~chrom.isin({"X", "Y", "MT", "M", "0", "XY"})]
    return d["genotype"].astype(str).str.upper().str.replace(" ", "", regex=False)


def king_kinship(df_a: pd.DataFrame, df_b: pd.DataFrame) -> dict:
    """Compute KING-robust kinship + IBS0 over shared autosomal biallelic SNPs.

    Fully vectorised: deduplicates non-unique indices, joins on shared rsIDs,
    and counts hetero/homozygote states with boolean masks — scales to
    whole-chip (600k+) comparisons in well under a second."""
    ga = _prep(df_a)
    gb = _prep(df_b)
    shared = ga.index.intersection(gb.index)
    ga = ga.loc[shared]
    gb = gb.loc[shared]

    # Keep only biallelic ACGT 2-char genotypes in BOTH.
    valid = ga.str.fullmatch(r"[ACGT]{2}") & gb.str.fullmatch(r"[ACGT]{2}")
    ga = ga[valid]
    gb = gb[valid]
    n_shared = len(ga)
    if n_shared == 0:
        return {"n_shared_snps": 0, "n_hethet": 0, "n_ibs0": 0,
                "ibs0_rate": None, "n_het_a": 0, "n_het_b": 0,
                "kinship": None, "concordance": None}

    a0 = ga.str[0].to_numpy()
    a1 = ga.str[1].to_numpy()
    b0 = gb.str[0].to_numpy()
    b1 = gb.str[1].to_numpy()

    het_a = a0 != a1
    het_b = b0 != b1
    hom_a = ~het_a
    hom_b = ~het_b

    n_het_a = int(het_a.sum())
    n_het_b = int(het_b.sum())
    n_hethet = int((het_a & het_b).sum())
    # opposite homozygotes → IBS0 (both homozygous, different allele)
    n_ibs0 = int((hom_a & hom_b & (a0 != b0)).sum())

    # concordance on strand-normalised genotypes ("AG" == "GA")
    norm_a = np.where(a0 <= a1, np.char.add(a0.astype(str), a1.astype(str)),
                      np.char.add(a1.astype(str), a0.astype(str)))
    norm_b = np.where(b0 <= b1, np.char.add(b0.astype(str), b1.astype(str)),
                      np.char.add(b1.astype(str), b0.astype(str)))
    n_concordant = int((norm_a == norm_b).sum())

    denom = n_het_a + n_het_b
    kinship = ((n_hethet - 2 * n_ibs0) / denom) if denom > 0 else None
    ibs0_rate = (n_ibs0 / n_shared) if n_shared else None
    concordance = (n_concordant / n_shared) if n_shared else None

    return {
        "n_shared_snps": n_shared,
        "n_hethet": n_hethet,
        "n_ibs0": n_ibs0,
        "ibs0_rate": ibs0_rate,
        "n_het_a": n_het_a,
        "n_het_b": n_het_b,
        "kinship": kinship,
        "concordance": concordance,
    }


def classify_relationship(kinship: float | None, ibs0_rate: float | None,
                          n_shared: int) -> dict:
    """Map KING kinship + IBS0 to a relationship-degree estimate with a
    confidence caveat."""
    if kinship is None or n_shared < 200:
        return {
            "degree": "Indeterminate",
            "label": "Not enough overlapping markers to estimate relationship",
            "confidence": "none",
            "detail": (f"Only {n_shared} shared biallelic autosomal SNPs — too "
                       "few for a reliable KING estimate (usually happens across "
                       "different testing vendors). Concordance is still shown."),
        }

    if kinship > 0.354:
        degree, label = "Duplicate / MZ twin", "Identical samples or monozygotic twins"
    elif kinship > 0.177:
        # 1st-degree — refine PO vs FS using IBS0
        if ibs0_rate is not None and ibs0_rate < 0.0025:
            label = ("1st-degree — most consistent with **parent-child** "
                     "(they share an allele at essentially every locus; IBS0≈0)")
        else:
            label = ("1st-degree — most consistent with **full siblings** "
                     "(non-zero opposite-homozygote rate)")
        degree = "1st-degree"
    elif kinship > 0.0884:
        degree, label = "2nd-degree", ("2nd-degree — grandparent/grandchild, "
                                       "aunt/uncle-niece/nephew, half-sibling, or "
                                       "double first cousin")
    elif kinship > 0.0442:
        degree, label = "3rd-degree", "3rd-degree — most consistent with first cousins"
    else:
        degree, label = "Unrelated", ("No detectable close relationship "
                                      "(unrelated or more distant than 3rd-degree)")

    conf = "high" if n_shared >= 5000 else ("moderate" if n_shared >= 1000 else "low")
    return {
        "degree": degree, "label": label, "confidence": conf,
        "detail": (f"KING-robust kinship φ = {kinship:.4f} over {n_shared:,} shared "
                   f"autosomal SNPs (IBS0 rate {ibs0_rate:.4f})." if ibs0_rate is not None
                   else f"KING-robust kinship φ = {kinship:.4f} over {n_shared:,} SNPs."),
    }


def couple_carrier_risk(df_a: pd.DataFrame, df_b: pd.DataFrame) -> dict | None:
    """If both individuals are carriers of the same recessive condition, a child
    has a 25% risk of being affected. Uses the project carrier module if present."""
    try:
        from carrier import analyze_carriers
    except Exception:
        return None
    try:
        ca = analyze_carriers(df_a)
        cb = analyze_carriers(df_b)
    except Exception:
        return None

    def _recessive_carrier_diseases(res):
        out = {}
        for e in (res.get("carriers", []) + res.get("affected", [])):
            if "recessive" in str(e.get("inheritance", "")).lower():
                out[e["disease"]] = e
        return out

    a_dis = _recessive_carrier_diseases(ca)
    b_dis = _recessive_carrier_diseases(cb)
    shared = sorted(set(a_dis) & set(b_dis))
    items = []
    for d in shared:
        items.append({
            "disease": d,
            "gene": a_dis[d].get("gene"),
            "variant": a_dis[d].get("variant"),
            "child_affected_risk": 0.25,
            "note": ("Both partners carry a recessive variant for this condition "
                     "→ each child has a 25% chance of being affected, 50% of "
                     "being a carrier. Genetic counselling strongly advised."),
        })
    return {
        "shared_recessive_conditions": items,
        "n_shared": len(items),
        "note": ("Couple carrier overlap is chip-limited; a clinical expanded "
                 "carrier screen (250+ conditions) is the standard before "
                 "conception if both partners want certainty."),
    }


def analyze_multi_person(df_a: pd.DataFrame, df_b: pd.DataFrame,
                         label_a: str = "Person A",
                         label_b: str = "Person B") -> dict:
    """Full two-genome comparison: kinship/relationship, concordance, and
    couple carrier risk."""
    king = king_kinship(df_a, df_b)
    relationship = classify_relationship(
        king["kinship"], king["ibs0_rate"], king["n_shared_snps"])
    couple = couple_carrier_risk(df_a, df_b)

    return {
        "available": king["n_shared_snps"] > 0,
        "label_a": label_a,
        "label_b": label_b,
        "king": king,
        "relationship": relationship,
        "couple_carrier": couple,
        "concordance_note": (
            "Genotype concordance is descriptive only — two unrelated people of "
            "the same ancestry share most common-allele genotypes by chance. "
            "Relationship is estimated from KING kinship, not from concordance."),
        "disclaimer": (
            "Research-grade relationship estimate, not a legal or clinical "
            "paternity/relationship test. Accuracy depends on shared-marker "
            "count (best when both files are from the same vendor) and is "
            "reduced by genotyping error."),
    }


# ── Standalone HTML page (never threaded into the single-person report) ────────

def _esc(s) -> str:
    return ("" if s is None else str(s)
            .replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def render_multi_person_html(result: dict, version: str = "v6.19.0") -> str:
    """Render the two-genome comparison as a standalone HTML page."""
    import datetime as _dt
    if not result or not result.get("available"):
        return "<html><body><p>No overlapping markers to compare.</p></body></html>"

    king = result["king"]
    rel = result["relationship"]
    la, lb = _esc(result["label_a"]), _esc(result["label_b"])

    kin = king["kinship"]
    kin_txt = f"{kin:.4f}" if kin is not None else "—"
    conc = king["concordance"]
    conc_txt = f"{conc*100:.1f}%" if conc is not None else "—"

    couple_html = ""
    cc = result.get("couple_carrier")
    if cc and cc.get("shared_recessive_conditions"):
        rows = "".join(
            f"<tr><td><strong>{_esc(i['disease'])}</strong></td>"
            f"<td>{_esc(i['gene'])} {_esc(i['variant'])}</td>"
            # Use the computed value. Hardcoding 25% here meant a future
            # inheritance model that yields anything else would be silently
            # misreported by the display.
            f"<td style='color:#b3261e;font-weight:700'>"
            f"{float(i.get('child_affected_risk', 0.25)):.0%} per child</td>"
            f"<td>{_esc(i['note'])}</td></tr>"
            for i in cc["shared_recessive_conditions"])
        couple_html = f"""
<h2 style="color:#b3261e">⚠️ Shared recessive carrier status</h2>
<p>Both individuals carry a recessive variant for the following — each future
child would have a 25% chance of being affected. Genetic counselling advised.</p>
<table><thead><tr><th>Condition</th><th>Variant</th><th>Child risk</th><th>Notes</th></tr></thead>
<tbody>{rows}</tbody></table>
<p style="font-size:.9em;color:#5b6673">{_esc(cc.get("note", ""))}</p>"""
    elif cc:
        couple_html = ("<h2>Couple carrier overlap</h2><p>No shared recessive "
                       "carrier conditions detected among chip-testable variants. "
                       + _esc(cc.get("note", "")) + "</p>")

    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Genome Comparison — {la} vs {lb}</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;
  max-width:860px;margin:0 auto;padding:32px 24px 80px;color:#1f2328;line-height:1.6}}
h1{{font-size:24px;margin-bottom:4px}} h2{{font-size:18px;margin:24px 0 8px;
  border-bottom:1px solid #d0d7de;padding-bottom:6px}}
.sub{{color:#8b949e;font-size:13px;font-style:italic}}
.hero{{display:flex;gap:24px;flex-wrap:wrap;background:linear-gradient(135deg,#f4f8fc,#eef2f7);
  border:1.5px solid #12467a;border-radius:14px;padding:20px 24px;margin:18px 0}}
.metric{{text-align:center;min-width:150px}}
.metric .v{{font-size:2.2em;font-weight:800;color:#12467a;line-height:1}}
.metric .l{{font-size:.72em;color:#5b6673;text-transform:uppercase;letter-spacing:.05em}}
.rel{{flex:1;min-width:240px}}
table{{border-collapse:collapse;width:100%;font-size:13px;margin-top:8px}}
th,td{{border:1px solid #d0d7de;padding:7px 9px;text-align:left;vertical-align:top}}
th{{background:#f6f8fa}}
.disc{{background:rgba(210,153,34,.08);border-left:4px solid #d29922;border-radius:6px;
  padding:14px 18px;margin-top:24px;font-size:13px}}
</style></head><body>
<h1>Genome Comparison</h1>
<div class="sub">{la} vs {lb} · generated {_dt.datetime.now().strftime("%B %d, %Y")} · {_esc(version)}</div>

<div class="hero">
  <div class="metric"><div class="v">{kin_txt}</div><div class="l">KING kinship φ</div></div>
  <div class="metric"><div class="v">{conc_txt}</div><div class="l">Genotype concordance</div></div>
  <div class="rel">
    <div style="font-size:1.3em;font-weight:800;color:#12467a">{_esc(rel['degree'])}</div>
    <div style="margin:4px 0">{_esc(rel['label'])}</div>
    <div style="color:#5b6673;font-size:.9em">{_esc(rel['detail'])}</div>
    <div style="color:#8a94a3;font-size:.82em;margin-top:4px">Confidence: {_esc(rel['confidence'])} · {king['n_shared_snps']:,} shared autosomal SNPs</div>
  </div>
</div>

<p style="color:#5b6673;font-size:.9em">{_esc(result['concordance_note'])}</p>

<h2>How the relationship was estimated</h2>
<p style="font-size:.92em">KING-robust kinship φ = (N_HetHet − 2·N_IBS0) /
(N_Het<sub>A</sub> + N_Het<sub>B</sub>) over shared autosomal biallelic SNPs —
robust to population structure, needs no allele frequencies. Thresholds:
φ&gt;0.354 duplicate/twin · 0.177–0.354 1st-degree · 0.088–0.177 2nd-degree ·
0.044–0.088 3rd-degree · &lt;0.044 unrelated. Parent-child vs full-sibling is
refined with IBS0 (opposite-homozygote rate ≈0 for parent-child).</p>

{couple_html}

<div class="disc"><strong>Important:</strong> {_esc(result['disclaimer'])}</div>
</body></html>"""
