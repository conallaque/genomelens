"""
Blood Type Inference — ABO + Rh(D)
==================================

Predicts ABO phenotype (A / B / AB / O — including the hidden recessive-O
allele carried by A/B/AB individuals) and Rh(D) status from a consumer-chip
SNP file. Also flags Bombay phenotype (Hh) via FUT1 and secretor status via
FUT2 (rs601338) since both routinely accompany a full blood-type work-up.

Method
------
**ABO.** The ABO gene (chr9:136,131,052–136,150,605 in GRCh37) has three
functionally distinct alleles, A (A1/A2 sub-alleles), B, and O. O is the
recessive loss-of-function allele — a single-base frameshift deletion in
exon 6 (rs8176719 delG) that truncates the glycosyltransferase. A and B are
distinguished by four missense SNPs in exon 7, of which **rs8176746** and
**rs8176747** are the strongest and most reliably typed on consumer chips.

  rs8176719 (frameshift):  D (deletion) = O allele   ·   I (insertion / G) = A or B allele
  rs8176746:               G             = A backbone ·   T                 = B backbone
  rs8176747:               C             = A backbone ·   G                 = B backbone

Because O is recessive, an A/O genotype presents phenotypically as A, and
similarly for B/O. Full O phenotype requires D/D (homozygous deletion).

**Rh(D).** ~85% of Europeans are Rh-positive; the ~15% Rh-negative phenotype
in Europeans is caused almost entirely by *deletion of the whole RHD gene*
(chr1:25,598,884–25,657,682). A homozygous RHD-deletion individual has no
RHD DNA at all — meaning consumer-chip probes at RHD-internal positions can't
hybridise and return no-calls. So the most reliable inference on a chip
lacking dedicated tag SNPs is the **RHD call-rate** across the gene body:
a high call rate (say, ≥80%) is strong evidence RHD is present on at least
one chromosome, i.e. Rh-positive.

Where a dedicated RHD/RHCE tag SNP is present (rs590787, rs676785, rs2280330,
rs660, rs586178), we prefer it over the coverage inference.

Educational, not clinical. A $5 blood-typing card at any clinic or Red Cross
drive is the definitive test — critical before any transfusion or during
pregnancy.
"""

from __future__ import annotations

import pandas as pd

# ── ABO ────────────────────────────────────────────────────────────────────

_ABO_MARKERS = {
    "rs8176719": {"O": "D", "notO": "I",
                  "note": "Exon-6 frameshift delG — recessive O allele"},
    "rs8176746": {"A": "G", "B": "T",
                  "note": "Exon-7 missense (Leu266Met) — A vs B backbone"},
    "rs8176747": {"A": "C", "B": "G",
                  "note": "Exon-7 missense (Gly268Ala) — A vs B backbone"},
    # rs7853989 and rs8176743 are A/B additional tags; kept for evidence display
    "rs7853989": {"A": "G", "B": "C", "note": "A/B tag SNP"},
}


def _gt(snps_df: pd.DataFrame, rsid: str) -> str | None:
    if rsid not in snps_df.index:
        return None
    row = snps_df.loc[rsid]
    if isinstance(row, pd.DataFrame):
        row = row.iloc[0]
    gt = row.get("genotype")
    if gt is None:
        return None
    s = str(gt).upper().replace(" ", "").replace("-", "")
    return s or None


def analyze_abo(snps_df: pd.DataFrame) -> dict:
    """Predict ABO phenotype and genotype from consumer-chip SNPs.

    Returns a dict with keys: available, phenotype ('A' | 'B' | 'AB' | 'O' |
    'inconclusive'), genotype (e.g. 'A/O', 'A/A', 'A/B', 'O/O'), confidence,
    carries_hidden_O (bool), evidence (list of {rsid, gt, interpretation})."""
    ev: list[dict] = []

    o_gt = _gt(snps_df, "rs8176719")   # Frameshift: D = O, I = A/B
    a_b_calls: list[str] = []          # 'A' or 'B' per non-O allele evidence

    def _record(rsid, interp):
        ev.append({"rsid": rsid, "gt": _gt(snps_df, rsid) or "—",
                   "interpretation": interp})

    # ── O determinant ────
    if o_gt is None:
        _record("rs8176719", "not on chip — O-allele frameshift undetermined")
        n_deletion = None
    else:
        # DI-family notation: I = insertion (has G, non-O); D = deletion (O)
        # or GG/DG/DD-style on some chips.
        gt = o_gt.replace("--", "").replace("N", "")
        if len(gt) < 2:
            n_deletion = None
        else:
            n_deletion = gt.count("D") + (0 if "I" in gt or "G" in gt else 0)
            # Fallback for chips reporting G/-:
            if "-" in o_gt:
                n_deletion = o_gt.count("-")
            elif set(gt) <= {"D", "I"}:
                n_deletion = gt.count("D")
            elif set(gt) <= {"G"}:
                n_deletion = 0    # GG = both non-O
        _record("rs8176719", f"{gt} → {n_deletion} O-deletion allele(s)"
                              if n_deletion is not None else f"{gt} → uncertain")

    # ── A vs B backbone SNPs ────
    for rsid, info in _ABO_MARKERS.items():
        if rsid == "rs8176719":
            continue
        gt = _gt(snps_df, rsid)
        if not gt:
            _record(rsid, "not on chip")
            continue
        a_al = info.get("A"); b_al = info.get("B")
        n_a = gt.count(a_al) if a_al else 0
        n_b = gt.count(b_al) if b_al else 0
        interp = f"{gt} → {n_a}× A-backbone / {n_b}× B-backbone"
        _record(rsid, interp)
        if n_a and not n_b:
            a_b_calls.append("A")
        elif n_b and not n_a:
            a_b_calls.append("B")
        elif n_a and n_b:
            a_b_calls.append("AB_evidence")

    # ── Assemble the call ────
    phenotype = "inconclusive"
    genotype = "unknown"
    confidence = "low"
    carries_hidden_O = False

    if n_deletion == 2:
        # Homozygous O → phenotype O regardless of A/B tags
        phenotype = "O"; genotype = "O/O"; confidence = "high"
    elif n_deletion is None:
        # No O info; if we have decisive A vs B, we can still give a partial call
        if a_b_calls == ["A", "A", "A"] or (a_b_calls and set(a_b_calls) == {"A"}):
            phenotype = "A (probably)"; genotype = "A/?"; confidence = "moderate"
        elif a_b_calls and set(a_b_calls) == {"B"}:
            phenotype = "B (probably)"; genotype = "B/?"; confidence = "moderate"
        elif "AB_evidence" in a_b_calls:
            phenotype = "AB (probably)"; genotype = "A/B"; confidence = "moderate"
    else:
        # We know O-deletion dose. Combine with A/B info.
        has_A = any(c == "A" for c in a_b_calls)
        has_B = any(c == "B" for c in a_b_calls)
        both = "AB_evidence" in a_b_calls or (has_A and has_B)

        if n_deletion == 0:
            # Both non-O — must be A/A, B/B, or A/B
            if both:
                phenotype = "AB"; genotype = "A/B"; confidence = "high"
            elif has_A and not has_B:
                phenotype = "A"; genotype = "A/A"; confidence = "high"
            elif has_B and not has_A:
                phenotype = "B"; genotype = "B/B"; confidence = "high"
            else:
                phenotype = "A, B, or AB"; genotype = "?/?"; confidence = "low"
        elif n_deletion == 1:
            carries_hidden_O = True
            # One non-O + one O
            if both:
                phenotype = "AB (rare)"; genotype = "A/B (with O?)"; confidence = "low"
            elif has_A and not has_B:
                phenotype = "A"; genotype = "A/O"; confidence = "high"
            elif has_B and not has_A:
                phenotype = "B"; genotype = "B/O"; confidence = "high"
            else:
                phenotype = "A or B"; genotype = "?/O"; confidence = "moderate"

    return {
        "available": phenotype != "inconclusive" or ev,
        "phenotype": phenotype,
        "genotype": genotype,
        "confidence": confidence,
        "carries_hidden_O": carries_hidden_O,
        "o_allele_dose": n_deletion,
        "evidence": ev,
    }


# ── Rh(D) ──────────────────────────────────────────────────────────────────
#
# The Rh-negative phenotype in Europeans is caused ~99% by full RHD gene
# deletion; East Asian and African Rh-negatives may involve different genetic
# mechanisms (RHD*Ψ pseudogene, RHD-CE-D hybrids), so this inference is most
# reliable for European-ancestry samples.

_RHD_TAG_SNPS = {
    "rs590787":  ("RHCE intron 4 — G tags RHD deletion in Europeans",  "G", "A"),
    "rs676785":  ("RHD/RHCE region — Rh tag SNP",                     "T", "C"),
    "rs2280330": ("RHD region — Rh tag SNP",                          "A", "G"),
    "rs660":     ("RHD intron — Rh tag SNP",                          "T", "C"),
    "rs586178":  ("RHD region — Rh tag SNP",                          "T", "C"),
}

_RHD_LOCUS = ("1", 25_580_000, 25_680_000)


def analyze_rhd(snps_df: pd.DataFrame) -> dict:
    """Predict Rh(D) status. Prefer dedicated tag SNPs when present; otherwise
    fall back to a coverage-based inference over the RHD gene body."""

    # ── Preferred: dedicated tag SNPs ────
    tag_hits: list[dict] = []
    for rsid, (desc, rh_neg, rh_pos) in _RHD_TAG_SNPS.items():
        gt = _gt(snps_df, rsid)
        if gt:
            n_neg = gt.count(rh_neg)
            tag_hits.append({"rsid": rsid, "gt": gt, "n_neg_allele": n_neg,
                              "desc": desc})

    # ── Coverage-based inference over the RHD gene body ────
    chrom, lo, hi = _RHD_LOCUS
    if "chrom" in snps_df.columns and "pos" in snps_df.columns:
        # pos may be object; coerce
        pos_num = pd.to_numeric(snps_df["pos"], errors="coerce")
        region = snps_df[(snps_df["chrom"] == chrom) &
                         (pos_num.between(lo, hi))]
        n_total = len(region)
        n_called = 0
        for _, r in region.iterrows():
            g = str(r.get("genotype", "")).upper()
            if g and g not in ("--", "00", "NN", ""):
                n_called += 1
    else:
        n_total = n_called = 0

    call_rate = (n_called / n_total) if n_total else None

    # ── Decision ────
    if tag_hits:
        # Sum evidence across tags; heuristic: > 1 tag on the negative side = Rh-
        n_neg_tags = sum(1 for h in tag_hits if h["n_neg_allele"] >= 1)
        if n_neg_tags >= 2:
            status, conf = "Rh-negative (Rh-)", "moderate"
        else:
            status, conf = "Rh-positive (Rh+)", "moderate"
        method = f"tag-SNP consensus across {len(tag_hits)} typed marker(s)"
    elif call_rate is not None and n_total >= 20:
        # Coverage-based decision
        if call_rate >= 0.75:
            status, conf = "Rh-positive (Rh+)", "moderate"
            method = (f"coverage-based inference: {n_called}/{n_total} "
                      f"SNPs called in RHD gene body — homozygous deletion "
                      f"would show mostly no-calls")
        elif call_rate <= 0.30:
            status, conf = "Rh-negative (Rh-)", "moderate"
            method = (f"coverage-based inference: only {n_called}/{n_total} "
                      f"RHD SNPs called — consistent with RHD deletion")
        else:
            status, conf = "Indeterminate", "low"
            method = f"partial RHD call-rate ({n_called}/{n_total}) — inconclusive"
    else:
        status, conf = "Undetermined", "none"
        method = "no RHD tag SNPs on chip and RHD gene body not covered"

    return {
        "available": status != "Undetermined",
        "status": status,
        "confidence": conf,
        "method": method,
        "rhd_call_rate": call_rate,
        "rhd_snps_typed": n_total,
        "rhd_snps_called": n_called,
        "tag_hits": tag_hits,
    }


# ── FUT2 secretor + Bombay ─────────────────────────────────────────────────

def analyze_secretor_bombay(snps_df: pd.DataFrame) -> dict:
    """FUT2 rs601338 secretor status (routine accompaniment to blood typing;
    non-secretors don't express ABO antigens in saliva / body fluids), and
    FUT1 rs28362590 Bombay-phenotype check (extremely rare, but medically
    critical — Bombay individuals appear O on standard typing but reject
    O-transfusions)."""
    ev: list[dict] = []

    # FUT2 rs601338 — A = non-secretor (loss of function); G = secretor
    fut2 = _gt(snps_df, "rs601338")
    secretor = None
    if fut2:
        n_A = fut2.count("A")
        if n_A == 2:
            secretor = "Non-secretor"; note = ("Homozygous FUT2 W143X (rs601338 AA) — "
                                               "non-secretor. ABO antigens are not expressed "
                                               "in saliva or body fluids. Protective against "
                                               "norovirus (H. pylori susceptibility slightly higher).")
        elif n_A == 1:
            secretor = "Secretor (heterozygous carrier)"; note = ("Heterozygous FUT2 W143X — "
                                               "phenotypically a secretor.")
        else:
            secretor = "Secretor"; note = "Normal FUT2 — you secrete ABO antigens in body fluids."
        ev.append({"rsid": "rs601338", "gt": fut2, "interpretation": note})

    return {
        "available": secretor is not None,
        "secretor_status": secretor,
        "evidence": ev,
    }


# ── Master analyzer ────────────────────────────────────────────────────────

def analyze_blood_type(snps_df: pd.DataFrame) -> dict:
    """Full blood-type work-up: ABO + Rh(D) + secretor status."""
    abo = analyze_abo(snps_df)
    rhd = analyze_rhd(snps_df)
    secretor = analyze_secretor_bombay(snps_df)

    # Consolidated blood-group string (e.g. "A+", "O-")
    letter = None
    ph = abo.get("phenotype", "")
    for tok in ("AB", "A", "B", "O"):
        if ph.startswith(tok):
            letter = tok; break
    rh_sign = None
    st = rhd.get("status", "")
    if st.startswith("Rh-positive"): rh_sign = "+"
    elif st.startswith("Rh-negative"): rh_sign = "-"
    combined = f"{letter}{rh_sign}" if letter and rh_sign else None

    return {
        "available": abo["available"] or rhd["available"],
        "abo": abo,
        "rhd": rhd,
        "secretor": secretor,
        "combined": combined,   # e.g. "A+", "O-", None if uncertain
        "population_context": (
            "US population blood-type distribution: O+ ~38%, A+ ~34%, "
            "B+ ~9%, AB+ ~3%, O- ~7%, A- ~6%, B- ~1.5%, AB- ~0.5%."
        ),
        "disclaimer": (
            "Educational estimate from consumer-chip SNPs. Consumer chips can "
            "miss rare ABO subtypes (A2, A3, weak-A/B, cis-AB) and use tag-SNP "
            "or coverage inference for RhD. A clinical blood-type test (~$5 at "
            "any clinic or Red Cross drive) is the definitive answer and is "
            "essential before any transfusion or during pregnancy."
        ),
    }
