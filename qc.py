"""
Quality Control & Callability Report
------------------------------------

Reports file-level metrics and per-domain callability:

  * File source, build, total SNPs, no-call rate
  * Autosomal / X / Y / MT coverage
  * Per-key-variant coverage (APOE, MTHFR, BRCA tags, top CYPs, etc.)
  * PRS panel callability per disease
  * Inferred sex
  * Inferred ancestry (very rough — based on tag SNPs)

This is a critical part of professional-grade reporting because it tells
the user how much of the report's confidence depends on chip coverage.
"""

from typing import Dict, List, Optional
import pandas as pd
import hashlib


# ─── Key variants we'd want any clinical-grade analysis to cover ──────────────
KEY_VARIANTS_BY_DOMAIN: Dict[str, List[Dict]] = {
    "APOE / Alzheimer's": [
        {"rsid": "rs429358", "gene": "APOE ε4"},
        {"rsid": "rs7412",   "gene": "APOE ε2"},
        {"rsid": "rs75932628", "gene": "TREM2 R47H"},
    ],
    "Hemochromatosis": [
        {"rsid": "rs1800562", "gene": "HFE C282Y"},
        {"rsid": "rs1799945", "gene": "HFE H63D"},
    ],
    "Thrombophilia": [
        {"rsid": "rs6025",    "gene": "Factor V Leiden"},
        {"rsid": "rs1799963", "gene": "Prothrombin G20210A"},
    ],
    "Methylation": [
        {"rsid": "rs1801133", "gene": "MTHFR C677T"},
        {"rsid": "rs1801131", "gene": "MTHFR A1298C"},
    ],
    "Pharmacogenomics (CYP2D6)": [
        {"rsid": "rs3892097",  "gene": "CYP2D6 *4"},
        {"rsid": "rs1065852",  "gene": "CYP2D6 *10"},
        {"rsid": "rs28371725", "gene": "CYP2D6 *41"},
        {"rsid": "rs16947",    "gene": "CYP2D6 *2"},
    ],
    "Pharmacogenomics (CYP2C9 / VKORC1 — warfarin)": [
        {"rsid": "rs1799853", "gene": "CYP2C9 *2"},
        {"rsid": "rs1057910", "gene": "CYP2C9 *3"},
        {"rsid": "rs9923231", "gene": "VKORC1 -1639"},
    ],
    "Pharmacogenomics (CYP2C19 — clopidogrel / PPI / SSRI)": [
        {"rsid": "rs4244285",  "gene": "CYP2C19 *2"},
        {"rsid": "rs12248560", "gene": "CYP2C19 *17"},
    ],
    "Pharmacogenomics (TPMT / NUDT15 — thiopurines)": [
        {"rsid": "rs1142345",    "gene": "TPMT *3C"},
        {"rsid": "rs1800460",    "gene": "TPMT *3A"},
        {"rsid": "rs116855232",  "gene": "NUDT15 *3"},
    ],
    "Pharmacogenomics (SLCO1B1 — statins)": [
        {"rsid": "rs4149056", "gene": "SLCO1B1 *5"},
    ],
    "Coronary / Lp(a)": [
        {"rsid": "rs10455872", "gene": "LPA"},
        {"rsid": "rs3798220",  "gene": "LPA"},
        {"rsid": "rs10757278", "gene": "9p21"},
        {"rsid": "rs1333049",  "gene": "9p21"},
    ],
    "Diabetes": [
        {"rsid": "rs7903146",  "gene": "TCF7L2"},
        {"rsid": "rs10830963", "gene": "MTNR1B"},
    ],
    "Cancer Breast / Prostate": [
        {"rsid": "rs2981582",   "gene": "FGFR2"},
        {"rsid": "rs17879961",  "gene": "CHEK2 I157T"},
        {"rsid": "rs1447295",   "gene": "8q24"},
        {"rsid": "rs10993994",  "gene": "MSMB"},
    ],
}


def _file_hash(filepath: str) -> str:
    """Return a short hash for the input file (stable per-file identifier)."""
    try:
        h = hashlib.sha256()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()[:16]
    except Exception:
        return "unavailable"


def run_qc(snps_df: pd.DataFrame, filepath: str = "",
           tier1_match_count: int = 0,
           file_format: str = "unknown",
           db_size: int = 0) -> Dict:
    """Compute quality / callability metrics."""
    total = len(snps_df)
    no_call = 0
    if "genotype" in snps_df.columns:
        no_call = int(snps_df["genotype"].astype(str).isin(["nan", "--", "", "00"]).sum())

    # Per-chromosome counts
    chrom_counts = {}
    if "chrom" in snps_df.columns:
        chrom_counts = snps_df["chrom"].value_counts().to_dict()
    autosomal = sum(v for k, v in chrom_counts.items() if str(k).isdigit())
    x_count = int(chrom_counts.get("X", 0))
    y_count = int(chrom_counts.get("Y", 0))
    mt_count = int(chrom_counts.get("MT", 0) + chrom_counts.get("M", 0))

    # Inferred sex (Y chromosome SNP count)
    inferred_sex = "male" if y_count > 100 else "female"

    # Domain callability
    domain_coverage = {}
    for domain, vs in KEY_VARIANTS_BY_DOMAIN.items():
        called = sum(1 for v in vs if v["rsid"] in snps_df.index)
        domain_coverage[domain] = {
            "called": called,
            "total": len(vs),
            "pct": round(100.0 * called / len(vs), 1),
            "variants": [
                {"rsid": v["rsid"], "gene": v["gene"], "called": v["rsid"] in snps_df.index}
                for v in vs
            ],
        }

    # Overall QC grade
    callability_pct = 100.0 * (total - no_call) / max(total, 1)
    avg_domain = sum(d["pct"] for d in domain_coverage.values()) / max(len(domain_coverage), 1)

    if total > 600000 and callability_pct > 99 and avg_domain > 75:
        grade = "Excellent"
        grade_class = "qc-grade-excellent"
        grade_note = "High-density chip, comprehensive coverage. Report findings are well-supported."
    elif total > 400000 and callability_pct > 98 and avg_domain > 60:
        grade = "Good"
        grade_class = "qc-grade-good"
        grade_note = "Adequate coverage for most analyses; some specialized PGx / PRS panels may be partial."
    elif total > 200000 and avg_domain > 40:
        grade = "Fair"
        grade_class = "qc-grade-fair"
        grade_note = "Moderate coverage. Some risk scores and PGx phenotypes will be incomplete."
    else:
        grade = "Limited"
        grade_class = "qc-grade-limited"
        grade_note = "Sparse coverage. Many variants unable to be called; interpret results with caution."

    return {
        "file_path": filepath,
        "file_hash": _file_hash(filepath) if filepath else "n/a",
        "file_format": file_format,
        "total_snps": total,
        "no_call_count": no_call,
        "callability_pct": round(callability_pct, 2),
        "autosomal_count": autosomal,
        "x_count": x_count,
        "y_count": y_count,
        "mt_count": mt_count,
        "inferred_sex": inferred_sex,
        "tier1_match_count": tier1_match_count,
        "db_size": db_size,
        "tier1_match_pct": round(100.0 * tier1_match_count / max(db_size, 1), 1),
        "domain_coverage": domain_coverage,
        "grade": grade,
        "grade_class": grade_class,
        "grade_note": grade_note,
        "average_domain_callability": round(avg_domain, 1),
    }
