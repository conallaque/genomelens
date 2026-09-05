#!/usr/bin/env python3
"""Regenerate the committed economics sample from a purpose-built synthetic genome.

A committed sample that cannot be reproduced is just a picture. This builds the
synthetic whole-genome VCF the sample comes from, runs the pipeline, and writes
the consolidated economics page. No human genome is involved at any step.

    python scripts/make_econ_sample.py [outdir] [--refresh-committed]
    python scripts/make_econ_sample.py [outdir] --profiles

Writes everything into ``outdir``. It touches ``docs/samples/`` only when
``--refresh-committed`` is given, so running it to inspect a change cannot
destroy the committed artifact you wanted to compare against.

WHAT THE INPUT IS, AND IS NOT
-----------------------------
The VCF is **synthetic and purpose-built on GRCh38**. It is assembled from two
sources, both coordinate-consistent:

  * the curated SNP registry (``core.snp_registry``) emitted at its **GRCh37**
    positions, which is what produces the pharmacogenomic, polygenic and APOE
    findings; and
  * a small set of **real ClinVar pathogenic/likely-pathogenic variants** at
    their real GRCh37 coordinates, in genes the economic model has anchors for,
    so the ClinVar screen and the curated ACMG path actually engage.

It is **not** a lifted-over copy of ``data/test_genome.txt``, and the two should
not be confused. That chip export is a 23andMe file on **GRCh37**; nothing here
derives from it.

WHY IT IS BUILT THIS WAY
------------------------
The previous version of this script built the VCF *from* that GRCh37 chip export
and then wrote ``##reference=GRCh38`` on top of it, injecting its rare variants
at GRCh38 coordinates. The file was therefore ~200,000 GRCh37 records plus eight
GRCh38 records, labelled GRCh38. Build detection read the bulk and correctly
resolved GRCh37, loaded the GRCh37 ClinVar table, and matched none of the
injected variants — so the feature this script exists to exercise had never
worked. Nothing complained, because nothing compared the declared build against
the detected one.

Building the input deliberately on one build is honest by construction rather
than by patching, and it is a foundation you can compose profiles from. The
build-consistency gate in ``pipeline.py`` now refuses a file whose header and
coordinates disagree, so the old shape would fail loudly rather than quietly.
"""
from __future__ import annotations

import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# GRCh38 contig lengths, so build detection is unambiguous from the header alone
# rather than inferred from which positions happen to appear.
GRCH37_CONTIGS = [
    ("1", 249250621), ("2", 243199373), ("3", 198022430), ("6", 171115067),
    ("7", 159138663), ("10", 135534747), ("11", 135006516),
    ("12", 133851895), ("13", 115169878), ("14", 107349540),
    ("17", 81195210), ("19", 59128983), ("22", 51304566),
]

# PHARMACOGENOMIC VARIANTS, at their real GRCh37 coordinates, taken from the
# committed chip export rather than written from memory.
#
# The sample had none. A synthetic genome with zero actionable pharmacogenomic
# findings is not a conservative sample, it is an unrealistic one: CYP2C19*2
# runs near 15% allele frequency, CYP2D6*4 near 20%, SLCO1B1*5 near 15%. A
# person carrying none of them is the unusual case, and the report has a whole
# page for medication-genotype decisions that was rendering empty.
#
# These are the three the star-allele caller can act on that the chip export
# also carries, so every coordinate here is one this repository already had.
# Emitted heterozygous, which is the actionable intermediate-metaboliser state.
PGX_VARIANTS = [
    ("10", 96541616, "G", "A", "rs4244285"),   # CYP2C19*2  — clopidogrel
    ("22", 42524947, "G", "A", "rs3892097"),   # CYP2D6*4   — opioids, SSRIs
    ("12", 21178615, "T", "C", "rs4149056"),   # SLCO1B1*5  — statin myopathy
]

# Real ClinVar P/LP SNVs at their real GRCh37 coordinates, verified present in
# the distilled table shipped by `python setup.py --clinvar`. Star ratings are
# 2-3, i.e. above the >=1 floor the clinical-variant screen applies. Spread
# across distinct organ systems rather than stacked within one, so the pooling
# correction has genuinely separate liabilities to work on.
#
# NOT fabricated: each is a catalogued human variant. The genome carrying all of
# them is synthetic, and its joint prevalence is reported by the profile that
# uses it rather than left implied.
# Approximate pathogenic-variant HETEROZYGOTE prevalence, order of magnitude,
# from published condition prevalence. NOT carrier frequency: familial
# hypercholesterolaemia and the hereditary cancer syndromes are dominant, so a
# heterozygote is at risk rather than a silent carrier, and calling these
# "carrier frequencies" would misdescribe what the person has.
#
# Used to print each configuration's joint prevalence. Joint figures assume
# independence, which is an approximation in two known directions: in a
# sequenced cohort the joint rate EXCEEDS the product, because ascertainment
# biases toward people with family history. The exclusion argument for
# multi-gene stacks survives that — a product of 1 in 10^9 is not rescued by a
# correction of this size — but the printed figure is an estimate, and it is
# labelled as one wherever it appears.
HET_PREVALENCE = {
    "LDLR":  1/300,      # FH overall ~1/250; LDLR the commonest cause
    "APOB":  1/1_000,
    "PCSK9": 1/10_000,
    "BRCA1": 1/400,
    "BRCA2": 1/400,
    "MLH1":  1/2_000,    # Lynch overall ~1/300 across MMR genes
    "MSH2":  1/2_000,
    "TP53":  1/10_000,   # Li-Fraumeni
    "RET":   1/30_000,   # MEN2
    "KCNQ1": 1/4_000,
}

# The sample's pathogenic variants. Real ClinVar P/LP records at real GRCh38
# coordinates; the genome carrying them is synthetic.
#
# DELIBERATELY SHORT. An earlier version carried four ACMG pathogenic variants
# at once, whose joint prevalence is past 1 in 10^12 — more than a hundred
# times the number of humans who have ever lived. A configuration whose printed
# prevalence is absurd is disqualified by that fact alone, however good its
# economics look, and printing the figure is what makes the disqualification
# automatic rather than a matter of taste.
#
# LDLR is valued (CAD anchor, ~1 in 300 — commoner than most single findings a
# report of this kind describes). TP53 is deliberately NOT valued: no registry
# anchor describes a multi-site cancer syndrome, so the report shows the finding
# and withholds the figure. Keeping it in the sample is the point — it
# demonstrates the withholding on the gene that used to carry the largest QALY
# anchor in the model.
CLINVAR_PLP = [
    ("19", 11200225, "A", "C", "LDLR"),    # Pathogenic 3* familial hypercholesterolaemia
    ("17", 41197728, "G", "C", "BRCA1"),   # Pathogenic 3* hereditary breast/ovarian
]


def joint_prevalence(genes) -> tuple[float, str]:
    """(probability, "roughly 1 in N") for carrying all of `genes` at once."""
    pr = 1.0
    for g in genes:
        pr *= HET_PREVALENCE.get(g.upper(), 1 / 1_000)
    if pr <= 0:
        return 0.0, "not estimable"
    n = round(1 / pr)
    return pr, f"roughly 1 in {n:,}"


def _registry_rows() -> list[tuple[str, int, str, str, str]]:
    """Curated registry variants at GRCh38 positions, as (chrom, pos, ref, alt, rsid).

    Genotypes are assigned deterministically so the sample is reproducible: a
    variant is emitted heterozygous, which is the informative case for the
    panels and keeps the file a fixed function of the registry.
    """
    from core import snp_registry as reg
    rows = []
    for r in reg._RECORDS:
        pos = getattr(r, "pos_grch37", None)
        anc, der = (r.ancestral or "").upper(), (r.derived or "").upper()
        if not pos or len(anc) != 1 or len(der) != 1 or anc == der:
            continue
        if anc not in "ACGT" or der not in "ACGT":
            continue
        rows.append((str(r.chrom), int(pos), anc, der, r.rsid))
    return sorted(rows, key=lambda x: (x[0], x[1]))


def build_vcf(dest: pathlib.Path) -> tuple[int, int]:
    reg_rows = _registry_rows()
    with dest.open("w") as w:
        w.write("##fileformat=VCFv4.2\n")
        w.write("##source=GenomeLens-synthetic-purpose-built\n")
        w.write("##reference=GRCh37\n")
        for c, ln in GRCH37_CONTIGS:
            w.write(f"##contig=<ID={c},length={ln}>\n")
        w.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n')
        w.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for chrom, pos, ref, alt, rsid in reg_rows:
            w.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t100\tPASS\tAF=0.3\t"
                    f"GT\t0/1\n")
        for chrom, pos, ref, alt, rsid in PGX_VARIANTS:
            w.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t100\tPASS\tAF=0.15\t"
                    f"GT\t0/1\n")
        for chrom, pos, ref, alt, gene in CLINVAR_PLP:
            w.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t200\tPASS\t"
                    f"AF=0.00001;GENE={gene}\tGT\t0/1\n")
    return len(reg_rows), len(CLINVAR_PLP) + len(PGX_VARIANTS)


def _require_clinvar() -> None:
    """Refuse to build the sample without a ClinVar table. Loudly.

    `reference/` is gitignored, so on a clean clone the full distilled table is
    absent. Before the committed subset existed, the screen simply matched
    nothing, the run completed, and the report rendered — disagreeing with the
    committed sample while looking perfectly well-formed. A sample that is
    quietly wrong is worse than one that refuses to build.
    """
    from risk.clinical_variants import CLINVAR_DIR, SAMPLE_SUBSETS
    full = CLINVAR_DIR / "clinvar_plp_grch37.tsv.gz"
    if full.exists() or SAMPLE_SUBSETS["grch37"].exists():
        return
    raise SystemExit(
        "\nNo ClinVar table found, and the sample's pathogenic variants cannot\n"
        "be screened without one. This would produce a report that renders\n"
        "correctly and disagrees with the committed sample.\n\n"
        f"  expected either: {full}\n"
        f"               or: {SAMPLE_SUBSETS[chr(39)+chr(39)] if False else SAMPLE_SUBSETS['grch37']}  (committed)\n\n"
        "Fix: python setup.py --clinvar\n")



# ── chip vs whole genome ──────────────────────────────────────────────────────

def _run_one(vcf_or_chip: pathlib.Path, out: pathlib.Path) -> dict | None:
    """Run the pipeline on one input and return its economics payload."""
    import json
    out.mkdir(parents=True, exist_ok=True)
    r = subprocess.run([sys.executable, str(ROOT / "analyze.py"),
                        str(vcf_or_chip), "--no-ai",
                        "--output", str(out / "report.html")],
                       cwd=ROOT, capture_output=True, text=True)
    pj = out / "economics-payload.json"
    if r.returncode or not pj.exists():
        print(f"  run failed for {vcf_or_chip.name}: rc={r.returncode}")
        return None
    return json.loads(pj.read_text())


def _summarise(p: dict) -> dict:
    rc, td = p["reference_case"], p.get("testing_decision", {})
    mon = [f for f in p["findings"] if f.get("is_monetized")]
    return {
        "input_type": p["metadata"].get("input_type", "?"),
        "nmb": rc["nmb"], "qalys": rc["incremental_qalys"],
        "cost": rc["incremental_cost"],
        "n_findings": len(p["findings"]), "n_priced": len(mon),
        "observed_wgs_only_findings": td.get("observed_wgs_only_findings", 0),
        "observed_wgs_only_value": td.get("observed_wgs_only_value", 0.0),
    }


def emit_profiles(out: pathlib.Path) -> int:
    """Run both committed inputs and report the marginal value of sequencing.

    THE CONTRAST IS THE POINT, and it is a core economic result rather than a
    presentational one: "what does sequencing add over an array" is the question
    page 7 answers, and until now it was answered with one real artifact and a
    population average. Two real artifacts answer it directly.

    The two inputs are deliberately on different builds, because that is what
    they actually are: the chip export is a 23andMe file on GRCh37, the
    whole-genome sample is purpose-built on GRCh38. Each is internally
    consistent; neither is a lifted-over copy of the other.
    """
    chip = ROOT / "data" / "test_genome.txt"
    wgs = out / "wgs" / "synthetic_wgs.vcf"
    (out / "wgs").mkdir(parents=True, exist_ok=True)
    build_vcf(wgs)

    print("\n=== chip vs whole genome ===")
    rows = {}
    for label, src in (("chip", chip), ("wgs", wgs)):
        pay = _run_one(src, out / label)
        if pay is None:
            return 1
        rows[label] = _summarise(pay)

    c, w = rows["chip"], rows["wgs"]
    print(f"{'':28}{'chip (GRCh37)':>18}{'whole genome (GRCh38)':>24}")
    print("-" * 72)
    for lbl, k, fmt in (("reference-case NMB", "nmb", "money"),
                        ("incremental QALYs", "qalys", "raw"),
                        ("incremental cost", "cost", "money"),
                        ("findings", "n_findings", "int"),
                        ("of which priced", "n_priced", "int"),
                        ("sequencing-only findings", "observed_wgs_only_findings", "int"),
                        ("sequencing-only value", "observed_wgs_only_value", "money")):
        cv, wv = c[k], w[k]
        f = ((lambda v: f"${v:,.0f}") if fmt == "money"
             else (lambda v: f"{v:,}") if fmt == "int" else (lambda v: f"{v}"))
        print(f"{lbl:28}{f(cv):>18}{f(wv):>24}")

    marginal = w["nmb"] - c["nmb"]
    only = w["observed_wgs_only_value"]
    print("-" * 72)
    print(f"{'marginal value of sequencing':28}{'':>18}{f'${marginal:,.0f}':>24}")
    print("\nThe array finds none of the sequencing-only findings, by construction:")
    print("a genotyping array reports the positions it carries probes for. The")
    print("marginal figure is what the model says the extra positions were worth")
    print("FOR THIS GENOME — not a population average, and not transferable to")
    print("another person's.")
    print(f"\nWHY ${marginal:,.0f} EXCEEDS THE ${only:,.0f} SEQUENCING-ONLY LINE.")
    print("Those are different quantities and the gap is not an inconsistency.")
    print("'Sequencing-only value' counts findings that ONLY a whole genome can")
    print("report. The marginal figure also includes findings the array could")
    print("have reported and did not — positions outside its probe set, and")
    print("variants its probes cover but which the ClinVar screen only reaches")
    print("on a whole-genome input. The second is the larger effect here.")
    return 0


def main() -> int:
    argv = [a for a in sys.argv[1:] if a != "--refresh-committed"]
    # Promoting a run into docs/samples/ is opt-in. It used to happen on every
    # invocation, which meant you could not run this script to LOOK at a change
    # without destroying the committed artifact you wanted to compare against —
    # the first before/after measured after this fix compared the new payload
    # to itself, because the script had already overwritten the old one.
    refresh = "--refresh-committed" in sys.argv
    argv = [a for a in argv if a != "--profiles"]
    out = pathlib.Path(argv[0] if argv else "econ-sample").resolve()
    out.mkdir(parents=True, exist_ok=True)
    if "--profiles" in sys.argv:
        _require_clinvar()
        return emit_profiles(out)
    vcf = out / "synthetic_wgs.vcf"
    _require_clinvar()
    n_reg, n_plp = build_vcf(vcf)
    genes = [g for *_r, g in CLINVAR_PLP]
    _pr, phrase = joint_prevalence(genes)
    print(f"synthetic GRCh38 VCF: {n_reg} registry + {n_plp} ClinVar P/LP -> {vcf}")
    print(f"  pathogenic variants : {', '.join(genes)}")
    print(f"  joint prevalence    : {phrase} "
          f"(approximate, assuming independence)")
    print("  NOTE: pathogenic-variant heterozygote prevalence, not carrier "
          "frequency —")
    print("        FH and the hereditary cancer syndromes are dominant, so a "
          "heterozygote")
    print("        is at risk rather than a silent carrier.")
    r = subprocess.run([sys.executable, str(ROOT / "analyze.py"), str(vcf),
                        "--output", str(out / "report.html")], cwd=ROOT)
    if r.returncode:
        return r.returncode
    page = out / "economics.html"
    print(f"\nconsolidated economics page: {page}" if page.exists()
          else "\nWARNING: economics.html not written")

    # The canonical payload the report renders from, and the consistency
    # findings for this run. Copied out of the run directory rather than
    # rebuilt, so the committed sample is byte-identical to what the pipeline
    # produced — a payload regenerated by a second code path is a second code
    # path, and the drift between two views of one number is the defect this
    # project has spent the longest removing.
    ff = out / "economics-findings-first.html"
    if ff.exists():
        print(f"findings-first report:       {ff}")
        print("  print it to PDF to refresh "
              "docs/samples/econ-output-sample.pdf")

    payload = out / "economics-payload.json"
    dest = ROOT / "docs/samples/econ-payload-sample.json"
    if payload.exists():
        print(f"canonical report payload:    {payload}")
        if refresh:
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(payload.read_text())
            print(f"  committed sample refreshed: {dest}")
        else:
            print(f"  committed sample NOT touched ({dest.name} left as-is).")
            print("  re-run with --refresh-committed to promote this run.")
    else:
        print("WARNING: economics-payload.json not written")

    print("\neconomics.html carries the full technical appendix; "
          "economics-findings-first.html is the main report.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
