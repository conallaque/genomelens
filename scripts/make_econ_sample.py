#!/usr/bin/env python3
"""Regenerate the committed economics sample from a purpose-built synthetic genome.

A committed sample that cannot be reproduced is just a picture. This builds the
synthetic whole-genome VCF the sample comes from, runs the pipeline, and writes
the consolidated economics page. No human genome is involved at any step.

    python scripts/make_econ_sample.py [outdir] [--refresh-committed]

Writes everything into ``outdir``. It touches ``docs/samples/`` only when
``--refresh-committed`` is given, so running it to inspect a change cannot
destroy the committed artifact you wanted to compare against.

WHAT THE INPUT IS, AND IS NOT
-----------------------------
The VCF is **synthetic and purpose-built on GRCh38**. It is assembled from two
sources, both coordinate-consistent:

  * the curated SNP registry (``core.snp_registry``) emitted at its **GRCh38**
    positions, which is what produces the pharmacogenomic, polygenic and APOE
    findings; and
  * a small set of **real ClinVar pathogenic/likely-pathogenic variants** at
    their real GRCh38 coordinates, in genes the economic model has anchors for,
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
GRCH38_CONTIGS = [
    ("1", 248956422), ("2", 242193529), ("3", 198295559), ("7", 159345973),
    ("10", 133797422), ("11", 135086622), ("13", 114364328),
    ("17", 83257441), ("19", 58617616), ("22", 50818468),
]

# Real ClinVar P/LP SNVs at their real GRCh38 coordinates, verified present in
# the distilled table shipped by `python setup.py --clinvar`. Star ratings are
# 2-3, i.e. above the >=1 floor the clinical-variant screen applies. Spread
# across distinct organ systems rather than stacked within one, so the pooling
# correction has genuinely separate liabilities to work on.
#
# NOT fabricated: each is a catalogued human variant. The genome carrying all of
# them is synthetic, and its joint prevalence is reported by the profile that
# uses it rather than left implied.
CLINVAR_PLP = [
    ("17", 7670669,  "G", "T", "TP53"),    # Pathogenic 3* Li-Fraumeni
    ("19", 11089549, "A", "C", "LDLR"),    # Pathogenic 3* familial hypercholesterolaemia
    ("3",  36993548, "A", "G", "MLH1"),    # Pathogenic 3* Lynch syndrome
]


def _registry_rows() -> list[tuple[str, int, str, str, str]]:
    """Curated registry variants at GRCh38 positions, as (chrom, pos, ref, alt, rsid).

    Genotypes are assigned deterministically so the sample is reproducible: a
    variant is emitted heterozygous, which is the informative case for the
    panels and keeps the file a fixed function of the registry.
    """
    from core import snp_registry as reg
    rows = []
    for r in reg._RECORDS:
        pos = getattr(r, "pos_grch38", None)
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
        w.write("##reference=GRCh38\n")
        for c, ln in GRCH38_CONTIGS:
            w.write(f"##contig=<ID={c},length={ln}>\n")
        w.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n')
        w.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for chrom, pos, ref, alt, rsid in reg_rows:
            w.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t100\tPASS\tAF=0.3\t"
                    f"GT\t0/1\n")
        for chrom, pos, ref, alt, gene in CLINVAR_PLP:
            w.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t200\tPASS\t"
                    f"AF=0.00001;GENE={gene}\tGT\t0/1\n")
    return len(reg_rows), len(CLINVAR_PLP)


def main() -> int:
    argv = [a for a in sys.argv[1:] if a != "--refresh-committed"]
    # Promoting a run into docs/samples/ is opt-in. It used to happen on every
    # invocation, which meant you could not run this script to LOOK at a change
    # without destroying the committed artifact you wanted to compare against —
    # the first before/after measured after this fix compared the new payload
    # to itself, because the script had already overwritten the old one.
    refresh = "--refresh-committed" in sys.argv
    out = pathlib.Path(argv[0] if argv else "econ-sample").resolve()
    out.mkdir(parents=True, exist_ok=True)
    vcf = out / "synthetic_wgs.vcf"
    n_reg, n_plp = build_vcf(vcf)
    print(f"synthetic GRCh38 VCF: {n_reg} registry + {n_plp} ClinVar P/LP -> {vcf}")
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
              "docs/samples/econ-output-sample.pdf (8 pages)")

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
          "economics-findings-first.html is the eight-page main report.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
