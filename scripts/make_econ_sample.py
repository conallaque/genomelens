#!/usr/bin/env python3
"""Regenerate the committed economics sample from a synthetic whole genome.

A committed sample that cannot be reproduced is just a picture. This builds the
synthetic VCF the sample came from, runs the pipeline, and writes the
consolidated economics page. No human genome is involved at any step.

    python scripts/make_econ_sample.py [outdir]

The VCF carries the committed synthetic chip genome's variants plus a small set
of rare missense variants in genes the economic model has cost anchors for, so
the whole-genome-only path (offline AlphaMissense scoring -> predicted
pathogenic rare -> the Pathogenic condition anchor) actually produces something.
On a chip-shaped input the sequencing-value question is answerable only in the
abstract, which is exactly what the sample is meant to show it no longer is.
"""
from __future__ import annotations

import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent

# Rare missense variants in anchored genes, GRCh38 coordinates.
RARE = [
    ("17", 43093464, "BRCA1", "G", "A"),
    ("13", 32340301, "BRCA2", "T", "C"),
    ("3",  37020409, "MLH1",  "C", "T"),
    ("17", 7674220,  "TP53",  "C", "T"),
    ("19", 11113489, "LDLR",  "G", "A"),
    ("11", 2782777,  "KCNQ1", "C", "G"),
    ("2",  47799469, "MSH2",  "G", "A"),
    ("1",  55039974, "PCSK9", "C", "T"),
]


def build_vcf(dest: pathlib.Path) -> int:
    n = 0
    with dest.open("w") as w:
        w.write("##fileformat=VCFv4.2\n##source=GenomeLens-synthetic-WGS\n")
        w.write("##reference=GRCh38\n")
        w.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n')
        w.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for line in (ROOT / "data/test_genome.txt").open(errors="replace"):
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            rsid, chrom, pos, gt = p[0], p[1], p[2], p[3].strip()
            if len(gt) != 2 or any(c not in "ACGT" for c in gt):
                continue
            ref = gt[0]
            hom = gt[1] == gt[0]
            alt = ("A" if ref != "A" else "G") if hom else gt[1]
            w.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t100\tPASS\tAF=0.3\t"
                    f"GT\t{'0/0' if hom else '0/1'}\n")
            n += 1
        for chrom, pos, gene, ref, alt in RARE:
            w.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t200\tPASS\t"
                    f"AF=0.00004;GENE={gene}\tGT\t0/1\n")
    return n


def main() -> int:
    out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "econ-sample").resolve()
    out.mkdir(parents=True, exist_ok=True)
    vcf = out / "synthetic_wgs.vcf"
    n = build_vcf(vcf)
    print(f"synthetic VCF: {n:,} chip-equivalent + {len(RARE)} rare missense -> {vcf}")
    r = subprocess.run([sys.executable, str(ROOT / "analyze.py"), str(vcf),
                        "--output", str(out / "report.html")], cwd=ROOT)
    if r.returncode:
        return r.returncode
    page = out / "economics.html"
    print(f"\nconsolidated economics page: {page}" if page.exists()
          else "\nWARNING: economics.html not written")
    print("Print it to PDF to reproduce docs/samples/econ-output-sample.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
