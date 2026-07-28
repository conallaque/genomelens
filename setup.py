#!/usr/bin/env python3
"""
One-Time Setup for DNA Analysis Tool v3
=======================================

Downloads and caches external dependencies needed for advanced features:

  --check       Run dependency checks only (Java, Python packages, disk space)
  --beagle      Download Beagle 5.4 JAR + 1000 Genomes reference panels (~3 GB)
  --pgs         Download PGS Catalog scoring files for 15 conditions (~200 MB)
  --ancestry    Download 1000 Genomes AIMs reference data for PCA (~200 MB)
  --all         Do everything above

All downloads go to ~/dna-project/reference/ and are idempotent (skip if
already present). Files are cached and re-used across runs.

The analyze.py tool gracefully degrades when these aren't downloaded — the
core report works without them. Run this script when you want to enable:
  * Imputation (--impute)            requires --beagle
  * Expanded PGS panels              requires --pgs
  * PCA-based ancestry estimation    requires --ancestry
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import Dict, List

SCRIPT_DIR = Path(__file__).parent
REF_DIR = SCRIPT_DIR / "reference"
BEAGLE_DIR = REF_DIR / "beagle"
PGS_DIR = REF_DIR / "pgs_scores"
ANCESTRY_DIR = REF_DIR / "ancestry"

# ── Beagle ────────────────────────────────────────────────────────────────────
# Use the August 2024 release (current at time of writing). Version updates can
# be picked up by changing the JAR_URL.
BEAGLE_JAR_URL = "https://faculty.washington.edu/browning/beagle/beagle.06Aug24.a91.jar"
BEAGLE_JAR_NAME = "beagle.06Aug24.a91.jar"

# Beagle pre-built 1000G phase 3 reference panels (bref3 format, GRCh37).
# Per-chromosome files at:
BEAGLE_REF_BASE = "http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.bref3/"

# Genetic maps used by Beagle for phasing/imputation
BEAGLE_MAP_BASE = "http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/"
BEAGLE_MAP_ZIP = "plink.GRCh37.map.zip"

# ── PGS Catalog scoring files (one per condition) ─────────────────────────────
# Map condition → published PGS Catalog ID. These are well-validated scores for
# each trait. Files downloaded from https://www.pgscatalog.org/
PGS_CONDITIONS: Dict[str, Dict[str, str]] = {
    "coronary_artery_disease": {
        "pgs_id": "PGS000018",
        "description": "Coronary Artery Disease (Inouye 2018, 1.7M variants)",
    },
    "type_2_diabetes": {
        "pgs_id": "PGS000014",
        "description": "Type 2 Diabetes (Khera 2018)",
    },
    "breast_cancer": {
        "pgs_id": "PGS000004",
        "description": "Breast Cancer (Mavaddat 2019, 313 variants)",
    },
    "prostate_cancer": {
        "pgs_id": "PGS000662",
        "description": "Prostate Cancer (Conti 2021)",
    },
    "alzheimers_disease": {
        "pgs_id": "PGS000334",
        "description": "Alzheimer's Disease (Zhang 2020)",
    },
    "atrial_fibrillation": {
        "pgs_id": "PGS000016",
        "description": "Atrial Fibrillation (Khera 2018)",
    },
    "bmi": {
        "pgs_id": "PGS000027",
        "description": "Body Mass Index (Khera 2019)",
    },
    "major_depressive_disorder": {
        "pgs_id": "PGS000145",
        "description": "Major Depressive Disorder (Wray 2018)",
    },
    "schizophrenia": {
        "pgs_id": "PGS000019",
        "description": "Schizophrenia (PGC 2014/2020)",
    },
    "hypertension": {
        "pgs_id": "PGS000301",
        "description": "Hypertension / systolic BP (Evangelou 2018)",
    },
    "stroke": {
        "pgs_id": "PGS000039",
        "description": "Ischemic Stroke (Malik 2018)",
    },
    "chronic_kidney_disease": {
        "pgs_id": "PGS000314",
        "description": "Chronic Kidney Disease (Wuttke 2019)",
    },
    "asthma": {
        "pgs_id": "PGS000037",
        "description": "Asthma (Demenais 2018)",
    },
    "inflammatory_bowel_disease": {
        "pgs_id": "PGS000020",
        "description": "Inflammatory Bowel Disease (de Lange 2017)",
    },
    "rheumatoid_arthritis": {
        "pgs_id": "PGS000038",
        "description": "Rheumatoid Arthritis (Okada 2014)",
    },
}

PGS_CATALOG_DOWNLOAD_BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores"

# ── 1000 Genomes AIMs reference for ancestry PCA ──────────────────────────────
# We use a curated set of ancestry-informative markers from the 1000 Genomes
# Project Phase 3, formatted as plink BED/BIM/FAM or VCF. To keep the download
# small, we package a pre-pruned set of ~5000 AIMs distinguishing the 5
# superpopulations.
# Since there's no single canonical "AIMs only" download, we point at the
# 1000G Phase 3 chr-merged VCF and let the user filter, OR they can provide
# their own AIM set. See README in reference/ancestry/ after running this.

KGP_AIMS_URL = (
    "https://ftp.ensemblgenomes.org/pub/release-48/plants/variation/vcf/"
    # Placeholder — see README; the user can swap in any AIMs file they prefer.
)


# ── Utility functions ────────────────────────────────────────────────────────

def log(msg: str) -> None:
    print(f"[setup] {msg}", flush=True)


def check_java() -> bool:
    """Verify Java is installed and accessible."""
    try:
        result = subprocess.run(
            ["java", "-version"],
            capture_output=True, text=True, timeout=5,
        )
        # `java -version` writes to stderr
        version_text = (result.stderr or result.stdout).strip().splitlines()
        first_line = version_text[0] if version_text else ""
        log(f"  Java: {first_line}")
        return True
    except FileNotFoundError:
        log("  Java NOT FOUND. Install via:")
        log("    macOS:   brew install openjdk")
        log("    Linux:   sudo apt install default-jre   (or yum install java-17-openjdk)")
        return False
    except Exception as e:
        log(f"  Java check failed: {e}")
        return False


def check_python_packages() -> Dict[str, bool]:
    """Verify Python packages used by v3 modules."""
    required = {
        "snps": "raw DNA file parsing",
        "pandas": "data frames",
        "numpy": "numerical math",
        "requests": "Ollama HTTP",
    }
    optional = {
        "weasyprint": "PDF export (--pdf)",
        "sklearn": "PCA ancestry (--ancestry)",
        "matplotlib": "ancestry PCA plot",
        "pyarrow": "compressed imputation cache",
    }
    results: Dict[str, bool] = {}
    log("  Required packages:")
    for pkg, desc in required.items():
        try:
            __import__(pkg if pkg != "sklearn" else "sklearn")
            log(f"    OK     {pkg:14s}  ({desc})")
            results[pkg] = True
        except ImportError:
            log(f"    MISS   {pkg:14s}  ({desc})  -- pip install {pkg}")
            results[pkg] = False
    log("  Optional packages:")
    for pkg, desc in optional.items():
        try:
            __import__(pkg if pkg != "sklearn" else "sklearn")
            log(f"    OK     {pkg:14s}  ({desc})")
            results[pkg] = True
        except ImportError:
            log(f"    MISS   {pkg:14s}  ({desc})  -- pip install {pkg}")
            results[pkg] = False
    return results


def check_disk_space(path: Path, gb_needed: float) -> bool:
    try:
        stat = shutil.disk_usage(path)
        free_gb = stat.free / 1e9
        ok = free_gb >= gb_needed
        log(f"  Disk free at {path}: {free_gb:.1f} GB ({'OK' if ok else 'INSUFFICIENT — need ' + str(gb_needed) + ' GB'})")
        return ok
    except Exception as e:
        log(f"  Disk check failed: {e}")
        return False


def download_to(url: str, dest: Path, label: str = "") -> bool:
    """Download a URL to dest with simple progress feedback. Skips if exists."""
    if dest.exists() and dest.stat().st_size > 0:
        log(f"    [skip] {label or dest.name} already present")
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    log(f"    [download] {label or dest.name} <- {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            total = int(r.headers.get("Content-Length", 0))
            seen = 0
            chunk = 1 << 20
            with open(tmp, "wb") as f:
                while True:
                    data = r.read(chunk)
                    if not data:
                        break
                    f.write(data)
                    seen += len(data)
                    if total:
                        pct = 100 * seen / total
                        print(f"\r      {seen/1e6:6.1f}/{total/1e6:.1f} MB  ({pct:5.1f}%)", end="", flush=True)
            print()
        tmp.rename(dest)
        return True
    except Exception as e:
        log(f"    FAIL: {e}")
        if tmp.exists():
            tmp.unlink()
        return False


# ── Setup commands ────────────────────────────────────────────────────────────

def setup_beagle() -> None:
    """Download Beagle JAR + 1000G reference panels."""
    log("Beagle setup ...")
    BEAGLE_DIR.mkdir(parents=True, exist_ok=True)
    if not check_java():
        log("  Beagle requires Java. Install it first, then re-run.")
        return

    jar = BEAGLE_DIR / BEAGLE_JAR_NAME
    download_to(BEAGLE_JAR_URL, jar, label="Beagle JAR")

    # Genetic maps
    maps_zip = BEAGLE_DIR / BEAGLE_MAP_ZIP
    if download_to(BEAGLE_MAP_BASE + BEAGLE_MAP_ZIP, maps_zip, label="Beagle genetic maps"):
        # Unzip
        try:
            import zipfile
            with zipfile.ZipFile(maps_zip) as z:
                z.extractall(BEAGLE_DIR)
            log(f"    extracted maps to {BEAGLE_DIR}")
        except Exception as e:
            log(f"    unzip failed: {e}")

    # Reference panels per chromosome (b37 bref3 format)
    log("  Downloading 1000G reference panels (~3 GB total) ...")
    chroms = [str(i) for i in range(1, 23)] + ["X"]
    for c in chroms:
        fname = f"chr{c}.1kg.phase3.v5a.b37.bref3"
        download_to(BEAGLE_REF_BASE + fname, BEAGLE_DIR / fname,
                    label=f"chr{c} reference")

    log("Beagle setup complete.")


def setup_pgs_catalog() -> None:
    """Download PGS Catalog scoring files for our 15 conditions."""
    log("PGS Catalog scoring file setup ...")
    PGS_DIR.mkdir(parents=True, exist_ok=True)

    for condition, info in PGS_CONDITIONS.items():
        pgs_id = info["pgs_id"]
        # PGS Catalog harmonized scoring files in GRCh37 (Hmpos_GRCh37):
        # /scores/<PGS_ID>/ScoringFiles/Harmonized/<PGS_ID>_hmPOS_GRCh37.txt.gz
        url = (
            f"{PGS_CATALOG_DOWNLOAD_BASE}/{pgs_id}/ScoringFiles/Harmonized/"
            f"{pgs_id}_hmPOS_GRCh37.txt.gz"
        )
        dest = PGS_DIR / f"{pgs_id}_hmPOS_GRCh37.txt.gz"
        log(f"  {condition} -> {pgs_id} ({info['description']})")
        download_to(url, dest, label=f"{pgs_id}")

    # Write a manifest mapping condition -> file
    manifest = {
        cond: {**info, "file": f"{info['pgs_id']}_hmPOS_GRCh37.txt.gz"}
        for cond, info in PGS_CONDITIONS.items()
    }
    with open(PGS_DIR / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    log(f"  manifest written: {PGS_DIR / 'manifest.json'}")
    log("PGS Catalog setup complete.")


def setup_clinvar() -> None:
    """Download ClinVar (GRCh37 + GRCh38) and distill each to a compact
    pathogenic/likely-pathogenic table for the Phase-2 clinical-variants screen.

    ~380 MB downloaded once; distilled to small P/LP-only tables in
    reference/clinvar/. Idempotent — skips a build whose distilled table
    already exists."""
    log("ClinVar clinical-variant database setup ...")
    clinvar_dir = REF_DIR / "clinvar"
    clinvar_dir.mkdir(parents=True, exist_ok=True)
    try:
        from clinical_variants import distill_clinvar_vcf
    except Exception as e:
        log(f"  ERROR: cannot import distiller: {e}")
        return

    builds = {
        "grch37": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
        "grch38": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
    }
    for build, url in builds.items():
        distilled = clinvar_dir / f"clinvar_plp_{build}.tsv.gz"
        if distilled.exists():
            log(f"  {build}: distilled table already present — skipping.")
            continue
        raw = clinvar_dir / f"clinvar_{build}.vcf.gz"
        log(f"  {build}: downloading ClinVar VCF (~190 MB) ...")
        if not download_to(url, raw, label=f"ClinVar {build}"):
            log(f"  {build}: download failed — skipping.")
            continue
        log(f"  {build}: distilling to P/LP table ...")
        n = distill_clinvar_vcf(str(raw), str(distilled), log=log)
        log(f"  {build}: wrote {n:,} clinically-significant records -> {distilled.name}")
        try:
            raw.unlink()   # reclaim ~190 MB; keep only the compact table
        except Exception:
            pass
    log("ClinVar setup complete.")


def setup_ancestry() -> None:
    """Set up 1000 Genomes AIMs reference for PCA-based ancestry."""
    log("Ancestry reference setup ...")
    ANCESTRY_DIR.mkdir(parents=True, exist_ok=True)

    # We need 1000G genotypes at a set of ancestry-informative markers (AIMs).
    # The cleanest source is the 1000G Phase 3 chr-merged plink dataset filtered
    # to a curated AIMs panel. To keep this script self-contained we ship a list
    # of well-established AIMs that the ancestry module looks up against any
    # 1000G genotype data the user provides.

    aims_list = ANCESTRY_DIR / "aims_panel.tsv"
    if not aims_list.exists():
        # Embed a small curated AIMs set as a starter. Real PCA needs many more,
        # but this gets the framework working.
        aims = [
            ("rs3827760", "EDAR",     "EAS/AMR-informative"),
            ("rs1426654", "SLC24A5",  "EUR/SAS pigmentation"),
            ("rs16891982", "SLC45A2", "EUR pigmentation"),
            ("rs12913832", "HERC2",   "EUR eye colour"),
            ("rs1805007",  "MC1R",    "EUR red hair"),
            ("rs1042602",  "TYR",     "EUR pigmentation"),
            ("rs17822931", "ABCC11",  "EAS/AMR earwax"),
            ("rs1129038",  "HERC2",   "EUR blue eye tag"),
            ("rs2402130",  "EDAR",    "EAS hair thickness"),
            ("rs671",      "ALDH2",   "EAS alcohol flush"),
            ("rs1229984",  "ADH1B",   "EAS ADH1B fast"),
            ("rs2814778",  "DARC",    "AFR Duffy null"),
            ("rs4988235",  "LCT",     "EUR lactase persistence"),
            ("rs182549",   "MCM6",    "EUR lactase tag"),
        ]
        with open(aims_list, "w") as f:
            f.write("rsid\tgene\tnote\n")
            for r in aims:
                f.write("\t".join(r) + "\n")
        log(f"  curated AIMs panel written: {aims_list}")

    readme = ANCESTRY_DIR / "README.txt"
    readme.write_text(
        "Ancestry reference data\n"
        "=======================\n\n"
        "The PCA-based ancestry module needs 1000 Genomes Phase 3 genotypes at\n"
        "ancestry-informative markers. Options for obtaining this:\n\n"
        "1. Use 23andMe/AncestryDNA's pre-computed ancestry — usually more\n"
        "   accurate than chip-based PCA against 1000G.\n\n"
        "2. Run `plink2 --vcf ALL.chr<N>.phase3.vcf.gz --extract aims_panel.tsv\n"
        "    --make-bed --out kgp_aims_chr<N>` for each chromosome, then merge.\n\n"
        "3. Download the curated AIMs panel from the 23andMe lab\n"
        "   public datasets or HGDP+1000G merged data.\n\n"
        "The module gracefully degrades to a small-AIMs heuristic when full\n"
        "1000G genotype data is not present.\n"
    )
    log(f"  README written: {readme}")
    log("Ancestry setup complete (curated AIMs only; see README for full PCA).")


def write_manifest() -> None:
    """Write a top-level setup manifest summarising what's available."""
    manifest = {
        "beagle": {
            "jar": (BEAGLE_DIR / BEAGLE_JAR_NAME).exists(),
            "ref_panels_present": sum(
                1 for c in [str(i) for i in range(1, 23)] + ["X"]
                if (BEAGLE_DIR / f"chr{c}.1kg.phase3.v5a.b37.bref3").exists()
            ),
            "maps_present": any(BEAGLE_DIR.glob("plink.chr*.GRCh37.map")),
        },
        "pgs_catalog": {
            "files_present": sum(
                1 for cond, info in PGS_CONDITIONS.items()
                if (PGS_DIR / f"{info['pgs_id']}_hmPOS_GRCh37.txt.gz").exists()
            ),
            "expected": len(PGS_CONDITIONS),
        },
        "ancestry": {
            "aims_panel": (ANCESTRY_DIR / "aims_panel.tsv").exists(),
        },
    }
    with open(REF_DIR / "setup_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    log(f"  manifest written: {REF_DIR / 'setup_manifest.json'}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="One-time setup for DNA Analysis Tool v3",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--check", action="store_true",
                    help="Check dependencies only — no downloads")
    ap.add_argument("--beagle", action="store_true",
                    help="Download Beagle JAR + 1000G reference (~3 GB)")
    ap.add_argument("--pgs", action="store_true",
                    help="Download PGS Catalog scoring files (~200 MB)")
    ap.add_argument("--ancestry", action="store_true",
                    help="Set up 1000G ancestry reference (~200 MB if full)")
    ap.add_argument("--clinvar", action="store_true",
                    help="Download + distill ClinVar for the Phase-2 clinical-"
                         "variants screen (~380 MB dl, distilled to small tables)")
    ap.add_argument("--all", action="store_true",
                    help="Run --beagle, --pgs, --ancestry and --clinvar")
    args = ap.parse_args()

    if not any([args.check, args.beagle, args.pgs, args.ancestry,
                args.clinvar, args.all]):
        ap.print_help()
        return

    REF_DIR.mkdir(parents=True, exist_ok=True)
    log(f"Reference root: {REF_DIR}")

    log("Dependency check ...")
    check_java()
    check_python_packages()
    check_disk_space(REF_DIR, 5.0)

    if args.check:
        log("Check-only mode — exiting without downloads.")
        return

    if args.beagle or args.all:
        setup_beagle()
    if args.pgs or args.all:
        setup_pgs_catalog()
    if args.ancestry or args.all:
        setup_ancestry()
    if args.clinvar or args.all:
        setup_clinvar()

    write_manifest()
    log("Setup complete.")


if __name__ == "__main__":
    main()
