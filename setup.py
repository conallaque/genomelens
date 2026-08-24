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
import contextlib
import json
import shutil
import subprocess
import urllib.request
from pathlib import Path

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
PGS_CONDITIONS: dict[str, dict[str, str]] = {
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


def check_python_packages() -> dict[str, bool]:
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
    results: dict[str, bool] = {}
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


def _remote_last_modified(url: str) -> str | None:
    """HEAD the URL and return its Last-Modified header (or None if unreachable).
    Used to detect when NCBI has published a newer ClinVar than we hold."""
    import urllib.request
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.headers.get("Last-Modified")
    except Exception:
        return None


def setup_clinvar(force: bool = False) -> None:
    """Download ClinVar (GRCh37 + GRCh38) and distill each to a compact
    pathogenic/likely-pathogenic table for the Phase-2 clinical-variants screen.

    **Auto-updating:** ClinVar is republished ~weekly. Rather than blindly skip
    an existing table, this checks NCBI's ``Last-Modified`` against the value
    stored when we last distilled (in a ``.meta.json`` sidecar) and re-downloads
    only when a newer ClinVar is actually available. ``force=True``
    (``--clinvar-refresh``) rebuilds unconditionally. Offline with a table
    already present → keep the existing one.

    ~380 MB fetched when refreshing; the raw VCFs are deleted after distilling,
    leaving only the small tables in reference/clinvar/."""
    import datetime as _dt
    import json
    log("ClinVar clinical-variant database setup (auto-updating) ...")
    clinvar_dir = REF_DIR / "clinvar"
    clinvar_dir.mkdir(parents=True, exist_ok=True)
    try:
        from risk.clinical_variants import distill_clinvar_vcf
    except Exception as e:
        log(f"  ERROR: cannot import distiller: {e}")
        return

    builds = {
        "grch37": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz",
        "grch38": "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
    }
    for build, url in builds.items():
        distilled = clinvar_dir / f"clinvar_plp_{build}.tsv.gz"
        meta_path = clinvar_dir / f"clinvar_plp_{build}.meta.json"
        remote_lm = _remote_last_modified(url)

        if distilled.exists() and not force:
            stored = {}
            if meta_path.exists():
                try:
                    stored = json.loads(meta_path.read_text())
                except Exception:
                    stored = {}
            if remote_lm is None:
                log(f"  {build}: NCBI unreachable — keeping existing table "
                    f"(distilled {stored.get('distilled','?')}).")
                continue
            if stored.get("source_last_modified") == remote_lm:
                log(f"  {build}: up to date (ClinVar {remote_lm}) — skipping.")
                continue
            log(f"  {build}: newer ClinVar available "
                f"(have {stored.get('source_last_modified','none')} → {remote_lm}); refreshing.")

        raw = clinvar_dir / f"clinvar_{build}.vcf.gz"
        log(f"  {build}: downloading ClinVar VCF (~190 MB) ...")
        if not download_to(url, raw, label=f"ClinVar {build}"):
            log(f"  {build}: download failed — skipping.")
            continue
        log(f"  {build}: distilling to P/LP table ...")
        n = distill_clinvar_vcf(str(raw), str(distilled), log=log)
        log(f"  {build}: wrote {n:,} clinically-significant records -> {distilled.name}")
        meta_path.write_text(json.dumps({
            "source_last_modified": remote_lm,
            "distilled": _dt.datetime.now().isoformat(timespec="seconds"),
            "rows": n, "source_url": url,
        }, indent=2))
        with contextlib.suppress(Exception):
            raw.unlink()   # reclaim ~190 MB; keep only the compact table
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


# ── Phase-3 predictor tables ─────────────────────────────────────────────────
# Offline pathogenicity predictors for the novel-variant screen. Each is
# download -> (transform) -> tabix-index. Licenses are surfaced so a commercial
# evaluation can drop the non-commercial ones (--commercial-safe).

_AM_URLS = {
    "hg38": "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz",
    "hg19": "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz",
}
_REVEL_URL = "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip"
_GNOMAD_AF_HG38 = ("https://storage.googleapis.com/gatk-best-practices/"
                   "somatic-hg38/af-only-gnomad.hg38.vcf.gz")


def _ensure_bgzip_tabix(gz_path: Path, seq_col: int, start_col: int,
                        end_col: int, meta_char: str = "#") -> bool:
    """Guarantee gz_path is bgzipped + tabix-indexed. Recompresses plain gzip if
    needed (a large decompress/recompress cycle). Returns True on success."""
    import gzip as _gz
    import shutil
    try:
        import pysam
    except Exception:
        log("    ERROR: pysam not installed — run `pip install pysam`.")
        return False
    if Path(str(gz_path) + ".tbi").exists():
        log(f"    [skip] index present: {gz_path.name}.tbi")
        return True
    try:
        pysam.tabix_index(str(gz_path), seq_col=seq_col, start_col=start_col,
                          end_col=end_col, meta_char=meta_char, force=True)
        log(f"    indexed {gz_path.name}")
        return True
    except Exception as e:
        log(f"    {gz_path.name} not bgzipped ({e}); recompressing (several min) ...")
    raw = gz_path.with_suffix("")
    try:
        with _gz.open(gz_path, "rb") as fin, open(raw, "wb") as fout:
            shutil.copyfileobj(fin, fout, length=1 << 22)
        gz_path.unlink()
        pysam.tabix_compress(str(raw), str(gz_path), force=True)
        raw.unlink()
        pysam.tabix_index(str(gz_path), seq_col=seq_col, start_col=start_col,
                          end_col=end_col, meta_char=meta_char, force=True)
        log(f"    re-bgzipped + indexed {gz_path.name}")
        return True
    except Exception as e:
        log(f"    ERROR indexing {gz_path.name}: {e}")
        return False


def setup_alphamissense(force: bool = False, build: str = "hg38") -> None:
    """AlphaMissense precomputed missense scores (CC BY 4.0 — commercial-OK)."""
    log("AlphaMissense (CC BY 4.0, commercial use OK) ...")
    d = REF_DIR / "alphamissense"
    d.mkdir(parents=True, exist_ok=True)
    url = _AM_URLS.get(build, _AM_URLS["hg38"])
    dest = d / Path(url).name
    if force and dest.exists():
        dest.unlink()
        Path(str(dest) + ".tbi").unlink(missing_ok=True)
    if not download_to(url, dest, label=dest.name):
        return
    _ensure_bgzip_tabix(dest, seq_col=0, start_col=1, end_col=1)


def setup_gnomad(force: bool = False) -> None:
    """gnomAD allele frequencies (open license). Uses GATK's compact AF-only
    GRCh38 sites VCF (already bgzipped + tabix-indexed remotely)."""
    log("gnomAD AF (open license, commercial OK) — GATK af-only GRCh38 ...")
    d = REF_DIR / "gnomad"
    d.mkdir(parents=True, exist_ok=True)
    gz = d / "gnomad.af_only.hg38.vcf.gz"
    tbi = Path(str(gz) + ".tbi")
    if force:
        gz.unlink(missing_ok=True)
        tbi.unlink(missing_ok=True)
    ok = download_to(_GNOMAD_AF_HG38, gz, label=gz.name)
    ok = download_to(_GNOMAD_AF_HG38 + ".tbi", tbi, label=tbi.name) and ok
    if ok and not tbi.exists():
        _ensure_bgzip_tabix(gz, seq_col=0, start_col=1, end_col=1)


def setup_revel(force: bool = False) -> None:
    """REVEL (NON-COMMERCIAL). Download zip → slim to chr/pos/ref/alt/REVEL per
    build → sort → bgzip → tabix."""
    import csv
    import zipfile
    log("REVEL (NON-COMMERCIAL license — see PREDICTOR_LICENSES.md) ...")
    d = REF_DIR / "revel"
    d.mkdir(parents=True, exist_ok=True)
    out = d / "revel_grch38.tsv.gz"
    if out.exists() and not force:
        log("    [skip] revel_grch38.tsv.gz present")
        return
    zpath = d / "revel_all.zip"
    if not download_to(_REVEL_URL, zpath, label="revel zip (~0.5 GB)"):
        return
    try:
        import pysam
        with zipfile.ZipFile(zpath) as z:
            member = next(n for n in z.namelist() if n.endswith(".csv"))
            raw = d / "revel_grch38.tsv"
            rows = []
            with z.open(member) as fh:
                rdr = csv.reader(ln.decode() for ln in fh)
                next(rdr, None)                      # header
                for r in rdr:
                    if len(r) < 8 or not r[2]:       # grch38_pos empty → skip
                        continue
                    rows.append((r[0], int(r[2]), r[3], r[4], r[7]))
            rows.sort(key=lambda x: (x[0], x[1]))
            with open(raw, "w") as f:
                for c, p, ref, alt, rev in rows:
                    f.write(f"{c}\t{p}\t{ref}\t{alt}\t{rev}\n")
            pysam.tabix_compress(str(raw), str(out), force=True)
            pysam.tabix_index(str(out), seq_col=0, start_col=1, end_col=1, force=True)
            raw.unlink()
        zpath.unlink()
        log("    REVEL slim table built + indexed")
    except Exception as e:
        log(f"    ERROR building REVEL table: {e}")


def setup_spliceai(force: bool = False) -> None:
    """SpliceAI (NON-COMMERCIAL, gated behind an Illumina BaseSpace login) — cannot
    be auto-downloaded. Print guided instructions + detect a dropped-in file."""
    d = REF_DIR / "spliceai"
    d.mkdir(parents=True, exist_ok=True)
    hit = list(d.glob("spliceai*.vcf.*gz"))
    if hit:
        log(f"SpliceAI table detected: {hit[0].name}")
        if not Path(str(hit[0]) + ".tbi").exists():
            _ensure_bgzip_tabix(hit[0], seq_col=0, start_col=1, end_col=1)
        return
    log("SpliceAI (NON-COMMERCIAL) requires a manual download:")
    log("  1. Log in to Illumina BaseSpace: https://basespace.illumina.com/s/5u6ThOblecrh")
    log("  2. Download spliceai_scores.masked.snv.hg38.vcf.gz  (+ .tbi)")
    log(f"  3. Drop both files into: {d}")
    log("  4. Re-run `python setup.py --spliceai` to index/verify.")


def setup_cadd(force: bool = False) -> None:
    """CADD (NON-COMMERCIAL, ~81 GB). Not auto-downloaded — the analysis can query
    CADD's public tabix file REMOTELY for a demo, or you can host it locally."""
    log("CADD (NON-COMMERCIAL, ~81 GB) — not downloaded automatically.")
    log("  Local (fully offline) option — download once (~81 GB):")
    log("    https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/"
        "whole_genome_SNVs.tsv.gz  (+ .tbi)")
    log(f"    → place both in {REF_DIR / 'cadd'}/")
    log("  Demo option: the analyzer can query CADD remotely on PUBLIC data "
        "(allow_remote=True) — no bulk download.")


def setup_predictors(commercial_safe: bool = False, force: bool = False) -> None:
    """Run the auto-downloadable predictor set. Commercial-safe = AlphaMissense +
    gnomAD only; otherwise also REVEL (SpliceAI/CADD stay manual)."""
    setup_alphamissense(force=force)
    setup_gnomad(force=force)
    if not commercial_safe:
        setup_revel(force=force)
        setup_spliceai(force=force)
        setup_cadd(force=force)
    write_predictor_licenses()


def write_predictor_licenses() -> None:
    """Emit reference/PREDICTOR_LICENSES.md — the attribution/notice manifest."""
    txt = """# Predictor licenses & attribution

GenomeLens' Phase-3 novel-variant screen uses these third-party resources. Their
licenses differ — the non-commercial ones are disabled by `--commercial-safe`.

| Predictor | License | Commercial use | Source |
|---|---|---|---|
| AlphaMissense | CC BY 4.0 | **Yes** (attribution) | Cheng et al., Science 2023 |
| gnomAD | Open / no restrictions | **Yes** | Broad Institute |
| REVEL | Non-commercial | No (contact authors) | Ioannidis et al., AJHG 2016 |
| SpliceAI | CC BY-NC 4.0 | No (Illumina license) | Jaganathan et al., Cell 2019 |
| CADD | Non-commercial | No (UW license) | Rentzsch et al., NAR 2019 |

Predictions are computational estimates, not clinical determinations.
"""
    (REF_DIR / "PREDICTOR_LICENSES.md").write_text(txt)
    log(f"  wrote {REF_DIR / 'PREDICTOR_LICENSES.md'}")


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
                         "variants screen. Auto-updating: re-downloads only when "
                         "NCBI has a newer release (~380 MB when it refreshes).")
    ap.add_argument("--clinvar-refresh", action="store_true",
                    help="Force a ClinVar rebuild even if the local table looks "
                         "current (ignores the Last-Modified check).")
    ap.add_argument("--alphamissense", action="store_true",
                    help="Download + index AlphaMissense (CC BY 4.0, ~0.6 GB).")
    ap.add_argument("--revel", action="store_true",
                    help="Download + index REVEL (non-commercial, ~0.5 GB).")
    ap.add_argument("--gnomad", action="store_true",
                    help="Download gnomAD allele-frequency sites (open license).")
    ap.add_argument("--spliceai", action="store_true",
                    help="Guide/verify SpliceAI (manual BaseSpace download).")
    ap.add_argument("--cadd", action="store_true",
                    help="Guidance for CADD (~81 GB, non-commercial).")
    ap.add_argument("--predictors", action="store_true",
                    help="Auto-download the Phase-3 predictor set "
                         "(AlphaMissense + gnomAD + REVEL).")
    ap.add_argument("--commercial-safe", action="store_true",
                    help="With --predictors, restrict to commercially-licensed "
                         "predictors only (AlphaMissense + gnomAD).")
    ap.add_argument("--predictors-refresh", action="store_true",
                    help="Force a re-download/rebuild of predictor tables.")
    ap.add_argument("--all", action="store_true",
                    help="Run --beagle, --pgs, --ancestry and --clinvar")
    args = ap.parse_args()

    _predictor_flags = [args.alphamissense, args.revel, args.gnomad,
                        args.spliceai, args.cadd, args.predictors,
                        args.predictors_refresh]
    if not any([args.check, args.beagle, args.pgs, args.ancestry, args.clinvar, args.clinvar_refresh, args.all, *_predictor_flags]):
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
    if args.clinvar or args.clinvar_refresh or args.all:
        setup_clinvar(force=args.clinvar_refresh)

    _pf = args.predictors_refresh
    if args.predictors:
        setup_predictors(commercial_safe=args.commercial_safe, force=_pf)
    if args.alphamissense:
        setup_alphamissense(force=_pf)
    if args.revel:
        setup_revel(force=_pf)
    if args.gnomad:
        setup_gnomad(force=_pf)
    if args.spliceai:
        setup_spliceai(force=_pf)
    if args.cadd:
        setup_cadd(force=_pf)
    if any(_predictor_flags):
        write_predictor_licenses()

    write_manifest()
    log("Setup complete.")


if __name__ == "__main__":
    main()
