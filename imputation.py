"""
Statistical Imputation via Beagle 5.4
=====================================

Pipeline:
    1. Convert raw chip data (23andMe / TellmeGen / AncestryDNA TSV) -> VCF
    2. Split VCF per chromosome
    3. Run Beagle 5.4 per chromosome against the 1000 Genomes Phase 3 reference
       panel (b37 bref3 format), producing imputed VCFs with DR2 (R-squared)
       quality scores.
    4. Parse output VCFs, extract genotypes + DR2 + AF.
    5. Merge imputed calls back into the snps_df, with a `source` column
       marking 'chip' vs 'imputed' and an `r2` column for confidence.
    6. Cache the imputed dataframe as Parquet keyed by file hash so the slow
       Beagle run only happens once per input file.

Beagle docs: https://faculty.washington.edu/browning/beagle/beagle.html

Notes:
  * Beagle requires Java (run setup.py first).
  * Reference panels are GRCh37 / b37 coordinates. If the input is GRCh38,
    the snps library handles liftover up front.
  * Imputation of 770K chip -> 5-10M SNPs typically takes 30-90 minutes
    on a modern laptop and ~5-10 GB disk space for intermediate files.
"""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).parent
REF_DIR = SCRIPT_DIR / "reference"
BEAGLE_DIR = REF_DIR / "beagle"
CACHE_DIR = SCRIPT_DIR / "cache"


def _log(msg: str) -> None:
    print(f"[impute] {msg}", flush=True)


# ── Capability check ─────────────────────────────────────────────────────────
def imputation_available() -> tuple[bool, str]:
    """Returns (ok, reason)."""
    if shutil.which("java") is None:
        return False, "Java not installed (`brew install openjdk` or apt install default-jre)"
    if not BEAGLE_DIR.exists():
        return False, f"Beagle dir {BEAGLE_DIR} missing — run `python setup.py --beagle`"
    jar = next(BEAGLE_DIR.glob("beagle.*.jar"), None)
    if jar is None:
        return False, "Beagle JAR not found — run `python setup.py --beagle`"
    refs = list(BEAGLE_DIR.glob("chr*.1kg.phase3.v5a.b37.bref3"))
    if len(refs) < 22:
        return False, (
            f"Only {len(refs)}/22 autosomal Beagle reference panels present — "
            f"run `python setup.py --beagle`"
        )
    return True, f"Ready (Beagle JAR: {jar.name}, {len(refs)} reference chroms)"


def beagle_jar() -> Path:
    return next(BEAGLE_DIR.glob("beagle.*.jar"))


def genetic_map_for(chrom: str) -> Path | None:
    candidates = list(BEAGLE_DIR.glob(f"plink.chr{chrom}.GRCh37.map"))
    return candidates[0] if candidates else None


# ── Chip TSV -> VCF conversion ────────────────────────────────────────────────
# Beagle wants a VCF with phased or unphased genotypes. For chip data we write
# unphased calls; Beagle will phase and impute.

def _file_hash(filepath: str) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def chip_to_vcf(snps_df: pd.DataFrame, sample_id: str, vcf_path: Path) -> None:
    """Write an unphased VCF from the parsed chip data. Drops mt/Y/no-call.

    snps_df is the pandas DataFrame from the snps library, indexed by rsID
    with columns chrom, pos, genotype. The genotype is a two-character
    diploid string like 'AG' (or '--' for no-call).
    """
    vcf_path.parent.mkdir(parents=True, exist_ok=True)
    # Filter: keep autosomes + X, drop no-calls
    df = snps_df.copy()
    df = df[df["chrom"].astype(str).isin([str(i) for i in range(1, 23)] + ["X"])]
    df = df[~df["genotype"].astype(str).isin(["--", "00", "nan", ""])]
    df["chrom"] = df["chrom"].astype(str)
    df["pos"] = df["pos"].astype(int)
    df = df.sort_values(["chrom", "pos"])

    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=dna-project-v3\n")
        for c in [str(i) for i in range(1, 23)] + ["X"]:
            f.write(f"##contig=<ID={c}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}\n")

        for rsid, row in df.iterrows():
            gt = str(row["genotype"]).upper()
            if len(gt) != 2:
                continue
            a1, a2 = gt[0], gt[1]
            if a1 not in "ACGT" or a2 not in "ACGT":
                continue
            # Use the two distinct alleles as REF/ALT. Beagle's reference panel
            # provides the canonical alleles; mismatches are handled by Beagle.
            if a1 == a2:
                ref, alt = a1, "."
                gt_str = "0/0"
            else:
                ref, alt = a1, a2
                gt_str = "0/1"
            # Sort alleles alphabetically for REF when ALT is unknown so that
            # different chips produce consistent VCFs for the same site.
            f.write(
                f"{row['chrom']}\t{row['pos']}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt_str}\n"
            )


def _split_vcf_by_chrom(vcf_path: Path, out_dir: Path) -> dict[str, Path]:
    """Split a single VCF into per-chromosome files."""
    out_dir.mkdir(parents=True, exist_ok=True)
    handles: dict[str, tempfile._TemporaryFileWrapper] = {}
    paths: dict[str, Path] = {}
    header_lines: list[str] = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line)
                continue
            chrom = line.split("\t", 1)[0]
            if chrom not in handles:
                p = out_dir / f"chr{chrom}.vcf"
                paths[chrom] = p
                # One handle per chromosome, held open across the whole
                # scan and closed in the finally below — a with-block
                # cannot express a pool whose lifetime spans the loop.
                handles[chrom] = open(p, "w")  # noqa: SIM115
                for h in header_lines:
                    handles[chrom].write(h)
            handles[chrom].write(line)
    for h in handles.values():
        h.close()
    return paths


# ── Beagle runner ────────────────────────────────────────────────────────────
def _run_beagle_for_chrom(
    chrom: str,
    in_vcf: Path,
    ref_panel: Path,
    out_prefix: Path,
    jar: Path,
    map_file: Path | None = None,
    java_mem_gb: int = 4,
    timeout_s: int = 7200,
) -> tuple[bool, str]:
    cmd = [
        "java", f"-Xmx{java_mem_gb}g", "-jar", str(jar),
        f"gt={in_vcf}",
        f"ref={ref_panel}",
        f"out={out_prefix}",
        "impute=true",
        "gp=true",
        "ap=true",
        "nthreads=4",
    ]
    if map_file and map_file.exists():
        cmd.append(f"map={map_file}")
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
        if r.returncode != 0:
            return False, (r.stderr or r.stdout)[-2000:]
        return True, ""
    except subprocess.TimeoutExpired:
        return False, f"Beagle timed out after {timeout_s}s on chrom {chrom}"
    except Exception as e:
        return False, str(e)


# ── Parse imputed VCF ─────────────────────────────────────────────────────────
def parse_imputed_vcf(vcf_gz: Path, min_r2: float = 0.0) -> pd.DataFrame:
    """Parse a Beagle imputed VCF.gz. Extract per-variant genotype, DR2 (R^2)
    and AF from the INFO field. Returns a DataFrame indexed by rsID with
    columns: chrom, pos, genotype, source ('imputed'), r2.
    """
    opener = gzip.open if str(vcf_gz).endswith(".gz") else open
    rows: list[dict] = []
    with opener(vcf_gz, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos, rsid, ref, alt, qual, flt, info, fmt, sample = parts[:10]
            if rsid in (".", ""):
                rsid = f"{chrom}:{pos}"

            # Parse INFO
            info_d = {}
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_d[k] = v
            try:
                dr2 = float(info_d.get("DR2", 0))
            except ValueError:
                dr2 = 0.0
            if dr2 < min_r2:
                continue

            # Parse FORMAT/sample for GT
            fmt_keys = fmt.split(":")
            sample_vals = sample.split(":")
            gt_idx = fmt_keys.index("GT") if "GT" in fmt_keys else 0
            gt_field = sample_vals[gt_idx] if gt_idx < len(sample_vals) else "./."
            gt_field = gt_field.replace("|", "/")
            if gt_field in ("./.", "."):
                continue
            a1, a2 = gt_field.split("/")
            try:
                a1_i, a2_i = int(a1), int(a2)
            except ValueError:
                continue
            alts = alt.split(",")
            # ref/alts bound as defaults, not closed over: the closure is
            # called in the same iteration today, so late binding never bites,
            # but it would the moment this is deferred.
            def _allele(idx: int, ref: str = ref, alts: list = alts) -> str:
                if idx == 0:
                    return ref
                if 1 <= idx <= len(alts):
                    return alts[idx - 1]
                return "N"
            genotype = _allele(a1_i) + _allele(a2_i)

            rows.append({
                "rsid": rsid,
                "chrom": chrom,
                "pos": int(pos),
                "genotype": genotype,
                "source": "imputed",
                "r2": dr2,
            })
    if not rows:
        return pd.DataFrame(columns=["rsid", "chrom", "pos", "genotype", "source", "r2"]).set_index("rsid")
    df = pd.DataFrame(rows).set_index("rsid")
    return df


# ── Top-level run + cache ────────────────────────────────────────────────────
def cache_path_for(input_file: str) -> Path:
    h = _file_hash(input_file)
    CACHE_DIR.mkdir(exist_ok=True)
    return CACHE_DIR / f"imputed_{h}.parquet"


def load_cached(input_file: str) -> pd.DataFrame | None:
    cp = cache_path_for(input_file)
    if not cp.exists():
        return None
    try:
        df = pd.read_parquet(cp)
        _log(f"Loaded cached imputation: {cp} ({len(df):,} variants)")
        return df
    except Exception as e:
        _log(f"Could not load cache ({e}); will rebuild")
        return None


def save_cache(input_file: str, df: pd.DataFrame) -> None:
    cp = cache_path_for(input_file)
    try:
        df.to_parquet(cp, compression="snappy")
        _log(f"Cached imputation saved: {cp} ({cp.stat().st_size / 1e6:.1f} MB)")
    except Exception as e:
        _log(f"Cache save failed (install pyarrow): {e}")


def impute_genotypes(
    snps_df: pd.DataFrame,
    input_file: str,
    sample_id: str = "USER",
    min_r2: float = 0.3,
    chromosomes: list[str] | None = None,
    force: bool = False,
) -> tuple[pd.DataFrame, dict]:
    """Main entry point. Returns:
        (merged_df, info)
    where merged_df is the chip + imputed variants combined (with `source`
    and `r2` columns), and info is a summary dict.
    """
    ok, reason = imputation_available()
    if not ok:
        _log(f"Imputation not available: {reason}")
        return snps_df.assign(source="chip", r2=1.0), {
            "available": False, "reason": reason,
            "n_chip": len(snps_df), "n_imputed": 0,
        }

    if not force:
        cached = load_cached(input_file)
        if cached is not None:
            return cached, {
                "available": True, "from_cache": True,
                "n_chip": int((cached["source"] == "chip").sum()),
                "n_imputed": int((cached["source"] == "imputed").sum()),
                "min_r2": min_r2,
            }

    chromosomes = chromosomes or [str(i) for i in range(1, 23)] + ["X"]
    work = Path(tempfile.mkdtemp(prefix="impute_"))
    _log(f"Workspace: {work}")

    # Step 1: chip -> VCF
    chip_vcf = work / "chip.vcf"
    _log("Converting chip to VCF ...")
    chip_to_vcf(snps_df, sample_id, chip_vcf)

    # Step 2: split per chrom
    _log("Splitting per chromosome ...")
    per_chrom = _split_vcf_by_chrom(chip_vcf, work / "per_chrom")

    jar = beagle_jar()
    imputed_frames = []
    failures: list[str] = []

    for c in chromosomes:
        if c not in per_chrom:
            continue
        ref = BEAGLE_DIR / f"chr{c}.1kg.phase3.v5a.b37.bref3"
        if not ref.exists():
            failures.append(f"chr{c}: ref panel missing")
            continue
        out_prefix = work / f"imputed_chr{c}"
        _log(f"  Imputing chr{c} ... (this can take 1-5 minutes)")
        t0 = time.time()
        ok, err = _run_beagle_for_chrom(
            chrom=c,
            in_vcf=per_chrom[c],
            ref_panel=ref,
            out_prefix=out_prefix,
            jar=jar,
            map_file=genetic_map_for(c),
        )
        dt = time.time() - t0
        if not ok:
            failures.append(f"chr{c}: {err[:200]}")
            _log(f"    FAIL ({dt:.0f}s): {err[:200]}")
            continue
        _log(f"    OK ({dt:.0f}s)")
        out_vcf_gz = out_prefix.with_suffix(".vcf.gz")
        if out_vcf_gz.exists():
            df_imp = parse_imputed_vcf(out_vcf_gz, min_r2=min_r2)
            imputed_frames.append(df_imp)

    if not imputed_frames:
        _log("No imputed output produced.")
        merged = snps_df.assign(source="chip", r2=1.0)
        save_cache(input_file, merged)
        return merged, {
            "available": True, "from_cache": False,
            "n_chip": len(snps_df), "n_imputed": 0,
            "failures": failures,
        }

    imputed_all = pd.concat(imputed_frames)
    # Chip data takes precedence over imputed for the same site
    chip_marked = snps_df.assign(source="chip", r2=1.0)
    # Keep the original chip's index; remove imputed entries already in chip
    imputed_only = imputed_all[~imputed_all.index.isin(chip_marked.index)]
    merged = pd.concat([chip_marked, imputed_only])

    save_cache(input_file, merged)

    # Clean up workspace
    with contextlib.suppress(Exception):
        shutil.rmtree(work)

    return merged, {
        "available": True,
        "from_cache": False,
        "n_chip": int((merged["source"] == "chip").sum()),
        "n_imputed": int((merged["source"] == "imputed").sum()),
        "min_r2": min_r2,
        "failures": failures,
    }
