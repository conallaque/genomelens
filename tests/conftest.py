"""
Shared pytest fixtures.

Every V6 module is designed to be a pure function of its inputs (snps_df,
upstream module results). That makes fixtures simple — we build the inputs
deterministically and reuse them across tests.

Two scopes:
  • `synthetic_snps_df` — small in-memory DataFrame for fast unit tests
  • `fixture_chip_path` — same data on disk in 23andMe-style TSV format, used
                          when a test needs the full parse_dna_file path
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

# Make the project root importable so tests can `import supplements` etc.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


# ── Synthetic genotypes ───────────────────────────────────────────────────────
# Hand-picked to exercise every V6 rule + a representative subset of the core
# modules. Genotypes chosen to produce predictable downstream results:
#   • Heterozygous risk for the headline supplement rules
#   • Mixed power/endurance profile (ACTN3 CT, ACE-D)
#   • Slow caffeine metaboliser, ALDH2 wild-type
#   • Below-average vitamin D genotype
#   • One R1b-typical Y marker so y_haplogroup has something to chew on
#
# The keys MUST stay stable — golden-snapshot tests are byte-exact against
# the resulting JSON. If you change a genotype here, you change every snapshot.

_SYNTHETIC_GENOTYPES: dict[str, dict[str, object]] = {
    # rsid:        (chrom, pos_GRCh37,  genotype)
    # ── V6 supplement rule SNPs ────────────────────────────────────────────
    "rs1801133":  {"chrom": "1",  "pos": 11856378, "genotype": "CT"},   # MTHFR C677T het
    "rs1801131":  {"chrom": "1",  "pos": 11854476, "genotype": "AC"},   # MTHFR A1298C het
    "rs4680":     {"chrom": "22", "pos": 19951271, "genotype": "AG"},   # COMT Val/Met
    "rs10741657": {"chrom": "11", "pos": 14914878, "genotype": "GG"},   # CYP2R1 (low vit D)
    "rs2282679":  {"chrom": "4",  "pos": 72618323, "genotype": "CC"},   # GC/VDBP (low vit D)
    "rs2228570":  {"chrom": "12", "pos": 48272895, "genotype": "TT"},   # VDR FokI active
    "rs602662":   {"chrom": "19", "pos": 49206985, "genotype": "GA"},   # FUT2 (B12 absorption)
    "rs1801198":  {"chrom": "22", "pos": 31019043, "genotype": "GC"},   # TCN2
    "rs855791":   {"chrom": "22", "pos": 37462936, "genotype": "GG"},   # TMPRSS6
    "rs174547":   {"chrom": "11", "pos": 61570783, "genotype": "TC"},   # FADS1 het
    "rs1799941":  {"chrom": "17", "pos":  7533423, "genotype": "GA"},   # SHBG
    "rs1800795":  {"chrom": "7",  "pos": 22766645, "genotype": "GG"},   # IL6 high-CRP het
    "rs1695":     {"chrom": "11", "pos": 67352689, "genotype": "AG"},   # GSTP1
    "rs6721961":  {"chrom": "2",  "pos": 178098964, "genotype": "GT"},  # NRF2
    "rs762551":   {"chrom": "15", "pos": 75041917, "genotype": "AC"},   # CYP1A2 slow-het
    "rs1815739":  {"chrom": "11", "pos": 66328095, "genotype": "CT"},   # ACTN3 R/X
    "rs6265":     {"chrom": "11", "pos": 27679916, "genotype": "GG"},   # BDNF Val/Val
    "rs1801260":  {"chrom": "4",  "pos": 56412708, "genotype": "TT"},   # CLOCK 3111 evening
    "rs9939609":  {"chrom": "16", "pos": 53820527, "genotype": "AT"},   # FTO
    "rs7903146":  {"chrom": "10", "pos": 114758349, "genotype": "CT"},  # TCF7L2
    "rs429358":   {"chrom": "19", "pos": 45411941, "genotype": "TC"},   # APOE ε4 SNP
    "rs7412":     {"chrom": "19", "pos": 45412079, "genotype": "CC"},   # APOE ε4 SNP
    "rs4988235":  {"chrom": "2",  "pos": 136608646, "genotype": "GA"},  # LCT lactase persistence
    "rs671":      {"chrom": "12", "pos": 112241766, "genotype": "GG"},  # ALDH2 wild-type
    "rs1229984":  {"chrom": "4",  "pos": 100239319, "genotype": "GG"},  # ADH1B wild-type
    "rs699":      {"chrom": "1",  "pos": 230845794, "genotype": "AG"},  # AGT M235T
    "rs2187668":  {"chrom": "6",  "pos": 32605884, "genotype": "CC"},   # HLA-DQ2.5 tag (neg)
    "rs7454108":  {"chrom": "6",  "pos": 32814869, "genotype": "AA"},   # HLA-DQ8 tag (neg)
    "rs12722":    {"chrom": "9",  "pos": 137721567, "genotype": "TT"},  # COL5A1 protective
    "rs1800012":  {"chrom": "17", "pos": 48275363, "genotype": "GG"},   # COL1A1 wild-type
    "rs5751876":  {"chrom": "22", "pos": 24827015, "genotype": "TT"},   # ADORA2A
    "rs2032597":  {"chrom": "Y",  "pos": 22719028, "genotype": "T"},    # M9 (haplogroup K)
    "rs9786153":  {"chrom": "Y",  "pos": 13470467, "genotype": "T"},    # R/M207 derived
    "rs9786184":  {"chrom": "Y",  "pos":  2887824, "genotype": "A"},    # R1b-M343 derived
    # ── Carrier-panel anchors ─────────────────────────────────────────────
    "rs1800562":  {"chrom": "6",  "pos": 26093141, "genotype": "GG"},   # HFE C282Y wild-type
    "rs1799945":  {"chrom": "6",  "pos": 26091179, "genotype": "CG"},   # HFE H63D carrier
    "rs6025":     {"chrom": "1",  "pos": 169519049, "genotype": "GG"},  # F5 Leiden wild-type
    "rs1799963":  {"chrom": "11", "pos": 46761055, "genotype": "GG"},   # F2 prothrombin wild-type
    # ── PGx (CYP2C19 *2 / VKORC1 / TPMT — minimal set) ────────────────────
    "rs4244285":  {"chrom": "10", "pos": 96541616, "genotype": "GA"},   # CYP2C19 *2 het
    "rs9923231":  {"chrom": "16", "pos": 31107689, "genotype": "CT"},   # VKORC1 warfarin
    "rs1142345":  {"chrom": "6",  "pos": 18130918, "genotype": "TT"},   # TPMT *3C wild-type
}


@pytest.fixture(scope="session")
def synthetic_snps_df() -> pd.DataFrame:
    """In-memory DataFrame mirroring the structure of `parse_dna_file` output."""
    rows = []
    for rsid, info in _SYNTHETIC_GENOTYPES.items():
        rows.append({
            "rsid": rsid,
            "chrom": info["chrom"],
            "pos": info["pos"],
            "genotype": info["genotype"],
        })
    df = pd.DataFrame(rows).set_index("rsid")
    return df


@pytest.fixture(scope="session")
def fixture_chip_path(tmp_path_factory) -> Path:
    """
    Synthetic 23andMe-style tab-separated chip file written once per session.
    Use when a test needs to exercise the full parse_dna_file path.
    """
    p = tmp_path_factory.mktemp("chip") / "synthetic_genome.txt"
    lines = ["# This data file generated by synthetic test fixture",
             "# rsid\tchromosome\tposition\tgenotype"]
    for rsid, info in _SYNTHETIC_GENOTYPES.items():
        gt = info["genotype"]
        # Pad single-allele Y entries to 2 chars so consumer parsers don't choke
        if info["chrom"] == "Y" and len(gt) == 1:
            gt = gt + gt
        lines.append(f"{rsid}\t{info['chrom']}\t{info['pos']}\t{gt}")
    p.write_text("\n".join(lines) + "\n")
    return p


# ── Golden-snapshot helpers ──────────────────────────────────────────────────

SNAPSHOT_DIR = Path(__file__).parent / "snapshots"


import re as _re

# Pattern matches FHIR cross-resource references like "Observation/a7c39…"
# (where the trailing component is a generated UUID hex). The resource type is
# stable and meaningful; the ID after the slash is volatile.
_FHIR_REF_RE = _re.compile(
    r"^(Patient|Observation|DiagnosticReport|Provenance|Bundle|Practitioner|"
    r"Specimen|MolecularSequence|Organization|RiskAssessment)/[A-Za-z0-9-]+$"
)


def _normalise(obj):
    """Make dict/list comparable across runs: drop volatile keys (timestamps,
    UUIDs) and normalise FHIR resource references to their type only."""
    if isinstance(obj, dict):
        out = {}
        for k, v in obj.items():
            if k in {"id", "fullUrl", "timestamp", "recorded",
                     "occurredDateTime", "effectiveDateTime",
                     "issued", "lastUpdated", "generated"}:
                continue
            # Reference fields: keep the resource type, drop the generated ID
            if k == "reference" and isinstance(v, str) and _FHIR_REF_RE.match(v):
                out[k] = v.split("/")[0] + "/<id>"
                continue
            out[k] = _normalise(v)
        return out
    if isinstance(obj, list):
        return [_normalise(v) for v in obj]
    if isinstance(obj, str) and obj.startswith("urn:uuid:"):
        return "urn:uuid:<volatile>"
    return obj


@pytest.fixture
def assert_snapshot(request):
    """
    Compare a result dict to a JSON snapshot on disk.

    Usage:
        def test_foo(assert_snapshot):
            result = some_function()
            assert_snapshot("test_foo", result)

    Run with `pytest --snapshot-update` to regenerate snapshots intentionally.
    """
    update = request.config.getoption("--snapshot-update", default=False)

    def _compare(name: str, actual: dict) -> None:
        SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
        path = SNAPSHOT_DIR / f"{name}.json"
        actual_norm = _normalise(actual)
        if update or not path.exists():
            path.write_text(json.dumps(actual_norm, indent=2, sort_keys=True,
                                       default=str))
            if update:
                pytest.skip(f"Snapshot updated: {path}")
            return
        expected = json.loads(path.read_text())
        assert actual_norm == expected, (
            f"Snapshot mismatch for {name}.\n"
            f"  Expected: {path}\n"
            f"  Run `pytest --snapshot-update -k {request.node.name}` to accept "
            f"the new output (and review the diff before committing).\n"
        )

    return _compare


def pytest_addoption(parser):
    parser.addoption(
        "--snapshot-update",
        action="store_true",
        default=False,
        help="Regenerate golden snapshots instead of comparing against them.",
    )
