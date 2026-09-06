"""Unit tests for novel_variants.py (Phase 3 — computational predictor screen).

Builds a real synthetic tabix-indexed AlphaMissense table in tmp_path and runs the
full annotate → classify path, plus the graceful-degradation and build-gate paths.
"""
from __future__ import annotations

import pytest

from risk import novel_variants as nv

pysam = pytest.importorskip("pysam")


# ── fixtures ──────────────────────────────────────────────────────────────────

def _make_am_table(root, rows):
    """rows: (chrom, pos, ref, alt, am_pathogenicity, am_class). Writes a sorted,
    bgzipped, tabix-indexed AlphaMissense-shaped table under root/alphamissense/."""
    d = root / "alphamissense"
    d.mkdir(parents=True, exist_ok=True)
    raw = d / "AlphaMissense_hg38.tsv"
    rows = sorted(rows, key=lambda r: (r[0], int(r[1])))
    with open(raw, "w") as f:
        for c, p, ref, alt, am, cls in rows:
            f.write(f"{c}\t{p}\t{ref}\t{alt}\thg38\tU1\tT1\tp.X\t{am}\t{cls}\n")
    gz = str(raw) + ".gz"
    pysam.tabix_compress(str(raw), gz, force=True)
    pysam.tabix_index(gz, seq_col=0, start_col=1, end_col=1, force=True)
    return d


def _write_vcf(root, variants, build_contig="248956422"):
    """variants: (chrom, pos, ref, alt, gt). Writes a minimal single-sample VCF."""
    p = root / "user.vcf"
    with open(p, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID=1,length={build_contig}>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        for c, pos, ref, alt, gt in variants:
            f.write(f"{c}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")
    return str(p)


# ── tests ─────────────────────────────────────────────────────────────────────

def test_predicted_pathogenic_rare(tmp_path, monkeypatch):
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [
        ("1", "1000", "C", "T", "0.98", "likely_pathogenic"),
        ("1", "2000", "G", "A", "0.10", "likely_benign"),
        ("1", "3000", "A", "G", "0.55", "ambiguous"),
    ])
    vcf = _write_vcf(tmp_path, [
        ("1", "1000", "C", "T", "0/1"),   # damaging → flagged
        ("1", "2000", "G", "A", "1/1"),   # benign → dropped
        ("1", "3000", "A", "G", "0/1"),   # ambiguous → ambiguous bucket
        ("1", "9999", "T", "C", "0/1"),   # not in table → no annotation
    ])
    r = nv.analyze_novel_variants(vcf, "grch38")
    assert r["available"] is True
    assert r["n_predicted_pathogenic"] == 1        # the ambiguous one is not "pathogenic"
    assert r["n_ambiguous"] == 1
    dmg = r["buckets"]["predicted_pathogenic_rare"]
    assert len(dmg) == 1 and dmg[0]["pos"] == 1000
    assert dmg[0]["am_class"] == "likely_pathogenic"
    assert dmg[0]["confidence"] == "higher"        # am ≥ 0.9
    assert dmg[0]["rarity"] == "unknown"           # no gnomAD table present
    assert "AlphaMissense" in dmg[0]["evidence"]


def test_dedup_against_clinvar(tmp_path, monkeypatch):
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [("1", "1000", "C", "T", "0.97", "likely_pathogenic")])
    vcf = _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/1")])
    clinvar = {"findings": [{"chrom": "1", "pos": 1000, "ref": "C", "alt": "T"}]}
    r = nv.analyze_novel_variants(vcf, "grch38", clinvar_result=clinvar)
    assert r["available"] is True
    assert r["n_predicted_pathogenic"] == 0        # already a ClinVar hit → skipped


def test_dedup_handles_chip_result_without_findings_key(tmp_path, monkeypatch):
    # On chip input the clinvar dict is {"available": False, ...} with NO findings key.
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [("1", "1000", "C", "T", "0.97", "likely_pathogenic")])
    vcf = _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/1")])
    r = nv.analyze_novel_variants(vcf, "grch38",
                                  clinvar_result={"available": False, "reason": "x"})
    assert r["available"] is True
    assert r["n_predicted_pathogenic"] == 1        # no crash on missing 'findings'


def test_graceful_degrade_no_tables(tmp_path, monkeypatch):
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)  # empty → no predictor tables
    vcf = _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/1")])
    r = nv.analyze_novel_variants(vcf, "grch38")
    assert r["available"] is False
    assert "setup.py" in r["reason"]
    assert r["negative_disclaimer"] and r["disclaimer"]


def test_build_gate(tmp_path, monkeypatch):
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [("1", "1000", "C", "T", "0.97", "likely_pathogenic")])
    vcf = _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/1")])
    r = nv.analyze_novel_variants(vcf, "unknown")
    assert r["available"] is False
    assert "assume-build" in r["reason"] or "Build" in r["reason"]


# ── predictor-table build guard ───────────────────────────────────────────────
#
# THE DEFECT THESE PIN. `resolve()` took the first glob hit, and the `build`
# argument was checked for membership in ("grch37","grch38") and then never used
# again — so a GRCh37 genome was scored against whichever table happened to be on
# disk, and the only table setup.py downloaded by default was hg38. A lookup
# against the wrong build is not a degraded answer: it returns None where the
# coordinate is unused, and some OTHER variant's score where the coordinate is
# also valid in the other build. Nothing in the result said so.

def test_table_build_reads_the_build_from_a_filename():
    assert nv.table_build("AlphaMissense_hg38.tsv.gz") == "grch38"
    assert nv.table_build("AlphaMissense_hg19.tsv.gz") == "grch37"
    assert nv.table_build("gnomad.af_only.hg38.vcf.gz") == "grch38"
    # Declaring no build is "unverifiable", deliberately distinguished from
    # declaring the wrong one — the first stays eligible, the second is fatal.
    assert nv.table_build("revel_all_chromosomes.tsv.gz") is None


def test_predictor_refuses_a_table_keyed_on_the_other_build(tmp_path, monkeypatch):
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [("1", "1000", "C", "T", "0.97", "likely_pathogenic")])
    vcf = _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/1")])

    # hg38 table on disk, GRCh37 genome in hand: refuse rather than score.
    r = nv.analyze_novel_variants(vcf, "grch37")
    assert r["available"] is False, "a wrong-build table must not be scored"
    assert "build" in r["reason"].lower()
    assert "grch38" in r["reason"], "the reason must name the table's own build"
    assert "setup.py" in r["reason"], "and say how to get the right table"

    # Same file, same table, on the build they belong to: still works. So the
    # guard rejects the mismatch, not the table.
    ok = nv.analyze_novel_variants(vcf, "grch38")
    assert ok["available"] is True
    assert ok["n_predicted_pathogenic"] == 1


def test_matching_build_table_is_selected_over_a_mismatched_one(tmp_path, monkeypatch):
    """With both tables present, selection must be by build, not sort order.

    'hg19' sorts before 'hg38', so the previous first-glob-hit behavior would
    have handed a GRCh38 genome the hg19 table.
    """
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [("1", "1000", "C", "T", "0.97", "likely_pathogenic")])
    d = tmp_path / "alphamissense"
    raw = d / "AlphaMissense_hg19.tsv"
    raw.write_text("1\t1000\tC\tT\thg19\tU1\tT1\tp.X\t0.99\tlikely_pathogenic\n")
    gz = str(raw) + ".gz"
    pysam.tabix_compress(str(raw), gz, force=True)
    pysam.tabix_index(gz, seq_col=0, start_col=1, end_col=1, force=True)

    for build in ("grch38", "grch37"):
        r = nv.analyze_novel_variants(
            _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/1")]), build)
        assert r["available"] is True, f"{build}: {r.get('reason')}"
        used = [p for p in r["predictors_used"] if p["name"] == "AlphaMissense"]
        assert used, f"{build}: AlphaMissense not used"
        assert used[0]["table_build"] == build, (
            f"input {build} must select the {build} table, got "
            f"{used[0]['table_build']}")


def test_input_build_is_reported_on_the_result(tmp_path, monkeypatch):
    # The screen has to say which build it ran on, so a reader can tell whether
    # the coordinates it reports mean what they appear to mean.
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [("1", "1000", "C", "T", "0.97", "likely_pathogenic")])
    vcf = _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/1")])
    r = nv.analyze_novel_variants(vcf, "grch38")
    assert r["input_build"] == "grch38"
    assert r["dropped_build_mismatch"] == []


def test_commercial_safe_drops_noncommercial(tmp_path, monkeypatch):
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [("1", "1000", "C", "T", "0.97", "likely_pathogenic")])
    vcf = _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/1")])
    r = nv.analyze_novel_variants(vcf, "grch38", commercial_safe=True)
    assert r["available"] is True
    used = {p["name"] for p in r["predictors_used"]}
    assert used == {"AlphaMissense"}               # only commercial-OK survives
    assert set(r["dropped_noncommercial"]) >= {"REVEL", "SpliceAI", "CADD"}


def test_hemizygous_and_homozygous_carried(tmp_path, monkeypatch):
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [("1", "1000", "C", "T", "0.95", "likely_pathogenic")])
    # not carried (0/0) → must be skipped
    vcf = _write_vcf(tmp_path, [("1", "1000", "C", "T", "0/0")])
    r = nv.analyze_novel_variants(vcf, "grch38")
    assert r["n_predicted_pathogenic"] == 0


# ── gene assignment for predicted variants ────────────────────────────────────
#
# A predicted variant used to arrive at the economics with no gene, so every one
# routed to the same generic bucket regardless of where it landed. The predictor
# knew which protein it had scored and nothing asked.

def test_uniprot_accession_resolves_to_a_gene():
    assert nv.gene_for_uniprot("P38398") == "BRCA1"
    assert nv.gene_for_uniprot("p04637") == "TP53"        # case-insensitive
    # Canonical, not an isoform. First-hit derivation picked Q7L775 for MLH1;
    # majority vote over ~60 supporting coordinates gives the canonical entry.
    assert nv.gene_for_uniprot("P40692") == "MLH1"
    # NAMING IS NOT ANCHORING. This used to assert P05091 -> "", on the
    # reasoning that a gene with no economic anchor should resolve to nothing.
    # That conflated two jobs. ALDH2 is simply the correct name for P05091, and
    # withholding it does not protect anything — it just produced unnamed rows.
    # On a whole genome that was the whole story: 40 of 40 predictions came back
    # unnamed and collapsed into one "predicted-pathogenic ?" row. The property
    # that actually matters is that an unanchored gene is not PRICED, and that
    # is enforced downstream by _gene_to_econ, asserted just below.
    assert nv.gene_for_uniprot("P05091") == "ALDH2"
    assert nv.gene_for_uniprot("P38398-2") == "BRCA1"     # isoform suffix
    assert nv.gene_for_uniprot("") == ""
    # Nonsense accessions still resolve to nothing rather than a near neighbor.
    assert nv.gene_for_uniprot("NOT_AN_ACCESSION") == ""

    from econ.value_of_information import _gene_to_econ
    assert _gene_to_econ("ALDH2")[0] is None or not _gene_to_econ("ALDH2")[0]
    assert _gene_to_econ("BRCA1")[0]                      # anchored, still priced


def test_every_anchor_gene_has_an_accession():
    """The map must cover exactly the genes that can be priced or withheld.

    A gene with an anchor but no accession can never be recognized from a
    predictor row; a gene with an accession but no anchor adds a label and no
    economics. Both are worth knowing about, so both are asserted.
    """
    from econ import gene_anchors as ga
    anchored = set(ga.GENE_ANCHORS) | set(ga.NOT_VALUED_GENES)
    mapped = set(nv.UNIPROT_TO_GENE.values())
    assert not (anchored - mapped), (
        f"anchor genes with no UniProt accession: {sorted(anchored - mapped)}")
    assert len(nv.UNIPROT_TO_GENE) == len(mapped), "duplicate gene in the map"


def test_gene_basis_records_how_the_gene_was_assigned(tmp_path, monkeypatch):
    """A gene from the predictor and a gene from an interval are different
    claims, and the record says which it is.

    A coordinate inside a gene's span is a GUESS: spans are not exon models, so
    an intronic position or one inside an overlapping gene is assigned
    confidently and wrongly. Recording the basis is what lets a reader discount
    it appropriately instead of reading both as equally certain.
    """
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    # BRCA1 GRCh38 span is 43,044,295-43,125,483; this sits inside it.
    _make_am_table(tmp_path, [("17", "43050000", "C", "T", "0.97",
                               "likely_pathogenic")])
    vcf = _write_vcf(tmp_path, [("17", "43050000", "C", "T", "0/1")])
    r = nv.analyze_novel_variants(vcf, "grch38")
    assert r["available"] is True
    dmg = r["buckets"]["predicted_pathogenic_rare"]
    assert dmg, "expected a predicted-pathogenic finding"
    f = dmg[0]
    # The synthetic table carries a placeholder accession, so this falls through
    # to the window path — which is exactly the case that must be labeled.
    assert f["gene"] == "BRCA1"
    assert f["gene_basis"] == "coordinate_window"


def test_unassignable_variant_is_labeled_unassigned(tmp_path, monkeypatch):
    monkeypatch.setattr(nv, "REFERENCE_DIR", tmp_path)
    _make_am_table(tmp_path, [("1", "999999", "C", "T", "0.97",
                               "likely_pathogenic")])
    vcf = _write_vcf(tmp_path, [("1", "999999", "C", "T", "0/1")])
    r = nv.analyze_novel_variants(vcf, "grch38")
    f = r["buckets"]["predicted_pathogenic_rare"][0]
    assert f["gene"] == ""
    assert f["gene_basis"] == "unassigned"
