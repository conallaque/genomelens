"""Checks on the v1.11 demonstration reports and the reproduction kit. They read only files committed to this repository."""
import pathlib
import re
import shutil
import subprocess
import sys

import fitz
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TABLE = ROOT / "reproduction" / "kastrinos2021_table2.tsv"
NEEDS_TABLE = pytest.mark.skipif(not TABLE.exists(), reason="the published table is not shipped (publisher terms); supply it as described in reproduction/README.md to run this check")
DEMO = ROOT / "demonstrations" / "v1.11"
README = (ROOT / "README.md").read_text(encoding="utf-8")
PDFS = sorted(DEMO.glob("GenomeLens_public_demonstration_*.pdf"))
ALLOWED_BASIS = {"Published model", "Hypothetical, assumption-based scenario", "No basis label (a drug-response decision assembled from published inputs)", "n/a (no modeled value)",
                 "Measured null: the evidence finds no incremental benefit"}
EXPECTED_REPORTS = 16     # the 15 demonstration genomes at the standard stated age, plus the published-model demonstration at its source's age


def test_the_demonstration_reports_are_published():
    assert len(PDFS) == EXPECTED_REPORTS


@pytest.mark.parametrize("pdf", PDFS, ids=lambda p: p.stem.split("demonstration_")[-1])
def test_a_demonstration_pdf_carries_only_its_title_and_no_local_path_or_draft_text(pdf):
    d = fitz.open(pdf)
    meta = {k: v for k, v in (d.metadata or {}).items() if v and k not in ("format", "title", "encryption")}
    assert not meta, meta
    assert not d.get_xml_metadata()
    text = "\n".join(p.get_text() for p in d)
    links = " ".join(l.get("uri", "") for p in d for l in p.get_links())
    for bad in ("file:" + "///", "/Us" + "ers/", "/private/" + "tmp", "Draft " + "notice", "owner's " + "wording", "{" + "{"):
        assert bad not in text and bad not in links, bad
    assert "labeled UNRESOLVED" in text or "stated for this run" in text


def test_the_readme_demonstration_table_matches_the_published_files():
    rows = [l for l in README.splitlines() if l.startswith("| ") and "demonstrations/v1.11/GenomeLens_public_demonstration_" in l]
    assert len(rows) == EXPECTED_REPORTS
    for r in rows:
        cells = [c.strip() for c in r.strip("|").split("|")]
        assert (ROOT / re.search(r"\((demonstrations/v1\.11/[^)]+\.pdf)\)", r).group(1)).exists()
        assert cells[4] in ALLOWED_BASIS or cells[4] == "", cells[4]


def test_the_four_primary_demonstrations_are_stated_in_order_and_link_to_published_reports():
    block = README[README.index("## Four demonstrations"):README.index("## What a GenomeLens report looks like")]
    heads = re.findall(r"^\*\*([1-4])\. ([^*]+)\*\*", block, flags=re.M)
    assert [h[0] for h in heads] == ["1", "2", "3", "4"], heads
    assert heads[0][1].startswith("A published-model decision") and heads[1][1].startswith("A hypothetical, assumption-based scenario")
    assert heads[2][1].startswith("A correctly classified measured null") and heads[3][1].startswith("A clinically meaningful finding whose monetization is withheld")
    links = re.findall(r"\[Report\]\((demonstrations/v1\.11/[^)]+\.pdf)\)", block)
    assert len(links) == 4 and all((ROOT / l).exists() for l in links), links
    assert "Class: published model" in block and "Class: hypothetical, assumption-based scenario" in block
    assert "withholds a dollar figure and prints the reason; the result is not zero" in block


def test_the_published_model_demonstration_is_the_source_age_report_and_the_hypothetical_one_is_the_standard_age_report():
    rows = {re.match(r"\| ([^|]+) \|", l).group(1).strip(): l for l in README.splitlines()
            if l.startswith("| hu4AE0F2") and "demonstrations/v1.11/" in l}
    old, new = rows["hu4AE0F2 (stated age 30)"], rows["hu4AE0F2"]
    assert "| Published model |" in old and "| Hypothetical, assumption-based scenario |" in new


def test_the_readme_keeps_published_results_apart_from_assumption_based_scenarios_and_claims_no_external_validation():
    assert "{{" not in README
    for phrase in ("Published-model result", "Hypothetical, assumption-based scenario", "Measured null", "Not monetized",
                   "no model has been externally validated"):
        assert phrase.lower() in README.lower(), phrase


@NEEDS_TABLE
def test_the_published_table_reproduction_matches_its_recorded_output():
    out = subprocess.run([sys.executable, str(ROOT / "reproduction" / "reproduce_table2.py")], capture_output=True, text=True, check=True).stdout
    assert out == (ROOT / "reproduction" / "EXPECTED_OUTPUT.txt").read_text(encoding="utf-8")
    # the table compares schedules among carriers who are already known; its optimal strategy has a different net monetary benefit against
    # each of two comparators, and neither is the value of identifying a carrier
    assert "MSH2 optimal: Q2Y, Start age 25" in out
    assert "against the next efficient strategy (Q3Y, Start age 25): QALYs +0.028, cost +908.26, NMB 1,891.74" in out
    assert "against the least intensive strategy listed (Q5Y, Start age 40): QALYs +0.664, cost -25,130.91, NMB 91,530.91" in out
    # one printed ratio does not reproduce from the printed totals, and the script says so rather than hiding it
    assert out.count("does not reproduce") == 1 and "MSH2 Q2Y, Start age 25 vs Q3Y, Start age 25: printed 29,298.71, recomputed 32,437.86" in out
    assert out.count("dominated (printed Dominated)") == 12 and out.count("frontier (printed") == 8      # 20 printed strategies


@NEEDS_TABLE
@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is not installed")
def test_the_r_reproduction_prints_the_same_recorded_output():
    out = subprocess.run(["Rscript", str(ROOT / "reproduction" / "reproduce_table2.R")], capture_output=True, text=True, check=True).stdout
    assert out == (ROOT / "reproduction" / "EXPECTED_OUTPUT.txt").read_text(encoding="utf-8")


def _png_chunk_types(path):
    data = path.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    i, out = 8, []
    while i < len(data):
        n = int.from_bytes(data[i:i + 4], "big")
        out.append(data[i + 4:i + 8].decode("latin-1"))
        i += 12 + n
    return out


def test_the_preview_image_carries_no_metadata():
    pngs = list((DEMO / "previews").glob("*.png"))
    assert pngs
    for p in pngs:
        assert not [c for c in _png_chunk_types(p) if c in ("tEXt", "iTXt", "zTXt", "eXIf", "tIME")], p.name


# ── the second reading of the release (three independent reviewers, 2026-10-23 evening) ──────────────────────────────────────────
def _block(start, end):
    return README[README.index(start):README.index(end)]


def test_demonstration_two_names_the_published_model_it_does_not_count_and_says_the_two_values_are_not_comparable():
    block = _block("## Four demonstrations", "## What a GenomeLens report looks like")
    two = block[block.index("**2. "):block.index("**3. ")]
    assert "Bennette 2015" in two and "does not count it" in two and "not corrected for ascertainment" in two
    assert "not comparable" in two and "not an age effect" in two
    assert "No published model reports" not in README


def test_demonstration_one_says_the_source_finds_a_lynch_only_test_not_cost_effective_and_the_value_is_restated():
    one = _block("## Four demonstrations", "**2. ")
    assert "restated to a per-carrier value by stated arithmetic" in one and "not cost-effective" in one and "per QALY gained" in one


def test_the_measured_null_demonstration_states_its_evidence_and_is_not_called_a_missing_value():
    three = _block("**3. ", "**4. ")
    assert "0.64 to 1.55" in three and "37.4% vs 35.6%" in three and "not a missing value" in three and "no total is formed" not in three


def test_the_withheld_demonstration_gives_its_whole_reason():
    four = _block("**4. ", "## What a GenomeLens report looks like")
    assert "Fuerte 2025" in four and "Roberts 2025" in four and "TPMT" in four and "…" not in four


def test_no_public_report_prints_a_placeholder_figure_for_a_row_whose_genotype_call_is_not_established():
    for pdf in PDFS:
        text = "\n".join(p.get_text() for p in fitz.open(pdf))
        for bad in ("conditional on an unconfirmed call", "SOURCE-GENRE", "zero by rule", "The evidence measures no effect"):
            assert bad.lower() not in text.lower(), (pdf.name, bad)


def test_the_readme_discloses_ai_authorship_by_name_and_carries_no_drafting_note():
    assert "Anthropic's Claude" in README and "ChatGPT" in README
    for bad in ("provisional wording", "pending the owner", "recomputed independently", "have the work checked", "Reproducible public-genome"):
        assert bad.lower() not in README.lower(), bad


def test_the_reproduction_kit_prints_the_discounting_check_and_says_its_convention_is_unresolved():
    out = (ROOT / "reproduction" / "EXPECTED_OUTPUT.txt").read_text(encoding="utf-8")
    assert "4. A check of the stated discounting" in out and "all above 26.50" in out and "unresolved" in out
    if TABLE.exists():
        assert "PMC9330543_efetch" not in TABLE.read_text(encoding="utf-8")


def test_the_public_scanner_covers_the_reproduction_kit_and_no_text_file_names_the_build_machine():
    from scan_engine import TEXT_SUFFIXES
    assert {".R", ".tsv"} <= set(TEXT_SUFFIXES)
    bad = ("claude-" + "501", "Desktop-" + "Claude", "/Us" + "ers/")
    for p in sorted(ROOT.rglob("*")):
        if p.is_file() and p.suffix in TEXT_SUFFIXES and ".git" not in p.parts and p.name != "test_release_security.py":
            t = p.read_text(encoding="utf-8", errors="replace")
            assert not any(b in t for b in bad), (str(p.relative_to(ROOT)), [b for b in bad if b in t])


def test_no_stray_backup_file_ships():
    stray = [str(p.relative_to(ROOT)) for p in ROOT.rglob("*") if p.is_file() and ".git" not in p.parts
             and any(k in p.name for k in ("backup", ".bak", ".orig"))]
    assert not stray, stray


@pytest.mark.skipif(TABLE.exists(), reason="the table has been supplied")
def test_without_the_published_table_each_script_says_how_to_supply_it_and_runs_nothing():
    for cmd in ([sys.executable, str(ROOT / "reproduction" / "reproduce_table2.py")],
                ["Rscript", str(ROOT / "reproduction" / "reproduce_table2.R")] if shutil.which("Rscript") else None):
        if cmd is None:
            continue
        r = subprocess.run(cmd, capture_output=True, text=True)
        assert r.returncode == 2 and r.stdout == "" and "reproduction/README.md" in r.stderr, (cmd[0], r.returncode, r.stderr[:200])


@pytest.mark.skipif(TABLE.exists(), reason="the table has been supplied")
def test_the_kit_ships_the_scripts_a_template_and_a_note_but_not_the_publishers_table():
    kit = ROOT / "reproduction"
    assert {p.name for p in kit.iterdir()} == {"EXPECTED_OUTPUT.txt", "README.md", "kastrinos2021_table2.TEMPLATE.tsv", "reproduce_table2.R", "reproduce_table2.py"}
    rows = [l for l in (kit / "kastrinos2021_table2.TEMPLATE.tsv").read_text(encoding="utf-8").splitlines() if l and not l.startswith("#")]
    assert rows == ["gene\tstrategy_as_printed\tQALYs\tCost\tICER"]          # a header only: no printed value is copied
