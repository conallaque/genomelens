"""Correctness guarantees for maternal-lineage calling.

Two classes of error matter here and neither shows up as a crash: a marker
recorded with its derived and ancestral alleles the wrong way round, and a
clade attached to the wrong parent. Both produce a confident, plausible,
wrong haplogroup. Both existed in this module before these tests.
"""

import pandas as pd

import mt_haplogroup as mt


def _df(rows):
    """rows: list of (rsid, pos, genotype) on chrMT."""
    return pd.DataFrame([(r, "MT", p, g) for r, p, g in rows],
                        columns=["rsid", "chrom", "pos", "genotype"])


def _nodes(tree, out=None):
    out = [] if out is None else out
    out.append(tree)
    for c in tree.get("children", []):
        _nodes(c, out)
    return out


def _find(hg):
    return next(n for n in _nodes(mt.MT_TREE) if n["haplogroup"] == hg)


def _parent(hg):
    for n in _nodes(mt.MT_TREE):
        if any(c["haplogroup"] == hg for c in n.get("children", [])):
            return n["haplogroup"]
    return None


# ══════════════════════════════════════════════════════════════════════════
# rCRS orientation — the marker-inversion trap
# ══════════════════════════════════════════════════════════════════════════

def test_R_derived_allele_matches_the_reference_sequence():
    # THE BUG THIS CATCHES. The revised Cambridge Reference Sequence is itself
    # a haplogroup H2a2a1 sequence, and H sits inside R. So R's derived allele
    # at 12705 must be the rCRS allele, C. The flat marker list recorded T,
    # which claims R carries the non-R state — an inversion at a backbone node
    # that misroutes every sample passing through it.
    m = next(x for x in mt.MTDNA_MARKERS if x["pos"] == 12705)
    assert m["derived"] == "C", (
        "12705 derived allele must be C: rCRS is a haplogroup H sequence, H is "
        "inside R, and rCRS carries C here")
    assert m["ancestral"] == "T"


def test_the_two_R_defining_markers_agree_with_each_other():
    # 12705 and 16223 both discriminate R from non-R. Before the fix they
    # disagreed: 16223 treated T as the non-R ancestral state while 12705
    # treated T as R-defining. Two markers for one branch point cannot point
    # in opposite directions.
    r = _find("R")
    by_pos = {pos: (der, anc) for (_n, _r, pos, der, anc) in r["markers"]}
    assert by_pos[12705][0] == "C"
    assert by_pos[16223][0] == "C"


def test_tree_and_flat_marker_list_do_not_contradict_each_other():
    # The flat list is still consumed by the report's marker table while the
    # tree drives the lineage chain. If they disagree, the two halves of the
    # section describe different people.
    flat = {m["pos"]: (m["derived"], m["ancestral"]) for m in mt.MTDNA_MARKERS}
    for node in _nodes(mt.MT_TREE):
        for (_name, _rsids, pos, der, anc) in node.get("markers", []):
            if pos in flat:
                assert flat[pos] == (der, anc), (
                    f"position {pos} disagrees between the tree "
                    f"({der}/{anc}) and the flat marker list {flat[pos]}")


def test_clades_containing_H_use_the_reference_allele_as_derived():
    # Generalises the 12705 rule. rCRS is H, so every clade H belongs to must
    # have its derived state equal to the rCRS allele at that position. These
    # are the known rCRS bases at the backbone positions H passes through.
    rcrs = {12705: "C", 16223: "C", 14766: "C", 7028: "C"}
    for hg in ("R", "HV", "H"):
        for (_n, _r, pos, der, _anc) in _find(hg)["markers"]:
            if pos in rcrs:
                assert der == rcrs[pos], (
                    f"{hg} marker at {pos} has derived={der}; H lies inside "
                    f"{hg} and rCRS carries {rcrs[pos]} there")


# ══════════════════════════════════════════════════════════════════════════
# Topology
# ══════════════════════════════════════════════════════════════════════════

def test_I_W_and_X_hang_off_N_not_R():
    # These are N lineages, not R lineages. The old priority-list classifier
    # had no way to express parentage at all, so someone carrying both an R
    # marker and an X marker got whichever the priority list happened to
    # rank first.
    for hg in ("I", "W", "X"):
        assert _parent(hg) == "N", f"{hg} must descend from N, not {_parent(hg)}"


def test_european_backbone_parentage_is_correct():
    for child, parent in (("R", "N"), ("HV", "R"), ("H", "HV"), ("V", "HV"),
                          ("JT", "R"), ("J", "JT"), ("T", "JT"),
                          ("U", "R"), ("K", "U"), ("M", "mt-MRCA"),
                          ("N", "mt-MRCA")):
        assert _parent(child) == parent, (
            f"{child} should descend from {parent}, got {_parent(child)}")


def test_every_node_except_the_root_and_N_carries_a_marker():
    # N has no marker on consumer arrays and is inferred; everything else must
    # be observable, or it can never be confirmed.
    for n in _nodes(mt.MT_TREE):
        if n["haplogroup"] in ("mt-MRCA", "N"):
            continue
        assert n["markers"], f"{n['haplogroup']} has no defining marker"


def test_no_haplogroup_appears_twice_in_the_tree():
    names = [n["haplogroup"] for n in _nodes(mt.MT_TREE)]
    assert len(names) == len(set(names)), "duplicate clade in the phylogeny"


# ══════════════════════════════════════════════════════════════════════════
# The lineage chain
# ══════════════════════════════════════════════════════════════════════════

def test_lineage_path_is_reported_not_just_the_endpoint():
    # THE POINT OF THE CHANGE. Y-DNA has always reported the chain of branch
    # points; mtDNA reported a bare label.
    r = mt.analyze_mt_haplogroup(_df([
        ("rs28358580", 12705, "C"), ("rs193302980", 14766, "C"),
        ("rs2854122", 7028, "C")]))
    assert r["haplogroup"] == "H"
    assert r["haplogroup_path"] == "mt-MRCA > N > R > HV > H"
    assert [n["haplogroup"] for n in r["path"]] == \
        ["mt-MRCA", "N", "R", "HV", "H"]


def test_each_link_carries_its_own_evidence():
    r = mt.analyze_mt_haplogroup(_df([
        ("rs28358580", 12705, "C"), ("rs193302980", 14766, "C"),
        ("rs2854122", 7028, "C")]))
    confirmed = [n for n in r["path"] if n["snp_status"] == "confirmed"]
    assert {n["haplogroup"] for n in confirmed} == {"R", "HV", "H"}
    for n in confirmed:
        assert n["evidence"], f"{n['haplogroup']} confirmed with no evidence"


def test_unobserved_links_are_marked_inferred_not_confirmed():
    # N is never typed on a consumer array. Claiming it as confirmed would
    # overstate the data; omitting it would break the chain.
    r = mt.analyze_mt_haplogroup(_df([
        ("rs28358580", 12705, "C"), ("rs2854122", 7028, "C")]))
    n_node = next(n for n in r["path"] if n["haplogroup"] == "N")
    assert n_node["snp_status"] == "inferred"
    assert "N" in r["chip_gaps"]


def test_the_chain_never_runs_deeper_than_the_evidence():
    # With only an R marker, the call must stop at R rather than guessing at
    # the far more common H below it.
    r = mt.analyze_mt_haplogroup(_df([("rs28358580", 12705, "C")]))
    assert r["haplogroup_path"].endswith("R")
    assert r["haplogroup"] == "R"


def test_a_U_lineage_routes_through_R_and_reaches_K():
    r = mt.analyze_mt_haplogroup(_df([
        ("rs28358580", 12705, "C"), ("rs28359175", 12308, "G"),
        ("rs28358571", 9055, "A")]))
    assert [n["haplogroup"] for n in r["path"]] == \
        ["mt-MRCA", "N", "R", "U", "K"]


def test_an_X_lineage_does_not_route_through_R():
    r = mt.analyze_mt_haplogroup(_df([("rs_x", 6221, "C")]))
    # Compare node names, not the joined string — "R" is a substring of
    # "mt-MRCA" and a naive containment check passes for the wrong reason.
    assert [n["haplogroup"] for n in r["path"]] == ["mt-MRCA", "N", "X"], (
        "X is an N lineage; routing it through R is the parentage bug")
    assert r["haplogroup"] == "X"


def test_conflicting_sibling_branches_are_flagged_not_silently_resolved():
    # One person has one maternal line. Two confirmed mutually exclusive
    # branches means a mis-called genotype, and hiding it would present a
    # data problem as a result.
    r = mt.analyze_mt_haplogroup(_df([
        ("rs28358580", 12705, "C"), ("rs193302980", 14766, "C"),
        ("rs3088309", 4216, "G")]))
    assert r["contradictions"], "conflicting HV and JT calls should be flagged"


def test_branch_point_counts_are_consistent_with_the_path():
    r = mt.analyze_mt_haplogroup(_df([
        ("rs28358580", 12705, "C"), ("rs193302980", 14766, "C"),
        ("rs2854122", 7028, "C")]))
    assert r["n_branch_points"] == len(
        [n for n in r["path"] if n["haplogroup"] != "mt-MRCA"])
    assert r["n_confirmed_branch_points"] == len(
        [n for n in r["path"] if n["snp_status"] == "confirmed"])


# ══════════════════════════════════════════════════════════════════════════
# Degradation and backward compatibility
# ══════════════════════════════════════════════════════════════════════════

def test_no_mtdna_data_degrades_with_every_key_present():
    # Consumers index these keys directly; a missing one is a crash in the
    # report rather than a graceful "unknown".
    r = mt.analyze_mt_haplogroup(
        pd.DataFrame(columns=["rsid", "chrom", "pos", "genotype"]))
    assert r["status"] == "no_data"
    for key in ("haplogroup", "path", "haplogroup_path", "terminal_haplogroup",
                "chip_gaps", "contradictions", "n_branch_points",
                "n_confirmed_branch_points", "matched_markers"):
        assert key in r, f"missing {key} on the no-data path"
    assert r["path"] == [] and r["haplogroup"] == "Unknown"


def test_existing_result_keys_are_preserved():
    # The report, narrative, deep-ancestry and family-planning modules all
    # read this dict. Adding the chain must not remove anything.
    r = mt.analyze_mt_haplogroup(_df([("rs2854122", 7028, "C")]))
    for key in ("status", "haplogroup", "confidence", "confidence_note",
                "matched_markers", "n_markers_matched", "n_markers_derived",
                "n_markers_expected", "mt_snp_count", "migration",
                "ancient_dna", "further_testing", "message"):
        assert key in r, f"dropped previously-present key {key}"


def test_ancestral_genotype_excludes_a_branch():
    # A confirmed ancestral call at H must not yield an H lineage.
    r = mt.analyze_mt_haplogroup(_df([
        ("rs28358580", 12705, "C"), ("rs2854122", 7028, "T")]))
    assert not r["haplogroup_path"].endswith("H")
