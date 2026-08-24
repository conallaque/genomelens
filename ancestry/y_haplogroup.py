"""
Y-DNA Haplogroup Analysis Module
Walks a decision tree of Y-chromosome haplogroup-defining SNPs covering the
major backbone of the human Y phylogeny (CT → CF → F → {G, IJK → {IJ → {I, J},
K → {LT → {L, T}, K2 → {NO → {N, O}, P → {Q, R → R1 → {R1a, R1b}}}}}}).

Marker data (rsID, GRCh37 position, ancestral→derived) is taken from the ISOGG
Y-SNP index. Each node carries SEVERAL co-defining markers and the call is by
majority vote across whichever of them the chip actually typed. This makes the
result robust to a single mistyped/missing SNP.

Tree-walking algorithm:
  1. For each node, look up every defining marker by rsID then by chrY position.
     DERIVED votes vs ANCESTRAL votes decide the node's status.
  2. ANCESTRAL  → this branch is ruled out; do not descend.
  3. DERIVED    → confirmed; descend into children.
  4. NOT_FOUND  → chip gap; descend ONLY toward a branch that contains a
                  genuinely-confirmed downstream marker (otherwise stop — never
                  guess a deeper subclade than the data supports).
  5. The reported terminal haplogroup is the DEEPEST CONFIRMED node. Inferred
     (chip-gap) ancestors between two confirmed markers are kept but labelled.

Strand safety: markers whose ancestral/derived alleles are a complementary pair
(A/T or C/G) cannot be oriented from genotype alone and are skipped — clades are
defined by non-ambiguous co-markers instead.
"""


import contextlib

import pandas as pd

# ── Complement helper ──────────────────────────────────────────────────────────
_COMP = str.maketrans("ACGT", "TGCA")


def complement(base: str) -> str:
    return base.upper().translate(_COMP)


def _comp(base: str) -> str:
    return base.upper().translate(_COMP)


# ── Node builder ────────────────────────────────────────────────────────────────
#
# markers: list of (name, rsid, GRCh37_pos, ancestral, derived). rsid may be ""
# (position-only lookup). Strand-ambiguous pairs (A/T, C/G) are skipped at
# lookup time, so a clade should always carry at least one non-ambiguous marker.

def _node(haplogroup: str, markers: list[tuple], description: str = "",
          migration: str = "", further: str = "", children: list | None = None) -> dict:
    rsids = [m[1] for m in markers if m[1]]
    return {
        "haplogroup": haplogroup,
        "snp_name": markers[0][0] if markers else haplogroup,
        "rsids": rsids,
        "pos": markers[0][2] if markers else None,
        "markers": markers,
        "description": description,
        "migration": migration,
        "further": further,
        "children": children or [],
    }


# ── Haplogroup decision tree (backbone) ─────────────────────────────────────────

HAPLOGROUP_TREE: dict = _node(
    "CT", [("M168", "rs2032595", 14813991, "C", "T"),
           ("M5576", "", 2744386, "G", "T"),
           ("M5577", "", 2757670, "C", "T")],
    description="All non-African lineages plus most African ones (everyone except A and B).",
    migration=(
        "Haplogroup CT (M168) marks the common patrilineal ancestor of the vast majority of "
        "men alive today, living in Africa roughly 70,000 years ago, just before the "
        "out-of-Africa expansion."
    ),
    further="Any confirmed downstream marker implicitly confirms CT.",
    children=[
        _node(
            "E", [("M66", "rs2032627", 21881573, "A", "C"),
                  ("M155", "", 21736331, "G", "A"),
                  ("M156", "", 21717227, "T", "C")],
            description="Common in Africa, the Middle East and southern Europe (approximates DE/E).",
            migration=(
                "Haplogroup E arose in Africa or the Near East ~50,000 years ago and is today "
                "the most common lineage in Africa, with significant frequencies in the Levant "
                "and Mediterranean Europe."
            ),
        ),
        _node(
            "CF", [("M3690", "", 15203676, "A", "G"),
                   ("P143", "rs4141886", 14197867, "G", "A")],
            description="Ancestor of haplogroups C and F — nearly all non-African men.",
            children=[
                _node(
                    "F", [("M89", "rs2032652", 21917313, "C", "T"),
                          ("M235", "rs7067496", 14832620, "T", "G"),
                          ("P135", "rs9786502", 21618856, "C", "T")],
                    description="Ancestor of >90% of all men outside Africa.",
                    migration=(
                        "Haplogroup F (M89) appeared ~48,000 years ago in South Asia or the Near "
                        "East. Its descendants — G, H, I, J, K and everything below K — account "
                        "for the overwhelming majority of non-African paternal lineages."
                    ),
                    children=[
                        _node(
                            "G", [("M201", "rs2032636", 15027529, "G", "T"),
                                  ("M3242", "", 7145960, "C", "T"),
                                  ("M3248", "", 7565637, "G", "A")],
                            description="Caucasus, Anatolia and Mediterranean Europe; early Neolithic farmers.",
                            migration=(
                                "Haplogroup G (M201) is associated with the spread of early "
                                "farming from the Near East and Caucasus into Neolithic Europe."
                            ),
                        ),
                        _node(
                            "IJK", [("M522", "rs9786714", 7173143, "G", "A"),
                                    ("M523", "rs9786139", 6753519, "A", "G")],
                            description="Ancestor of haplogroups I, J and K.",
                            children=[
                                _node(
                                    "IJ", [("P123", "rs17315821", 19166861, "T", "C"),
                                           ("P127", "rs7892893", 8590752, "C", "T"),
                                           ("P129", "rs17306699", 14144593, "A", "G")],
                                    description="Ancestor of haplogroups I and J.",
                                    children=[
                                        _node(
                                            "I", [("M170", "rs2032597", 14847792, "A", "C"),
                                                  ("M161", "", 21717515, "G", "T")],
                                            description="Indigenous European lineage (e.g. I1 Nordic, I2 Balkan/Sardinian).",
                                            migration=(
                                                "Haplogroup I (M170) is the oldest European-specific "
                                                "lineage, present among Palaeolithic hunter-gatherers "
                                                "and still common across Europe."
                                            ),
                                            children=[
                                                _node("I1", [("M253", "rs9341296", 15022707, "C", "T"),
                                                             ("P30", "rs112707890", 14496753, "G", "A"),
                                                             ("P40", "rs113686221", 14484394, "C", "T")],
                                                      description="Scandinavia / Northwest Europe."),
                                                _node("I2", [("M438", "rs17307294", 16638804, "A", "G"),
                                                             ("PF3781", "rs35547782", 18700150, "C", "T")],
                                                      description="Southeastern / Central Europe, Sardinia."),
                                            ],
                                        ),
                                        _node(
                                            "J", [("M304", "rs13447352", 22749853, "A", "C"),
                                                  ("M280", "rs13447367", 21878762, "G", "A"),
                                                  ("M289", "rs13447368", 21878708, "G", "A")],
                                            description="Near East, Arabia, Caucasus, Mediterranean.",
                                            migration=(
                                                "Haplogroup J (M304) originated in the Near East and "
                                                "spread with Neolithic farmers and later Semitic and "
                                                "Mediterranean populations."
                                            ),
                                            children=[
                                                _node("J1", [("M267", "rs9341313", 22741818, "T", "G")],
                                                      description="Arabia, Levant, Caucasus."),
                                                _node("J2", [("M172", "rs2032604", 14969634, "T", "G"),
                                                             ("L228", "", 7771358, "C", "T")],
                                                      description="Anatolia, Levant, Mediterranean."),
                                            ],
                                        ),
                                    ],
                                ),
                                _node(
                                    "K", [("P128", "rs17250121", 20837553, "C", "T"),
                                          ("P131", "rs9786043", 15472863, "C", "T"),
                                          ("P132", "rs3853054", 8679843, "G", "T")],
                                    description="Macro-haplogroup K — ancestor of L, T, N, O, Q, R and more.",
                                    migration=(
                                        "Macro-haplogroup K (M9 and equivalents) emerged ~45,000 years "
                                        "ago in South or Central Asia. Its descendants dominate Europe, "
                                        "Asia, Oceania and the Americas."
                                    ),
                                    further=(
                                        "Defined here by P128/P131/P132 rather than M9, which is not on "
                                        "most consumer chips."
                                    ),
                                    children=[
                                        _node(
                                            "LT", [("P326", "", 8467290, "T", "C"),
                                                   ("PF5525", "", 6994764, "G", "A"),
                                                   ("PF5531", "", 8628308, "C", "T")],
                                            description="Ancestor of haplogroups L and T (also called K1).",
                                            children=[
                                                _node("L", [("M20", "rs3911", 21733454, "A", "G"),
                                                            ("M11", "rs3902", 21730647, "A", "G"),
                                                            ("M185", "rs2032607", 14904859, "C", "T")],
                                                      description="South Asia, with branches in the Near East."),
                                                _node(
                                                    "T", [("M272", "rs9341308", 22738775, "A", "G"),
                                                          ("M320", "rs13447374", 15030767, "T", "G"),
                                                          ("PF5597", "", 6794129, "G", "A")],
                                                    description="A relatively rare lineage of the Near East, East Africa and the Mediterranean.",
                                                    migration=(
                                                        "Haplogroup T (M184/M272) is an old and "
                                                        "geographically scattered lineage found at low "
                                                        "frequency around the Near East, East Africa, the "
                                                        "Mediterranean and the Horn of Africa, with deep "
                                                        "roots tracing back to the early diversification "
                                                        "of macro-haplogroup K."
                                                    ),
                                                    children=[
                                                        _node(
                                                            "T1a", [("M70", "rs2032672", 21893881, "A", "C"),
                                                                    ("PF7472.1", "", 16212441, "G", "A")],
                                                            description="The dominant and most widespread branch of haplogroup T.",
                                                            children=[
                                                                _node(
                                                                    "T1a1", [("L454", "", 14577272, "C", "T"),
                                                                             ("FGC3945.2", "", 8032311, "G", "A"),
                                                                             ("CTS5542", "", 16350661, "A", "C")],
                                                                    description="Major sub-branch of T1a, common in Europe and the Near East.",
                                                                    children=[
                                                                        _node(
                                                                            "T1a1a", [("CTS2611", "", 14389760, "G", "A"),
                                                                                      ("L905", "", 6659212, "A", "C")],
                                                                            description="A widespread T1a1 subclade found across Europe, the Near East and North Africa.",
                                                                            further=(
                                                                                "FTDNA Big Y-700 or equivalent Y-sequencing "
                                                                                "would resolve subclades below T1a1a."
                                                                            ),
                                                                        ),
                                                                    ],
                                                                ),
                                                            ],
                                                        ),
                                                    ],
                                                ),
                                            ],
                                        ),
                                        _node(
                                            "K2", [("M526", "rs2033003", 23550924, "A", "C")],
                                            description="Ancestor of N, O, P, Q and R.",
                                            children=[
                                                _node(
                                                    "NO", [("M2313", "", 8674808, "C", "T"),
                                                           ("CTS11667", "", 23208284, "G", "A")],
                                                    description="Ancestor of haplogroups N and O.",
                                                    children=[
                                                        _node("N", [("M231", "rs9341278", 15469724, "G", "A"),
                                                                    ("M232", "", 15437152, "C", "T")],
                                                              description="Northern Eurasia, Siberia, Finland, Baltic.",
                                                              migration=(
                                                                  "Haplogroup N (M231) spread across northern "
                                                                  "Eurasia and is common among Uralic-speaking "
                                                                  "and Siberian populations.")),
                                                        _node("O", [("M297", "rs13447345", 22746689, "A", "G"),
                                                                    ("M1530", "", 7060243, "G", "A")],
                                                              description="East and Southeast Asia — the most common lineage there.",
                                                              migration=(
                                                                  "Haplogroup O (M175) is the dominant paternal "
                                                                  "lineage of East and Southeast Asia.")),
                                                    ],
                                                ),
                                                _node(
                                                    "P", [("P295", "rs895530", 7963031, "T", "G"),
                                                          ("PF5862", "rs7892927", 7628900, "G", "A")],
                                                    description="Ancestor of haplogroups Q and R (also called K2b2).",
                                                    children=[
                                                        _node("Q", [("M242", "rs8179021", 15018582, "C", "T"),
                                                                    ("M1064", "", 6778043, "G", "A")],
                                                              description="Siberia and the indigenous Americas.",
                                                              migration=(
                                                                  "Haplogroup Q (M242) crossed Beringia and is "
                                                                  "the predominant lineage of indigenous "
                                                                  "peoples of the Americas.")),
                                                        _node(
                                                            "R", [("M207", "rs2032658", 15581983, "A", "G"),
                                                                  ("P224", "rs17307398", 17285993, "C", "T")],
                                                            description="Ancestor of R1a and R1b — most common in Europe and South Asia.",
                                                            migration=(
                                                                "Haplogroup R (M207) is the most common lineage "
                                                                "in Europe and is widespread across South and "
                                                                "Central Asia."
                                                            ),
                                                            children=[
                                                                _node(
                                                                    "R1", [("M173", "rs2032624", 15026424, "A", "C"),
                                                                           ("M306", "rs1558843", 22750583, "C", "A"),
                                                                           ("P225", "rs17307070", 15590342, "G", "T")],
                                                                    description="Ancestor of R1a and R1b.",
                                                                    children=[
                                                                        _node("R1a", [("M420", "rs17250535", 23473201, "T", "A"),
                                                                                      ("M64.2", "rs2032626", 21903383, "A", "G"),
                                                                                      ("M87", "rs2032644", 21906109, "T", "C")],
                                                                              description="Eastern Europe, Central and South Asia."),
                                                                        _node(
                                                                            "R1b", [("M343", "rs9786184", 2887824, "C", "A"),
                                                                                    ("M228.1", "rs9341273", 15591445, "T", "C")],
                                                                            description="Western Europe — the most common European lineage.",
                                                                            children=[
                                                                                _node("R1b1a2", [("M269", "rs9786153", 22739367, "T", "C"),
                                                                                                 ("PF6399", "rs2058276", 2668456, "C", "T")],
                                                                                      description="The dominant Western-European R1b branch."),
                                                                            ],
                                                                        ),
                                                                    ],
                                                                ),
                                                            ],
                                                        ),
                                                    ],
                                                ),
                                            ],
                                        ),
                                    ],
                                ),
                            ],
                        ),
                    ],
                ),
            ],
        ),
    ],
)


# ── Lookup tables ────────────────────────────────────────────────────────────────

def _build_lookup(snps_df: pd.DataFrame) -> tuple[dict[str, str], dict[int, str]]:
    """
    Build two lookup tables from the parsed DataFrame:
      rsid_map  : rsid (str) → genotype (str, upper-case)
      pos_map   : chr-Y integer position → genotype
    Handles both hemizygous single-char ('A') and doubled ('AA') genotypes.
    Filters to chromosome Y rows only.
    """
    rsid_map: dict[str, str] = {}
    pos_map: dict[int, str] = {}

    def _is_y(chrom_val) -> bool:
        return str(chrom_val).strip().upper() in ("Y", "24")

    chrom_col = "chrom" if "chrom" in snps_df.columns else None

    for rsid_raw, row in snps_df.iterrows():
        if chrom_col and not _is_y(row[chrom_col]):
            continue

        gt = str(row.get("genotype", "")).strip().upper()
        if not gt or gt in ("NAN", "--", "0", ""):
            continue

        if len(gt) == 2 and gt[0] == gt[1]:
            gt = gt[0]

        rsid_str = str(rsid_raw).strip()
        rsid_map[rsid_str] = gt

        pos = row.get("pos", None)
        if pos is not None:
            with contextlib.suppress(ValueError, TypeError):
                pos_map[int(pos)] = gt

    return rsid_map, pos_map


# ── Node lookup (multi-marker majority vote) ────────────────────────────────────

def _lookup(node: dict, rsid_map: dict, pos_map: dict) -> tuple[str, int, int, list[str]]:
    """
    Vote across all of a node's defining markers.

    Returns (status, n_derived, n_ancestral, evidence):
      status ∈ {'derived', 'ancestral', 'conflict', 'not_found'}
      evidence : human-readable per-marker calls, e.g. ['M89:T✓', 'P135:C–']

    Strand-ambiguous markers (ancestral/derived a complementary pair) are skipped
    because the observed base cannot be oriented.
    """
    n_der = n_anc = 0
    evidence: list[str] = []

    for marker in node.get("markers", []):
        name, rsid, pos, anc, der = marker
        anc = anc.upper()
        der = der.upper()
        if _comp(anc) == der:        # strand-ambiguous → uninformative
            continue

        gt = ""
        if rsid:
            gt = rsid_map.get(rsid, "")
        if not gt and pos:
            gt = pos_map.get(int(pos), "")
        if not gt:
            continue

        g = gt[0]
        if g == der or g == _comp(der):
            n_der += 1
            evidence.append(f"{name}:{g}✓")
        elif g == anc or g == _comp(anc):
            n_anc += 1
            evidence.append(f"{name}:{g}–")

    if n_der and n_der > n_anc:
        status = "derived"
    elif n_anc and n_anc > n_der:
        status = "ancestral"
    elif n_der and n_anc:
        status = "conflict"          # equal split — treat as a gap, do not assert
    else:
        status = "not_found"
    return status, n_der, n_anc, evidence


# ── Tree walker ────────────────────────────────────────────────────────────────

def _make_entry(node: dict, status: str, n_der: int, n_anc: int, evidence: list[str]) -> dict:
    return {
        "haplogroup": node["haplogroup"],
        "snp_name": node["snp_name"],
        "rsids": node.get("rsids", []),
        "pos": node.get("pos"),
        "description": node.get("description", ""),
        "migration": node.get("migration", ""),
        "further": node.get("further", ""),
        "children": node.get("children", []),
        "snp_status": "confirmed" if status == "derived" else "chip_gap",
        "found_genotype": ", ".join(evidence) if evidence else "–",
        "n_derived": n_der,
        "n_ancestral": n_anc,
        "contradiction": [],
    }


def _walk(node: dict, rsid_map: dict, pos_map: dict, prefix: list[dict], depth: int):
    """
    Recursive walk. Returns (path, status) or None if this node is ruled out
    (ancestral). Never descends past a node unless a downstream marker is
    genuinely confirmed, so the reported lineage cannot run deeper than the data.
    """
    if depth > 30:
        return prefix, "max_depth"

    status, n_der, n_anc, evidence = _lookup(node, rsid_map, pos_map)
    if status == "ancestral":
        return None                                  # branch excluded

    entry = _make_entry(node, status, n_der, n_anc, evidence)
    path = [*prefix, entry]
    children = node.get("children", [])

    if not children:
        return path, ("resolved" if status == "derived" else "chip_gap")

    confirmed: list[tuple[dict, int]] = []
    gaps: list[dict] = []
    for child in children:
        cs, cder, _canc, _ev = _lookup(child, rsid_map, pos_map)
        if cs == "derived":
            confirmed.append((child, cder))
        elif cs in ("not_found", "conflict"):
            gaps.append(child)
        # ancestral children are excluded

    if confirmed:
        confirmed.sort(key=lambda x: -x[1])
        if len(confirmed) > 1:
            # Biologically impossible: a man belongs to ONE lineage. Flag it.
            entry["contradiction"] = [c["haplogroup"] for c, _ in confirmed]
        res = _walk(confirmed[0][0], rsid_map, pos_map, path, depth + 1)
        if res is None:
            return path, ("resolved" if status == "derived" else "chip_gap")
        return res

    # Only chip-gap children: descend only toward a confirmed descendant.
    best = None
    best_confirmed = 0
    branches_with_support: list[str] = []
    for child in gaps:
        res = _walk(child, rsid_map, pos_map, path, depth + 1)
        if res is None:
            continue
        rpath, _rstatus = res
        # count confirmations strictly BELOW the current node
        below = sum(1 for n in rpath[len(path):] if n["snp_status"] == "confirmed")
        if below > 0:
            branches_with_support.append(child["haplogroup"])
        if below > best_confirmed:
            best_confirmed = below
            best = res
    if best and best_confirmed > 0:
        if len(branches_with_support) > 1:
            # Confirmed markers in two mutually-exclusive branches — impossible.
            entry["contradiction"] = sorted(set(branches_with_support))
        return best

    return path, ("resolved" if status == "derived" else "chip_gap")


# ── Public API ─────────────────────────────────────────────────────────────────

def _empty(status: str, y_count: int, message: str, **extra) -> dict:
    base = {
        "status": status,
        "confidence": "none",
        "y_snp_count": y_count,
        "haplogroup_path": "Unknown",
        "terminal_haplogroup": "Unknown",
        "terminal_description": "",
        "terminal_migration": None,
        "path": [],
        "chip_gaps": [],
        "not_tested_branches": [],
        "further_testing": "",
        "message": message,
    }
    base.update(extra)
    return base


def analyze_y_haplogroup(snps_df: pd.DataFrame) -> dict:
    """
    Analyse Y-chromosome SNPs and return a result dict. The terminal haplogroup
    reported is always the DEEPEST CONFIRMED node — markers that are merely
    inferred across chip gaps never drive the call.
    """
    rsid_map, pos_map = _build_lookup(snps_df)
    y_count = len(rsid_map)

    if y_count == 0:
        return _empty(
            "no_y_data", 0,
            "No Y-chromosome SNPs found in this file. Possible reasons: the person is "
            "female, the chip does not genotype Y markers, or Y-SNPs are stored under "
            "an unrecognised chromosome label.",
        )

    if y_count < 10:
        return _empty(
            "insufficient_y_data", y_count,
            f"Insufficient Y-chromosome SNPs ({y_count} found, minimum 10 required) "
            "to make a reliable haplogroup call.",
            haplogroup_path="Insufficient Y-chromosome SNPs",
            terminal_haplogroup="Insufficient Y-chromosome SNPs",
        )

    res = _walk(HAPLOGROUP_TREE, rsid_map, pos_map, [], 0)

    if res is None:
        # Root (CT) is ancestral → haplogroup A or B (pre-CT), or no usable markers.
        return _empty(
            "pre_ct", y_count,
            "Y-chromosome does not carry the CT (M168) derived allele — the lineage "
            "appears to belong to haplogroup A or B, the deepest branches of the human "
            "Y tree. Specialist testing is needed to confirm and place it.",
            haplogroup_path="A or B (pre-CT)",
            terminal_haplogroup="A or B",
            further_testing="Contact a Y-chromosome phylogenetics specialist or FTDNA.",
        )

    full_path, walk_status = res

    confirmed_idx = [i for i, n in enumerate(full_path) if n["snp_status"] == "confirmed"]
    if not confirmed_idx:
        # No marker confirmed anywhere — cannot assert anything beyond the root.
        gap_names = [n["snp_name"] for n in full_path]
        return _empty(
            "partial", y_count,
            "Could not confirm any haplogroup-defining marker on this chip.",
            haplogroup_path="Unresolved (defining markers not on chip)",
            chip_gaps=gap_names,
            further_testing=HAPLOGROUP_TREE.get("further", ""),
        )

    # Terminal = deepest CONFIRMED node. Trailing inferred guesses are dropped.
    deepest = confirmed_idx[-1]
    path = full_path[: deepest + 1]
    terminal = path[-1]

    confirmed_count = sum(1 for n in path if n["snp_status"] == "confirmed")
    inferred = [n["snp_name"] for n in path if n["snp_status"] == "chip_gap"]
    contradictions = [c for n in path for c in n.get("contradiction", [])]

    # Untested branches immediately below the terminal node.
    not_tested: list[dict] = []
    for child in terminal.get("children", []):
        cs, _d, _a, _ev = _lookup(child, rsid_map, pos_map)
        if cs in ("not_found", "conflict"):
            not_tested.append({
                "haplogroup": child["haplogroup"],
                "snp_name": child["snp_name"],
                "description": child.get("description", ""),
                "further": child.get("further", ""),
            })

    is_leaf = not terminal.get("children")
    resolved = is_leaf or (walk_status == "resolved")

    if contradictions:
        confidence = "low"
    elif confirmed_count >= 3:
        confidence = "high"
    elif confirmed_count >= 2:
        confidence = "moderate"
    else:
        confidence = "low"

    confidence_note = (
        f"{confirmed_count} of {len(path)} markers on the assigned path confirmed "
        f"(derived by majority vote of co-defining SNPs); the remainder inferred across "
        f"chip gaps. {y_count:,} Y-SNPs typed total."
    )
    if contradictions:
        confidence_note += (
            "  WARNING: conflicting derived calls were seen for branches that are "
            f"mutually exclusive ({', '.join(sorted(set(contradictions)))}) — the "
            "underlying genotypes or marker definitions may be unreliable here."
        )

    return {
        "status": "resolved" if resolved else "partial",
        "confidence": confidence,
        "confidence_note": confidence_note,
        "y_snp_count": y_count,
        "n_markers_confirmed": confirmed_count,
        "n_markers_on_path": len(path),
        "haplogroup_path": " > ".join(n["haplogroup"] for n in path),
        "terminal_haplogroup": terminal["haplogroup"],
        "terminal_description": terminal.get("description", ""),
        "terminal_migration": terminal.get("migration", ""),
        "path": path,
        "chip_gaps": inferred,
        "not_tested_branches": not_tested,
        "contradictions": sorted(set(contradictions)),
        "further_testing": terminal.get("further", "") or HAPLOGROUP_TREE.get("further", ""),
        "message": (
            f"Y-DNA: {terminal['haplogroup']} "
            f"({'resolved' if resolved else 'partial — chip gaps'}, "
            f"{confirmed_count} confirmed marker{'s' if confirmed_count != 1 else ''})"
        ),
    }
