"""
Y-DNA Haplogroup Analysis Module
Walks a decision tree of Y-chromosome haplogroup-defining SNPs, focused on
Macro-haplogroup K and all downstream subclades (K1, K2, K2a, K2b → N, O, Q,
P → R1a, R1b and its major European branches).

Tree-walking algorithm:
  1. Look up each node's defining SNP by rsID, then by chr Y position.
  2. DERIVED   → confirmed in this haplogroup; recurse into children.
  3. ANCESTRAL → this branch is ruled out; stop.
  4. NOT FOUND → chip gap; note it, still attempt children so downstream
                  confirmed markers can be reported (they implicitly confirm
                  every ancestor).
  5. After walking, any gap-walked ancestors are labelled "inferred".
"""

from typing import Optional, Dict, List, Tuple
import pandas as pd

# ── Complement helper ──────────────────────────────────────────────────────────
_COMP = str.maketrans("ACGT", "TGCA")


def complement(base: str) -> str:
    return base.upper().translate(_COMP)


# ── Haplogroup decision tree ───────────────────────────────────────────────────
#
# Each node:
#   haplogroup  – name shown in report (e.g. "R1b")
#   snp_name    – ISOGG marker name (e.g. "M343")
#   rsids       – list of rsIDs to try, in priority order (chip-version-aware)
#   pos         – GRCh37 chr Y position for position-based fallback (int | None)
#   derived     – derived allele character (uppercase, forward strand)
#   ancestral   – ancestral allele character (uppercase, forward strand)
#   description – one-liner shown next to breadcrumb node
#   migration   – migration narrative paragraph shown in the report section
#   further     – what Big Y / dedicated panel would add
#   children    – mutually exclusive downstream branches (list of nodes)
#
# Position values are approximate (GRCh37/hg19) from ISOGG Y-DNA SNP Index.
# rsID values are those most commonly cited in consumer-chip contexts; the same
# physical mutation sometimes appears under different rsIDs in different chip
# versions — listing several maximises hit rate.

HAPLOGROUP_TREE: Dict = {
    "haplogroup": "K",
    "snp_name": "M9",
    "rsids": ["rs2032597"],
    "pos": 22_719_028,
    "derived": "T",
    "ancestral": "G",
    "description": "Macro-haplogroup K — ancestor of most non-African men",
    "migration": (
        "Macro-haplogroup K (M9) emerged approximately 45,000–47,000 years ago, most "
        "likely in South or Central Asia, shortly after the 'Out of Africa' dispersal. "
        "A single man carried the M9 mutation, and today his Y-chromosome descendants "
        "account for the vast majority of men outside sub-Saharan Africa — including "
        "virtually all men of European, East Asian, South Asian, Oceanian, and Native "
        "American descent. Haplogroups N, O, Q, R, S, M, and T all trace back to this "
        "common K ancestor."
    ),
    "further": (
        "FTDNA Big Y-700 or equivalent long-read Y-chromosome sequencing would confirm "
        "K and immediately resolve which downstream branch you belong to. If your chip "
        "does not include M9, any confirmed downstream K marker (M207, M175, M231, "
        "M242, M343, M269 …) implicitly confirms K ancestry."
    ),
    "children": [
        # ── T (M70) ────────────────────────────────────────────────────────────
        {
            "haplogroup": "T",
            "snp_name": "M70",
            "rsids": ["rs9786474", "rs2032605"],
            "pos": 21_569_730,
            "derived": "C",
            "ancestral": "A",
            "description": "Haplogroup T — Middle East, East Africa, Mediterranean",
            "migration": (
                "Haplogroup T (M70) is an ancient non-African lineage found at moderate "
                "frequencies in East Africa (particularly Ethiopia and the Horn), the Middle "
                "East, Anatolia, and Mediterranean Europe. It is especially notable among the "
                "Lemba people of southern Africa and certain Ethiopian populations. T reached "
                "southern Europe via Neolithic expansions from the Fertile Crescent ~8,000–10,000 "
                "years ago, and remains a minor haplogroup in Europe — more common in Sardinia, "
                "the Canary Islands, and parts of Italy. Famous proposed T carriers include "
                "Thomas Jefferson (T1a1)."
            ),
            "further": (
                "Big Y-700 would identify your T subclade: T1a1 (Middle Eastern / Lemba), "
                "T1a2 (South Asian), or T1b (rare). FTDNA's Haplogroup T Project has thousands "
                "of members with detailed geographic correlations."
            ),
            "children": [],
        },
        # ── K1 (P226) ──────────────────────────────────────────────────────────
        {
            "haplogroup": "K1",
            "snp_name": "P226",
            "rsids": [],
            "pos": None,
            "derived": "T",
            "ancestral": "C",
            "description": "Haplogroup K1 — rare Oceanian/Australian branch",
            "migration": (
                "Haplogroup K1 (P226, P228, P230) is one of the rarest Y-chromosome lineages "
                "in the world, found almost exclusively among Aboriginal Australians, some "
                "Melanesian populations, and very rarely in Southeast Asia. It represents an "
                "extremely ancient split from macro-haplogroup K that accompanied the earliest "
                "human migrations into the Sahul landmass (Australia and New Guinea) roughly "
                "50,000–65,000 years ago. K1 is essentially absent in all European, East Asian, "
                "and African populations — finding K1 outside Oceania would be remarkable."
            ),
            "further": (
                "K1 is too rare and region-specific to be resolved further with consumer chips. "
                "Whole Y-chromosome sequencing and consultation with specialists in Australian "
                "Aboriginal genetics would be needed."
            ),
            "children": [],
        },
        # ── K2a (P295) — ancestor of S and M ───────────────────────────────────
        {
            "haplogroup": "K2a",
            "snp_name": "P295",
            "rsids": [],
            "pos": 6_826_764,
            "derived": "C",
            "ancestral": "T",
            "description": "Haplogroup K2a — ancestor of S and M (Melanesia)",
            "migration": (
                "Haplogroup K2a is the common ancestor of haplogroups S and M, both found "
                "almost exclusively in Papua New Guinea, Melanesia, and nearby island groups. "
                "K2a lineages represent an ancient migration into the Sahul landmass and are "
                "essentially absent outside Oceania."
            ),
            "further": (
                "K2a subclades S and M are highly region-specific. Specialized Y-chromosome "
                "sequencing and collaboration with Melanesian genetic research groups would "
                "be needed to resolve your specific subclade."
            ),
            "children": [
                {
                    "haplogroup": "S",
                    "snp_name": "M230",
                    "rsids": [],
                    "pos": 14_843_410,
                    "derived": "T",
                    "ancestral": "G",
                    "description": "Haplogroup S — Highland Papua New Guinea",
                    "migration": (
                        "Haplogroup S (B254/M230) is the dominant Y-chromosome lineage in the "
                        "highlands of Papua New Guinea, where it reaches frequencies exceeding "
                        "50% in some populations. It traces back to the founding populations of "
                        "Sahul and has diversified extensively over the past 40,000+ years of "
                        "relative isolation. S subclades correlate closely with linguistic "
                        "groups and geographic regions within Papua New Guinea."
                    ),
                    "further": (
                        "Full Y-chromosome sequencing is the only way to resolve S subclades "
                        "meaningfully. The S tree is still being built by research groups "
                        "studying Papuan genetics."
                    ),
                    "children": [],
                },
                {
                    "haplogroup": "M",
                    "snp_name": "M4",
                    "rsids": [],
                    "pos": 8_602_381,
                    "derived": "C",
                    "ancestral": "T",
                    "description": "Haplogroup M — Melanesia and eastern Indonesia",
                    "migration": (
                        "Haplogroup M (M4/P256) is found primarily in Papua New Guinea, "
                        "Melanesia, and the Maluku Islands of eastern Indonesia. Like S, it "
                        "represents a deep founding Sahul lineage. Its subclades (M1–M5) track "
                        "geographic and linguistic boundaries across Melanesian islands."
                    ),
                    "further": (
                        "Specialized Y-chromosome sequencing and FTDNA's Polynesian/Melanesian "
                        "projects would provide subclade resolution for M."
                    ),
                    "children": [],
                },
            ],
        },
        # ── K2b (M526) — the big branch: N, O, Q, P → R ────────────────────────
        {
            "haplogroup": "K2b",
            "snp_name": "M526",
            "rsids": [],
            "pos": 15_017_066,
            "derived": "T",
            "ancestral": "C",
            "description": "Haplogroup K2b — ancestor of N, O, Q, and R (most Eurasian men)",
            "migration": (
                "Haplogroup K2b (M526) is the common ancestor of haplogroups N, O, Q, and R — "
                "the four lineages that dominate Eurasian, East Asian, Oceanian, and Native "
                "American Y-chromosomes. K2b likely arose in Central or South Asia approximately "
                "40,000–45,000 years ago. From there, its descendants dispersed in all "
                "directions: N and O moved eastward into Asia, Q eventually crossed into the "
                "Americas, and the P/R branch spread across Eurasia with the later "
                "Indo-European expansions."
            ),
            "further": (
                "M526 is rarely included on consumer chips. Any confirmed downstream marker "
                "(N-M231, O-M175, Q-M242, or R-M207) implicitly confirms K2b. Big Y-700 "
                "would definitively resolve which K2b branch you are in."
            ),
            "children": [
                # ── N (M231) ──────────────────────────────────────────────────
                {
                    "haplogroup": "N",
                    "snp_name": "M231",
                    "rsids": ["rs9341278", "rs9785941", "rs2032630"],
                    "pos": 9_388_483,
                    "derived": "T",
                    "ancestral": "C",
                    "description": "Haplogroup N — North Eurasia, Finland, Siberia, Turkic peoples",
                    "migration": (
                        "Haplogroup N (M231) spans a vast arc from northeastern Europe across "
                        "Siberia to East Asia and is the defining Y-chromosome lineage of "
                        "Uralic-speaking peoples. Today ~60% of Finnish men, ~67% of Estonian "
                        "men, and high fractions of Sami, Nenets, Selkup, and Yakut men carry N. "
                        "N likely originated in East or Central Asia ~20,000–30,000 years ago and "
                        "expanded westward into northeastern Europe during the Mesolithic or "
                        "Neolithic, probably carried by ancestral Uralic speakers. Its European "
                        "expansion correlates with the spread of Finnish, Estonian, and "
                        "historically related languages. Major subclades include N1a1a (formerly "
                        "N3, dominant in Finnic and Baltic peoples) and N1b (Siberian)."
                    ),
                    "further": (
                        "Big Y-700 + FTDNA's Haplogroup N Project would resolve N1a vs N1b and "
                        "identify geographically specific branches (Finnish N1a1a1, Siberian N1b, "
                        "etc.). N1a1a is well-studied with hundreds of subclades mapped to "
                        "specific Uralic language groups."
                    ),
                    "children": [
                        # ── N1 (CTS11726/L735) ────────────────────────────────
                        {
                            "haplogroup": "N1",
                            "snp_name": "CTS11726/L735",
                            "rsids": [],
                            "pos": 18_380_393,
                            "derived": "G",
                            "ancestral": "A",
                            "description": "Haplogroup N1 — main branch covering nearly all living N men",
                            "migration": (
                                "N1 (CTS11726, equivalent SNP L735) contains essentially every "
                                "living N-bearing man outside of a handful of rare N2 lineages. "
                                "The split between N1 and the sibling N2 branch occurred deep in "
                                "Asia ~20,000+ years ago. From N1, two major subclades emerged: "
                                "N1a (which became the dominant northern Eurasian / Uralic "
                                "lineage) and N1b (Siberian)."
                            ),
                            "further": (
                                "CTS11726 and L735 are typically only on Big Y-700 or dedicated "
                                "Y-SNP panels — not on consumer autosomal+Y chips."
                            ),
                            "children": [
                                # ── N1a (M2291/F1206) ─────────────────────────
                                {
                                    "haplogroup": "N1a",
                                    "snp_name": "M2291/F1206",
                                    "rsids": [],
                                    "pos": 8_430_571,
                                    "derived": "G",
                                    "ancestral": "A",
                                    "description": "Haplogroup N1a — ancestor of the Finnic/Uralic N1a1 expansion",
                                    "migration": (
                                        "N1a (M2291, equivalent F1206) is the parent of N1a1 "
                                        "(M46/Tat), the lineage that came to dominate northeastern "
                                        "European and Uralic-speaking populations. N1a likely "
                                        "diverged in Central or East Asia and spread westward "
                                        "across Siberia before the major N1a1 expansion."
                                    ),
                                    "further": (
                                        "M2291/F1206 are Big Y-discovered SNPs without dbSNP IDs "
                                        "and almost never appear on consumer chips."
                                    ),
                                    "children": [
                                        # ── N1a1 (M46/Tat) ────────────────────
                                        {
                                            "haplogroup": "N1a1",
                                            "snp_name": "M46/Tat",
                                            "rsids": ["rs2032673", "rs9785945"],
                                            "pos": 14_179_811,
                                            "derived": "T",
                                            "ancestral": "C",
                                            "description": "Haplogroup N1a1 — Tat-positive, classical 'northern N'",
                                            "migration": (
                                                "N1a1 (M46/Tat) is the famous Tat-C lineage of "
                                                "Finnic, Baltic, Sami, and northern Russian "
                                                "populations. It expanded westward across Siberia "
                                                "into northeastern Europe during the Bronze Age, "
                                                "carried by ancestral Uralic speakers, and is the "
                                                "single most common Y-haplogroup among Finns "
                                                "(~58%), Estonians (~34%), and Sami (~40%)."
                                            ),
                                            "further": (
                                                "Tat-positive testing is widely available, but "
                                                "downstream M178/L708/CTS9976/L1026/Z1936 "
                                                "resolution requires Big Y-700."
                                            ),
                                            "children": [
                                                # ── N1a1a (M178) ──────────────
                                                {
                                                    "haplogroup": "N1a1a",
                                                    "snp_name": "M178",
                                                    "rsids": ["rs367573274"],
                                                    "pos": 21_717_307,
                                                    "derived": "A",
                                                    "ancestral": "T",
                                                    "description": "N1a1a — the dominant European N branch (formerly 'N3')",
                                                    "migration": (
                                                        "N1a1a (M178) is the European face of N: "
                                                        "almost all Finnish, Estonian, Lithuanian, "
                                                        "Latvian, and northern Russian N men belong "
                                                        "here. M178 marks the founder lineage that "
                                                        "swept westward into the eastern Baltic "
                                                        "~2,500–3,500 years ago, very likely with "
                                                        "early Finno-Ugric speakers."
                                                    ),
                                                    "further": (
                                                        "Big Y-700 resolves M178 into L708 → "
                                                        "CTS9976 → L1026 (Finnish/Baltic) or Z1936 "
                                                        "(Siberian/Uralic) and onward into terminal "
                                                        "branches with surname-level resolution."
                                                    ),
                                                    "children": [
                                                        # ── N1a1a1 (L708) ─────
                                                        {
                                                            "haplogroup": "N1a1a1",
                                                            "snp_name": "L708",
                                                            "rsids": [],
                                                            "pos": 13_954_389,
                                                            "derived": "A",
                                                            "ancestral": "G",
                                                            "description": "N1a1a1 — phylo-equivalent layer below M178",
                                                            "migration": (
                                                                "N1a1a1 (L708) sits just below M178 "
                                                                "and contains essentially all living "
                                                                "M178-positive men. It is one of "
                                                                "several phylo-equivalent SNPs in "
                                                                "this region."
                                                            ),
                                                            "further": (
                                                                "L708 was discovered by FTDNA's Walk "
                                                                "Through the Y / Big Y. It has no "
                                                                "dbSNP rsID and is not on any consumer "
                                                                "chip."
                                                            ),
                                                            "children": [
                                                                # ── N1a1a1a (CTS9976) ───
                                                                {
                                                                    "haplogroup": "N1a1a1a",
                                                                    "snp_name": "CTS9976",
                                                                    "rsids": [],
                                                                    "pos": 17_713_080,
                                                                    "derived": "T",
                                                                    "ancestral": "C",
                                                                    "description": "N1a1a1a — parent of the L1026 vs Z1936 split",
                                                                    "migration": (
                                                                        "CTS9976 defines the node "
                                                                        "above the major L1026 "
                                                                        "(Finnic/Baltic) vs Z1936 "
                                                                        "(Siberian/Uralic) split that "
                                                                        "occurred ~4,000 years ago in "
                                                                        "western Siberia or the Urals."
                                                                    ),
                                                                    "further": (
                                                                        "CTS-series SNPs are Big "
                                                                        "Y-only markers."
                                                                    ),
                                                                    "children": [
                                                                        {
                                                                            "haplogroup": "N1a1a1a1",
                                                                            "snp_name": "L1026",
                                                                            "rsids": [],
                                                                            "pos": 16_887_278,
                                                                            "derived": "G",
                                                                            "ancestral": "A",
                                                                            "description": "N1a1a1a1 — Finnish/Baltic core lineage",
                                                                            "migration": (
                                                                                "L1026 is the defining "
                                                                                "SNP of the Finnish "
                                                                                "and Baltic N branch. "
                                                                                "~50% of Finnish men "
                                                                                "and very high "
                                                                                "fractions of Estonian, "
                                                                                "Latvian, and "
                                                                                "Lithuanian men are "
                                                                                "L1026-positive. Its "
                                                                                "subclades (VL29, "
                                                                                "Z1934, CTS1737) map "
                                                                                "cleanly onto Finnic, "
                                                                                "Baltic, and Russian-"
                                                                                "speaking populations."
                                                                            ),
                                                                            "further": (
                                                                                "Big Y-700 + FTDNA's "
                                                                                "N-L1026 Project place "
                                                                                "you in a terminal "
                                                                                "branch (e.g. N-VL29, "
                                                                                "N-Z1934, N-CTS1737) "
                                                                                "linked to specific "
                                                                                "Finnic or Baltic "
                                                                                "regional lineages."
                                                                            ),
                                                                            "children": [],
                                                                        },
                                                                        {
                                                                            "haplogroup": "N1a1a1a2",
                                                                            "snp_name": "Z1936",
                                                                            "rsids": [],
                                                                            "pos": 16_439_111,
                                                                            "derived": "C",
                                                                            "ancestral": "T",
                                                                            "description": "N1a1a1a2 — Siberian/Uralic branch (Sami, Volga-Uralic)",
                                                                            "migration": (
                                                                                "Z1936 is the parallel "
                                                                                "branch to L1026, "
                                                                                "dominant among Sami "
                                                                                "(~40%+), Khanty, "
                                                                                "Mansi, and Volga-"
                                                                                "Uralic peoples "
                                                                                "(Mari, Udmurt, Komi). "
                                                                                "It traces a more "
                                                                                "eastern path through "
                                                                                "the Urals into "
                                                                                "northern Fennoscandia."
                                                                            ),
                                                                            "further": (
                                                                                "Big Y-700 + FTDNA's "
                                                                                "N-Z1936 Project "
                                                                                "resolves Sami-specific "
                                                                                "vs Volga-Uralic "
                                                                                "branches."
                                                                            ),
                                                                            "children": [],
                                                                        },
                                                                    ],
                                                                },
                                                            ],
                                                        },
                                                    ],
                                                },
                                            ],
                                        },
                                    ],
                                },
                                # ── N1b (P43) ─────────────────────────────────
                                {
                                    "haplogroup": "N1b",
                                    "snp_name": "P43",
                                    "rsids": [],
                                    "pos": 21_466_748,
                                    "derived": "C",
                                    "ancestral": "T",
                                    "description": "Haplogroup N1b — Siberian / Samoyedic N branch",
                                    "migration": (
                                        "N1b (P43) is the sister branch to N1a and is concentrated "
                                        "in Samoyedic peoples (Nenets, Nganasan, Selkup) and "
                                        "northern Siberian populations, with notable frequencies "
                                        "among Northern Khanty and some Turkic groups. It is "
                                        "virtually absent in Europe. N1b diverged from N1a deep "
                                        "in Siberia and stayed east of the Urals."
                                    ),
                                    "further": (
                                        "P43 is on a few dedicated Y-SNP panels but rarely on "
                                        "consumer chips. Big Y-700 resolves P43 into B187 vs B188 "
                                        "and onward."
                                    ),
                                    "children": [
                                        {
                                            "haplogroup": "N1b1",
                                            "snp_name": "B187",
                                            "rsids": [],
                                            "pos": None,
                                            "derived": "T",
                                            "ancestral": "C",
                                            "description": "N1b1 — one of two major P43 sub-branches",
                                            "migration": (
                                                "B187 defines one of the two main P43 lineages, "
                                                "found in Nganasan and some Selkup populations."
                                            ),
                                            "further": (
                                                "B-series SNPs were discovered through Big Y / "
                                                "academic WGS and have no dbSNP rsIDs. Only Big "
                                                "Y-700 (or equivalent WGS) can call them."
                                            ),
                                            "children": [],
                                        },
                                        {
                                            "haplogroup": "N1b2",
                                            "snp_name": "B188",
                                            "rsids": [],
                                            "pos": None,
                                            "derived": "T",
                                            "ancestral": "C",
                                            "description": "N1b2 — second major P43 sub-branch",
                                            "migration": (
                                                "B188 is the parallel P43 sub-branch to B187, "
                                                "concentrated in Nenets and northern Samoyedic "
                                                "groups."
                                            ),
                                            "further": (
                                                "Big Y-700 is required — B188 is not on any "
                                                "consumer chip."
                                            ),
                                            "children": [],
                                        },
                                    ],
                                },
                            ],
                        },
                        # ── N2 (Y6503) ────────────────────────────────────────
                        {
                            "haplogroup": "N2",
                            "snp_name": "Y6503",
                            "rsids": [],
                            "pos": None,
                            "derived": "T",
                            "ancestral": "C",
                            "description": "Haplogroup N2 — rare deep-rooted sister of N1",
                            "migration": (
                                "N2 (Y6503) is an extremely rare deep branch of N, sister to "
                                "the entire N1 clade. It has been reported in a handful of "
                                "Vietnamese, Han Chinese, and other East Asian samples and is "
                                "essentially absent in Europe and Siberia. Y6503 was identified "
                                "by YFull from full Y-sequencing of academic samples."
                            ),
                            "further": (
                                "Y6503 has no dbSNP rsID and is only callable via Big Y-700 or "
                                "full Y-chromosome sequencing. Anyone landing here would be a "
                                "valuable contribution to FTDNA's Haplogroup N Project."
                            ),
                            "children": [],
                        },
                    ],
                },
                # ── O (M175) ──────────────────────────────────────────────────
                {
                    "haplogroup": "O",
                    "snp_name": "M175",
                    "rsids": ["rs2032658"],
                    "pos": 13_907_985,
                    "derived": "A",
                    "ancestral": "G",
                    "description": "Haplogroup O — dominant East and Southeast Asian lineage",
                    "migration": (
                        "Haplogroup O (M175) accounts for 80–90% of Y-chromosomes in East and "
                        "Southeast Asia — Chinese, Japanese, Korean, Vietnamese, Thai, Malay, and "
                        "Filipino populations all show very high O frequencies. O likely emerged "
                        "in Southeast Asia or southern China ~30,000–40,000 years ago. Its major "
                        "expansion correlates with the rise of rice agriculture in the Yangtze "
                        "River valley and subsequent Han Chinese demographic growth. Subclades: "
                        "O1 (Southeast Asia, Austronesian speakers, Taiwan Aboriginal), O2 (Han "
                        "Chinese, Korean, Japanese), O3 (widespread across East and Southeast "
                        "Asia, associated with Sino-Tibetan expansions). O is essentially absent "
                        "in Europe, Africa, and the Americas."
                    ),
                    "further": (
                        "Big Y-700 or dedicated East Asian Y-DNA panels would identify your "
                        "O subclade (O1a, O1b, O2, O3) with geographic precision, distinguishing "
                        "between different East and Southeast Asian ethnic groups."
                    ),
                    "children": [],
                },
                # ── Q (M242) ──────────────────────────────────────────────────
                {
                    "haplogroup": "Q",
                    "snp_name": "M242",
                    "rsids": ["rs2032646"],
                    "pos": 8_028_953,
                    "derived": "T",
                    "ancestral": "C",
                    "description": "Haplogroup Q — founding haplogroup of most Native Americans",
                    "migration": (
                        "Haplogroup Q (M242) is most famous as the founding Y-chromosome lineage "
                        "of the majority of indigenous Americans. Q-bearing populations crossed "
                        "the Bering land bridge approximately 15,000–20,000 years ago, and in "
                        "the Americas, Q-M3 (a subclade) reaches 90%+ frequency in many "
                        "indigenous groups from Alaska to Tierra del Fuego. Outside the Americas, "
                        "Q is found in Siberia (Ket, Selkup, Yeniseian speakers), Central Asia, "
                        "and at low frequency in South Asia, the Middle East, and Europe. "
                        "Notably, ~5% of Ashkenazi Jewish men carry Q1b2, a Middle Eastern "
                        "branch with a distinct history from the Native American Q1a branch."
                    ),
                    "further": (
                        "Big Y-700 would distinguish Q1a (Siberian/Native American) from Q1b "
                        "(Middle Eastern/Ashkenazi) and resolve specific subclade with geographic "
                        "and population significance. The Q-M3 branch is essentially "
                        "diagnostic of Native American ancestry."
                    ),
                    "children": [],
                },
                # ── P (P331) — ancestor of R ────────────────────────────────
                {
                    "haplogroup": "P",
                    "snp_name": "P331",
                    "rsids": ["rs2032652"],
                    "pos": 22_738_626,
                    "derived": "C",
                    "ancestral": "A",
                    "description": "Haplogroup P — ancestor of R (major Eurasian lineage)",
                    "migration": (
                        "Haplogroup P (P331/M45/P226) is the direct ancestor of haplogroup R — "
                        "the dominant Y-chromosome lineage in Europe and Central Asia. P itself "
                        "likely arose in Central or South Asia ~35,000–40,000 years ago, and "
                        "virtually all living P-bearing men today belong to downstream R subclades. "
                        "P* (without R) is extremely rare."
                    ),
                    "further": (
                        "Finding P without R confirmed downstream would be unusual. Big Y-700 "
                        "would immediately resolve R1a vs R1b placement."
                    ),
                    "children": [
                        # ── R (M207) ──────────────────────────────────────────
                        {
                            "haplogroup": "R",
                            "snp_name": "M207",
                            "rsids": ["rs9785952", "rs2032658"],
                            "pos": 8_019_819,
                            "derived": "T",
                            "ancestral": "C",
                            "description": "Haplogroup R — dominant European and Central/South Asian lineage",
                            "migration": (
                                "Haplogroup R (M207) is one of the most widespread Y-chromosome "
                                "haplogroups in Eurasia. It arose ~27,000–35,000 years ago, probably "
                                "in Central Asia or South Asia, before the Last Glacial Maximum. "
                                "After the LGM, R diversified rapidly into R1a and R1b. R1a spread "
                                "eastward and westward with the Indo-European expansions from the "
                                "Pontic-Caspian steppe, while R1b came to dominate Western Europe "
                                "following the Yamnaya expansion ~5,000 years ago. Together, R1a "
                                "and R1b account for roughly 60–70% of all European Y-chromosomes."
                            ),
                            "further": (
                                "Big Y-700 would immediately resolve R1a vs R1b, which have very "
                                "different geographic distributions and migration histories. "
                                "Both branches have been deeply studied and have dedicated FTDNA "
                                "project pages with tens of thousands of members."
                            ),
                            "children": [
                                # ── R1 (M173) ─────────────────────────────────
                                {
                                    "haplogroup": "R1",
                                    "snp_name": "M173",
                                    "rsids": ["rs9786153"],
                                    "pos": 13_470_467,
                                    "derived": "T",
                                    "ancestral": "C",
                                    "description": "Haplogroup R1 — European and Central/South Asian R",
                                    "migration": (
                                        "Haplogroup R1 contains nearly all R-bearing men in Europe, "
                                        "Central Asia, and South Asia, split into R1a (eastward "
                                        "expansions, Eastern European and South Asian) and R1b "
                                        "(westward expansions, Western European)."
                                    ),
                                    "further": (
                                        "Distinguishing R1a (M420) from R1b (M343) is the critical "
                                        "next step. Many consumer chips include markers for this "
                                        "split. M420 and M343 are the defining SNPs to look for."
                                    ),
                                    "children": [
                                        # ── R1a (M420) ────────────────────────
                                        {
                                            "haplogroup": "R1a",
                                            "snp_name": "M420",
                                            "rsids": ["rs2032655", "rs3908938"],
                                            "pos": 14_734_005,
                                            "derived": "T",
                                            "ancestral": "C",
                                            "description": "Haplogroup R1a — Indo-European steppe expansion (Eastern Europe, South Asia)",
                                            "migration": (
                                                "Haplogroup R1a (M420/M17) is the signature Y-DNA lineage "
                                                "of the Indo-European expansion from the Pontic-Caspian "
                                                "steppe. Today it is found at highest frequencies in Eastern "
                                                "Europe: ~55–65% in Poland, ~40–50% in Russia and Ukraine, "
                                                "~45% in the Czech Republic, and up to 70%+ in certain "
                                                "South Asian Brahmin caste populations. R1a entered Europe "
                                                "~4,900 years ago with the Corded Ware culture, and reached "
                                                "South Asia with Indo-Aryan migrations ~3,500 years ago. "
                                                "R1a-Z282 (Z283, M458, Z280) defines the Eastern European "
                                                "branch; R1a-Z93 (Z94, M780, L342) defines the South/Central "
                                                "Asian branch. R1a men are descendants of the Yamnaya steppe "
                                                "herders who revolutionised Eurasian prehistory with horses, "
                                                "wheeled vehicles, and Proto-Indo-European language."
                                            ),
                                            "further": (
                                                "Big Y-700 would distinguish R1a-Z282 (Eastern European) "
                                                "from R1a-Z93 (South/Central Asian), and within Z282 "
                                                "identify branches like M458 (Polish/Czech) or Z280 "
                                                "(Russian/Baltic). FTDNA's R1a Project has >30,000 members "
                                                "with richly annotated geographic and surname data."
                                            ),
                                            "children": [
                                                {
                                                    "haplogroup": "R1a1",
                                                    "snp_name": "SRY10831.2",
                                                    "rsids": ["rs35284194"],
                                                    "pos": 2_711_327,
                                                    "derived": "A",
                                                    "ancestral": "G",
                                                    "description": "R1a1 — near-universal R1a subclade (contains Z282 and Z93)",
                                                    "migration": (
                                                        "R1a1 (SRY10831.2/M17) contains almost all R1a men. "
                                                        "The primary split is between Z282 (European) and "
                                                        "Z93 (South/Central Asian). These are not typically "
                                                        "on consumer chips."
                                                    ),
                                                    "further": (
                                                        "Z282 (European) vs Z93 (South Asian) is the next "
                                                        "critical branching point. Both are well-studied and "
                                                        "are the first markers tested in targeted R1a panels."
                                                    ),
                                                    "children": [],
                                                },
                                            ],
                                        },
                                        # ── R1b (M343) ────────────────────────
                                        {
                                            "haplogroup": "R1b",
                                            "snp_name": "M343",
                                            "rsids": ["rs9786153", "rs2032655"],
                                            "pos": 2_887_824,
                                            "derived": "A",
                                            "ancestral": "C",
                                            "description": "Haplogroup R1b — dominant Western European lineage",
                                            "migration": (
                                                "Haplogroup R1b (M343) is the most common Y-chromosome "
                                                "lineage in Western Europe, reaching 80–95% in Ireland, "
                                                "Wales, the Basque Country, and Atlantic coastal regions. "
                                                "R1b arrived in Europe via the Yamnaya steppe expansion "
                                                "~5,000 years ago, sweeping through the continent in "
                                                "association with the Bell Beaker culture and rapidly "
                                                "displacing earlier Neolithic and Mesolithic Y-lineages. "
                                                "This expansion is archaeogenetically one of the most "
                                                "dramatic demographic replacements in prehistory: within "
                                                "~500 years, Yamnaya-related ancestry replaced ~90% of "
                                                "the previous male lineages in much of Western Europe. "
                                                "R1b is also found at significant frequencies in Central "
                                                "Asia (R1b-M73, Bashkirs ~50%), the Middle East, and "
                                                "North Africa."
                                            ),
                                            "further": (
                                                "R1b is extremely well studied. Big Y-700 + FTDNA projects "
                                                "for R1b-U106 (Germanic), R1b-P312 (Celtic/Latin), or "
                                                "R1b-L21 (Irish/Scottish/Welsh/Breton) would place you "
                                                "in a terminal subclade connecting to genetic families "
                                                "with common ancestors 500–1,500 years ago — often "
                                                "traceable to specific surnames, regions, or clan lineages."
                                            ),
                                            "children": [
                                                {
                                                    "haplogroup": "R1b1",
                                                    "snp_name": "L278",
                                                    "rsids": [],
                                                    "pos": 22_748_026,
                                                    "derived": "A",
                                                    "ancestral": "G",
                                                    "description": "R1b1 — primary R1b subclade",
                                                    "migration": (
                                                        "R1b1 contains essentially all R1b men. Its main "
                                                        "branch R1b1a (P297) separates Western European R1b "
                                                        "(M269) from Central Asian R1b (M73)."
                                                    ),
                                                    "further": "M269 is the key marker separating Western European from Central Asian R1b.",
                                                    "children": [
                                                        {
                                                            "haplogroup": "R1b1a",
                                                            "snp_name": "P297",
                                                            "rsids": [],
                                                            "pos": 22_752_468,
                                                            "derived": "T",
                                                            "ancestral": "C",
                                                            "description": "R1b1a — ancestor of M269 (European) and M73 (Central Asian) R1b",
                                                            "migration": (
                                                                "R1b1a (P297) is found in Europe, Central "
                                                                "Asia, and the Middle East. It predates the "
                                                                "split between the dominant Western European "
                                                                "M269 branch and the Central Asian M73 branch."
                                                            ),
                                                            "further": "M269 vs M73 testing distinguishes Western European from Central Asian R1b1a.",
                                                            "children": [
                                                                # ── R1b-M269 ──────────────────────
                                                                {
                                                                    "haplogroup": "R1b1a2",
                                                                    "snp_name": "M269",
                                                                    "rsids": ["rs9786153", "rs9384893"],
                                                                    "pos": 17_231_092,
                                                                    "derived": "C",
                                                                    "ancestral": "T",
                                                                    "description": "R1b-M269 — dominant Western European R1b",
                                                                    "migration": (
                                                                        "R1b-M269 is the most common Y-chromosome "
                                                                        "haplogroup in Western Europe. All modern "
                                                                        "Irish, Welsh, English, French, Spanish, "
                                                                        "Portuguese, and Italian R1b men are "
                                                                        "M269-positive. This lineage arrived in "
                                                                        "Europe during the Bell Beaker expansion "
                                                                        "~4,800–5,000 years ago, sweeping through "
                                                                        "the continent from the Atlantic coast "
                                                                        "inward. M269 frequency: Ireland/Wales ~90%, "
                                                                        "England ~70%, France ~60%, Spain ~55%, "
                                                                        "Germany ~45%, Italy ~40%. It is essentially "
                                                                        "absent east of the Carpathians, where R1a "
                                                                        "dominates."
                                                                    ),
                                                                    "further": (
                                                                        "Downstream M269 branches are geographically "
                                                                        "precise and rarely on consumer chips. Big "
                                                                        "Y-700 or dedicated R1b panels are needed:\n"
                                                                        "• U106/Z381 — Germanic (German, Dutch, "
                                                                        "English, Danish, ~25–35% of those regions)\n"
                                                                        "• P312/S116 — Western/Southern European "
                                                                        "(British Isles, France, Iberia, Italy)\n"
                                                                        "  – L21/M529 — Irish/Scottish/Welsh/Breton\n"
                                                                        "  – DF27/S250 — Iberian / S. French\n"
                                                                        "  – U152/S28 — Italian/Alpine/Gaulish\n"
                                                                        "Once in a terminal subclade, FTDNA project "
                                                                        "members with matching haplogroups often "
                                                                        "share common ancestors 500–1,500 years ago."
                                                                    ),
                                                                    "children": [
                                                                        # ── U106 ──────────────────────────
                                                                        {
                                                                            "haplogroup": "R1b-U106",
                                                                            "snp_name": "U106",
                                                                            "rsids": [],
                                                                            "pos": None,
                                                                            "derived": "T",
                                                                            "ancestral": "C",
                                                                            "description": "R1b-U106 — Germanic / North-Central European",
                                                                            "migration": (
                                                                                "R1b-U106 (Z381) is concentrated in Northern "
                                                                                "and Central Europe — Germany, Netherlands, "
                                                                                "England, Denmark, Belgium — where it reaches "
                                                                                "25–35%. Its spread is associated with Germanic "
                                                                                "tribal expansions of the 1st millennium CE. "
                                                                                "It is also found at lower levels in Scotland, "
                                                                                "Ireland, Poland, and Scandinavia. Subclades "
                                                                                "Z156, Z18, and Z301 further partition U106 "
                                                                                "by regional Germanic ancestry."
                                                                            ),
                                                                            "further": (
                                                                                "FTDNA Big Y-700 + R1b-U106 Project (one of "
                                                                                "the largest Y-DNA projects, >10,000 members) "
                                                                                "would identify your terminal subclade and "
                                                                                "connect you to surname lineage groups with "
                                                                                "common ancestors in the medieval period."
                                                                            ),
                                                                            "children": [],
                                                                        },
                                                                        # ── P312 ──────────────────────────
                                                                        {
                                                                            "haplogroup": "R1b-P312",
                                                                            "snp_name": "P312",
                                                                            "rsids": [],
                                                                            "pos": None,
                                                                            "derived": "C",
                                                                            "ancestral": "T",
                                                                            "description": "R1b-P312 — Western/Southern European (Celtic, Iberian, Latin)",
                                                                            "migration": (
                                                                                "R1b-P312 (S116) is the dominant R1b subclade "
                                                                                "in Western and Southern Europe — British Isles, "
                                                                                "France, Spain, Portugal, and Italy. Its three "
                                                                                "main branches align with major European "
                                                                                "historical populations: L21 (Celtic Northwest "
                                                                                "Europe), DF27 (Iberian/Basque), and U152 "
                                                                                "(Italian/Alpine/Gaulish). P312 is tightly "
                                                                                "associated with the Bell Beaker culture's "
                                                                                "Atlantic expansion."
                                                                            ),
                                                                            "further": (
                                                                                "FTDNA's R1b-P312 Project and Big Y-700 would "
                                                                                "identify your branch: L21 (Irish/British), "
                                                                                "DF27 (Iberian), U152 (Italian), or DF19 (rare)."
                                                                            ),
                                                                            "children": [
                                                                                {
                                                                                    "haplogroup": "R1b-L21",
                                                                                    "snp_name": "L21",
                                                                                    "rsids": [],
                                                                                    "pos": None,
                                                                                    "derived": "T",
                                                                                    "ancestral": "C",
                                                                                    "description": "R1b-L21 — Irish, Scottish, Welsh, Breton",
                                                                                    "migration": (
                                                                                        "R1b-L21 (M529/S145) is found at very high "
                                                                                        "frequency in Ireland (>80%), Scotland (~65%), "
                                                                                        "Wales (~75%), Brittany (~70%), and is common "
                                                                                        "in England (~35%) and NW France. It is "
                                                                                        "strongly associated with Celtic-speaking "
                                                                                        "peoples and their ancestors from the "
                                                                                        "Bronze Age Atlantic facade. L21 subclades "
                                                                                        "DF21 (Scottish), DF13 (Irish, widespread), "
                                                                                        "Z253 (Irish), and L513 (linked to medieval "
                                                                                        "Irish royal dynasties) allow direct "
                                                                                        "connection to named historical lineages."
                                                                                    ),
                                                                                    "further": (
                                                                                        "Big Y-700 + FTDNA R1b-L21 Project is the "
                                                                                        "gold standard for Celtic genetic genealogy. "
                                                                                        "Terminal subclades often connect to named "
                                                                                        "Irish septs, Scottish clans, Welsh families, "
                                                                                        "or Breton lineages with common ancestors "
                                                                                        "500–1,500 years ago."
                                                                                    ),
                                                                                    "children": [],
                                                                                },
                                                                                {
                                                                                    "haplogroup": "R1b-DF27",
                                                                                    "snp_name": "DF27",
                                                                                    "rsids": [],
                                                                                    "pos": None,
                                                                                    "derived": "T",
                                                                                    "ancestral": "C",
                                                                                    "description": "R1b-DF27 — Iberian Peninsula and southern France",
                                                                                    "migration": (
                                                                                        "R1b-DF27 (S250/Z209) is the dominant R1b "
                                                                                        "subclade in Spain (~70%), Portugal (~65%), "
                                                                                        "and the Basque Country (~85%), and is common "
                                                                                        "in southern France. It is strongly associated "
                                                                                        "with pre-Roman Iberian and Gaulish populations. "
                                                                                        "The Basque association is striking — despite "
                                                                                        "Basques speaking a language isolate, their "
                                                                                        "Y-chromosomes show heavy Bell Beaker/R1b input."
                                                                                    ),
                                                                                    "further": (
                                                                                        "FTDNA's Iberian R1b Project and Big Y-700 "
                                                                                        "resolve DF27 into ZZ12 (widespread) and Z195 "
                                                                                        "(Iberian-specific) with regional geographic "
                                                                                        "correlations."
                                                                                    ),
                                                                                    "children": [],
                                                                                },
                                                                                {
                                                                                    "haplogroup": "R1b-U152",
                                                                                    "snp_name": "U152",
                                                                                    "rsids": [],
                                                                                    "pos": None,
                                                                                    "derived": "T",
                                                                                    "ancestral": "C",
                                                                                    "description": "R1b-U152 — Italian, Swiss, Alpine",
                                                                                    "migration": (
                                                                                        "R1b-U152 (S28/Z36) is the dominant R1b "
                                                                                        "subclade in Italy, Switzerland, and the "
                                                                                        "Alpine regions, reaching 40–55% in northern "
                                                                                        "Italy and parts of France. It is associated "
                                                                                        "with Gaulish, Celtic, and early Roman-period "
                                                                                        "populations of the Western Alps and Po Valley."
                                                                                    ),
                                                                                    "further": (
                                                                                        "FTDNA's R1b-U152 Project and Big Y-700 "
                                                                                        "resolve to Italian regional subclades and "
                                                                                        "Alpine genetic families with exceptional "
                                                                                        "geographic precision."
                                                                                    ),
                                                                                    "children": [],
                                                                                },
                                                                            ],
                                                                        },
                                                                    ],
                                                                },
                                                                # ── R1b-M73 (Central Asian) ────────
                                                                {
                                                                    "haplogroup": "R1b1a1",
                                                                    "snp_name": "M73",
                                                                    "rsids": [],
                                                                    "pos": None,
                                                                    "derived": "T",
                                                                    "ancestral": "C",
                                                                    "description": "R1b-M73 — Central Asian R1b (Bashkir, Kazakh, Mongol)",
                                                                    "migration": (
                                                                        "R1b-M73 is the primary Central Asian R1b lineage, "
                                                                        "found at high frequencies among Bashkirs (~50%), "
                                                                        "Mongolians, and some Turkic and Caucasian populations. "
                                                                        "It represents an eastern branch of R1b1a that "
                                                                        "remained on the steppe rather than participating "
                                                                        "in the Yamnaya expansion into Western Europe."
                                                                    ),
                                                                    "further": "Big Y-700 would definitively confirm M73 vs M269 placement.",
                                                                    "children": [],
                                                                },
                                                            ],
                                                        },
                                                    ],
                                                },
                                            ],
                                        },
                                        # ── R2 (M479) ─────────────────────────
                                        {
                                            "haplogroup": "R2",
                                            "snp_name": "M479",
                                            "rsids": [],
                                            "pos": 7_763_788,
                                            "derived": "A",
                                            "ancestral": "G",
                                            "description": "Haplogroup R2 — South Asian R lineage",
                                            "migration": (
                                                "Haplogroup R2 (M479) is found primarily in South "
                                                "Asia — India, Pakistan, Afghanistan — and at low "
                                                "frequency in the Middle East and Central Asia. It "
                                                "is the second most common R subclade in India after "
                                                "R1a, reaching ~10–15% in Pakistan and NW India, "
                                                "and is associated with ancient South Asian lineages "
                                                "predating the Indo-Aryan migrations."
                                            ),
                                            "further": (
                                                "Big Y-700 + FTDNA's South Asian projects would "
                                                "identify R2 subclades with population specificity "
                                                "within South Asian groups."
                                            ),
                                            "children": [],
                                        },
                                    ],
                                },
                            ],
                        },
                    ],
                },
            ],
        },
    ],
}


# ── Lookup helpers ─────────────────────────────────────────────────────────────

def _build_lookup(snps_df: pd.DataFrame) -> Tuple[Dict[str, str], Dict[int, str]]:
    """
    Build two lookup tables from the parsed DataFrame:
      rsid_map  : rsid (str) → genotype (str, upper-case)
      pos_map   : chr-Y integer position → genotype
    Handles both hemizygous single-char ('A') and doubled ('AA') genotypes.
    Filters to chromosome Y rows only.
    """
    rsid_map: Dict[str, str] = {}
    pos_map: Dict[int, str] = {}

    # normalise chromosome label to Y across formats
    def _is_y(chrom_val) -> bool:
        return str(chrom_val).strip().upper() in ("Y", "24")

    chrom_col = "chrom" if "chrom" in snps_df.columns else None

    for rsid_raw, row in snps_df.iterrows():
        if chrom_col and not _is_y(row[chrom_col]):
            continue

        gt = str(row.get("genotype", "")).strip().upper()
        if not gt or gt in ("NAN", "--", "0", ""):
            continue

        # collapse 'AA' → 'A', keep 'AC' as-is (shouldn't occur on Y but be safe)
        if len(gt) == 2 and gt[0] == gt[1]:
            gt = gt[0]

        rsid_str = str(rsid_raw).strip()
        rsid_map[rsid_str] = gt

        pos = row.get("pos", None)
        if pos is not None:
            try:
                pos_map[int(pos)] = gt
            except (ValueError, TypeError):
                pass

    return rsid_map, pos_map


def _lookup(node: Dict, rsid_map: Dict, pos_map: Dict) -> Tuple[str, str]:
    """
    Determine whether this node's SNP is DERIVED / ANCESTRAL / NOT_FOUND.
    Checks rsIDs first, then position, then strand-complement of each.

    Returns (status, genotype_found).
    """
    derived = node.get("derived", "").upper()
    ancestral = node.get("ancestral", "").upper()
    comp_derived = complement(derived) if derived else ""
    comp_ancestral = complement(ancestral) if ancestral else ""

    # By rsID
    for rsid in node.get("rsids", []):
        gt = rsid_map.get(rsid, "")
        if not gt:
            continue
        g = gt[0] if len(gt) >= 1 else ""
        if g == derived or g == comp_derived:
            return "derived", gt
        if g == ancestral or g == comp_ancestral:
            return "ancestral", gt

    # By position
    pos = node.get("pos")
    if pos:
        gt = pos_map.get(int(pos), "")
        if gt:
            g = gt[0] if gt else ""
            if g == derived or g == comp_derived:
                return "derived", gt
            if g == ancestral or g == comp_ancestral:
                return "ancestral", gt

    return "not_found", ""


# ── Tree walker ────────────────────────────────────────────────────────────────

def _walk(
    node: Dict,
    rsid_map: Dict,
    pos_map: Dict,
    path: List[Dict],
    depth: int,
) -> Tuple[List[Dict], str, List[str]]:
    """
    Recursive tree walk. Returns (path, final_status, gap_snp_names).

    Strategy for NOT_FOUND (chip gap):
      - Still descend into children so that a confirmed downstream marker
        (e.g. R-M207 present even though K-M9 is absent) can be reported.
      - All gap ancestors are marked snp_status='chip_gap' / 'inferred'.
    """
    if depth > 25:
        return path, "max_depth", []

    status, gt = _lookup(node, rsid_map, pos_map)

    if status == "ancestral":
        return path, "ruled_out", []   # this branch definitively excluded

    node_entry = {k: v for k, v in node.items() if k != "children"}
    node_entry["found_genotype"] = gt if gt else "–"
    node_entry["snp_status"] = "confirmed" if status == "derived" else "chip_gap"

    new_path = path + [node_entry]
    children = node.get("children", [])

    if not children:
        final_st = "resolved" if status == "derived" else "chip_gap"
        return new_path, final_st, ([] if status == "derived" else [node["snp_name"]])

    # Check children to find which branch (if any) is confirmed or gap-traversable
    confirmed_child = None
    gap_children: List[str] = []
    ruled_out: List[str] = []

    for child in children:
        cs, cgt = _lookup(child, rsid_map, pos_map)
        if cs == "derived":
            confirmed_child = child
            break
        elif cs == "not_found":
            gap_children.append(child["snp_name"])
        else:
            ruled_out.append(child["snp_name"])

    if confirmed_child:
        return _walk(confirmed_child, rsid_map, pos_map, new_path, depth + 1)

    if gap_children:
        unruled = [c for c in children if c["snp_name"] not in ruled_out]
        if len(unruled) == 1:
            return _walk(unruled[0], rsid_map, pos_map, new_path, depth + 1)
        # Multiple gap children: lookahead to find branch with most confirmed descendants
        best_result = None
        best_confirmed = 0
        for child in unruled:
            r_path, r_status, r_gaps = _walk(child, rsid_map, pos_map, new_path, depth + 1)
            confirmed_in = sum(1 for n in r_path if n.get("snp_status") == "confirmed")
            if confirmed_in > best_confirmed:
                best_confirmed = confirmed_in
                best_result = (r_path, r_status, r_gaps)
        if best_result and best_confirmed > 0:
            return best_result
        return new_path, "chip_gap", gap_children

    # All children were ancestral — current node is the terminal haplogroup
    return new_path, "resolved", []


# ── Public API ─────────────────────────────────────────────────────────────────

def analyze_y_haplogroup(snps_df: pd.DataFrame) -> Dict:
    """
    Analyse Y-chromosome SNPs and return a result dict with keys:
      status              – 'resolved' | 'partial' | 'no_y_data' | 'not_k'
      y_snp_count         – number of Y-SNPs found in the file
      haplogroup_path     – human-readable path string, e.g. "K > K2b > P > R > R1b"
      terminal_haplogroup – name of deepest confirmed/inferred haplogroup
      path                – list of node dicts (each has snp_status: confirmed/chip_gap)
      chip_gaps           – SNP names of markers not on chip that block further resolution
      further_testing     – text describing what FTDNA Big Y would add
      message             – one-line summary for the terminal log
    """
    rsid_map, pos_map = _build_lookup(snps_df)
    y_count = len(rsid_map)

    if y_count == 0:
        return {
            "status": "no_y_data",
            "confidence": "none",
            "y_snp_count": 0,
            "haplogroup_path": "Unknown",
            "terminal_haplogroup": "Unknown",
            "terminal_description": "",
            "terminal_migration": None,
            "path": [],
            "chip_gaps": [],
            "not_tested_branches": [],
            "further_testing": "",
            "message": (
                "No Y-chromosome SNPs found in this file. Possible reasons: "
                "(1) the person is female and has no Y chromosome, "
                "(2) the chip or lab does not genotype Y-chromosome markers, or "
                "(3) Y-SNPs are stored under a chromosome label not yet recognised."
            ),
        }

    if y_count < 10:
        return {
            "status": "insufficient_y_data",
            "confidence": "none",
            "y_snp_count": y_count,
            "haplogroup_path": "Insufficient Y-chromosome SNPs",
            "terminal_haplogroup": "Insufficient Y-chromosome SNPs",
            "terminal_description": "",
            "terminal_migration": None,
            "path": [],
            "chip_gaps": [],
            "not_tested_branches": [],
            "further_testing": "",
            "message": (
                f"Insufficient Y-chromosome SNPs ({y_count} found, minimum 10 required) "
                "to make a reliable haplogroup call."
            ),
        }

    walked_path, walk_status, chip_gaps = _walk(
        HAPLOGROUP_TREE, rsid_map, pos_map, [], 0
    )

    if not walked_path:
        if walk_status == "ruled_out":
            return {
                "status": "not_k",
                "confidence": "low",
                "y_snp_count": y_count,
                "haplogroup_path": "Pre-K (A, B, C, D, E, F, G, H, I, or J)",
                "terminal_haplogroup": "Pre-K",
                "terminal_description": "Y-chromosome branches that diverged before Haplogroup K arose",
                "terminal_migration": (
                    "Your Y-chromosome lineage appears to belong to a non-K haplogroup — "
                    "one of A, B, C, D, E, F, G, H, I, or J. These haplogroups represent "
                    "deep branches of the human Y-chromosome tree that diverged from the K "
                    "ancestor more than 45,000 years ago. Haplogroup E is the most common "
                    "non-K lineage outside Africa, found widely in Africa, and at lower "
                    "levels in the Middle East and Southern Europe. Haplogroups G, I, J "
                    "are common in the Middle East, Caucasus, and parts of Europe."
                ),
                "path": [],
                "chip_gaps": [],
                "not_tested_branches": [],
                "further_testing": (
                    "Contact FTDNA or a specialist in human Y-chromosome phylogenetics "
                    "to confirm your haplogroup assignment and identify your specific subclade."
                ),
                "message": "Y-chromosome does not appear to be in Haplogroup K.",
            }
        # chip gap at root or other issue
        return {
            "status": "partial",
            "confidence": "none",
            "y_snp_count": y_count,
            "haplogroup_path": "Unresolved (key markers not on chip)",
            "terminal_haplogroup": "Unknown",
            "terminal_description": "",
            "terminal_migration": None,
            "path": [],
            "chip_gaps": chip_gaps or [HAPLOGROUP_TREE["snp_name"]],
            "not_tested_branches": [],
            "further_testing": HAPLOGROUP_TREE.get("further", ""),
            "message": "Could not resolve haplogroup — defining markers not on chip.",
        }

    terminal = walked_path[-1]

    # Collect untested branches downstream of the terminal node
    not_tested: List[Dict] = []
    for child in terminal.get("children", []):
        cs, _ = _lookup(child, rsid_map, pos_map)
        if cs == "not_found":
            not_tested.append({
                "haplogroup": child["haplogroup"],
                "snp_name": child["snp_name"],
                "description": child.get("description", ""),
                "further": child.get("further", ""),
            })

    resolved = walk_status == "resolved"
    confirmed_count = sum(1 for n in walked_path if n["snp_status"] == "confirmed")
    n_path = len(walked_path)

    # Confidence in the terminal call is driven by how many path markers were
    # actually confirmed (derived) vs. inferred across chip gaps.
    if resolved and confirmed_count >= 3:
        confidence = "high"
    elif confirmed_count >= 2:
        confidence = "moderate"
    else:
        confidence = "low"
    confidence_note = (
        f"{confirmed_count} of {n_path} markers on the assigned path confirmed "
        f"(derived); the rest inferred across chip gaps. "
        f"{y_count:,} Y-SNPs typed total. Terminal-branch resolution from chip "
        "data is limited — targeted Y-SNP or Y-sequencing refines it."
    )

    return {
        "status": "resolved" if resolved else "partial",
        "confidence": confidence,
        "confidence_note": confidence_note,
        "y_snp_count": y_count,
        "n_markers_confirmed": confirmed_count,
        "n_markers_on_path": n_path,
        "haplogroup_path": " > ".join(n["haplogroup"] for n in walked_path),
        "terminal_haplogroup": terminal["haplogroup"],
        "terminal_description": terminal.get("description", ""),
        "terminal_migration": terminal.get("migration", ""),
        "path": walked_path,
        "chip_gaps": chip_gaps,
        "not_tested_branches": not_tested,
        "further_testing": terminal.get("further", ""),
        "message": (
            f"Y-DNA: {terminal['haplogroup']} "
            f"({'resolved' if resolved else 'partial — chip gaps'}, "
            f"{confirmed_count} confirmed marker{'s' if confirmed_count != 1 else ''})"
        ),
    }
