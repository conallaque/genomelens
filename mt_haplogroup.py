"""
mtDNA Haplogroup Analysis Module
--------------------------------

Maternal-lineage classification from autosomal-chip mtDNA SNPs.

This is *approximate*. Consumer autosomal arrays type only a handful of
mitochondrial positions, which is sufficient to place most samples into one
of the macro-haplogroups (L, M, N, R, HV, H, V, J, T, U, K, I, W, X) but
rarely sufficient to resolve modern fine-grained subclades (e.g. H1a vs
H1b vs H1c). For high-resolution subclade calling, FTDNA's mtFull (full
mitochondrial sequence) is the appropriate test.

The module:
  1. Filters chrMT SNPs from the parsed DNA file.
  2. Looks up each defining marker by rsID first, then by mtDNA position.
  3. Walks a decision tree of mtDNA haplogroup-defining variants based on
     rCRS coordinates.
  4. Returns the best-guess haplogroup with confidence, matched markers,
     a migration narrative, and ancient-DNA comparisons.
"""


import pandas as pd

# ─── Defining markers ─────────────────────────────────────────────────────────
#
# Each marker is identified by:
#   name      — conventional notation (e.g. "T7028C")
#   pos       — mtDNA position in rCRS coordinates
#   rsids     — list of rsIDs sometimes used for this position (chip-dependent)
#   derived   — the derived allele (the variant that signals downstream lineage)
#   ancestral — the rCRS allele
#   defines   — the haplogroup the derived state defines
#   level     — broad-vs-fine: "macro" (L/M/N/R), "major" (H/V/J/T/U/K/I/W/X),
#               or "subclade" (e.g. H1, J1c)
#   description — short description shown in the report

MTDNA_MARKERS: list[dict] = [
    # ── Macro-haplogroups ─────────────────────────────────────────────────────
    {
        # Recorded from the non-R side until 2026-08-23 ("derived" T meaning
        # "retains the ancestral state, therefore not R"). That framing is
        # defensible in isolation but put this row in literal contradiction
        # with the 12705 row, which discriminates the same branch point. Both
        # now describe R in the same direction, so the marker table and the
        # lineage chain cannot disagree.
        "name": "T16223C", "pos": 16223, "rsids": ["rs41323649"],
        "derived": "C", "ancestral": "T",
        "defines": "R (and descendants)", "level": "macro",
        "description": "Position 16223 ancestral T is associated with non-R lineages "
                       "(L, M, N); derived C is the R-and-descendant state.",
    },
    {
        "name": "C7028T", "pos": 7028, "rsids": ["rs2854122", "rs3937033"],
        "derived": "C", "ancestral": "T",
        "defines": "H", "level": "major",
        "description": "Defining marker for haplogroup H — the most common European "
                       "haplogroup (~40% frequency). H descends from HV.",
    },
    {
        "name": "T14766C", "pos": 14766, "rsids": ["rs193302980"],
        "derived": "C", "ancestral": "T",
        "defines": "HV (and H/V descendants)", "level": "major",
        "description": "Defines haplogroup HV — the parent of H and V.",
    },
    {
        "name": "C4580T", "pos": 4580, "rsids": ["rs2032658"],
        "derived": "T", "ancestral": "C",
        "defines": "V", "level": "major",
        "description": "Defines haplogroup V — found at high frequency in Saami "
                       "and Iberian populations.",
    },
    {
        "name": "A4216G", "pos": 4216, "rsids": ["rs3088309"],
        "derived": "G", "ancestral": "A",
        "defines": "JT (and J/T descendants)", "level": "major",
        "description": "Defining marker for JT clade (parent of J and T).",
    },
    {
        "name": "A10398G", "pos": 10398, "rsids": ["rs2853826"],
        "derived": "G", "ancestral": "A",
        "defines": "J, I, N1a", "level": "subclade",
        "description": "Shared marker between J, I, and N1a — needs context "
                       "(e.g. 4216 status) to disambiguate.",
    },
    {
        "name": "A11251G", "pos": 11251, "rsids": ["rs28359178"],
        "derived": "G", "ancestral": "A",
        "defines": "T (within JT)", "level": "major",
        "description": "Defines haplogroup T (when combined with JT marker 4216).",
    },
    {
        # INVERTED UNTIL 2026-08-23. This read derived="T", which claims
        # haplogroup R carries T at 12705. The rCRS is itself a haplogroup H
        # sequence — H sits inside R — and rCRS has C here, so R carries C and
        # the entry contradicted both reality and the 16223 marker four rows
        # up, which correctly treats T as the non-R ancestral state.
        "name": "C12705T", "pos": 12705, "rsids": ["rs28358580"],
        "derived": "C", "ancestral": "T",
        "defines": "R (and descendants)", "level": "macro",
        "description": "Defines macro-haplogroup R — ancestor of HV, JT, U, B, F.",
    },
    {
        "name": "A12308G", "pos": 12308, "rsids": ["rs28359175"],
        "derived": "G", "ancestral": "A",
        "defines": "U (and K descendants)", "level": "major",
        "description": "Defining marker for haplogroup U — diverse Eurasian "
                       "lineage with deep subclades (U2–U8, K).",
    },
    {
        "name": "T9055A", "pos": 9055, "rsids": ["rs28358571"],
        "derived": "A", "ancestral": "T",
        "defines": "K (within U8)", "level": "subclade",
        "description": "Defines haplogroup K — common in Ashkenazi Jewish "
                       "and broader European populations.",
    },
    {
        "name": "T4529A", "pos": 4529, "rsids": ["rs3928306"],
        "derived": "A", "ancestral": "T",
        "defines": "I (within N1a-c)", "level": "subclade",
        "description": "Defines haplogroup I — northern European lineage.",
    },
    {
        "name": "G8994A", "pos": 8994, "rsids": ["rs41544112"],
        "derived": "A", "ancestral": "G",
        "defines": "W (within N2)", "level": "subclade",
        "description": "Defines haplogroup W — Central/Eastern European and "
                       "South Asian frequencies.",
    },
    {
        "name": "T6221C", "pos": 6221, "rsids": [],
        "derived": "C", "ancestral": "T",
        "defines": "X (within N)", "level": "subclade",
        "description": "Defines haplogroup X — uniquely distributed in Europe, "
                       "the Caucasus, and Native North American populations.",
    },
    # Position 489 — defines M
    {
        "name": "T489C", "pos": 489, "rsids": [],
        "derived": "C", "ancestral": "T",
        "defines": "M (and descendants)", "level": "macro",
        "description": "Defines macro-haplogroup M — Asian, Australian, "
                       "Oceanian, and Native American lineages.",
    },
]


# ─── Haplogroup metadata: migration & ancient DNA ─────────────────────────────
HAPLOGROUP_INFO: dict[str, dict[str, str]] = {
    "H": {
        "migration": (
            "Haplogroup H is the most common mtDNA lineage in Europe (~40% of "
            "Europeans, peaking at ~50% in some Western European populations). "
            "H originated in Southwest Asia approximately 25,000–30,000 years ago "
            "and spread into Europe during the Last Glacial Maximum and the "
            "subsequent post-glacial recolonization from Iberian and Italian "
            "refugia. The Neolithic farmer expansion from Anatolia ~8,000 years "
            "ago, and the Bronze Age steppe migrations ~5,000 years ago, both "
            "carried H subclades throughout Europe. Today H is also found at "
            "moderate frequency in North Africa, the Middle East, and Central Asia."
        ),
        "ancient_dna": (
            "Haplogroup H is widely represented in ancient European DNA: H is "
            "common in Linear Pottery Neolithic farmers (LBK, ~7,500 years ago), "
            "Bell Beaker culture individuals (~4,500 years ago), and Bronze Age "
            "samples from across Europe. Cheddar Man (~10,000-year-old Mesolithic "
            "Briton) carried haplogroup U5b1 — but the post-Neolithic European "
            "population is dominated by H lineages."
        ),
        "further_testing": (
            "FTDNA mtFull (full mitochondrial sequence) would resolve your "
            "H subclade (H1, H2, H3, … H11, etc.), of which there are dozens "
            "with different geographic distributions and historical contexts."
        ),
    },
    "V": {
        "migration": (
            "Haplogroup V is a Western European mtDNA lineage that arose from HV "
            "approximately 15,000 years ago, likely in the Franco-Cantabrian refugium "
            "during the Last Glacial Maximum. V reaches its highest frequencies in "
            "the Saami of northern Scandinavia (>50%), Basques (~10%), and certain "
            "Iberian populations. The post-glacial expansion northward carried V "
            "with the recolonization of Europe."
        ),
        "ancient_dna": (
            "V has been identified in Mesolithic and Neolithic European samples, "
            "particularly in northern populations. The Saami's distinctive V "
            "frequency reflects a founder effect and partial isolation."
        ),
        "further_testing": "FTDNA mtFull would resolve V subclades (V1, V2, V3, V7, etc.).",
    },
    "J": {
        "migration": (
            "Haplogroup J emerged in the Near East approximately 45,000 years ago "
            "and spread into Europe during the Neolithic transition. J is found at "
            "highest frequencies (~12%) in the Near East and Caucasus, with "
            "European frequencies of ~8–10%. The J1c subclade is strongly associated "
            "with the Neolithic farmer expansion from Anatolia ~8,000 years ago, "
            "and is found in many ancient farming sites across Europe."
        ),
        "ancient_dna": (
            "J haplogroups — particularly J1c — are common in Linear Pottery and "
            "Funnel Beaker Neolithic farmers. Ötzi the Iceman (~5,300 years ago, "
            "Italian Alps) carried haplogroup K1 (related to U8b'K), not J. The "
            "Neolithic farmer wave is well represented by J1c in central European "
            "Linear Pottery (LBK) sites."
        ),
        "further_testing": "FTDNA mtFull would resolve J subclades (J1a/b/c/d, J2, etc.).",
    },
    "T": {
        "migration": (
            "Haplogroup T also arose in the Near East ~30,000 years ago and "
            "expanded into Europe in two waves: one during the Upper Paleolithic "
            "and a major one during the Neolithic farming expansion ~8,000 years "
            "ago. T2 is the most common European T subclade. T is also notable "
            "for being the maternal lineage of Tsar Nicholas II of Russia and "
            "the rest of the Romanov family — used to identify their remains."
        ),
        "ancient_dna": (
            "T haplogroups are widely represented in European Neolithic farmer "
            "remains, including Linear Pottery and Funnel Beaker sites. The "
            "Romanov family carried haplogroup T (~T1)."
        ),
        "further_testing": "FTDNA mtFull would resolve T1 vs T2 vs T3 subclades.",
    },
    "U": {
        "migration": (
            "Haplogroup U is one of the oldest European mtDNA lineages, with "
            "U5 specifically being the dominant haplogroup among pre-Neolithic "
            "European hunter-gatherers (~50% of Mesolithic European mtDNA). "
            "U5 expanded across Europe during the late Paleolithic and Mesolithic "
            "(~30,000–10,000 years ago). U was substantially replaced in central "
            "and southern Europe by H during the Neolithic and Bronze Age but "
            "remains the dominant haplogroup in many northern European groups "
            "(Saami, Estonians, some northern Russians)."
        ),
        "ancient_dna": (
            "U5 is the haplogroup of Cheddar Man (~10,000-year-old Mesolithic "
            "Briton, U5b1) and many other Mesolithic European hunter-gatherers — "
            "including individuals from La Braña (Spain), Loschbour (Luxembourg), "
            "and Karelia (Russia). Pre-Neolithic European samples are dominated "
            "by U5 lineages."
        ),
        "further_testing": "FTDNA mtFull would resolve U2 vs U3 vs U4 vs U5 vs U6 vs U7 vs U8.",
    },
    "K": {
        "migration": (
            "Haplogroup K is a descendant of U8, originating in the Near East "
            "~30,000 years ago. K spread into Europe during the Neolithic and "
            "is found across Europe, the Near East, and North Africa. K reaches "
            "particularly high frequencies in Ashkenazi Jewish populations "
            "(~30%), reflecting both Near Eastern ancestry and significant "
            "founder effects during the medieval Ashkenazi diaspora."
        ),
        "ancient_dna": (
            "Ötzi the Iceman (~5,300 years ago, Alpine Copper Age) carried "
            "haplogroup K1f. K1a is common in Neolithic Linear Pottery (LBK) "
            "samples — particularly in Hungary and central Europe. K1a1b1a is "
            "one of four founder lineages of the modern Ashkenazi Jewish "
            "population."
        ),
        "further_testing": "FTDNA mtFull would resolve K1 vs K2 and many K subclades.",
    },
    "I": {
        "migration": (
            "Haplogroup I is a northern European and Caucasian mtDNA lineage "
            "within N1a-c, arising approximately 30,000 years ago. I is found "
            "at low frequencies (~2–3%) across Europe with somewhat higher "
            "frequencies in northern Europe and the Caucasus."
        ),
        "ancient_dna": (
            "I is found in some European Neolithic and Bronze Age samples but "
            "is not strongly associated with any single ancient migration."
        ),
        "further_testing": "FTDNA mtFull would resolve I1, I2, I3, I4, I5 subclades.",
    },
    "W": {
        "migration": (
            "Haplogroup W is found at low frequencies in Eastern Europe, the "
            "Caucasus, and South Asia. It arose within N2 approximately 25,000 "
            "years ago and likely spread with multiple Holocene migrations into "
            "Europe, including possible association with the Indo-European "
            "expansions."
        ),
        "ancient_dna": (
            "W has been identified in several Bronze Age and Iron Age European "
            "samples, with associations to steppe-origin populations."
        ),
        "further_testing": "FTDNA mtFull would resolve W1, W3, W6 subclades.",
    },
    "X": {
        "migration": (
            "Haplogroup X has an unusual distribution — it is found at low "
            "frequencies in Europe, the Caucasus, the Near East, and (notably) "
            "in Native North American populations (~3–5% in some Algonquin and "
            "Sioux groups). The Native American X2a subclade represents a "
            "founder lineage in the peopling of the Americas, complementing the "
            "more common A, B, C, D lineages. X arose within N approximately "
            "25,000 years ago."
        ),
        "ancient_dna": (
            "X has been identified in pre-Columbian Native American remains and "
            "in Neolithic European samples. The Kennewick Man controversy and "
            "later definitive Native American genetic studies confirmed that X2a "
            "in the Americas is part of the standard founder set, not evidence "
            "of an unusual migration."
        ),
        "further_testing": "FTDNA mtFull would resolve X1, X2 (and the Native American X2a).",
    },
    "HV": {
        "migration": (
            "Haplogroup HV is the parent of H and V. HV arose in the Near East / "
            "Caucasus region ~30,000 years ago. Outside Europe, HV* (without "
            "downstream H or V derivation) is found at low frequencies in the "
            "Near East and North Africa. Without further mtDNA testing, an HV* "
            "call may simply mean the test could not yet confirm derivation into "
            "H or V."
        ),
        "ancient_dna": (
            "HV* is rare in ancient remains; most ancient European samples are "
            "either H, V, or other macrohaplogroups."
        ),
        "further_testing": "Full mtDNA sequencing recommended to resolve H vs V vs HV*.",
    },
    "R": {
        "migration": (
            "Macro-haplogroup R is one of the major Eurasian mtDNA lineages, "
            "ancestor of H, V, HV, J, T, U, K, B, F, and many South Asian "
            "lineages. R arose in Asia ~55,000 years ago and seeded most of "
            "the maternal-line diversity of Eurasia and the Americas."
        ),
        "ancient_dna": (
            "R is the ancestor of the majority of Eurasian mtDNA in ancient remains."
        ),
        "further_testing": "Full mtDNA sequencing recommended to resolve specific subclade.",
    },
    "N": {
        "migration": (
            "Macro-haplogroup N is one of the two daughter clades of L3 that "
            "left Africa ~70,000 years ago — the other being M. N is the "
            "ancestor of most Western Eurasian mtDNA lineages (H, V, J, T, U, "
            "K, I, W, X). It originated in or near Arabia/the Near East."
        ),
        "ancient_dna": (
            "N is the macro-ancestor of most ancient European and Near Eastern "
            "samples."
        ),
        "further_testing": "Full mtDNA sequencing recommended to resolve specific subclade.",
    },
    "M": {
        "migration": (
            "Macro-haplogroup M arose in Asia from L3 shortly after the Out-of-"
            "Africa dispersal. M is the dominant maternal-line ancestor of "
            "South Asian, East Asian, Native American (haplogroups A, C, D), "
            "and Australian Aboriginal populations. M is essentially absent "
            "from autochthonous European populations."
        ),
        "ancient_dna": (
            "M is widely represented in ancient Asian and Native American DNA. "
            "M subclades C and D are major founder lineages of the peopling of "
            "the Americas."
        ),
        "further_testing": "Full mtDNA sequencing recommended to resolve subclade.",
    },
    "L": {
        "migration": (
            "Macro-haplogroup L is the root of the human mtDNA tree, found "
            "almost exclusively in sub-Saharan African populations. The "
            "deepest L lineages (L0, L1, L2) are found at highest frequency "
            "in Khoisan and Pygmy populations. L3 is the ancestor of all "
            "non-African mtDNA lineages."
        ),
        "ancient_dna": (
            "L is the ancestral state of all human mtDNA. Early human remains "
            "from Africa and the deepest non-African lineages all trace to L."
        ),
        "further_testing": "Full mtDNA sequencing recommended to resolve L subclade.",
    },
    "Unknown": {
        "migration": (
            "Insufficient mtDNA markers were available on this chip to confidently "
            "assign a haplogroup. Maternal-line testing via a dedicated mtDNA "
            "test (FTDNA mtFull) is needed for reliable classification."
        ),
        "ancient_dna": "",
        "further_testing": (
            "FTDNA mtFull provides full mitochondrial sequencing and high-"
            "resolution haplogroup placement. 23andMe and AncestryDNA report "
            "mtDNA haplogroups using curated chip-marker subsets if you also "
            "want to compare results."
        ),
    },
}


# ─── Analysis function ────────────────────────────────────────────────────────

def _normalise_genotype(gt: object) -> str | None:
    if gt is None:
        return None
    s = str(gt).upper().strip().replace(" ", "").replace("-", "")
    if s in ("", "NAN", "--"):
        return None
    return s


def _find_marker(snps_df: pd.DataFrame, marker: dict) -> dict | None:
    """Look up a marker by rsID list, then by chrMT+position. Returns dict
    with {found_genotype, rsid_hit, derived_status} or None if not found."""
    # rsID lookup
    for rsid in marker.get("rsids", []):
        if rsid in snps_df.index:
            row = snps_df.loc[rsid]
            gt = _normalise_genotype(row.get("genotype"))
            if gt:
                return {
                    "rsid": rsid,
                    "genotype": gt,
                    "is_derived": marker["derived"] in gt,
                }
    # Position lookup
    pos = marker.get("pos")
    if pos:
        mt_df = snps_df[snps_df["chrom"].isin(["MT", "M", "chrMT", "chrM"])]
        hits = mt_df[mt_df["pos"] == pos]
        if not hits.empty:
            row = hits.iloc[0]
            gt = _normalise_genotype(row.get("genotype"))
            if gt:
                return {
                    "rsid": str(row.name),
                    "genotype": gt,
                    "is_derived": marker["derived"] in gt,
                }
    return None


# ─── Phylogenetic tree ────────────────────────────────────────────────────────
#
# WHY THIS EXISTS. The flat marker list above answers "which haplogroup?" by
# scanning for derived markers and picking the most specific one off a fixed
# priority list. That returns a single label — "H" — and throws away the thing
# a maternal lineage actually is: a chain of branch points, each with its own
# marker, its own evidence, and its own date. The Y-DNA module has always
# reported that chain; mtDNA reported only the endpoint.
#
# The tree below is the same phylogeny the flat list encodes, written as the
# nested structure it always implicitly was, so it can be walked into a path:
#
#     mt-MRCA › N › R › HV › H
#
# Topology follows PhyloTree mtDNA build 17. Note that I, W and X hang off N
# directly, NOT off R — a detail the priority-list classifier could not
# express, and got wrong whenever someone carried both an R marker and an
# X marker.
#
# A NOTE ON READING rCRS COORDINATES. The revised Cambridge Reference Sequence
# is itself a haplogroup H2a2a1 sequence. So for any clade H belongs to
# (R, HV, H) the DERIVED allele is the rCRS allele, while for clades H does not
# belong to the derived allele differs from rCRS. Getting this backwards
# inverts a marker, and an inverted marker at a backbone node misroutes every
# sample that passes through it. ``test_mt_haplogroup.py`` pins the rule.


def _mt_node(haplogroup: str, markers: list[tuple], description: str = "",
             migration: str = "", further: str = "",
             children: list[dict] | None = None) -> dict:
    """One branch point. ``markers`` are (name, rsids, pos, derived, ancestral)."""
    return {
        "haplogroup": haplogroup,
        "snp_name": markers[0][0] if markers else haplogroup,
        "rsids": markers[0][1] if markers else [],
        "pos": markers[0][2] if markers else None,
        "markers": markers,
        "description": description,
        "migration": migration,
        "further": further,
        "children": children or [],
    }


MT_TREE: dict = _mt_node(
    "mt-MRCA", [],
    description="Mitochondrial Eve — the most recent common maternal ancestor "
                "of every living person.",
    migration="All maternal lineages alive today descend from one woman in "
              "Africa roughly 150,000-200,000 years ago. She was not the only "
              "woman alive; she is simply the one whose maternal line never "
              "died out.",
    children=[
        _mt_node(
            "M", [("T489C", [], 489, "C", "T"),
                  ("A10400G", ["rs28358279"], 10400, "G", "A")],
            description="Macro-haplogroup M — the great Asian maternal branch.",
            migration="M left Africa with the southern coastal migration around "
                      "60,000 years ago and became the dominant maternal lineage "
                      "of South, East and Southeast Asia, and of Native America "
                      "via C and D.",
            further="Full mtDNA sequencing resolves M into C, D, E, G, Q and Z."),
        _mt_node(
            "N", [],
            description="Macro-haplogroup N — the branch behind almost every "
                        "European and West Eurasian maternal line.",
            migration="N is the other great out-of-Africa maternal branch, "
                      "arising around 60,000 years ago. Nearly all indigenous "
                      "European maternal lineages sit inside it.",
            further="No single N-defining marker sits on consumer arrays, so N "
                    "is normally inferred from a confirmed descendant rather "
                    "than observed directly.",
            children=[
                _mt_node(
                    "I", [("T4529A", ["rs3928306"], 4529, "A", "T")],
                    description="Haplogroup I — a small, old northern European "
                                "lineage.",
                    migration="I is found at low frequency across northern and "
                              "eastern Europe and is associated with early "
                              "post-glacial resettlement."),
                _mt_node(
                    "W", [("G8994A", ["rs41544112"], 8994, "A", "G")],
                    description="Haplogroup W — eastern European and South Asian.",
                    migration="W spread from the Near East into eastern Europe "
                              "and South Asia after the last glacial maximum."),
                _mt_node(
                    "X", [("T6221C", [], 6221, "C", "T")],
                    description="Haplogroup X — unusually scattered.",
                    migration="X is one of the few maternal lineages found on "
                              "both sides of the Atlantic, present in Europe, "
                              "the Near East and among some Native American "
                              "populations — a distribution long argued over."),
                _mt_node(
                    "R", [("C12705T", ["rs28358580"], 12705, "C", "T"),
                          ("T16223C", ["rs41323649"], 16223, "C", "T")],
                    description="Haplogroup R — the branch of N containing most "
                                "European maternal lineages.",
                    migration="R arose within N around 55,000 years ago and gave "
                              "rise to H, V, J, T, U and K, which together cover "
                              "the great majority of European maternal lines.",
                    children=[
                        _mt_node(
                            "HV", [("T14766C", ["rs193302980"], 14766, "C", "T")],
                            description="Haplogroup HV — parent of H and V.",
                            migration="HV emerged in the Near East and moved into "
                                      "Europe before the last glacial maximum.",
                            children=[
                                _mt_node(
                                    "H", [("C7028T", ["rs2854122", "rs3937033"],
                                           7028, "C", "T")],
                                    description="Haplogroup H — the most common "
                                                "European maternal lineage, "
                                                "around 40% of Europeans.",
                                    migration="H expanded out of the Franco-"
                                              "Cantabrian refuge as the ice "
                                              "retreated roughly 15,000 years "
                                              "ago and repopulated Europe.",
                                    further="H is enormous and highly "
                                            "substructured. Full mtDNA "
                                            "sequencing is the only way to "
                                            "resolve H1 from H3 from H5."),
                                _mt_node(
                                    "V", [("C4580T", ["rs2032658"], 4580, "T", "C")],
                                    description="Haplogroup V — western European, "
                                                "notably Saami and Basque.",
                                    migration="V shares H's post-glacial "
                                              "expansion but stayed far rarer, "
                                              "peaking in northern Scandinavia "
                                              "and the Basque country."),
                            ]),
                        _mt_node(
                            "JT", [("A4216G", ["rs3088309"], 4216, "G", "A")],
                            description="Haplogroup JT — parent of J and T.",
                            migration="JT arose in the Near East and entered "
                                      "Europe with the spread of farming.",
                            children=[
                                _mt_node(
                                    "J", [("A10398G", ["rs2853826"], 10398, "G", "A")],
                                    description="Haplogroup J — Near Eastern origin, "
                                                "widespread in Europe.",
                                    migration="J is strongly associated with the "
                                              "Neolithic expansion of agriculture "
                                              "out of the Fertile Crescent."),
                                _mt_node(
                                    "T", [("A11251G", ["rs28359178"], 11251, "G", "A")],
                                    description="Haplogroup T — Near Eastern origin.",
                                    migration="T accompanied J into Europe with "
                                              "early farming communities."),
                            ]),
                        _mt_node(
                            "U", [("A12308G", ["rs28359175"], 12308, "G", "A")],
                            description="Haplogroup U — one of the oldest European "
                                        "maternal lineages.",
                            migration="U is old enough to predate the last glacial "
                                      "maximum in Europe and is found in some of "
                                      "the earliest European hunter-gatherer "
                                      "remains.",
                            children=[
                                _mt_node(
                                    "K", [("T9055A", ["rs28358571"], 9055, "A", "T")],
                                    description="Haplogroup K — a major subclade of "
                                                "U8.",
                                    migration="K is common in Europe and reaches "
                                              "high frequency in Ashkenazi Jewish "
                                              "populations, where a few founder "
                                              "lineages dominate."),
                            ]),
                    ]),
            ]),
    ])


def _mt_lookup(node: dict, snps_df) -> tuple:
    """Resolve one node's markers. Returns (status, n_derived, n_ancestral,
    evidence) using the same vocabulary as the Y module."""
    n_der = n_anc = 0
    evidence: list[str] = []
    for (name, rsids, pos, derived, ancestral) in node.get("markers", []):
        hit = _find_marker(snps_df, {"name": name, "rsids": rsids, "pos": pos,
                                     "derived": derived, "ancestral": ancestral})
        if hit is None:
            continue
        if hit["is_derived"]:
            n_der += 1
            evidence.append(f"{name}={hit['genotype']} (derived)")
        else:
            n_anc += 1
            evidence.append(f"{name}={hit['genotype']} (ancestral)")
    if n_der:
        return "derived", n_der, n_anc, evidence
    if n_anc:
        return "ancestral", n_der, n_anc, evidence
    return "not_found", 0, 0, evidence


def _mt_entry(node: dict, status: str, n_der: int, evidence: list[str]) -> dict:
    return {
        "haplogroup": node["haplogroup"],
        "snp_name": node["snp_name"],
        "rsids": node.get("rsids", []),
        "pos": node.get("pos"),
        "snp_status": "confirmed" if status == "derived" else "inferred",
        "n_derived": n_der,
        "found_genotype": (evidence[0].split("=", 1)[1].split(" ")[0]
                           if evidence else ""),
        "evidence": evidence,
        "description": node.get("description", ""),
        "migration": node.get("migration", ""),
        "further": node.get("further", ""),
    }


def _mt_walk(node: dict, snps_df, prefix: list[dict], depth: int = 0):
    """Walk the maternal tree, returning (path, status).

    Mirrors the Y-DNA walker deliberately, including its central rule: never
    descend past a node unless a marker below it is genuinely confirmed, so
    the reported lineage can never run deeper than the evidence. Returns None
    when a node is excluded by an ancestral call.
    """
    if depth > 20:
        return prefix, "max_depth"
    status, n_der, _n_anc, evidence = _mt_lookup(node, snps_df)
    if status == "ancestral" and node.get("markers"):
        return None
    entry = _mt_entry(node, status, n_der, evidence)
    path = [*prefix, entry]
    children = node.get("children", [])
    if not children:
        return path, ("resolved" if status == "derived" else "chip_gap")

    confirmed, gaps = [], []
    for child in children:
        cs, cder, _ca, _ce = _mt_lookup(child, snps_df)
        if cs == "derived":
            confirmed.append((child, cder))
        elif cs == "not_found":
            gaps.append(child)

    if confirmed:
        confirmed.sort(key=lambda x: -x[1])
        if len(confirmed) > 1:
            # One person, one maternal line. Two confirmed sibling branches
            # cannot both be true, and silently picking one would hide a real
            # data problem.
            entry["contradiction"] = [c["haplogroup"] for c, _ in confirmed]
        res = _mt_walk(confirmed[0][0], snps_df, path, depth + 1)
        return res if res else (path, "resolved")

    best, best_below = None, 0
    for child in gaps:
        res = _mt_walk(child, snps_df, path, depth + 1)
        if not res:
            continue
        below = sum(1 for n in res[0][len(path):] if n["snp_status"] == "confirmed")
        if below > best_below:
            best_below, best = below, res
    if best and best_below:
        return best
    return path, ("resolved" if status == "derived" else "chip_gap")


def _classify(matched: list[dict]) -> tuple:
    """From the set of confirmed-derived markers, pick the most-specific
    haplogroup call. Returns (haplogroup, confidence)."""
    derived_haplogroups: list[str] = []
    for m in matched:
        if m["status"] == "derived":
            # marker['defines'] may contain "/" — take the first clean token
            defines = m["defines"]
            token = defines.split()[0].split("/")[0]
            derived_haplogroups.append(token)

    if not derived_haplogroups:
        return "Unknown", "low"

    # Priority order from most-specific subclade to broadest macro
    priority = ["K", "V", "I", "W", "X", "T", "J", "H", "HV", "U", "R", "N", "M", "L"]
    for hg in priority:
        if hg in derived_haplogroups:
            # Confidence depends on how many supportive markers we matched
            support = sum(1 for m in matched if m["status"] in ("derived", "ancestral"))
            confidence = "high" if support >= 3 else ("moderate" if support >= 2 else "low")
            return hg, confidence

    return derived_haplogroups[0], "low"


def analyze_mt_haplogroup(snps_df: pd.DataFrame) -> dict:
    """Walk MTDNA_MARKERS against the parsed SNPs and return a result dict."""
    mt_df = snps_df[snps_df["chrom"].isin(["MT", "M", "chrMT", "chrM"])]
    mt_count = len(mt_df)

    matched_markers: list[dict] = []
    for marker in MTDNA_MARKERS:
        hit = _find_marker(snps_df, marker)
        if hit is None:
            continue
        matched_markers.append({
            "name": marker["name"],
            "haplogroup_marker": f"{marker['name']} → {marker['defines']}",
            "rsid": hit["rsid"],
            "genotype": hit["genotype"],
            "is_derived": hit["is_derived"],
            "status": "derived" if hit["is_derived"] else "ancestral",
            "defines": marker["defines"],
            "description": marker["description"],
        })

    if not matched_markers and mt_count == 0:
        return {
            "status": "no_data",
            "haplogroup": "Unknown",
            "confidence": "low",
            "matched_markers": [],
            "path": [],
            "haplogroup_path": "Unknown",
            "terminal_haplogroup": "Unknown",
            "walk_status": "no_data",
            "chip_gaps": [],
            "contradictions": [],
            "n_branch_points": 0,
            "n_confirmed_branch_points": 0,
            "mt_snp_count": 0,
            "migration": HAPLOGROUP_INFO["Unknown"]["migration"],
            "ancient_dna": "",
            "further_testing": HAPLOGROUP_INFO["Unknown"]["further_testing"],
            "message": (
                "No mtDNA SNPs found on this chip — most autosomal arrays include "
                "few or no mtDNA markers."
            ),
        }

    # Walk the phylogeny for the full maternal chain. The flat classifier
    # below still produces the terminal label, so every existing consumer of
    # this result keeps working; the tree adds the branch points that label
    # was hiding.
    walk = _mt_walk(MT_TREE, snps_df, [])
    path, walk_status = walk if walk else ([], "chip_gap")
    haplogroup_path = " > ".join(n["haplogroup"] for n in path) or "Unknown"
    terminal = path[-1]["haplogroup"] if path else "Unknown"
    chip_gaps = [n["haplogroup"] for n in path
                 if n["snp_status"] != "confirmed" and n["haplogroup"] != "mt-MRCA"]
    contradictions = [{"at": n["haplogroup"], "branches": n["contradiction"]}
                      for n in path if n.get("contradiction")]

    haplogroup, _ = _classify(matched_markers)
    # Prefer the tree's terminal call when the walk actually confirmed
    # something: it respects the phylogeny, whereas the flat classifier picks
    # off a fixed priority list and cannot tell an R lineage from an X one.
    if terminal not in ("Unknown", "mt-MRCA") and any(
            n["snp_status"] == "confirmed" for n in path):
        haplogroup = terminal
    info = HAPLOGROUP_INFO.get(haplogroup, HAPLOGROUP_INFO["Unknown"])

    # Coverage / confidence. Confidence in a *called* haplogroup is driven by
    # how many DERIVED markers support it (ancestral markers only rule clades
    # out). Autosomal chips carry very few mtDNA markers, so we never claim
    # "high" without at least 3 confirming derived markers.
    n_derived = sum(1 for m in matched_markers if m["status"] == "derived")
    n_expected = len(MTDNA_MARKERS)
    if haplogroup == "Unknown" or n_derived == 0:
        confidence = "low"
    elif n_derived >= 3:
        confidence = "high"
    elif n_derived == 2:
        confidence = "moderate"
    else:
        confidence = "low"

    return {
        "status": "called" if haplogroup != "Unknown" else "low_resolution",
        "haplogroup": haplogroup,
        "confidence": confidence,
        "confidence_note": (
            f"{n_derived} derived (haplogroup-defining) marker(s) of "
            f"{len(matched_markers)} matched, out of {n_expected} on this panel. "
            "mtDNA calls from autosomal-chip data are coarse — full mtDNA "
            "sequencing is the gold standard."
        ),
        # Full maternal lineage, root to terminal — the chain the Y-DNA
        # section has always shown and this one did not.
        "path": path,
        "haplogroup_path": haplogroup_path,
        "terminal_haplogroup": terminal,
        "walk_status": walk_status,
        "chip_gaps": chip_gaps,
        "contradictions": contradictions,
        "n_branch_points": len([n for n in path if n["haplogroup"] != "mt-MRCA"]),
        "n_confirmed_branch_points": len(
            [n for n in path if n["snp_status"] == "confirmed"]),
        "matched_markers": matched_markers,
        "n_markers_matched": len(matched_markers),
        "n_markers_derived": n_derived,
        "n_markers_expected": n_expected,
        "mt_snp_count": int(mt_count),
        "migration": info["migration"],
        "ancient_dna": info["ancient_dna"],
        "further_testing": info["further_testing"],
        "message": (
            f"mtDNA haplogroup called as {haplogroup} ({confidence} confidence) "
            f"from {n_derived} derived / {len(matched_markers)} matched markers."
        ),
    }
