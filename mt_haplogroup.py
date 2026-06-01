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

from typing import Dict, List, Optional
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

MTDNA_MARKERS: List[Dict] = [
    # ── Macro-haplogroups ─────────────────────────────────────────────────────
    {
        "name": "T16223C", "pos": 16223, "rsids": ["rs41323649"],
        "derived": "T", "ancestral": "C",
        "defines": "L/M/N (ancestral state)", "level": "macro",
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
        "name": "C12705T", "pos": 12705, "rsids": ["rs28358580"],
        "derived": "T", "ancestral": "C",
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
HAPLOGROUP_INFO: Dict[str, Dict[str, str]] = {
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

def _normalise_genotype(gt: object) -> Optional[str]:
    if gt is None:
        return None
    s = str(gt).upper().strip().replace(" ", "").replace("-", "")
    if s in ("", "NAN", "--"):
        return None
    return s


def _find_marker(snps_df: pd.DataFrame, marker: Dict) -> Optional[Dict]:
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


def _classify(matched: List[Dict]) -> tuple:
    """From the set of confirmed-derived markers, pick the most-specific
    haplogroup call. Returns (haplogroup, confidence)."""
    derived_haplogroups: List[str] = []
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


def analyze_mt_haplogroup(snps_df: pd.DataFrame) -> Dict:
    """Walk MTDNA_MARKERS against the parsed SNPs and return a result dict."""
    mt_df = snps_df[snps_df["chrom"].isin(["MT", "M", "chrMT", "chrM"])]
    mt_count = len(mt_df)

    matched_markers: List[Dict] = []
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
            "mt_snp_count": 0,
            "migration": HAPLOGROUP_INFO["Unknown"]["migration"],
            "ancient_dna": "",
            "further_testing": HAPLOGROUP_INFO["Unknown"]["further_testing"],
            "message": (
                "No mtDNA SNPs found on this chip — most autosomal arrays include "
                "few or no mtDNA markers."
            ),
        }

    haplogroup, _ = _classify(matched_markers)
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
