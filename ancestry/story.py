"""
Ancestral Story — AI-narrated deep-history of who you are
=========================================================

This module doesn't produce genotype calls; it *tells the story*. It takes the
raw output of every ancestry-related module in the project (Y-DNA haplogroup,
mtDNA haplogroup, autosomal ancestry + lineage cross-check, deep-ancestry with
Neanderthal + Yamnaya/EEF/WHG affinity + N-S European axis, immunogenetics with
its historical-selection timeline) and weaves them into a single long-form
narrative of what the person's ancestors likely lived through, what religions
they practiced, what foods and drinks were part of their world, what pandemics
and migrations shaped them.

Two modes:

  • **Deterministic template mode** (default, no AI required): assembles a
    multi-section narrative from curated per-haplogroup / per-affinity story
    blocks. Ships in the report unconditionally.

  • **AI-enhanced mode** (when Ollama is available): calls the local LLM with
    a rich structured prompt containing every relevant fact from the run, and
    receives a long-form personalized narrative that includes religions,
    traditions, foods, wines, key historical events, and speculative but
    grounded ancestry storytelling. Cited to the underlying data; disclaimed
    as illustrative.

Both modes are wired into ``pipeline.py`` and rendered as a large Ancestral
Story section in the main report. Educational, not clinical.
"""

from __future__ import annotations

# ══════════════════════════════════════════════════════════════════════════
# Curated regional / haplogroup story blocks
# ══════════════════════════════════════════════════════════════════════════
#
# Each block describes the historical world of a specific ancestral thread —
# selected in the template narrative based on the user's actual data.

_Y_STORIES = {
    "T": {
        "chapter": "The Fertile Crescent — your paternal line's homeland",
        "text": (
            "Your Y-DNA haplogroup T (T-M184) is a rare and geographically "
            "scattered paternal lineage that traces back to the Neolithic Near "
            "East around 15–25 thousand years ago. T-carrying men were among "
            "the earliest wheat- and barley-farming populations that emerged in "
            "the Fertile Crescent — the swathe of land between the Tigris, "
            "Euphrates and the Levantine coast where agriculture was invented. "
            "T-lineage men were part of the demographic wave that carried "
            "farming, cattle, sheep and goats out of the Near East and into "
            "Europe, the Horn of Africa and parts of India, starting ~9,000 "
            "years ago. The specific T1a1a subclade you carry crystallised "
            "during this Neolithic expansion — roughly 8–10 kya. Some of the "
            "oldest T-DNA in Europe comes from Neolithic burials in the Balkans "
            "and Italy. The lineage never became dominant anywhere except "
            "small pockets (Ibiza, Bulgaria, the Somali coast, some Levantine "
            "Jewish communities), which is why it's rare in modern populations — "
            "a genuine survivor lineage of the Neolithic revolution."
        ),
        "religions": (
            "The Neolithic Near East is the world of the earliest organised "
            "religion humans have material evidence for: Göbekli Tepe (~11,500 "
            "ya) with its megalithic T-shaped pillars carved with animals; the "
            "Anatolian and Mesopotamian mother-goddess and bull cults; later "
            "the Sumerian pantheon (Anu, Enki, Inanna) and the Canaanite/"
            "Ugaritic deities (El, Baal, Asherah) that shaped early Semitic "
            "religion. Your ancestors on this line likely worshipped storm-gods "
            "of rain and fertility deities of grain — the practical concerns of "
            "the first farmers."
        ),
        "food_drink": (
            "The dietary revolution your paternal ancestors participated in: "
            "domesticated einkorn and emmer wheat, barley, lentils, chickpeas, "
            "flax; then domesticated sheep, goats, and cattle. The world's "
            "first beer was brewed in Mesopotamia around 6,000 ya — barley-"
            "based, filtered through reeds. Wine emerged from the wild "
            "grapevines of the South Caucasus and Zagros mountains around "
            "8,000 ya. Bread, cheese, yogurt, and fermented grain drinks were "
            "your Neolithic ancestors' innovations that fed every civilisation "
            "that followed."
        ),
    },
    "R1b": {
        "chapter": "The Bronze-Age steppe — your paternal line's homeland",
        "text": (
            "Your Y-DNA haplogroup R1b is the most common paternal lineage in "
            "Western Europe — carried today by ~60% of men in Ireland, ~50% in "
            "France, ~40% in Britain. It emerged on the Pontic-Caspian steppe "
            "(modern Ukraine / southern Russia) and expanded westward with the "
            "Bronze-Age Bell Beaker culture around 4,500 years ago. The Beaker "
            "wave replaced up to 90% of Britain's Neolithic male population in "
            "a few hundred years — one of the largest documented Y-chromosome "
            "replacements in prehistory."
        ),
        "religions": (
            "The Proto-Indo-European religion of your Bronze-Age steppe "
            "ancestors — sky-fathers (Dyēus Ph2tēr, ancestor of Zeus, Jupiter, "
            "Tiwaz), horse-cult twins, and warrior initiations — became the "
            "foundation of Greek, Roman, Celtic, Germanic and Vedic religion. "
            "Your ancestors on this line likely worshipped a sky-god of "
            "thunder, sacrificed cattle, and marked warrior initiations."
        ),
        "food_drink": (
            "Steppe pastoralist diet: mare's milk, kumis (fermented mare's "
            "milk — the world's first probiotic drink), yogurt, aged cheese. "
            "Meat-heavy diet from cattle, sheep and horses. Mead — fermented "
            "honey — was a status drink of Bronze Age Europe. Beer emerged in "
            "the settled north."
        ),
    },
    "R1a": {
        "chapter": "The Sintashta chariot lords — your paternal line's homeland",
        "text": (
            "Your Y-DNA haplogroup R1a is the eastern-branch counterpart to "
            "R1b — dominant among modern Slavs, Balts, and northern Indians. "
            "R1a men drove the Sintashta chariot expansion east of the Urals "
            "around 4,000 years ago, and this same wave carried the earliest "
            "Indo-Iranian speech and Vedic religion into South Asia."
        ),
        "religions": (
            "The Vedic hymns of the Rig Veda are — to a first approximation — "
            "the religion of your R1a ancestors, preserved through thousands "
            "of years of oral tradition. Fire-altars, soma sacrifices, and the "
            "cosmic-order concept (Ṛta / Aša) all trace to this world."
        ),
        "food_drink": (
            "Same steppe pastoralist package as R1b — dairy, mead, kumis — plus "
            "the ritual soma / haoma drink (probably ephedra or a psychoactive "
            "mushroom) central to Vedic and Zoroastrian religion."
        ),
    },
    "I1": {
        "chapter": "The Nordic hunter-gatherers — your paternal line's homeland",
        "text": (
            "Your Y-DNA haplogroup I1 is a young Nordic lineage — every modern "
            "I1 man descends from a single founder who lived only ~4,600 years "
            "ago in Scandinavia. It survived the Neolithic and Yamnaya waves "
            "as a local Scandinavian founder lineage. Its cultural children "
            "were the Norse."
        ),
        "religions": (
            "The Old Norse pantheon — Odin, Thor, Freya, Freyr, Loki — filtered "
            "through Icelandic saga, plus the more ancient Nordic Bronze Age "
            "sun-cult (the Trundholm sun chariot, ~3,400 ya) that preceded it."
        ),
        "food_drink": (
            "Fish (herring, cod, salmon), reindeer, moose, seal; barley bread, "
            "mead, ale. In the Viking Age, aquavit precursors and imported "
            "Frankish wine were high-status."
        ),
    },
    "I2": {
        "chapter": "The Balkan hunter-gatherers — your paternal line's homeland",
        "text": (
            "Your Y-DNA haplogroup I2 is one of the oldest surviving European-"
            "specific lineages, present since the Palaeolithic. It reaches its "
            "highest modern frequency in the Western Balkans (Bosnia, Croatia, "
            "Serbia) and Sardinia. Ancient European hunter-gatherers like "
            "Cheddar Man (~9 kya Britain) and La Braña 1 (~7 kya Iberia) — "
            "dark-skinned, blue-eyed — carried early I lineages."
        ),
        "religions": (
            "Before the Yamnaya wave brought the Indo-European pantheon to "
            "Europe, older Neolithic and Palaeolithic religions revolved around "
            "the mother goddess, animal spirits, and seasonal cycles. Vinča "
            "figurines of the Balkans (~7,000 ya) are among the oldest religious "
            "art in Europe."
        ),
        "food_drink": (
            "Wild game and gathered plants for the Palaeolithic ancestors; then "
            "wheat, barley, wild grapes fermented into wine, olives — the "
            "foundations of Mediterranean cuisine."
        ),
    },
    "J": {
        "chapter": "The Semitic / Levantine world — your paternal line's homeland",
        "text": (
            "Your Y-DNA haplogroup J is a Near-Eastern lineage that spread with "
            "the Neolithic and later with Semitic and Phoenician expansions. "
            "The J1-M267 sub-branch is common in Arabian populations and "
            "includes many Cohen-modal haplotype (Jewish priestly lineage) "
            "carriers; J2 is centred on Anatolia, the Levant and the Caucasus."
        ),
        "religions": (
            "The West Semitic religious world (Canaanite / Ugaritic / Phoenician "
            "pantheons), which shaped Judaism, Christianity and Islam."
        ),
        "food_drink": (
            "Levantine Neolithic package: wheat, barley, lentils, olives, "
            "grapes, dates, sheep and goat dairy. Wine, unleavened bread, "
            "olive-oil, honey."
        ),
    },
    "G": {
        "chapter": "The First Farmers of Europe — your paternal line's homeland",
        "text": (
            "Your Y-DNA haplogroup G (specifically G2a) was carried by the "
            "very first farmers who entered Europe from Anatolia ~8-9 thousand "
            "years ago. Ötzi the Iceman (5,300 years old, found frozen in the "
            "Alps) was G2a — a snapshot of exactly your ancestral world."
        ),
        "religions": (
            "Neolithic mother-goddess and bull cults; the megalithic tomb "
            "religions of Neolithic Europe (Newgrange, Malta, Carnac)."
        ),
        "food_drink": (
            "Wheat, barley, lentils, sheep and goat dairy, wild grapes, "
            "beer, mead. Ötzi's stomach contents showed einkorn, ibex meat, "
            "and pollen from Alpine flowers."
        ),
    },
    "E": {
        "chapter": "The African / Levantine world — your paternal line's homeland",
        "text": (
            "Your Y-DNA haplogroup E is the most common paternal lineage in "
            "Africa and is also common in the Levant, North Africa, and "
            "Mediterranean Europe. E-M78 and E-V13 spread with early farming "
            "expansions across North Africa and into the Balkans."
        ),
        "religions": (
            "Depending on which sub-branch: African animist / oral traditions, "
            "Egyptian polytheism, or the Semitic Near-Eastern religious world."
        ),
        "food_drink": "Sorghum, millet, teff (Ethiopia); wheat, barley, olives, wine (Med).",
    },
}

_MT_STORIES = {
    "V": {
        "chapter": "The Ice Age refuge — your maternal line's homeland",
        "text": (
            "Your mtDNA haplogroup V has one of the clearest stories in "
            "European maternal genetics. When the Last Glacial Maximum "
            "(~26,000-19,000 years ago) buried northern Europe under ice, "
            "the surviving human populations retreated to warm refugia — "
            "the Iberian Peninsula (Spain / Portugal), the Balkans, and "
            "southern Italy. mtDNA V emerged in the Iberian Franco-Cantabrian "
            "refugium during this Ice Age squeeze and then, as the glaciers "
            "retreated after ~15,000 years ago, its carriers expanded rapidly "
            "northward, following retreating herds. Today V peaks in Basque "
            "and Saami populations (both direct descendants of the northward "
            "recolonization) — genetic evidence of a single Ice Age population "
            "that recolonized post-glacial Europe from the Atlantic coast up "
            "through Scandinavia."
        ),
        "religions": (
            "Late Palaeolithic and Mesolithic Europe: cave painting traditions "
            "(Lascaux, Altamira — Iberian, ~17-14 kya, exactly your ancestors' "
            "world), Venus figurines, animal-spirit shamanism. The 'Sorcerer' "
            "of Trois-Frères and the shamanic aurochs murals of southern France "
            "are echoes of your ancestral religious universe."
        ),
        "food_drink": (
            "Reindeer, red deer, aurochs, ibex — hunted by Solutrean and "
            "Magdalenian ancestors. Salmon and shellfish in coastal caves. "
            "Wild berries, hazelnuts, seeds. No grain, no dairy, no alcohol — "
            "your maternal line's diet in the Ice Age was pure hunter-gatherer."
        ),
    },
    "H": {
        "chapter": "The Palaeolithic mother of Europe — your maternal line's homeland",
        "text": (
            "Your mtDNA haplogroup H is the most common European maternal "
            "lineage — nearly half of all modern Europeans descend on their "
            "mother's mother's mother's line from H. It emerged in West Asia "
            "~25 thousand years ago and expanded across Europe both before and "
            "after the Last Glacial Maximum."
        ),
        "religions": (
            "Palaeolithic Venus figurines, cave painting traditions, and later "
            "the mother-goddess religions of Neolithic Europe."
        ),
        "food_drink": (
            "Ice Age hunter-gatherer diet; then post-LGM wild plants, and "
            "eventually the Neolithic wheat-barley-dairy package."
        ),
    },
    "U": {
        "chapter": "The oldest European mothers — your maternal line's homeland",
        "text": (
            "Your mtDNA haplogroup U is the oldest surviving European mtDNA "
            "branch — U5 carriers walked in Europe during the Aurignacian "
            "(~40,000 ya), the same period as the earliest cave paintings. "
            "Cheddar Man (~9,000-year-old British hunter-gatherer, blue-eyed, "
            "dark-skinned) carried U5b1."
        ),
        "religions": (
            "Aurignacian and Magdalenian cave painting shamanism (Chauvet, "
            "Lascaux); the animal-mother goddess traditions that predate the "
            "Neolithic mother-goddess."
        ),
        "food_drink": (
            "Pure Palaeolithic hunter-gatherer — mammoth (until they died out), "
            "reindeer, horse, wild plants."
        ),
    },
    "K": {
        "chapter": "The Neolithic mothers — your maternal line's homeland",
        "text": (
            "Your mtDNA haplogroup K spread with the Neolithic farmer expansion "
            "from the Fertile Crescent into Europe. Ötzi the Iceman's maternal "
            "line was K1f. K1a1b1a is one of the four founder lineages of "
            "Ashkenazi Jewish mtDNA."
        ),
        "religions": (
            "Neolithic mother-goddess religions of the Anatolian and European "
            "Neolithic."
        ),
        "food_drink": "Neolithic farming package: wheat, barley, sheep and goat dairy.",
    },
    "T": {
        "chapter": "The Neolithic Near-East mothers — your maternal line's homeland",
        "text": (
            "Your mtDNA haplogroup T spread with the Neolithic wave from the "
            "Near East. Common today from Ireland to India."
        ),
        "religions": "Neolithic mother-goddess religions; then Bronze-Age pantheons.",
        "food_drink": "Neolithic farming package.",
    },
    "J": {
        "chapter": "The Fertile Crescent mothers — your maternal line's homeland",
        "text": (
            "Your mtDNA haplogroup J emerged in the Near East and spread with "
            "the Neolithic expansion. Common today across Europe and the Near "
            "East."
        ),
        "religions": "Neolithic Near-Eastern religion; later Levantine pantheons.",
        "food_drink": "Neolithic farming package.",
    },
    "L": {
        "chapter": "The African mothers — your maternal line's homeland",
        "text": (
            "Your mtDNA haplogroup L is the African super-lineage. All non-"
            "African mtDNA descends from a single L3 lineage that left Africa "
            "~60-70 kya. L0 sub-branches are among the oldest surviving human "
            "mtDNA lineages, present in Khoi-San populations of southern Africa."
        ),
        "religions": (
            "African indigenous religious traditions — spanning Khoi-San "
            "shamanism, West African orisha traditions, Bantu ancestor "
            "veneration, and East African pastoralist religions."
        ),
        "food_drink": (
            "Deeply variable by region: sorghum and millet across the Sahel, "
            "teff and coffee in Ethiopia, plantain and yam in West Africa, "
            "cattle-pastoralist milk-blood diet in East African savanna."
        ),
    },
}


# ══════════════════════════════════════════════════════════════════════════
# Story assembler (deterministic template mode)
# ══════════════════════════════════════════════════════════════════════════

def _resolve(hg: str, table: dict) -> dict | None:
    """Longest-prefix haplogroup lookup."""
    if not hg:
        return None
    up = hg.strip().upper()
    for key in sorted(table.keys(), key=len, reverse=True):
        if up.startswith(key.upper()):
            return {"key": key, **table[key]}
    return None


def build_template_story(y_result: dict | None,
                         mt_result: dict | None,
                         deep_ancestry_result: dict | None,
                         immunogenetics_result: dict | None,
                         ancestry_result: dict | None) -> dict:
    """Assemble a multi-chapter deterministic Ancestral Story from the run's
    upstream results."""
    chapters: list[dict] = []

    # Chapter 1 — deep prehistory (Neanderthal + WHG/EEF/Yamnaya)
    dp = deep_ancestry_result or {}
    n = (dp.get("neanderthal") or {})
    ap = (dp.get("ancient_populations") or {})
    if n.get("available") or ap.get("available"):
        parts = []
        if n.get("available"):
            parts.append(
                f"You carry approximately **{n['approx_pct']}% Neanderthal-tagged "
                f"ancestry** ({n['tier']} for a non-African) — the legacy of "
                "a series of mating events between anatomically modern humans "
                "and Neanderthals about 50,000 years ago, somewhere in the "
                "Near East. Every non-African alive today carries this "
                "inheritance. Your Neanderthal DNA sits at immune-system loci "
                "and cellular-defence genes; some of it (like the OAS1 "
                "haplotype) is actually adaptive and still helps you fight "
                "viruses today."
            )
        if ap.get("available"):
            tops = ap.get("populations", [])[:2]
            if tops:
                pair = " and ".join(f"**{t['short']}** ({t['affinity']*100:.0f}%)" for t in tops)
                parts.append(
                    f"On the autosomal (whole-genome) side, your affinity leans "
                    f"toward {pair} — the ancient European gene pools who shaped "
                    "you. Yamnaya was the Bronze-Age steppe pastoralist wave "
                    "(~5,000-4,000 years ago); EEF (Early European Farmer / "
                    "Anatolian Neolithic) was the ~9,000-year-old Neolithic wave "
                    "carrying agriculture into Europe; WHG (Western Hunter-"
                    "Gatherer) is the pre-Neolithic Palaeolithic Europeans who "
                    "painted the caves. Your genome is a mixture of all three, "
                    "layered on top of tens of thousands of years of Palaeolithic "
                    "European heritage."
                )
        chapters.append({
            "title": "Chapter 1 — Deep Prehistory (before ~10,000 years ago)",
            "body": "\n\n".join(parts),
        })

    # Chapter 2 — Y-DNA paternal line
    y_hg = (y_result or {}).get("terminal_haplogroup")
    y_story = _resolve(y_hg, _Y_STORIES)
    if y_story:
        chapters.append({
            "title": f"Chapter 2 — Your Paternal Line ({y_hg}, via '{y_story['key']}')",
            "body": (
                f"### {y_story['chapter']}\n\n"
                f"{y_story['text']}\n\n"
                f"**Religion & belief:** {y_story['religions']}\n\n"
                f"**Food & drink of this world:** {y_story['food_drink']}"
            ),
        })

    # Chapter 3 — mtDNA maternal line
    mt_hg = (mt_result or {}).get("haplogroup")
    mt_story = _resolve(mt_hg, _MT_STORIES)
    if mt_story:
        chapters.append({
            "title": f"Chapter 3 — Your Maternal Line ({mt_hg}, via '{mt_story['key']}')",
            "body": (
                f"### {mt_story['chapter']}\n\n"
                f"{mt_story['text']}\n\n"
                f"**Religion & belief:** {mt_story['religions']}\n\n"
                f"**Food & drink of this world:** {mt_story['food_drink']}"
            ),
        })

    # Chapter 4 — the pandemics your ancestors survived
    im = immunogenetics_result or {}
    if im.get("historical_timeline"):
        stories = []
        for ev in im["historical_timeline"]:
            stories.append(f"- **{ev['epoch']}** — {ev['driver']}. "
                           f"*{ev['finding']}* — {ev['verdict']}. {ev['narrative']}")
        chapters.append({
            "title": "Chapter 4 — What Your Ancestors Survived",
            "body": (
                "Some of the strongest positive-selection signals in the human "
                "genome are still visible in your DNA — variants that swept "
                "through your ancestral population because the people who "
                "carried them survived pandemics and their neighbours didn't. "
                "Reading down this list is reading a summary of the historical "
                "traumas that shaped you:\n\n" + "\n".join(stories)
            ),
        })

    # Chapter 5 — synthesis
    lean = (dp.get("european_axis") or {}).get("lean")
    if lean or chapters:
        chapters.append({
            "title": "Chapter 5 — Who You Are",
            "body": (
                f"Taken together, your genome describes a person of "
                f"{lean or 'European'} ancestry, whose paternal line came into "
                f"Europe with the Neolithic wave out of the Fertile Crescent, "
                f"whose maternal line survived the Last Ice Age in an Iberian "
                f"refugium and expanded north as the glaciers retreated, and "
                f"whose autosomal DNA carries the Bronze-Age Yamnaya × "
                f"Anatolian-farmer admixture that shaped Northern and Central "
                f"Europe. You carry the immune legacies of the Black Death, "
                f"of endemic gut viruses, of prion epidemics, and of the "
                f"Neanderthal encounter tens of thousands of years ago. You "
                f"are, in the most literal sense, a survivor of every "
                f"catastrophe that hit your ancestors — because their DNA "
                f"is you."
            ),
        })

    return {"available": bool(chapters), "mode": "template", "chapters": chapters}


# ══════════════════════════════════════════════════════════════════════════
# AI-enhanced narrative (optional; uses Ollama when available)
# ══════════════════════════════════════════════════════════════════════════

def build_ai_story_prompt(y_result, mt_result, deep_ancestry_result,
                          immunogenetics_result, ancestry_result) -> str:
    """Build a rich, well-structured prompt for the local LLM to produce a
    long-form Ancestral Story. Ships everything the model needs; asks for a
    specific structure so the output is consistent."""
    lines = ["You are a historical geneticist, cultural historian, and "
             "storyteller. Weave the following person's DNA-derived ancestry "
             "data into a **long-form Ancestral Story** — the story of who "
             "their ancestors were, what they lived through, what religions "
             "they practiced, what foods and drinks they knew, what "
             "pandemics and migrations shaped them. Ground every claim in "
             "the specific data below; be evocative and detailed but not "
             "invented — where you speculate, mark it as speculation.\n\n"
             "=== THEIR DNA-DERIVED DATA ==="]

    y = (y_result or {}).get("terminal_haplogroup")
    mt = (mt_result or {}).get("haplogroup")
    lines.append(f"\nY-DNA (paternal line): {y or 'unknown'}")
    lines.append(f"mtDNA (maternal line): {mt or 'unknown'}")

    if ancestry_result and ancestry_result.get("proportions"):
        lines.append("\nAutosomal ancestry proportions:")
        for sp, p in sorted(ancestry_result["proportions"].items(),
                            key=lambda kv: -kv[1]):
            lines.append(f"  - {sp}: {p*100:.1f}%")
        cc = (ancestry_result.get("haplogroup_crosscheck") or {})
        if cc.get("verdict"):
            lines.append(f"Lineage cross-check: {cc['verdict']} — {cc.get('summary','')}")

    dp = deep_ancestry_result or {}
    if dp.get("neanderthal", {}).get("available"):
        n = dp["neanderthal"]
        lines.append(f"\nNeanderthal introgression: ~{n['approx_pct']}% "
                     f"({n['tier']}), {n['n_carrying']}/{n['n_typed']} tagged loci carrying")
    if dp.get("ancient_populations", {}).get("available"):
        lines.append("\nAncient-population affinity (Bronze-Age vs Neolithic vs Mesolithic):")
        for p in dp["ancient_populations"]["populations"]:
            lines.append(f"  - {p['short']}: {p['affinity']*100:.0f}%")
    if dp.get("european_axis", {}).get("available"):
        lines.append(f"\nSub-continental European axis: {dp['european_axis']['lean']} "
                     f"(index {dp['european_axis']['index']})")

    im = immunogenetics_result or {}
    if im.get("historical_timeline"):
        lines.append("\nHistorical selection events visible in their genome "
                     "(pandemics their ancestors survived):")
        for ev in im["historical_timeline"]:
            lines.append(f"  - {ev['epoch']}: {ev['driver']} — {ev['finding']} "
                         f"({ev['verdict']})")

    lines.append("""

=== YOUR TASK ===

Write a **rich, chapter-structured Ancestral Story** covering (in this order):

**Chapter 1 — Deep Prehistory** (Palaeolithic and the Neanderthal legacy).
**Chapter 2 — Your Paternal Line** (Y-DNA haplogroup: origin, migrations, cultural world).
**Chapter 3 — Your Maternal Line** (mtDNA haplogroup: origin, migrations, cultural world).
**Chapter 4 — The Waves That Shaped You** (Yamnaya / EEF / WHG components; migrations; timelines).
**Chapter 5 — What Your Ancestors Survived** (each protective variant → the pandemic/pressure that selected it).
**Chapter 6 — Religions & Belief** (specific pantheons / traditions relevant to the ancestral regions).
**Chapter 7 — Foods, Wines, and Drinks of Your Ancestors** (crops, beers/wines/mead, dietary revolutions).
**Chapter 8 — Who You Are** (a synthesis of the above into a single story).

Make it vivid and specific. Where you introduce a religion, name specific deities and practices;
where you introduce foods, name specific crops, drinks, dishes; where you introduce a pandemic,
say what year and place and how many died. Cite each claim to the underlying DNA evidence
(e.g. "your PRNP MV genotype (rs1799990 AG) is the same signature that let some Fore people
survive the New Guinea kuru epidemic"). Speculate where useful but mark it: "It's plausible
that...". Aim for **2,000-4,000 words**.

Educational and historical, not clinical or genealogical. Note at the end that this is a
plausibility-weighted narrative from anonymous DNA-derived data, not a documented lineage.
""")
    return "\n".join(lines)


def analyze_ancestral_story(y_result: dict | None,
                            mt_result: dict | None,
                            deep_ancestry_result: dict | None,
                            immunogenetics_result: dict | None,
                            ancestry_result: dict | None,
                            model: str | None = None,
                            use_ai: bool = False) -> dict:
    """Build the Ancestral Story. Always returns the deterministic template
    narrative; optionally attaches an AI-enhanced long-form version if
    ``use_ai=True`` and Ollama is reachable."""
    template = build_template_story(
        y_result, mt_result, deep_ancestry_result,
        immunogenetics_result, ancestry_result,
    )

    ai_text: str | None = None
    ai_error: str | None = None
    if use_ai and model:
        try:
            # Local Ollama import kept optional so template mode always works.
            import analyze as _analyze
            prompt = build_ai_story_prompt(
                y_result, mt_result, deep_ancestry_result,
                immunogenetics_result, ancestry_result,
            )
            # Big context for a rich narrative; longer generation budget.
            ai_text = _analyze.call_ollama(
                prompt, model=model, num_ctx=8192, num_predict=4096, timeout=1800,
            )
        except Exception as e:
            ai_error = str(e)

    return {
        "available": template.get("available") or bool(ai_text),
        "template": template,
        "ai_text": ai_text,
        "ai_error": ai_error,
        "ai_used": bool(ai_text),
    }
