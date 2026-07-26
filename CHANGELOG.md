# Changelog

All notable changes to this project are documented here. Format inspired by
[Keep a Changelog](https://keepachangelog.com/); versioning follows the
`{major}.{minor}.{tag}` scheme.

---

## [6.13.0-premium] — 2026-07-25 — Holistic Synthesis: cross-panel patterns + Genome Leverage Score

### Added

- **New `holistic_synthesis.py` module** — the cross-panel pattern-detection
  layer that catches insights only visible when multiple modules are combined:
  * **FUT2 non-secretor × elevated hs-CRP** — non-secretor microbiome
    baseline may inflate CRP; trend matters more than absolute values.
  * **Fasting glucose ≥100 × HbA1c <5.7** — flags acute-stress /
    poor-sleep morning-of-draw pattern rather than real dysglycemia.
  * **APOE ε4 × elevated LDL/ApoB** — flags the ~2× lipid amplification
    that ε4 carriers face.
  * **CHRNA5 A-carrier × non-smoker** — recognises an actively realised
    prevention success.
  * **Ancestral diet fit** — European + LCT + Yamnaya × EEF → Mediterranean
    + dairy diet is genetically appropriate, not just generically healthy.
  * **BDNF Val/Val + adaptive neurotype + young adult** — deliberate
    practice compounds materially more; window matters.
  * **HFE clear × ferritin > 300** — high ferritin without HFE points to
    inflammation, diet, alcohol, or fatty liver, not hemochromatosis.
  * **Coffee protocol synthesis** — COMT × MAOA × CYP1A2 → specific dose
    and timing recommendations.
- **Genome Leverage Score (0-100)** — composite of protective/adverse findings
  across APOE, PRS panels, longevity variants, immunogenetics headlines,
  neurochemistry composite, PhenoAge acceleration, and flagged clinical
  markers. Buckets into Very Favorable / Favorable / Balanced / Actionable Risk
  with a narrative on the "environment = trajectory" leverage implication.
- **Ranked priority actions** — top 6 cross-panel insights ranked by severity ×
  confidence × modifiability × actionability.
- New "Holistic Synthesis" report section rendered immediately after the
  Executive Summary (meta-view leads the report).
- `REPORT_VERSION` → 6.13.0-premium.

### Tests

- +6 unit tests (glucose/HbA1c discordance, FUT2×CRP, APOE×lipid,
  favorable-genome leverage, CHRNA5 prevention, empty-input); 266 passing.

## [6.12.0-premium] — 2026-07-25 — Neurochemistry module (COMT / MAOA / BDNF composite phenotype)

### Added

- **New `neurochemistry.py` module** packaging the "warrior vs worrier" COMT
  literature with MAOA, BDNF, DRD2/DRD4, 5-HTT, HTR2A, TPH2, CACNA1C, OPRM1,
  and CHRNA5 into a composite phenotype with concrete recommendations:
  * Stress-response profile (warrior/worrier/adaptive middle)
  * Plasticity tier (BDNF-driven learning capacity)
  * Stimulant response prediction, SSRI response prediction
  * Caffeine protocol (dose/timing tuned to genotype)
  * Meditation-style fit (focused-attention vs open-monitoring vs somatic)
  * Career neurotype signature
  * Addiction / substance flags (CHRNA5 smoking + OPRM1 opioid/naltrexone)
- 10 new SNPs registered (rs4633, rs6323, rs1800497, rs1800955, rs25531,
  rs1006737, rs6313, rs1799971, rs4570625, rs16969968).
- Wired into pipeline, report section (nav + body), chat context.
- `REPORT_VERSION` → 6.12.0-premium.

### Tests

- +8 unit tests (COMT classes, MAOA, BDNF, CHRNA5, composite completeness); 260 passing.

## [6.11.0-premium] — 2026-07-25 — Immunogenetics + Ancestral Story

Two big additions that connect the ancestry, viral-resistance, and historical-
selection threads into the story no consumer service tells.

### Added

- **New `immunogenetics.py` module** — comprehensive viral / bacterial /
  parasitic resistance + Historical Selection Timeline:
  * **Viral**: HIV/CCR5-Δ32 (rs333), norovirus & rotavirus/FUT2 (rs601338),
    hepatitis C spontaneous clearance/IL28B (rs12979860), hepatitis B/HLA-DPB1
    (rs9277535), COVID-19/OAS1 (rs2660) + 3p21.31 Neanderthal haplotype,
    influenza/IFITM3 (rs12252) + MX1, prion disease/PRNP codon 129 (rs1799990).
  * **Bacterial**: plague/ERAP2 rs2549794 (Black Death survivor allele — Klunk
    2022 Nature), sepsis & RSV/TLR4.
  * **Parasitic**: P. vivax malaria/Duffy (rs2814778), P. falciparum/HbS &
    G6PD deficiency.
  * **Autoimmune trade-offs**: PTPN22, STAT4, IRF5.
  * **Historical Selection Timeline** — every protective variant mapped to the
    historical pandemic that likely selected for it (Black Death, malaria
    endemic zones, kuru, endemic gut viruses, hepatotropic viruses).
- **New `ancestral_story.py` module** — long-form Ancestral Story narrative
  weaving haplogroups, deep-ancestry components (Yamnaya/EEF/WHG), and the
  immunogenetics selection timeline into a chapter-structured story. Two
  modes: deterministic template (always runs) and AI-enhanced (uses local
  Ollama when available — a rich structured prompt asks for religions,
  foods, wines, and historical events to be woven into a 2,000-4,000-word
  narrative).
- New report sections: **Immunogenetics** and **The Ancestral Story**, both
  wired into the pipeline, renderer, and chat context.
- `REPORT_VERSION` → 6.11.0-premium.

### Tests

- +8 unit tests (CCR5, FUT2, IL28B, PRNP, ERAP2, timeline, empty-input, story
  template mode); 252 passing.

## [6.10.0-premium] — 2026-07-25 — Blood type (ABO + RhD + FUT2 secretor)

### Added

- **New `blood_type.py` module** predicting ABO phenotype/genotype (A / B /
  AB / O — including hidden recessive-O carrier flag), Rh(D) status, and
  FUT2 secretor status. Wired into report + chat context.
  - **ABO**: uses rs8176719 (delG frameshift = O allele), rs8176746, rs8176747,
    rs7853989 for A-vs-B backbone discrimination. Correctly handles the
    recessive nature of O (A/O genotype → phenotype A).
  - **Rh(D)**: prefers dedicated tag SNPs (rs590787, rs676785, rs2280330, rs660,
    rs586178) when present. When absent, uses a **coverage-based inference**
    over the RHD gene body (chr1:25.58-25.68 Mb): homozygous RHD-deletion
    individuals have essentially no calls in this region.
  - **FUT2 rs601338 secretor status**: AA = non-secretor (protective against
    noroviruses), G-carriers = secretor.
  - Consolidated blood-group string (e.g. "A+") with population-frequency
    context and prominent "not a clinical typing" disclaimer.
- `REPORT_VERSION` → 6.10.0-premium.

### Tests

- +6 unit tests (Type-A/O, Type-O homozygous, Type-AB, RhD-negative via low
  coverage, secretor status, empty input); 244 passing.

## [6.9.0-premium] — 2026-07-25 — Deep ancestry: Neanderthal + ancient populations + N-S European axis + migration timelines

### Added

- **New `deep_ancestry.py` module** — a state-of-the-art upgrade of the ancestry
  section, comparable in depth to the blood-work engine. Four sub-analyses,
  every SNP grounded in published papers:
  1. **Neanderthal introgression estimate** using 10+ curated Neanderthal-
     derived / -tagged SNPs from Sankararaman 2014, Zeberg & Pääbo Nature 2020
     (3p21.31 COVID haplotype), Vernot & Akey 2014, and BNC2 / immune / X-linked
     adaptive-introgression loci. Reports affinity, an approximate percent
     bucket, and the specific carrying variants — with the honest caveat that
     consumer chips tag only a subset of the ~6,000 Neanderthal-introgressed
     SNPs (23andMe v5 typing).
  2. **Ancient-population affinity** — Yamnaya-Steppe / Anatolian-Neolithic
     Farmer (EEF) / Western Hunter-Gatherer (WHG) fingerprints from
     Mathieson 2015 & Allentoft 2015 (LCT lactase persistence = Yamnaya-derived,
     SLC24A5 = near-fixed in EEF, HERC2 blue eyes = WHG-first, etc.).
  3. **Sub-continental European axis** — a soft Northern-vs-Southern European
     index over LCT/HERC2/TYR/MC1R with a visual position marker.
  4. **Y-DNA & mtDNA migration timelines** — TMRCA (from ISOGG/YFull) and the
     specific migration narrative for the user's haplogroups (T1a1a → Near-East
     Neolithic ~8-15 kya; mtDNA V → Iberian post-LGM refugium ~15 kya; etc.).
- New report section **"Deep Ancestry"** rendered inline after the standard
  ancestry section, with a Neanderthal-affinity hero, ancient-population cards
  with progress bars, a visual N-S axis marker, and haplogroup timeline cards.
  Fed into the chat assistant context so questions about "Neanderthal /
  Yamnaya / Bronze Age / my haplogroups" ground in the user's actual variants.
- `REPORT_VERSION` → 6.9.0-premium.

### Tests

- +6 unit tests (Neanderthal scoring, ancient-pop ranking, N-S axis,
  haplogroup-timeline longest-prefix matching, end-to-end); 238 passing.

## [6.8.0-premium] — 2026-07-25 — Urologic & genitourinary panel

### Added

- **New Urologic module** (`urologic.py`) — a dedicated genotype screen for
  urologic conditions that the other panels didn't cover (specifically requested
  after the chat noted OAB was missing). 5 sub-panels, 15 registered SNPs:
  - **Bladder / OAB** — ADRB3 Trp64Arg (mirabegron target), CHRM3 M3 receptor
    (antimuscarinic OAB drugs), plus NAT2 acetylator × bladder-cancer risk
    (cross-referenced from detox).
  - **Prostate — BPH & cancer** — SRD5A2 V89L (5α-reductase / finasteride
    response) and A49T, **HOXB13 G84E** (the top hereditary prostate-cancer
    marker, flagged high-confidence), 8q24 loci (rs1447295, rs6983267), MSMB
    rs10993994.
  - **Kidney stones** — CLDN14 rs219780 (calcium-stone risk), SLC34A1
    rs4074995, CASR R990G, PKD1 rs2072499.
  - **Testicular germ-cell cancer** — KITLG rs995030, SPRY4 rs4324715.
  - **Androgen bioavailability** — SHBG rs1799941 (free-T interpretation).
- New report section **Urologic & Genitourinary Panel**, plus urologic context
  fed into the chat assistant so questions about OAB / BPH / prostate now
  ground in the user's actual variants.
- `REPORT_VERSION` → 6.8.0-premium.

### Tests

- +5 unit tests (registry consistency + panel behaviour); 232 passing.

## [6.7.2-premium] — 2026-07-25 — Fix cross-category synthesis 400 (per-call num_ctx)

### Fixed

- 6.7.1 lowered the shared `num_ctx` in `call_ollama` to 4096 for speed on
  per-category batches, which was correct for those calls but caused the whole-
  report **cross-category synthesis** (and to a lesser extent the executive
  summary) to be rejected by Ollama with a 400 Bad Request when the assembled
  prompt exceeded 4k tokens.
- `call_ollama` now takes optional `num_ctx` and `num_predict` overrides.
  Per-category batches keep 4096/1024. Executive summary uses 8192/1536.
  Cross-category synthesis uses **16384/1800** and its prompt is capped at
  48k characters as a safety net. `REPORT_VERSION` → 6.7.2-premium.

## [6.7.1-premium] — 2026-07-25 — Tier-2 AI: streaming + smaller batches + auto-halving retry

### Fixed

- **Tier-2 AI runs no longer time out on qwen3:14b in reasoning mode.**
  Three coordinated changes:
  * `call_ollama` now **streams** by default (`stream=True`), with a 90-second
    read-idle guard instead of a wall-clock timeout. Long-generating reasoning
    models no longer silently trip the socket timeout; a call only dies if the
    stream actually stalls.
  * **AI_BATCH_MAX 15 → 8** and `num_ctx 8192 → 4096` — each per-category call
    now reliably finishes in 2–4 min on consumer hardware.
  * **Auto-halving retry:** when a batch does time out, the runner splits it
    in half and retries (down to AI_BATCH_MIN=2) instead of dropping the whole
    batch. Only irreducible failures are recorded for `--retry-failed`.
- `REPORT_VERSION` → 6.7.1-premium.

## [6.7.0-premium] — 2026-07-25 — AI chat assistant: deep mode + streaming + rich context

### Changed

- **Rebuilt the `--chat` REPL** into a genuinely useful assistant:
  - **Deep-mode system prompt** with a structured answering rubric (direct
    answer → what the data shows → mechanism → personal fit → action plan →
    uncertainty → clinician trigger → follow-ups). Default length now 500–1200
    words; a `/brief` toggle switches to concise 3–6-sentence answers.
  - **Streaming responses** — tokens print live as the model generates them,
    with `<think>` blocks filtered and inline bold rendered on the fly.
  - **Vastly richer context**: now surfaces bloodwork (PhenoAge biological age,
    10-yr mortality risk, PhenoAge levers, PREVENT ASCVD, flagged biomarkers
    with their genotype-aware notes), detox / smoke-resilience tier, longevity
    variants (FOXO3, APOE, CETP, KLOTHO, IL6, TP53), personal 10-yr economics,
    plus ancestry lineage cross-check — in addition to the existing PRS/PGx/
    carrier/traits/interactions.
  - **REPL commands**: `/help`, `/deep`, `/brief`, `/model <name>`, `/topic
    <name>`, `/context`, `/save <path>`, `/reset`, `/suggest`.
  - **Follow-up suggestions** printed after every answer, tailored to the
    topic (PGx / cardio / longevity / diet / detox / diabetes / …).
  - `num_ctx` bumped to 16384; temperature 0.4 → 0.35; `num_predict` 2048 for
    deep mode.
- `REPORT_VERSION` → 6.7.0-premium.

## [6.6.1-premium] — 2026-07-25 — Biological-aging economics grounded in the mortality model

### Changed

- The economic sheet's **biological-aging** line is now **high confidence**,
  earned honestly: instead of a hand-wavy per-year cost, it is derived from the
  PhenoAge clock's own validated 10-year mortality output — the person's modeled
  mortality risk vs the baseline for their chronological age — and valued in
  QALYs. PhenoAge is mortality-calibrated (each year of acceleration ≈ 9% higher
  all-cause mortality, HR 1.09/yr; Levine, Aging 2018). `REPORT_VERSION` → 6.6.1.

### Tests

- +1 test locking the grounded, high-confidence biological-aging item (226).

## [6.6.0-premium] — 2026-07-25 — Standalone personal economic-impact sheet

### Added

- **`economic_analysis.html`** — a new standalone page written during a full
  run, modeling the individual's 10-year economic impact of acting on their
  results: expected medical-cost avoidance + monetised quality-of-life (QALY)
  gains, net of intervention cost, across pharmacogenomic/carrier/PRS findings,
  PREVENT cardiovascular event avoidance, prediabetes→T2D prevention, and
  biological-age. Headline net value, per-finding table, and ROI vs the ~$700
  one-time analysis cost. QALYs monetised at the standard $100k/QALY threshold;
  prominent "illustrative model, not financial/medical advice" disclaimer.
- `analyze_personal_economics()` + `render_economic_analysis_html()` in
  `health_economics.py`; wired into the pipeline and `.gitignore`d as a
  genotype-derived output. On a real profile: modeled net ~$115k, ROI ~165×.
- `REPORT_VERSION` → 6.6.0-premium.

### Tests

- +2 unit tests (personal-economics modeling + empty case); 225 passing.

---

## [6.5.0-premium] — 2026-07-25 — AHA PREVENT CVD risk + longitudinal tracking

### Added

- **AHA PREVENT 2023 10-year ASCVD risk** — the current American Heart
  Association guideline model for atherosclerotic-cardiovascular-disease risk.
  Implemented with the **exact sex-specific coefficients** from Khan et al.
  (*Circulation* 2024), sourced from the `preventr` R package (v0.11.0) and a
  cross-validated open-source implementation, with correct variable centering,
  SBP/eGFR splines, and interaction terms; validated numerically (healthy 50F →
  ~1.3%, high-risk 65F → ~17%). Uses TC/HDL/SBP/eGFR (eGFR derived from
  creatinine) with optional smoker/BP-med/statin/diabetes flags; flags the 7.5%
  statin-consideration threshold. (This was the item deliberately deferred in
  6.4.0 pending verified coefficients — now shipped.)
- **Longitudinal tracking.** A `labs.json` may now carry a `"history"` array of
  dated panels (or be a top-level list). The report computes biological age,
  health score, ASCVD risk and every key marker **per visit** and renders a
  Trajectory-Over-Time section with direction-aware SVG sparklines (green =
  improving, red = worsening). Single-panel JSON is unchanged.
- `REPORT_VERSION` → 6.5.0-premium.

### Tests

- +3 unit tests (PREVENT reference values, PREVENT via clinical path,
  longitudinal trajectory); 223 passing. Golden snapshot refreshed.

---

## [6.4.0-premium] — 2026-07-25 — Interactive bio-age simulator + genetics×aging

### Added

- **Interactive in-browser biological-age simulator.** Sliders for the nine
  PhenoAge markers, pre-filled with the user's values; the biological age,
  delta and colour recompute live in the browser using an embedded copy of the
  Levine formula (static fallback in PDF). A hands-on "what moves my age" tool.
- **Genetics × Aging tie-in.** Reads longevity-associated variants from the
  user's own genome — FOXO3 (rs2802292), APOE ε2/ε4, CETP (rs5882), KLOTHO
  KL-VS (rs9536314), IL6 (rs1800795), TP53 (rs1042522) — and shows the inherited
  longevity "lean" alongside the phenotypic PhenoAge clock, with the honest
  caveat that individual effect sizes are modest (meta-analysis of exceptional-
  longevity GWAS, Revelas 2018).
- `REPORT_VERSION` → 6.4.0-premium.

### Deliberately not shipped

- The AHA PREVENT 2023 10-year CVD risk equations were researched but **not**
  implemented: their exact coefficients live in journal supplementary tables
  that couldn't be verified end-to-end, and a CVD-risk calculator built on
  guessed coefficients would be clinically misleading. Left as a documented
  next step pending the verified coefficient set.

### Tests

- +2 unit tests (genetic-longevity reader, simulator render); 220 passing.
  Golden snapshot refreshed.

---

## [6.3.0-premium] — 2026-07-25 — Longevity simulator, visuals & expanded index panel

Builds on 6.2.0's biological-age clock with an interactive-style longevity
simulator, mortality estimation, data-visualisation, and four more validated
indices — 20 composite scores in total.

### Added

- **Longevity levers (counterfactual simulator).** For each of the nine
  PhenoAge biomarkers, computes the biological-age cost vs its optimal value
  (holding the others fixed) and the total years recoverable if all were
  optimal — a quantified, actionable "what-if" ranked by impact.
- **Estimated 10-year mortality risk** surfaced directly from the PhenoAge
  Gompertz mortality score.
- **Data-visualisation (inline SVG, no dependencies):** a biological-age gauge
  (younger↔older with chronological tick + biological marker) and a
  body-system radar/spider chart of the 12 panel scores.
- **Four more validated indices:** METS-IR (Bello-Chavolla 2018), Fatty Liver
  Index (Bedogni 2006), Aggregate Index of Systemic Inflammation (AISI), and
  the Prognostic Nutritional Index (Onodera 1984) — 20 composite scores total.
- `REPORT_VERSION` → 6.3.0-premium.

### Tests

- +4 unit tests (levers, mortality, METS-IR/FLI/PNI/AISI); 218 passing. Golden
  snapshot refreshed.

### Notes

- The longevity simulator is illustrative: it is a within-model counterfactual
  on the Levine PhenoAge equation, not a clinical outcome prediction.

---

## [6.2.0-premium] — 2026-07-25 — State-of-the-art blood-work: biological age + composite indices

Adds a research-grade layer of validated, literature-cited multi-marker indices
on top of the clinical engine — led by a mortality-calibrated biological-age
clock. Formulas were sourced from the primary papers (see citations in-report).

### Added

- **Biological Age (Levine PhenoAge, *Aging* 2018).** A 9-biomarker + age clock
  that predicts all-cause mortality better than chronological age, reported as
  biological age and age acceleration. **Inputs are SI-converted per the paper**
  (albumin g/L, creatinine µmol/L, glucose mmol/L, CRP mg/dL) — a step several
  popular online implementations skip, which silently corrupts the score.
  Validated against hand-computed vectors (healthy 40yo → ~32.7).
- **Cardiovascular:** Sampson–NIH LDL-C (current-generation calculation, better
  than Friedewald), Atherogenic Index of Plasma, Castelli Risk Index I.
- **Insulin resistance / metabolic:** TyG index, QUICKI, estimated average
  glucose (from HbA1c), and NCEP ATP III metabolic-syndrome criteria scoring.
- **Systemic inflammation:** SII, SIRI, and PLR composite blood-count indices.
- **Liver fibrosis:** FIB-4, APRI, and NAFLD Fibrosis Score.
- **Renal / acid–base:** albumin-corrected calcium and anion gap.
- New report section **"Advanced Risk & Aging Indices"** leading with the
  biological-age headline; every index carries its source-paper citation.
- `REPORT_VERSION` → 6.2.0-premium.

### Tests

- +8 unit tests validating each formula against hand-computed reference values
  (216 passing). Bloodwork golden snapshot refreshed to include the new
  ``advanced`` key.

---

## [6.1.0-premium] — 2026-07-24 — Comprehensive blood-work engine

Turns the bloodwork module from a genetics-only comparison into a full clinical
analysis. The original PheWAS genetic-prediction layer is preserved unchanged;
everything below is additive (attached under a new ``clinical`` key).

### Added

- **Clinical + optimal reference-range engine.** A ~50-biomarker catalogue
  (lipids, glycemic, inflammation, liver, kidney, thyroid, iron, CBC,
  electrolytes, hormones, vitamins, blood pressure) classifies each value
  against BOTH standard clinical ranges and tighter functional/optimal ranges
  (Optimal / Normal / Borderline / High / Low / Critical), with sex-specific
  ranges where relevant.
- **~12 calculated markers** derived from whatever inputs are present: non-HDL,
  Total:HDL, TG:HDL (insulin-resistance proxy), ApoB:ApoA1, remnant cholesterol,
  HOMA-IR, eGFR (CKD-EPI 2021), transferrin saturation, neutrophil:lymphocyte
  ratio, mean arterial pressure, BUN:creatinine, FIB-4, and a free-testosterone
  estimate.
- **Genotype-aware interpretation.** Flagged biomarkers pull the user's own
  variants for context — e.g. high ferritin + HFE C282Y → hemochromatosis flag;
  high LDL/ApoB + APOE ε4; high glucose + TCF7L2; high homocysteine + MTHFR;
  low vitamin D + GC/CYP2R1; B12 + FUT2 secretor status; uric acid + ABCG2;
  isolated high bilirubin + UGT1A1 (Gilbert's); Lp(a) as a genetic risk factor.
- **Per-system and overall health scores**, priority out-of-range flag list,
  and a visual reference/optimal range bar per marker in the report.
- `REPORT_VERSION` → 6.1.0-premium.

### Tests

- +6 unit tests for the clinical engine (status tiers, sex-specific ranges,
  derived markers, scoring/flags, HFE genotype note, clinical attachment).
  Bloodwork golden snapshot refreshed to include the new ``clinical`` key.

---

## [6.0.0-premium] — 2026-07-24 — Ancestry lineage cross-check + Detoxification panel

Two headline additions and one real bug fix. No synthesised data — every new
finding is computed from the sample's own genotype against curated,
registry-backed markers.

### Fixed

- **Autosomal ancestry mis-called European samples as South/East Asian.** Two
  compounding defects in `ancestry_pca.py`: (1) the `SLC45A2` AIM (rs16891982),
  a C/G *palindrome*, had its effect allele hard-coded to the wrong strand, so a
  European `GG` genotype scored as "zero European alleles"; (2) `_dosage` was not
  strand-aware, so markers reported on the opposite strand (e.g. `LCT` rs4988235
  `AA`) read as dosage 0. On a real T1a1a sample this flipped the call from
  European to South Asian. Dosage is now strand-aware, palindromic markers are
  displayed but excluded from scoring, and the AIM panel is expanded. The same
  sample now resolves European with a decisive evidence margin.

### Added

- **Uniparental (Y-DNA / mtDNA) geographic cross-check.** The ancestry section
  now maps the terminal Y-DNA and mtDNA haplogroups to geographic expectations
  over the five 1000G superpopulations and renders a concordant / plausible /
  discordant verdict against the autosomal call. A discordant deep-lineage check
  forces autosomal confidence to "low" — a small panel that contradicts the deep
  paternal/maternal lineage is treated as the suspect, not the haplogroup.
- **Detoxification & Environmental Resilience section (`detox.py`).** Phase I
  bioactivation (CYP1A1/1A2/1B1, AHR), Phase II conjugation (EPHX1, NAT2, NQO1,
  GSTM1/T1/P1), the NRF2 antioxidant axis (NFE2L2, SOD2, GPX1, CAT, HMOX1) and
  heavy-metal handling (ALAD, AS3MT, PON1, metallothioneins). Produces a
  composite wildfire-smoke resilience index (the "activate-but-don't-clear"
  Phase I/II mismatch) and a genotype-personalised protocol (sulforaphane/NRF2,
  glutathione/NAC, selenium, plus AQI/HEPA/N95 behavioural guidance framed for
  Great Lakes / Upper-Peninsula wildfire-smoke seasons). 14 new SNPs registered.
- **Metal Handling, Oxidative Defense & Neurodegeneration section.** Wires the
  previously-orphaned `metal_oxidative.py` (LRRK2/GBA, HFE, G6PD, ATP7B,
  metallothioneins, ZIP transporters, catalase) into the report — it computed
  findings but was never rendered.

### Tests

- 193 → 203 passing. New `tests/registry/test_detox_registry_consistency.py`
  (registry coverage + smoke-resilience scoring) and
  `tests/unit/test_ancestry_crosscheck.py` (strand-aware dosage, palindrome
  exclusion, European call, and lineage concordance/discordance).

---

## [Unreleased] — 2026-07-08 — Depth pass: surface computed-but-dropped data

No new genotype rules and no invented data. A systematic audit of all ~29
`build_*`/`render_*` functions found fields that analysis modules already
compute (and route to the renderer) but the renderer silently discarded. This
pass surfaces the high-value ones and fixes one live rendering bug. Every value
below is produced by an existing module from the sample's own genotype or is
author-curated in `snp_database.json` — nothing is synthesised.

### Fixed

- **Endurance training plans rendered blank weeks (live bug).** In
  `render_exercise_html`, the sport-plan week body used a fallback chain
  (`schedule → weekly_template → focus`) that matched none of the keys the
  marathon / road-cycling / triathlon plans actually use
  (`weekly_mileage_km`, `key_workouts`) — so every week of those plans printed
  its phase label followed by an empty body. The renderer now surfaces every
  descriptive field present per week, and additionally shows the powerlifting
  `accessory` prescription it had also been dropping.
- **Nutrition advanced layer dropped when no polygenic scores (live bug).**
  `_render_advanced_sections` early-returned `""` whenever `polygenic_scores`
  was falsy, silently discarding the cardiometabolic dashboard, inflammation,
  histamine, detox, glycemic, periodisation, food matrix, shopping list,
  recipes, **and** the entire protocols layer (30-day plan, glucose simulation,
  minerals, cycle-phase timing) — all computed and each already self-guarded.
  The PGS table is now the only block gated on `polygenic_scores`; every other
  section renders on its own data.

### Added — verified fields now shown

- **Curated variant catalogue** — each row gains an expandable "Recommendation
  & context" block: `recommendation` (all 612 entries; previously only the
  one-line `summary` showed), `chip_coverage_note` (27 entries — so a chip gap
  isn't read as "confirmed absent"), and `cross_references` (24 entries).
  `analyze.py` now threads `chip_coverage_note` into each tier-1 record.
- **Carrier status** — the affected/carrier tables now show the variant
  `rsid` (dbSNP-linked), the `pathogenic_allele` that defines carrier status,
  the ethnicity-specific `carrier_frequency` (9 entries, e.g. "European ~1:25"),
  and any indel `chip_caveat` — all previously computed and dropped.
- **Polygenic risk scores** — each panel card gains a "Contributing variants"
  table (`result.used[]`): rsID, gene, effect allele, the sample's genotype,
  effect-allele copies, per-allele odds ratio (exp of the published log-OR),
  and effect-allele frequency — plus a note listing panel variants not typed
  on this chip (`result.missing[]`). Previously collapsed to a bare count.
- **PGx dose simulation** — per-gene contributions (`clearance`,
  `dose_factor`, `ae_rr`) are shown inline per gene, so for multi-gene drugs
  (e.g. warfarin VKORC1×CYP2C9) the combined estimate is auditable rather than
  an opaque product.
- **Health economics** — the findings table gains a QALY-gain column with
  per-finding cost-per-QALY (the standard cost-effectiveness metric), replacing
  reliance on the single payer-level aggregate.
- **PheWAS** — each biomarker trait gains an expandable "Evidence & contributing
  variants" block: the GWAS reference, the standardised Z-score against the
  reference mean/SD, the driving variants (rsID / effect allele / your copies /
  per-allele β), and the untyped panel SNPs. `phewas._score_trait` now returns
  the `used_variants` / `missing_variants` it already computed (additive
  contract change; the enriched entries also carry effect allele + AF).
- **Ancestry** — surfaces the `evidence_margin_nats` (which the module calls
  "the honest measure of confidence" — the best-vs-runner-up log-likelihood
  gap) and an expandable table of the ancestry-informative markers used, with
  LD-redundant markers flagged so they aren't read as independent evidence.
- **Expanded PGS (PGS Catalog)** — per-panel low-r² imputed-variant count and
  not-covered count added to the coverage line; section-level Elevated/High
  headline summary added (previously only the curated-PRS renderer had one).
- **Genetic longevity** — each panel (longevity / telomere / skin) now shows
  its Z-score and variant coverage; especially important for the telomere proxy
  which can rest on ≤2 variants.
- **Bloodwork** — each predicted row shows the genetic coverage (% callable +
  SNP count) behind the prediction, i.e. how much genotype the estimate rests on.
- **HLA** — untested alleles now name the specific tag SNP(s) missing from the
  chip instead of an opaque "Untested".
- **References** — the curated ClinVar clinical-significance assertion is shown
  per reference (previously dropped).
- **QC** — no-call count and average per-domain callability added to the data-
  quality header.
- **Local ancestry** — each discordant segment now shows its per-superpopulation
  log-likelihoods, the confidence gap, and the number of AIMs in the window.

### Added — top prescribed-drugs pharmacogenomic screen

- **`top_prescribed_drugs.json` + `top_drugs_screen.py` + a new report section.**
  A curated menu of the ~630 most commonly prescribed U.S. medications (generic
  names + class + brand) is screened against the genome. For each drug the
  pharmacogenomic relevance is computed entirely from **real bundled data** —
  `drug_database.json` (CPIC/DPWG drug↔gene↔marker↔dosing) and
  `cpic_data/drugs.tsv` (ClinPGx/PharmGKB per-drug CPIC-pair & clinical-
  annotation levels) — cross-referenced with the user's per-gene metabolizer
  phenotypes from `pgx.py`. Drugs are tiered: (1) **genotype-actionable** —
  the drug's gene is one where the user is a non-normal metabolizer (open);
  (2) relevant gene typed & normal; (3) PGx-relevant but the gene is
  unresolved on this chip (→ `--impute` would resolve many); (4) no known PGx
  interaction in the bundled databases. Tiers 2–4 are collapsed.
  **Accuracy contract:** the prescribed-drugs file supplies *only* names and
  classes; every gene link and evidence level is from CPIC/PharmGKB, and drugs
  with no database entry are reported as "no known interaction," never as a
  fabricated finding. The section is explicitly not a dosing instruction — the
  star-allele `pgx.py` section remains authoritative. Wired into the pipeline
  behind a graceful `try/except ImportError`.

### Added — new data & provenance

- **ClinPGx / PharmGKB clinical-variant annotations** (`pharmgkb_clinical.py` +
  a new report section). The downloaded `cpic_data/clinicalVariants.tsv`
  (~4,300 rsID-keyed rows, 591 drugs) was previously read by no module. Every
  variant on the chip is now cross-referenced against it; positions the user
  carries are reported with gene, genotype, drug(s), evidence level, and
  phenotype. **Accuracy guardrails (deliberate):** (1) evidence tiers are
  labelled as *ClinPGx/PharmGKB clinical-annotation levels* (per the bundled
  README), explicitly **not** CPIC guideline strength; (2) high-evidence
  (Level 1A/1B/2A/2B) positions are shown openly, the weak/unreplicated
  Level 3/4 long tail behind a collapsed, labelled disclosure; (3) framed as
  "a genotype at an annotated position," **not** a direction-of-effect call,
  because the source table carries no risk allele — the dedicated `pgx.py`
  star-allele section remains authoritative for the actionable genes (flagged
  ↗ PGx). Wired into the pipeline behind a graceful `try/except ImportError`.
- **Imputation provenance now flows to results.** `analyze.tier1_lookup` reads
  the `source` (`chip`/`imputed`) and `r2` columns that `--impute` adds to
  `snps_df`, records them on each variant, and the catalogue flags imputed
  calls with an `imputed r²=…` badge. Audit of the impute path confirmed the
  merged chip+imputed frame is reassigned *before* tier1 matching and every
  analysis (so imputed variants already flow into all modules) and that no
  module filters them out; this closes the honesty gap so a large `--impute`
  run never presents a statistical call as a measured one.

### Notes

- Purely additive: no rows, fields, or sections were removed; existing columns
  are unchanged. Verified by the full suite (**193 passing**, +15 new tests:
  catalogue context block + imputed-provenance threading, exercise & nutrition
  regression tests, the PheWAS producer contract, the PharmGKB module, and the
  top-prescribed-drugs screen). Each renderer change was additionally exercised
  by calling the
  build function directly with a producer dict. Golden snapshots unaffected
  (they cover the V6 personalisation dicts, not the report HTML).
- The PRS "contributing variants" panel's note on untyped variants was
  corrected after review: untyped panel variants are **excluded** from the
  score (both raw score and expected mean/variance use typed variants only) —
  they are not treated as 0-copy reference, so there is no directional
  "downward bias"; lower coverage simply means a less complete estimate.
- **Registry mass-migration (58 → 612 records) was scoped and declined.**
  `snp_database.json` carries no GRCh37/38 coordinates or ancestral/derived
  alleles, so populating the coordinate-based registry from it would require
  unverifiable, fabricated coordinates. Deferred until a verified coordinate
  source (dbSNP/ClinVar dump) can be added offline.
- **Remaining audited opportunities (still deferred, lower value):** nutrition
  `polygenic_scores[].typed` per-locus genotypes + cycle-phase dosage lines;
  Y-DNA per-node migration narratives + structured `contradictions`; mtDNA
  per-marker derived/ancestral status; traits tested/not-tested counts;
  medications catalogued/uncatalogued drug counts. All follow the same
  fabrication-free "surface a computed-but-dropped field" pattern.

---

## [Unreleased] — 2026-06-16 — Y-DNA correctness & registry completion

### Fixed

- **Y-DNA haplogroup analysis reported the wrong haplogroup.** The decision
  tree in `y_haplogroup.py` had corrupted marker data — duplicate/mis-mapped
  rsIDs (e.g. `rs2032658`, actually M207/R, was attached to the O node), wrong
  positions and flipped ancestral/derived alleles, some copied from the
  synthetic test genome. Combined with first-derived-child-wins descent it
  produced confident wrong calls (e.g. terminal O for a haplogroup-T sample
  with only 1/3 path markers confirmed). The backbone was rebuilt from the
  ISOGG Y-SNP index (verified rsID + GRCh37 position + ancestral→derived),
  rooted at CT so non-K lineages (I, J, E, G) are representable. Lookup now
  votes across several co-defining markers per node, skips strand-ambiguous
  (A/T, C/G) SNPs, reports only the deepest **confirmed** node (no
  over-calling), still descends through chip gaps toward genuinely confirmed
  downstream markers, and flags biologically-impossible contradictory calls.

### Added

- **`snp_registry`: registered the gut-health and metals/oxidative panel
  SNPs** that `gut_health.py` and `metal_oxidative.py` referenced but the
  registry did not know about (`rs4988235`, `rs10156191`, `rs2066844`,
  `rs11209026`, `rs34637584`, `rs1001179`, `rs1061472`, `rs8052394`,
  `rs28366003`), each with verified GRCh37/GRCh38 coordinates and + strand
  ancestral/derived alleles.
- **`rs1061472` (ATP7B K832R) is a documented strand flip.** carrier.py labels
  the coding-strand `G`; the registry stores canonical + strand `T/C`
  (complement(C) == G), so it reconciles as a strand flip, not a biology
  disagreement. Added to the locked strand-flip set in the carrier
  consistency test.

---

## [Unreleased] — 2026-05-29 — Accuracy & honesty pass

No new features — existing outputs made correct and honest about uncertainty.

### Fixed

- **Ancestry AIM heuristic over-weighted AMR for Southern Europeans.** Two
  causes addressed in `ancestry_pca.py`:
  - LD-correlated AIMs (the HERC2 eye-colour pair, the EDAR pair, and the
    perfectly-correlated LCT/MCM6 pair) were treated as independent and
    double-/triple-counted the same signal. They are now grouped by `ld_block`
    and counted once.
  - Confidence is now driven by the **log-likelihood margin** between the best
    and runner-up population (not the softmax probability, which exaggerates
    small gaps on a tiny panel). Ambiguous calls (margin < 2.3 nats ≈ 10:1)
    are flagged and downgraded to *low*; the heuristic is capped at *moderate*.
    In a sweep of plausible Southern-European genotypes this downgrades ~91% of
    spurious AMR-primary calls to low confidence. The five-way numbers are now
    labelled *relative affinity*, not admixture proportions.

### Changed — explicit confidence + coverage on every score

Every score now reports an explicit confidence level and SNPs-covered-vs-expected,
and low-coverage results are downgraded or suppressed:

- **PGx** (`pgx.py`): genes with **0 defining variants typed** are now reported
  as *Indeterminate* instead of a confident "Normal Metabolizer", and an
  untyped HLA-B*57:01 tag is *Indeterminate* rather than a false "Negative".
  Per-gene confidence from defining-variant coverage.
- **Curated PRS** (`prs.py`): explicit confidence from callability; scores below
  3 typed variants suppressed; low-confidence scores excluded from headline
  findings.
- **Expanded PGS** (`pgs_catalog.py`): explicit confidence from coverage (a
  percentile from <40% of a scoring file is now flagged *low*).
- **mtDNA / Y haplogroups** (`mt_haplogroup.py`, `y_haplogroup.py`): explicit
  confidence and matched-vs-expected marker counts.

### Disclaimers

- Strengthened the global banner lead to "Informational and educational use
  only — not a clinical diagnostic" and added a "Reading the results" note on
  confidence/coverage labels (`renderers.py`).
- Added informational-use disclaimers to the Expanded-PGS and Y-DNA sections;
  relabelled the ancestry bars as relative affinity.

---

## [8.0.0] — 2026-05-17 — Decomposition + registry migration

V8 is a **foundation release**: zero new user-visible features, substantial
structural and correctness work under the hood. The pitch the V7 plan made
("the project is now contributable") is paid off here.

### Shipped

#### 1. analyze.py decomposition is complete

The monolithic `analyze.py` is split into focused modules:

| File | Lines | Role |
|---|---|---|
| `cli.py` | 196 | argparse parser + console-script entry point |
| `pipeline.py` | 818 | `run_pipeline(args)` — orchestration body |
| `renderers.py` | 3 097 | every `build_*_html` function + helpers |
| `analyze.py` | 1 055 | library functions (parse / tier1 / tier2) + lazy back-compat shim |

`analyze.py` dropped from **4 876 lines (pre-V7)** → **1 758 lines (post-V7)** →
**1 055 lines (post-V8)**: a 78 % reduction over the two foundation releases.

Backwards compatibility: `from analyze import build_html_report` still works
via a lazy `__getattr__` shim. The shim is *lazy* on purpose — an eager
re-export caused a hard circular-import deadlock when `analyze.py` is run as
a script (renderers needs `analyze.CATEGORY_ORDER` at module-load time,
which doesn't exist yet during the second-load of analyze.py as the
`analyze` module). Don't reinstate the eager shim.

#### 2. Registry migration — three modules

| Module | Status | rsIDs registered | Notes |
|---|---|---|---|
| `supplements.py` | ✅ full (V7 PoC) | all 16 rule SNPs | proof-of-concept; behaviour byte-preserved |
| `carrier.py` | ✅ partial | 13 direct + 3 strand-flipped (16 of 49) | 33 deferred (mostly indels) |
| `wellness.py` | ✅ partial | 25 of 34 | 9 obscure variants deferred |
| `traits.py` | ✅ partial | 16 of 48 | 32 phenotype SNPs deferred |

Every migrated module now ships an `audit_against_registry()` function that
runs at import time and **fails fast** if a rsID disagrees with the registry
on real biology (not strand). Strand flips (e.g. SERPINA1 PiZ, PAH R408W,
APOB R3500Q — historically reported on gene-coding strand vs registry's
+ strand convention) are explicitly tolerated and locked by tests.

Cross-module consistency tests live in `tests/registry/`:

- `test_carrier_registry_consistency.py` — 4 tests
- `test_wellness_registry_consistency.py` — 3 tests
- `test_traits_registry_consistency.py` — 3 tests

#### 3. Reconciliation log

**carrier.py — strand-flipped variants (kept on gene-coding strand for
clinical recognisability):**

- `rs28929474` SERPINA1 Z allele — carrier-side "A" is the gene-coding
  strand allele (Glu342→Lys). dbSNP reports the same biological change as
  "C→T" on the GRCh37 + strand of chr14 (SERPINA1 is on the − strand).
  No biology disagreement; registered as C→T on + strand.
- `rs5030858` PAH R408W — same situation. Carrier "A" = coding-strand
  variant; registry C→T = + strand of chr12 (PAH on − strand).
- `rs5742904` APOB R3500Q — same situation. Carrier "A" = coding-strand
  variant; registry C→T = + strand of chr2 (APOB on − strand).

These three were the only disagreements surfaced by the migration audit.
All resolved as strand conventions, not biology. The audit now treats
"local == complement(registry)" as `strand_flipped` rather than
`disagreed`, and tests lock that set against accidental future changes.

#### 4. Imputation provenance schema + reference-build auto-detection

New `provenance.py` module wires two annotations onto every parsed
DataFrame:

- **`snps_df["source"]`** — one of `"chip"`, `"imp_high_r2"`,
  `"imp_low_r2"`. `parse_dna_file` tags every row `"chip"`. The
  imputation hook (`tag_imputed_rows`) lets the imputation module
  overwrite rows with their Beagle DR2-based confidence tag. **Schema
  only — downstream consumers (PRS, PheWAS) do not yet read this
  column.** See "Deferred to V8.1" below.
- **`snps_df.attrs["build"]`** — `"grch37"`, `"grch38"`, `"mixed"`, or
  `"unknown"`. Auto-detected by sampling registry SNPs whose positions
  are known in both builds. The TellmeGen-style "mixed coordinates"
  case is now an explicit return value and a warning rather than
  silent breakage.

14 new tests in `tests/unit/test_provenance.py` cover both schema
operations + all five build-detection outcomes.

#### 5. Test infrastructure tightened

- **145 tests passing** in ~70 ms (was 121 at end of V7).
- 6 golden snapshots locked the V6 outputs byte-exactly through three
  major refactors (renderers extraction, pipeline extraction, registry
  migration). None drifted — the snapshots paid for themselves.
- New `tests/registry/` directory holds the cross-module consistency
  tests.

### Deferred to V8.1

**Confidence intervals on PRS and PheWAS.** Investigated and deferred:
`prs.py` stores per-variant `log_or` and `af` but **not** standard errors
or variance. Computing CIs without SE data would mean inventing them —
the exact "fake rigor" failure mode the advisor flagged. V8.1 should:

1. Extend `prs.py` panel entries with a `log_or_se` field, sourced from
   the cited GWAS summary statistics.
2. Add a `confidence_interval(score, vars_used)` helper that propagates
   per-variant SE through the score's variance.
3. Surface 95 % CIs alongside every percentile in the report (PRS,
   PheWAS, MR).

**Ancestry-stratified PRS.** Same problem: the codebase has only EUR
effect sizes. Per-population scaling would require either fabricating
scaling factors (rejected) or downloading per-population PGS Catalog
files and re-deriving (hours of data work). V8.1 should:

1. Use `ancestry_pca` output to infer a primary ancestry call
   (EUR / AFR / EAS / SAS / AMR).
2. Apply the cross-population PRS attenuation literature (typical
   ~30-50 % R² reduction in AFR vs EUR for EUR-derived scores) as
   either a global down-weight or per-panel calibrated factors.
3. Surface a "PRS calibration: derived for European populations,
   your inferred ancestry is X — interpret with caution" disclaimer
   on every non-EUR sample.

**Imputation provenance flowing downstream.** Schema landed; scoring
modules still ignore the `source` column. V8.1 should:

1. Update `prs._score_panel` to down-weight `imp_low_r2` variants
   (suggested: 0.5× weight) and skip them entirely if `--strict`.
2. Same for `phewas._score_trait`.
3. Surface per-panel provenance breakdown in the report cards
   ("LDL PRS: 18/22 variants typed directly, 3 imputed high-R², 1 low-R²").

**carrier.py indel schema extension.** 33 of carrier.py's 49 variants
are indels (CFTR ΔF508, HEXA 1278insTATC, BRCA1 185delAG, …). The
`SNPRecord` dataclass currently requires `ancestral` and `derived` to
be single-base A/T/C/G. V8.1 should:

1. Add an `IndelRecord` variant of `SNPRecord` (or extend `SNPRecord`
   with `kind: Literal["snv", "indel"]` and `inserted_seq` /
   `deleted_seq` fields).
2. Update `carrier.audit_against_registry` to handle indels.

**Remaining module migrations.** `pgx.py`, `phewas.py`, `prs.py`, `hla.py`
not yet migrated. These hold per-variant **effect sizes** alongside
metadata — they need the same registry schema extension that CIs need
(`log_or_se`, effect-size CI fields). One coherent V8.1 work-stream.

**`tests/registry/test_no_duplicate_rsids.py`.** The migration playbook
mentioned this — a meta-test that greps the codebase and asserts no
rsID metadata dict outside `snp_registry.py`. Worth adding once the
remaining modules are migrated (otherwise it would block them).

### Migration commands

```bash
git diff v7.0.0..v8.0.0 --stat | tail -20  # see what changed
pytest tests/                              # 145 tests, ~70 ms
pytest --snapshot-update tests/golden/     # only with deliberate review
```

### Reconciliation against published literature

Sources cited per registry entry in `snp_registry.py`. Notable sources for
the V8 additions:

- HFE C282Y, H63D — Feder 1996 (Nat Genet)
- F5 Leiden — Bertina 1994 (Nature); ClinVar 642
- F2 G20210A — Poort 1996; ClinVar 13310
- TREM2 R47H — Guerreiro 2013, Jonsson 2013 (NEJM)
- HBB sickle-cell — ClinVar 15333; Pauling 1949
- APOE ε4/ε2 — Corder 1993 (Science)
- SOD2 Ala16Val — Sutton 2003 (Pharmacogenetics)
- COL5A1 — Mokone 2006
- DIO2 Thr92Ala — Panicker 2009
- FOXO3 longevity — Willcox 2008 (PNAS)
- GPX1 Pro198Leu — Forsberg 2000
- TERT telomere length — Codd 2013

---

## [7.0.0] — 2026-05-17 — Consolidation

Foundation release: zero user-visible features, large durability gain.

- `snp_registry.py` — unified SNP metadata (single source of truth)
- `tests/` — first pytest suite, 121 tests, golden-snapshot regression
- `requirements.txt`, `pyproject.toml`, `.github/workflows/ci.yml`,
  `.pre-commit-config.yaml` — first time the project has any of these
- `cli.py` — CLI parser extracted from `analyze.py`
- `docs/MIGRATION_PLAYBOOK.md` — step-by-step migration recipe
- `supplements.py` migrated to the registry as a proof-of-concept

See `README.md` "What's new in V7" for the full V7 changelog.

---

## [6.0.0] — V6 Personalisation

- `bloodwork.py` — `--bloodwork labs.json` for predicted-vs-actual comparison
- `supplements.py`, `exercise.py`, `nutrition.py` — V6 personalisation
- `fhir_export.py` — `--fhir` for HL7 FHIR R4 clinical export
- `personalized_plan.py` — master dashboard synthesising the V6 outputs

---

## [5.0.0] and earlier

Pre-changelog. See git history for V1-V5 evolution. Major milestones:

- V1: curated SNP catalogue + APOE genotype
- V2: Y-DNA + mtDNA haplogroup analysis
- V3: imputation + medications + carrier-report + chat + compare
- V4: wellness predictions
- V5: HLA / ROH / local ancestry / PheWAS / MR / genetic age / PGx
      simulation / reproductive simulator / emergency card / narrative
