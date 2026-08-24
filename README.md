# GenomeLens

**An applied health-economics engine that answers one question: what is knowing your
genome actually worth, in decisions?** Cost-utility analysis, value of information and
budget impact — computed on real genomic data, entirely on your own machine.

### 📄 [See 26 pages of real output → `docs/samples/econ-output-sample.pdf`](docs/samples/econ-output-sample.pdf)

Faster than reading this page. **Synthetic input — no human genome and no personal health
data** were used to make it. Reproduce it with `python scripts/make_econ_sample.py`.

In it: cost and QALYs reported **separately** with the ICER *withheld* under dominance ·
a **double-counting correction** stating what the naive figure claimed and how much came
out · an **adherence discount** charged to benefit *and* ongoing cost · three statistical
corrections shown **before → after** · a footnote saying which figures rest on the sourced
registry and which do not.

![status](https://img.shields.io/badge/status-active-brightgreen)
![python](https://img.shields.io/badge/python-3.10%2B-blue)
![privacy](https://img.shields.io/badge/privacy-100%25%20local%20%C2%B7%20offline-purple)
[![CI](https://github.com/conallaque/genomelens/actions/workflows/ci.yml/badge.svg)](https://github.com/conallaque/genomelens/actions/workflows/ci.yml)
![tests](https://img.shields.io/badge/tests-780%20passing-brightgreen)
![input](https://img.shields.io/badge/input-chip%20%2B%20whole--genome%20VCF-blue)
![license](https://img.shields.io/badge/license-All%20Rights%20Reserved-red)
[![Buy Me a Coffee](https://img.shields.io/badge/buy%20me%20a%20coffee-support-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/caque)

### 🔁 [Independent check → `heor-model-replication`](https://github.com/conallaque/heor-model-replication)

Everything on this page is my own model assessing itself. That repo is not: it reproduces
**three peer-reviewed cohort state-transition models in Python**, matching every printed
cost, effect, ICER and dominance verdict exactly. If you want to know whether I can build
a cSTM that agrees with a published one before trusting the numbers here, start there.

**Methods and every equation:** [`docs/METHODS.md`](docs/METHODS.md) ·
**How to run it:** [`docs/USAGE.md`](docs/USAGE.md)

> **Not medical advice — educational and research use.** An illustrative
> decision-analytic model, not a formal economic evaluation.

---

## What this demonstrates as HEOR work

Anyone can run a cost-effectiveness analysis. These are the parts that are harder to
fake, and each is traceable to a named test.

| Competency | Evidence |
|---|---|
| **Cost–utility analysis done properly** | Cost, QALYs, ICER and INMB reported separately, never blended into one "value" figure. ICER suppressed in the dominance quadrants, because a negative ratio is ambiguous. `econ.engine.CEAResult` |
| **Finding my own errors** | Eight of twenty-one finding sources routed onto one cardiometabolic anchor and were **summed** — a 79% risk reduction, which is not a probability. Fixed by pooling on the risk scale; the report shows the size of its own correction. `test_stacked_findings_do_not_sum_their_risk_reductions` |
| **Uncertainty that is real** | An earlier version reported a strategy cost-saving in **100% of simulations** — the finding-level parameters were pinned outside the sampling loop. `test_psa_without_rebuild_understates_uncertainty` |
| **Parameter provenance, enforced** | Every figure carries a tier. `tier="assumption"` **may not** cite a source — the registry fails to load if it does. 47% of ~350 figures resolve to a PMID or DOI; the model reports its own coverage instead of claiming "sourced". `econ/params.py` |
| **Knowing what not to monetise** | Reproductive outcomes are never priced — attaching a figure to an affected birth prices a prospective child. Stated in code, enforced by a test, surfaced as a decision rather than an omission. `NOT_VALUED` |
| **Structural modelling** | Cohort state-transition model against US life-table mortality, Simpson's 1/3 within-cycle correction cross-checked against an independent implementation. `test_within_cycle_weights_match_the_published_implementation` |
| **Validated against published models** | Three peer-reviewed cohort state-transition models reproduced in Python, every printed cost, effect, ICER and dominance verdict matched exactly — the one claim here that is not self-assessed. [`heor-model-replication`](https://github.com/conallaque/heor-model-replication) |

The recurring theme: the model reports an unflattering answer as readily as a flattering
one, and several commits above exist because it did.

**Authorship, plainly.** The health-economics modelling, scientific decisions and product
direction are mine. The **software implementation was largely AI-generated** under my
direction and review — what is on offer here is the economics and the judgement.

---

## Why I built this

I've worked in health economics for five years — a BA focused on it, with a published
paper, and now a master's at Northeastern. One question keeps coming back: **what is
understanding your own health actually worth?**

Your genome is the best input to that question. Two things block it.

**Price.** Buying these analyses and the interpretation separately runs from a few
hundred to a few thousand dollars.

**Privacy.** To get them, you hand your genome to a company's cloud — the one piece of
data you can't change or take back — where it can be breached, sold, or reused.

That's a health-economics problem, not a technical one. The value of information in a
genome is real but individual: this model puts it anywhere from near zero to tens of
thousands, depending on what's actually in the file. It reports which one you got instead
of assuming the better answer. But access is gated by money and by a privacy cost most
people shouldn't accept, so the return genomics promises only reaches people who can pay
and are willing to give themselves up.

GenomeLens removes both. It runs locally, for free, on a laptop you already own. And it
doesn't just list findings — it prices them: what each one is worth, what *acting* on it
is worth, and whether acting is worth it at all.

I wasn't trying to build a genomics toy. The payoff of understanding your own biology
shouldn't cost a lot of money or your DNA. On a normal laptop, it doesn't have to.

---

## What else it tells you

The economics is the headline, but it is computed on top of a full analysis. Briefly,
what a run actually reports:

| Area | What you get |
|---|---|
| **Pharmacogenomics** | CPIC/DPWG star-allele calling across 217 drugs — which medications to avoid, dose-adjust, or expect a poor response from. HLA typing for drug-hypersensitivity risk. Drug–drug and gene–drug interactions. |
| **Clinical risk** | ClinVar pathogenic/likely-pathogenic screen, ACMG actionable genes, carrier status, and an offline computational screen (AlphaMissense) for rare variants beyond ClinVar. |
| **Polygenic risk** | PRS and PGS Catalog scores with coverage and an ancestry-transferability caveat, plus PheWAS biomarker predictions. |
| **Gut health** | Lactase persistence, FUT2 secretor status, coeliac HLA-DQ2.5/DQ8 tags, histamine/DAO clearance, NOD2 and IL23R. |
| **Fitness** | Power vs endurance lean, strength trainability, recovery speed, VO2max response, plus lift-level protocol cards and HR/FTP zones. |
| **Nutrition** | Macronutrient pressures, caffeine and alcohol handling, lactose and gluten, folate/B12 methylation, and a 30-day meal plan with a shopping list. |
| **Supplements** | A ranked stack with the genotype behind each entry, an evidence tier, what to avoid, and a monthly cost. |
| **Blood work** | Genotype-aware interpretation of a supplied lab panel — clinical *and* optimal ranges, biological age (Levine PhenoAge), AHA PREVENT 10-year ASCVD risk, and longitudinal tracking. |
| **Ancestry & lineage** | AIM-based population inference, Neanderthal and ancient-population affinity, chromosome painting, and Y-DNA / mtDNA lineage chains showing every branch point and where the chip runs out. |
| **Other panels** | Detoxification enzymes, metals and oxidative handling, immunogenetics, neurochemistry (COMT/MAOA/BDNF), urologic, addiction genetics, blood type, sleep and chronotype. |
| **Recommendations** | A master plan tying it together: decade-by-decade priorities, behavioural protocols, a printable emergency card for clinicians, and a plain-language summary of what to actually do first. |

Findings that have a **real signal but no proven intervention** carry no dollar value and
say so — a risk you cannot act on is reported as exactly that, not priced.

## Health economics — what the engine actually does

Full derivations, equations and citations: [`docs/METHODS.md`](docs/METHODS.md).
26 pages of output: [`docs/samples/econ-output-sample.pdf`](docs/samples/econ-output-sample.pdf).

**Cost-utility analysis.** Cost, QALYs, ICER and INMB reported separately. 3% discounting
on both costs and QALYs (0/3/5% in sensitivity). The ICER is **withheld** in the dominance
quadrants rather than reported as a negative number.

**Findings pooled, not summed.** Twenty-one finding sources feed the model and eight route
onto the same cardiometabolic anchor. Summing them claimed a 79% risk reduction — not a
probability. Pooling is complement-of-products on the risk scale with a correlated-signal
penalty, and each condition's cost of illness is charged once.

**Trial efficacy discounted to real-world effectiveness.** Adherence enters inside the
pooling product, keyed on what acting asks of the person rather than which organ is
involved. It is charged to the benefit *and* the ongoing intervention cost — but not to
the one-off test cost, which is why it moves cost-effectiveness at all.

**Uncertainty.** Probabilistic sensitivity analysis rebuilding the pools inside every
draw, a CEAC, a one-way tornado, EVPI, EVPPI and breakeven. An earlier version reported a
strategy cost-saving in 100% of simulations because the finding-level parameters were
pinned outside the loop.

**Cohort model.** Three-state state-transition model against US life-table mortality, with
`p = 1 − e^(−rΔt)` (never `r·Δt`), Simpson's 1/3 within-cycle correction cross-checked
against an independent implementation, and cohort-conservation validation.

**Budget impact**, on ISPOR conventions that deliberately differ from CEA: short horizon,
population-scaled, undiscounted, uptake phased in, reported PMPM.

**Provenance.** Every parameter carries a tier. `tier="assumption"` **may not** cite a
source — the registry fails to load if it does, which makes anti-assumption-laundering a
machine-checked rule rather than a promise. The report states which figures rest on the
registry and which rest on curated tables that do not.

**Reporting.** Second Panel dual perspective with an impact inventory, and a CHEERS 2022
checklist that names the items *not* addressed as well as those that are.

## Engineering notes

- **Unified, strand-aware SNP registry** — one source of truth for GRCh37/38 coordinates and ancestral/derived alleles; caught and fixed palindrome/strand bugs that silently mis-call ancestry.
- **KING-robust relationship inference** — a proper kinship estimator (not naive percent-identity), IBS0-refined for parent-child vs full-sibling.
- **No fabricated figures** — no invented polygenic percentiles; transmission ≠ disease penetrance; ClinVar review-star confidence; Phase-3 findings are labelled computational predictions, never clinical calls; and a **grounding guardrail rejects any AI-introduced figure absent from the deterministic data**, so the local LLM cannot invent a risk or a statistic.
- **Runs cleanly end-to-end on a public genome** — a *functional* end-to-end test (not a clinical accuracy validation): the full **GIAB HG001 (NA12878)** reference genome runs through build detection, Phase-2 ClinVar, the Phase-3 predictor screen, and the health-ROI engine with **no runtime errors** — 3 ClinVar pathogenic (incl. a carrier), 141 predicted-damaging rare variants, and a modelled ROI reported with a confidence interval. (Functional check — outputs are computational estimates, not accuracy-validated against a clinical truth set.)
- **780-test suite**, reference-build auto-detection (GRCh37/38 incl. rsID-less whole-genome VCFs), and graceful degradation when optional data or models are absent.

## Output files

One run writes a self-contained set of HTML pages beside the output path — no
server, no build step, no network.

| File | What is in it |
|---|---|
| `report.html` | The main report. Every module's section, including the pooled cost-effectiveness analysis, the plain-language summary, and the Y-DNA / mtDNA lineage chains. |
| `economic_analysis.html` | The individual's economic sheet — cost avoided, QALYs, net monetary benefit, and the cost-saving / cost-effective split kept explicitly apart. |
| `emergency_card.html` | One page, actionable findings only: drug hypersensitivities, clotting disorders, pharmacogenomic extremes. Designed to be printed and carried. |
| `supplements.html` | Ranked supplement stack with the genotype behind each entry and a tier for how well evidenced it is. |
| `nutrition.html` | Macronutrient pressures, food-specific guidance, and a 30-day meal plan. |
| `exercise.html` | Power/endurance bias, trainability tiers, and lift-level protocol cards. |
| `longevity.html` | Longevity composite and a quarterly year-long plan. |
| `personalized_plan.html` | The master dashboard tying the other pages together. |

## Privacy

- **No network calls** during analysis except optional Ollama (`localhost:11434`).
- DNA files never leave your machine.
- `tier1_results.json` and HTML reports contain genotypes and should be
  treated as PHI — the bundled `.gitignore` prevents them being committed.

---

## Support

GenomeLens is built and maintained by one person, in the open. If it saved you a
consult, taught you something about your own biology, or you just appreciate the
engineering, you can help fuel the next module:

[![Buy Me a Coffee](https://img.shields.io/badge/Buy%20Me%20a%20Coffee-caque-FFDD00?style=for-the-badge&logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/caque)

---

## Citation

If this tool informs research or teaching material, please cite it as:

```bibtex
@software{dna_analysis_tool,
  title  = {GenomeLens — a local, privacy-first health-economics engine for consumer genomics},
  author = {Aque, Conall R.},
  year   = {2026},
  url    = {https://github.com/conallaque/genomelens},
  note   = {Local chip/whole-genome analysis with HL7 FHIR R4 export (educational; not a certified clinical record)}
}
```

Underlying datasets and guidelines should be cited independently — CPIC for
pharmacogenomic recommendations, ClinVar for variant interpretation, PGS
Catalog for polygenic scores, ISOGG / YFull for Y-DNA phylogeny, the relevant
GWAS consortia (UK Biobank, GIANT, GLGC, MAGIC, PGC, Astle 2016, Yengo 2022,
…) for biomarker effect sizes, and the HL7 Clinical Genomics IG for the FHIR
output format.

---

## License

**Proprietary — © 2026 Conall Aque. All Rights Reserved.** This repository is
public for viewing and reference only. Use, copying, modification, or
distribution without the author's written permission is prohibited. No warranty.

---

## Acknowledgements

Built on top of the work of countless GWAS consortia (PGC, UK Biobank, GIANT,
GLGC, MAGIC, Astle 2016, …), CPIC, ClinVar, PharmGKB, ISOGG, YFull, and the
PGS Catalog. Reference genome positions are GRCh37/hg19 unless noted.
