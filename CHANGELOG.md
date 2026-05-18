# Changelog

All notable changes to this project are documented here. Format inspired by
[Keep a Changelog](https://keepachangelog.com/); versioning follows the
`{major}.{minor}.{tag}` scheme.

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
