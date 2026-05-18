# SNP Registry Migration Playbook

The V7 consolidation introduces `snp_registry.py` as the single source of
truth for rsID metadata. **Seven** legacy modules currently maintain their
own parallel registries (146 of 318 rsIDs in the project — 46 % — appear in
2+ files, often with contradictory metadata). This document is the playbook
for migrating each module, in priority order.

## Why this matters

When the same rsID lives in two files with different metadata, **silent
correctness bugs accumulate**:

- One module reads `rs1801133` on the + strand, another on the − strand →
  one of them is wrong about the carrier's MTHFR status.
- Position fields drift apart over time as one is updated and another isn't
  → position-based fallback lookups (used when a chip reports a SNP under a
  non-canonical rsID label) return inconsistent answers across modules.
- Adding a new rsID requires touching every consumer module, and forgetting
  one means the new variant is invisible to some pipelines.

`snp_registry.SNPRecord` is `@dataclass(frozen=True, slots=True)` — fast,
hashable, immutable, validates `ancestral != derived` and `letter ∈ {A,T,C,G}`
at construction time. Every entry carries dual-build coordinates (GRCh37 +
GRCh38), a literature `source`, and `last_verified` for staleness audit.

## Proof-of-concept migration (done)

`supplements.py` was migrated first because it had the cleanest entry point
(the strand-aware `_risk_dose` helper was the most isolated). The migration
preserved byte-exact behaviour — verified by the V6 golden snapshots in
`tests/snapshots/`. Total diff: one import, ~40 lines deleted from the
strand-aware helper (now a thin shim delegating to the registry).

## Migration order

Pick the next module by **rsID overlap with the registry seed × test
coverage**. Rough scoring:

| Module        | rsIDs | Existing tests | Risk  | Notes |
|---------------|-------|----------------|-------|-------|
| `carrier.py`  | ~12   | golden only    | low   | Small, well-defined panel; rsIDs already partially in registry. **Next.** |
| `wellness.py` | ~25   | none           | med   | Many overlaps with supplements. |
| `traits.py`   | ~80   | none           | med   | Largest; do after wellness establishes the pattern. |
| `phewas.py`   | ~120  | none           | high  | Effect sizes per variant — registry would need an extension. |
| `prs.py`      | ~150  | none           | high  | PGS weights live here; deeper schema change required. |
| `pgx.py`      | ~80   | none           | high  | CPIC pheno-bin tables; deeper schema change. |
| `hla.py`      | ~30   | none           | high  | Tag-SNP imputation has its own LD-quality field. |

## Step-by-step playbook (per module)

For each module the procedure is:

### 1. Inventory the local rsID dict

```bash
grep -nE '"rs[0-9]+"|"rsid":\s*"rs[0-9]+"' <module>.py
```

List every rsID, its locally-recorded ancestral, derived, position, and any
ad-hoc fields (effect size, GWAS source, LD quality, …).

### 2. Compare against `snp_registry.SNPS`

For each rsID present in both:

```python
from snp_registry import get
rec = get("rs1801133")
print(rec)
```

- Same metadata → drop the local copy, import via registry.
- **Different metadata → STOP.** This is the kind of latent disagreement
  the registry exists to surface. Reconcile against the literature `source`
  cited in the registry record; commit the resolution as a separate commit
  (`fix(registry): reconcile rs… per Frosst 1995`) BEFORE the rest of the
  migration. The fix may also require updating a golden snapshot — review
  the diff carefully.

### 3. Add missing rsIDs to the registry

For each rsID the local module has that the registry doesn't:

- Append a new `SNPRecord(...)` entry to `_RECORDS` in `snp_registry.py`.
- Include both GRCh37 and GRCh38 positions, cited `source`, and the current
  date as `last_verified`.
- Include any aliases (e.g. ``aliases=("C677T",)``) the module currently uses
  for display.

### 4. Refactor the consumer

- Replace local `_risk_dose` / `_dose` helpers with calls to
  `snp_registry.risk_dose_from_df(snps_df, rsid)` (alleles auto-resolved
  from the registry).
- Replace local `POSITIONS = {"rs…": …}` lookups with
  `snp_registry.get(rsid).pos_grch37`.
- Delete the local rsID dict entirely once every consumer-internal usage
  routes through the registry.

### 5. Run tests

```bash
pytest tests/ -o addopts="" -v
```

Golden snapshots are the safety net. If a behaviour change is intentional
(e.g. fixing a strand error), update the snapshot in a *separate* commit
with a clear message: `fix(<module>): correct strand for rs… (regen snapshot)`.

### 6. (Optional) Add unit tests for the migrated module

If the module had no unit tests before, add ≥1 per rule that touches a
registry SNP. Mirrors the `tests/unit/test_supplements.py` shape.

## Module-specific gotchas

### `carrier.py`
- Carrier rsIDs include the clinical-significance text ("Pathogenic when
  homozygous", etc.). Move that text into `SNPRecord.clinical_significance`
  during migration so the registry remains the single source.

### `prs.py` / `phewas.py` / `pgx.py`
- These store **per-variant effect sizes** alongside the rsID metadata.
  `SNPRecord` does NOT carry effect sizes (they're score-panel-specific, not
  variant-canonical). Keep the effect-size dicts in the score module, but
  strip the ancestral/derived/position fields from them and look those up
  via the registry instead.

### `hla.py`
- Tag-SNP records have an `ld_quality` field used by the imputer. Same
  pattern: keep the LD field local to `hla.py`, source ancestral/derived
  from the registry.

## Done criteria

The migration is complete when:

1. `snp_registry.SNPS` contains every rsID used anywhere in the project.
2. `grep -E '"rsid":\s*"rs|"rs[0-9]+":\s*\{' *.py` returns matches only in
   `snp_registry.py` itself (and possibly registry-extension modules for
   PRS / PGx / HLA effect tables).
3. `tests/registry/test_no_duplicate_rsids.py` (to be added) passes — it
   greps the codebase and asserts that no rsID appears in a local
   structure outside of `snp_registry`.
4. The `audit_registry()` summary printed at the end of CI shows
   `n_stale: 0` (every record verified within the last year).

## Status (as of V7-α)

- ✅ `snp_registry.py` created with 20 seed records covering 13 genes.
- ✅ `supplements.py` migrated (PoC, behaviour-preserved).
- ⬜ `carrier.py` — next.
- ⬜ `wellness.py`, `traits.py` — queue.
- ⬜ `phewas.py`, `prs.py`, `pgx.py`, `hla.py` — require small schema extension.

Track per-module migration commits with the conventional-commit prefix
`refactor(registry): migrate <module>` so the rollout history is greppable.
