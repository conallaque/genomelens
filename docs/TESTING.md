# Testing and validation

## Validation philosophy

GenomeLens separates **public verification** from **private production
testing**, because the two answer different questions.

The public suite in [`tests/public/`](../tests/public) verifies the claims this
repository makes: that the published benchmark numbers are what they say they
are, that the economic figures satisfy the identities defined in
[`ECONOMICS.md`](ECONOMICS.md), that the same number means the same thing in the
JSON, the report and the README, and that the release contains none of the
disclosure classes it is supposed to exclude.

The complete production suite stays private. Parts of it encode proprietary
implementation behavior — a test that asserts what the engine does under a given
input is a specification of that engine, and publishing enough of them would
describe the mechanism this repository deliberately withholds. The public suite
was therefore written from the published artifacts outward, not by exporting
private tests with assertions removed.

## Two suites, two different numbers

| | Private production suite | Public verification suite |
|---|---|---|
| Scale | **5,397 passing · 2 known residual failures · 24 skipped** | **55 test functions · 137 collected checks** |
| Tests | implementation behavior under given inputs | the claims and artifacts this repository publishes |
| Published | no | yes — [`tests/public/`](../tests/public) |

**The 137 public checks are not a sample of the 5,397.** They were written from
the published artifacts outward and verify different things. A private test
asserting what the engine does under a given input is a specification of that
engine; enough of them describe the mechanism this repository withholds. A
public test asserting that the benchmark figure in `summary.json` matches the
one in the report and the README verifies a claim without revealing how the
figure was produced.

## Current production-suite status

Release build of the private engine:

```
5,397 passed
    2 failed   (known, tracked)
   24 skipped
```

**This is not a fully passing suite, and it is not described as one.** Two
residual failures are carried deliberately:

- one concerns an architectural coverage question that predates this release and
  is tracked separately;
- one concerns an insertion-class variant that is **not present in any of the
  three published demonstrations** — verified directly rather than assumed.

Both are outside the HG002 / HG003 / HG004 paths published here. Neither is
described further, because the description would be the implementation detail.

## What the private suite covers

High-level categories only:

- genomic interpretation and callability semantics
- evidence handling and qualification
- health economics and estimand separation
- uncertainty and probabilistic sensitivity analysis
- reproducibility and determinism
- report generation
- regression safeguards against previously fixed defects

## What the public suite verifies

| Area | What it checks |
|---|---|
| **Benchmark integrity** | Per-sample and aggregate concordance, zero mismatches, pinned GIAB release, reference build and `chr1-22` scope, distinct sample identifiers, and that the trio's pedigree relationship stays visible next to the aggregate |
| **Cross-artifact consistency** | The same headline figure in `summary.json`, `report-public.html`, the rendered PDF and the README card — a report that silently disagrees with its own JSON is worse than one that publishes less |
| **Economic consistency** | Published NMB and canonical figures, ICER consistent with its own components, and dominance reported as a status rather than a ratio |
| **Estimand separation** | Reference-case and canonical results remain distinct fields, differ in value, are stated as never summed, and no collapsed "value of this genome" field exists |
| **Public schema** | A **closed** allowlist — every key at every level must be expected, so an unanticipated field cannot appear by being unanticipated |
| **Release security** | Private paths, credential shapes, commercial language and withheld figures, with positive and negative scanner controls, and extracted PDF text rather than raw bytes |
| **Claim discipline** | No unscoped accuracy, certification or clinical-validation claims; the scope qualifier must sit adjacent to the headline figure |
| **Missingness semantics** | An unresolved or not-assessed disposition is never rendered as benign, normal or absent; no non-computed value is published as zero |
| **Links and assets** | Every published link, report, summary and preview resolves |

Run them with:

```bash
pytest tests/public
```

Or verify the whole release in one command:

```bash
python tools/verify_public_release.py
```

Both read only files committed to this repository. No network access, no
credentials, no engine dependency.

## Two failure modes these tests were written against

The suite encodes two mistakes made during development, because both are easy to
repeat and neither is caught by ordinary testing.

**A scan that cannot fire reports a clean result forever.** An earlier scan
returned all-zeros while silently broken. Every scanning test here runs a
positive and a negative control first and fails if the controls misbehave.

**Unbounded substring matching produces confident nonsense.** An unbounded
search for `CAC` matched the letters inside "effi*cac*y" and reported commercial
language inside clinical prevention advice; a later unbounded scan flagged
`position` inside "dis*position*". The fix is word boundaries, not case
sensitivity — `\bCAC\b` does not match inside "efficacy" even case-insensitively,
because the letters there are surrounded by other letters.

Bounded matching is necessary but not sufficient. The genuinely ambiguous terms
— `CAC`, which is coronary artery calcium in clinical text and customer
acquisition cost in commercial text, and `payback`, which is cost-effectiveness
vocabulary as well as a commercial one — are adjudicated in context rather than
banned. Stripping them outright would remove legitimate clinical and
health-economic language from the published material.

## External benchmark

See [`VALIDATION.md`](VALIDATION.md), including the per-sample denominators and
the published selection chain. The benchmark driver used to produce those
figures reads subject identity from the VCF sample column and takes every truth
fact from a standard-library parser that does not import the production
resolver, so the comparison tests the resolver rather than agreeing with it by
construction.

## Reproducibility

Deterministic reruns were tested: the same engine, inputs, configuration and
seed reproduce the structured output exactly, with the sole difference being a
build identifier that embeds a timestamp.

Stochastic re-seeding was verified separately and is a different question.
Re-seeding moves the probabilistic sensitivity analysis — mean net monetary
benefit, cost-effectiveness probability, interval bounds and expected value of
perfect information all shift — while every deterministic reference-case figure
stays fixed. That is the expected signature, and observing it is what
establishes the published interpretation does not depend on a single Monte Carlo
realization. The observed shift in expected value of perfect information was a
small fraction of that quantity's own reported Monte Carlo standard error.

The default seed is not published: it is an implementation detail, and fixing a
public number to it would invite the figure to be read as more exact than it is.

## Scope and limitations

The public suite is **representative verification, not the complete engine test
suite**. It checks what this repository publishes. It does not, and cannot,
re-run the analysis — reproducing these results requires the production engine,
which is not published.

A passing public suite therefore establishes that the published artifacts are
internally consistent, correctly scoped and free of the disclosure classes
listed above. It does not establish clinical validity, and nothing here should
be read as claiming it does.
