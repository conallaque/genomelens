# Testing and validation

## Validation philosophy

GenomeLens separates **public verification** from **private production
testing**, because the two answer different questions.

The public suite in [`tests/public/`](../tests/public) verifies the claims this
repository makes: that the historical benchmark record is stated consistently and labeled as
historical, that every published report states its research and proprietary
status, that the dossier prints only numbers the public reports and documents
already print, and that the release contains none of the disclosure classes it
is supposed to exclude.

The complete production suite stays private. Parts of it encode proprietary
implementation behavior — a test that asserts what the engine does under a given
input is a specification of that engine, and publishing enough of them would
describe the mechanism this repository deliberately withholds. The public suite
was therefore written from the published artifacts outward, not by exporting
private tests with assertions removed.

## Two suites, two different numbers

| | Private production suite | Public verification suite |
|---|---|---|
| Scale | v1.11 release candidate (engine `e89f545c`): **11,165 passed, 68 skipped, 1 failed (an environmental disk-space failure in a test that starts a pipeline from sequencing reads), full private suite on a clean snapshot of the final tree**. Earlier builds were tested against their own private suites; those results are in the private release log and are not restated here. | **107 collected checks: 107 passed** |
| Tests | implementation behavior under given inputs | the claims and artifacts this repository publishes |
| Published | no | yes — [`tests/public/`](../tests/public) |

**The 107 public checks are not a sample of the private suite.** They were written from
the published artifacts outward and verify different things. A private test
asserting what the engine does under a given input is a specification of that
engine; enough of them describe the mechanism this repository withholds. A
public test asserting that the historical concordance figure in the README
matches the one in `VALIDATION.md` verifies a claim without revealing how the
figure was produced.

## Earlier builds

Earlier release builds were tested against their own private suites. Those results are recorded in the private release log and are not restated here, and no test figure describes any build other than the one named above.

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
| **Historical concordance record** | Per-sample and aggregate concordance as stated in the README and in `VALIDATION.md`, zero mismatches, pinned GIAB release, reference build and `chr1-22` scope, the trio's pedigree relationship stated next to the aggregate, the record labeled historical and not v1.11 validation, and no per-genome report of the earlier build left in the tree |
| **Release security** | Private paths, credential shapes and commercial language, with positive and negative scanner controls and extracted PDF text rather than raw bytes; a value-free scan for numeric value-of-information examples across text files and every PDF surface (text, link targets, properties, attachment names); synthetic-canary controls proving the scanner reports locations and counts but never matched text. It does **not** screen for exact withheld values: a list of those in public code would publish them, so that check belongs to a private release gate before publication |
| **Claim discipline** | No unscoped accuracy, certification or clinical-validation claims; the scope qualifier must sit adjacent to the headline figure; every demonstration PDF states its research status and that it is proprietary |
| **Links and assets** | Every published link and image resolves, and every PNG is a real image |
| **Presentation** | The repository serves its files verbatim and the README links nothing removed |
| **HEOR dossiers (public edition)** | The dossier PDF exists, is linked from the README, carries only its title as metadata and no link or annotation, says it is proprietary, and prints only decimals and whole numbers of 100 or more (apart from years and axis marks) that already appear in that genome's public reports, the public documents or the reproduction kit's recorded output |
| **Featured six-page report (public edition)** | The PDF exists, has exactly six pages, is linked and previewed in the README, carries only its title as metadata and no link but the PubMed citation links, states that it is hypothetical, proprietary and not externally validated, discloses that its comparator is an assumption that overstates the increment, and prints only numbers that already appear in the public reports and documents or are calculated from printed terms |
| **v1.11 demonstrations and reproduction** | The 16 demonstration PDFs carry only their title as metadata and contain no local path or draft text; the README's example table matches the published files and uses only the allowed model-basis labels; the README keeps published results apart from hypothetical scenarios; the published colorectal-surveillance recomputation prints exactly its recorded output in both Python and R when the published table is supplied (the table is not shipped, so those two checks skip, with that reason, until it is), and without it each script says how to supply it and computes nothing; the preview image carries no metadata |

Checks written for the earlier trio per-genome reports (their JSON schema, economic identities, estimand separation, missingness semantics and cross-artifact consistency) were retired on 2026-10-24 together with those reports. For the v1.11 reports the equivalent guarantees are the allowlist export that builds them, the v1.11 checks above and the dossier's number rule.

Run them with:

```bash
pytest tests/public
```

The current result is `105 passed, 2 skipped` out of 107 collected. The skips are the two reproduction checks that need the published table, which is not shipped (publisher terms); each reports that reason.

The verifier prints any skip reason rather than hiding it, because a skip whose
reason is unstated is indistinguishable from a check that was quietly dropped.

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

Deterministic reruns were tested: the same engine, inputs, configuration and seed reproduce the structured output exactly, apart from the generation timestamp printed on each document.

Stochastic re-seeding was verified separately and is a different question: reference-case figures are invariant to the random seed, and the probabilistic summaries move within their stated Monte Carlo error. The default seed is not published: fixing a public number to it would invite the figure to be read as more exact than it is.

## Scope and limitations

The public suite is **representative verification, not the complete engine test
suite**. It checks what this repository publishes. It does not, and cannot,
re-run the analysis — reproducing these results requires the production engine,
which is not published.

A passing public suite therefore establishes that the published artifacts are
internally consistent, correctly scoped and free of the disclosure classes
listed above. It does not establish clinical validity, and nothing here should
be read as claiming it does.
