# Changelog

Capability-level changes to the public record. Implementation detail,
parameter values and internal module structure are not published.

## 2026-09-15 — Canonical economic-credit aggregation validated

**One health consequence is now credited exactly once, and that is
proven rather than asserted.**

- A single credit-attribution authority governs every economic total.
  Competing mechanisms were traced, classified and reduced to one.
- Canonical economic identity is authored where the decision model is
  created. No downstream layer derives identity from a gene symbol, a
  display label, a condition name or an object address.
- Two pharmacogenes governing one drug dose now resolve to one economic
  decision. Previously they produced two.
- Quantities of different economic constructs can no longer be added.
  A consequence: there is deliberately no single scalar "genome value".
- Absent totals serialize as absent. A total of `$0` now means measured
  and found to be nothing.
- A refused pathway no longer names the clinical outcome its own trial
  failed to establish. Identity and adjudication are coupled at
  construction, so the two cannot drift apart.
- Permutation invariance, idempotence and byte-identical determinism are
  asserted on economic totals.
- Counting correction: distinct economic pathways are **15**, not the
  17 previously published. 4 distinct economic pathways currently have
  defensible numeric estimates within conditional/scenario frames; 0
  are currently eligible for a portfolio headline total. The earlier
  figures counted rows where the unit is pathways.

**Method note.** Four independent adversarial reviews were run on the
finished work, each told to assume it was wrong. They found sixteen real
defects, and most of the serious ones were in code written during that
same round — including a safety check that inspected source text rather
than behaviour, and three guards that survived their own mutation tests
and therefore proved nothing. All sixteen are fixed, and each is now
covered by a test that fails when the defect is reintroduced.

## 2026-09-14 — Lifetime survival engine

- National period life table vendored, checksummed and provenance-pinned;
  reconstructs published life expectancy to within rounding.
- Central death rate and annual death probability are typed distinctly,
  so the conversion for one cannot silently be applied to the other.
- Integration uses the quadrature convention the published table was
  built with. A legacy convention whose weights did not sum to the
  interval length is quarantined from production.
- Competing risks compose on the hazard scale.
- Unsourced economic defaults reduced from 23 to 12, reported rather
  than resolved.

## 2026-09-12 — Public reference-genome benchmark

- End-to-end run against a public, open-consent reference genome with a
  published truth set, cross-checked against an independent assay of the
  same individual.
- 100.00% concordance on explicit truth-set records; 99.68% against the
  independent assay, with the single discordance disclosed and resolved
  against the array rather than dropped.
- The economic layer published nothing: every result failed at least one
  evidence gate. That was the intended outcome.
