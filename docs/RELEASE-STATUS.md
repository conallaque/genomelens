# Release status

## Current development version

The current development version is v1.11. It is private and under active development, and it is not a public software release. A release candidate of v1.11 (engine commit `e89f545c`) is described in [`VALIDATION.md`](VALIDATION.md#v111-release-candidate), with its demonstration reports and the figures from its final logs. The public benchmark below describes the frozen v1 build and does not describe or validate v1.11.

## Deterministic Benchmark v1

Published: the [Public Evidence Brief](../benchmark/v1/GenomeLens_Deterministic_Benchmark_v1_Public_Evidence_Brief.pdf), the [Public Results](../benchmark/v1/GenomeLens_Deterministic_Benchmark_v1_Public_Results.pdf) record and a cover image, for the frozen public-genome benchmark (engine `4eaf729`, tooling `a9b3d3a`). Not published: internal benchmark packages, raw run outputs, the production engine, routing logic and parameter registries. The earlier trio per-genome reports were dated results from an earlier engine build that v1 supersedes; they were removed from this repository's current tree on 2026-10-24, and the sections below state what is and is not published.


What this release publishes, what it withholds, and why. Items elsewhere in the repository
that are withheld from this release link here.

Publishing this list is deliberate. A result that is not yet attributable to a
specific engine build is not published as evidence — that rule is stated in
[`VALIDATION.md`](VALIDATION.md#provenance), and applying it visibly is more
useful to a reader than quietly omitting the affected figures.

<a id="trio-reports"></a>
## Sanitized reports: the release candidate, and the historical trio snapshots

Each run produces a genomic and health-economic report. Reports are published
only through the allowlist export profile described in
[`PUBLIC-REPORT-SPEC.md`](PUBLIC-REPORT-SPEC.md), never by taking a production
report and removing fields. The sanitized reports of the v1.11 release
candidate (16 of them, one for HG003) are published under
[`demonstrations/v1.11/`](../demonstrations/v1.11/). Three sanitized reports from
an earlier engine build (HG002, HG003, HG004) were removed from the current tree
on 2026-10-24: their economic figures came from a build whose identifier is not
published, and they do not describe v1.11, v1 or any later version. They remain in
the repository's Git history. What is kept of that build is its technical
validation record (concordance with the Genome in a Bottle truth set), in
[`VALIDATION.md`](VALIDATION.md) and the README, labeled historical.

<a id="trio-economics"></a>
## Headline economics: the release candidate, and the historical trio snapshots

Headline economic results are regenerated against each engine build and are
published only alongside the build identifier that produced them. The headline
economics of the v1.11 release candidate are published in the demonstration
reports under [`demonstrations/v1.11/`](../demonstrations/v1.11/), with the build that
produced them named in [`VALIDATION.md`](VALIDATION.md#v111-release-candidate). The earlier trio snapshots'
economic figures came from an earlier engine build whose identifier is not
published with them, which does not meet the attribution standard in
[`VALIDATION.md`](VALIDATION.md#provenance). They were removed from the current
tree, remain only in the Git history, are not representative of v1.11
validation, and are not regenerated against v1 or any later version.

<a id="matrix-counts"></a>
## Economic-validation matrix counts — not in this release

The matrix requires every tracked pathway to terminate in a defined
disposition. The claim published here is one of coverage and accountability —
every tracked row reaches an explicit terminal disposition — and not a claim
that every row is empirically validated. The counts themselves are regenerated
per build and are not published in this release.

<a id="methods-hold"></a>
## Numerical value-of-information examples — withheld

Numerical value-of-information examples are withheld pending review. The
methods are published in [`ECONOMICS.md`](ECONOMICS.md#uncertainty); the numbers
are not. Publishing a number whose sign is not established would contradict the
evidence discipline the rest of this repository describes.

The public checks look for numeric value-of-information examples by shape and
context. They do not contain the withheld values and cannot confirm the absence
of any particular one: that exact screening is done before publication by a
private release gate that this repository does not include.

<a id="benchmark-manifest"></a>
## Benchmark manifest — partly published

The trio benchmark figure is published with its truth-set release (GIAB v4.2.1), reference build (GRCh38) and chromosome scope (chr1–22), in the README and in [`VALIDATION.md`](VALIDATION.md#per-sample). It is published without the confident-region file identity, the engine build identifier of the earlier build, the retrieval date, or the rule by which benchmarked positions were selected. [`VALIDATION.md`](VALIDATION.md#provenance) sets the standard that a result must be attributable to a specific build to count as evidence; the trio figure does not yet meet it. Publishing the remaining manifest fields is the next planned addition.

<a id="synthetic-cohort"></a>
## Synthetic stress cohort — designed, not yet run

Described in [`VALIDATION.md`](VALIDATION.md#synthetic-stress-testing--planned-not-yet-run)
as planned rather than existing. It answers a different question from the GIAB
benchmark: stability and coherence at scale, never real-world accuracy,
prevalence or clinical validity.

<a id="ai-assistance"></a>
## Creator role and AI-assisted development

How the project was directed and implemented is described in
[`PROJECT-SCOPE.md`](PROJECT-SCOPE.md#creator-role).

---

*This page describes one release. Items listed here are withheld, not
abandoned.*
