# Release status

What this release publishes, what it withholds, and why. Each item elsewhere in
the repository that is marked `OWNER REVIEW` or `METHODS HOLD` links here.

Publishing this list is deliberate. A result that is not yet attributable to a
specific engine build is not published as evidence — that rule is stated in
[`VALIDATION.md`](VALIDATION.md#provenance), and applying it visibly is more
useful to a reader than quietly omitting the affected figures.

<a id="trio-reports"></a>
## Sanitized trio reports — not in this release

Each run produces a genomic and health-economic report. Reports are published
only through the allowlist export profile described in
[`PUBLIC-REPORT-SPEC.md`](PUBLIC-REPORT-SPEC.md), never by taking a production
report and removing fields. No sanitized trio report is published in this
release.

<a id="trio-economics"></a>
## Headline trio economics — not in this release

Headline economic results are regenerated against each engine build and are
published only alongside the build identifier that produced them. None is
published in this release.

<a id="matrix-counts"></a>
## Economic-validation matrix counts — not in this release

The matrix requires every tracked pathway to terminate in a defined
disposition. The claim published here is one of coverage and accountability —
every tracked row reaches an explicit terminal disposition — and not a claim
that every row is empirically validated. The counts themselves are regenerated
per build and are not published in this release.

<a id="methods-hold"></a>
## Numerical value-of-information examples — withheld

`METHODS HOLD`. Numerical value-of-information examples are withheld pending
review of one treatment-effect parameter. The methods are published in
[`ECONOMICS.md`](ECONOMICS.md#uncertainty); the numbers are not.

The parameter under review sets a treatment threshold to which the results are
inversely proportional, so figures derived from it are not reported in either
direction until the review concludes. Publishing a number whose sign is not
established would contradict the evidence discipline the rest of this
repository describes.

<a id="synthetic-cohort"></a>
## Synthetic stress cohort — designed, not yet run

Described in [`VALIDATION.md`](VALIDATION.md#synthetic-stress-testing-planned-not-yet-run)
as planned rather than existing. It answers a different question from the GIAB
benchmark: stability and coherence at scale, never real-world accuracy,
prevalence or clinical validity.

<a id="ai-assistance"></a>
## Creator role and AI-assisted development

The disclosure in [`PROJECT-SCOPE.md`](PROJECT-SCOPE.md#creator-role) describes
how the project was directed and implemented. Its wording is under owner review.

---

*This page describes one release. Items listed here are withheld, not
abandoned.*
