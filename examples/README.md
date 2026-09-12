# `genomelens_partner_v1`

A small, versioned, machine-readable surface for a partner integrating against
GenomeLens. Deliberately narrower than the internal payload: it describes what
decisions a genome can change and how well evidenced each one is, and nothing
about how those conclusions were reached.

[`genomelens_partner_v1.example.json`](genomelens_partner_v1.example.json) is
generated from a synthetic genome.

## Shape

```
schema_version   "partner-1.0"
analysis         study_id, status, source_type, genome_build, status_reason
pathways[]       decision, mechanism, economic_state, beneficiary,
                 identity_status, expected_present_value, conditional_payoff,
                 genomic_attribution_resolved, blocked_by, unlocked_by
callability      per module: snp_callable_fraction, structurally_resolved,
                 and any structural limitation
```

## Economic states

| State | Meaning |
|---|---|
| `evidence_qualified` | bookable present value |
| `expected_future` | probability-weighted, discounted |
| `conditional` | a payoff if the triggering circumstance occurs |
| `scenario_only` | illustrative; never summed |
| `unpriced` | surfaced, deliberately not valued |

## Two fields that are never merged

`expected_present_value` is probability-weighted and discounted.
`conditional_payoff` is what the decision is worth **if** its trigger occurs.
They differ by the reciprocal of the exposure probability. A record carrying a
conditional payoff under the expected-value field is rejected by the schema
rather than shipped.

## Callability is not diplotype resolution

`snp_callable_fraction` of 1.0 does not mean a resolved diplotype. Where a
gene's phenotype depends on copy number, `structurally_resolved` is `false` and
the limitation is stated on that module — not in a global footnote.

## Versioning

A breaking field change requires a new `schema_version`. Validation fails
closed: a missing required field, an unknown version, an unknown economic
state, or an unknown beneficiary is rejected.
