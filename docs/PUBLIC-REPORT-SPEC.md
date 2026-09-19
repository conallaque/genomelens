# Public report specification

Two distinct outputs exist:

```
PRODUCTION_REPORT   internal, complete, not published
PUBLIC_REPORT       allowlist-built, published
```

`PUBLIC_REPORT` is **not** a production report with fields removed. It is constructed from an explicit allowlist. The distinction matters: a blacklist fails open — a field added to production tomorrow appears publicly by default — whereas an allowlist fails closed.

## Export rule

> **Unknown fields default to DO NOT EXPORT.**

A field is publishable only if it appears on the allowlist. Adding to the allowlist is a deliberate review decision.

## Allowlist

- public sample identifier
- public input provenance (source, build, version)
- public truth-set provenance and region definitions
- finding name and high-level interpretation
- public evidence citation
- disposition: supported · conditional · unresolved · refused
- high-level economic result
- incremental cost · incremental QALYs · NMB
- ICER, or dominance as a status
- evidence-basis category
- uncertainty summary
- limitations
- report and build provenance

## Prohibited

- internal registry keys and identity codes
- routing tokens and qualification-gate internals
- state-machine serialization
- internal beneficiary or arm identifiers
- raw intermediate payoff ledgers
- hidden or default parameter maps
- commercial pricing, margin, payback or acquisition-cost assumptions
- private filesystem paths and repository references
- production source
- debugging state

## Review gate

Before any public release: source identity · export mode · private-path scan · secret scan · internal-identifier scan · commercial-assumption scan · personal-data scan · raw-genomic-data scan · rendered PDF scan · HTML scan · JSON scan · diff review · independent fresh-clone verification.

Any failure blocks the release. Scans must carry positive controls — a scanner that silently matches nothing will report a clean result forever.
