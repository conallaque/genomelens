# Architecture

Conceptual. Module names, routing rules and reference assets are not published.

## The chain

```
┌──────────────────┐
│ genomic evidence │  variant calls + callable reference evidence
└────────┬─────────┘
         │  refuses if identity is not established
┌────────▼─────────┐
│ clinical context │  phenotype, medication, laboratory context where permitted
└────────┬─────────┘
         │  refuses if no intervention follows
┌────────▼─────────┐
│ decision pathway │  a named clinical decision the evidence could change
└────────┬─────────┘
         │  refuses if usual care would reach it anyway
┌────────▼─────────────┐
│ usual-care comparator│  what would have happened without sequencing
└────────┬─────────────┘
         │  refuses if no sourced effect estimate exists
┌────────▼──────────────────┐
│ health/economic consequence│
└────────┬──────────────────┘
         │  refuses if the trigger probability is not sourced
┌────────▼───────────────────┐
│ genomic-attributable value │  the share sequencing added
└────────┬───────────────────┘
         │
┌────────▼──────────────────┐
│ uncertainty and provenance │  every figure traces to a named parameter
└───────────────────────────┘
```

Each stage is a gate. A pathway that fails one is not discarded — it is
reported at the stage it stopped, with the evidence that would move it on.

## Two output layers, one computation

```
                  canonical economic payload
                   (all computation happens here)
                    │                      │
        partner executive view   detailed technical view
```

Renderers format. They do not compute. A figure on a page and the same figure
in the payload cannot disagree, because there is only one of them.

## Cohort layer

Per-genome results aggregate into cohort summaries without re-running or
re-deriving any economics. Aggregation carries its own typed semantics —
how a number was formed, over how many genomes, for which beneficiary, in
which realization state — so that quantities which are not commensurable
cannot be added.

## Assay handling

One genome may arrive as an array export, a block-compressed callset, or an
all-sites callset. These are three spellings of one biology and are required to
produce one answer. The resolution layer is canonical and shared; modules do
not implement their own.


## Boundaries

Four boundaries exist in the pipeline. Each is a place where something
is allowed to refuse, and nothing downstream may re-decide what an
upstream boundary settled.

```
Genome
  |
  |   canonical observation boundary
  |   one biological fact becomes one observation, whatever assay saw it
  v
Clinical interpretation
  |
  |   validation boundary
  |   nothing numeric passes unvalidated; a failing gate yields a
  |   refusal token, never a number
  v
Economic model
  |
  |   aggregation boundary
  |   one consequence is credited once; totals are produced here and
  |   nowhere else
  v
Payload
  |
  |   privacy boundary
  |   personal artifacts are not emitted unless the gate passes;
  |   an unknown privacy state refuses
  v
Renderer
```

**The renderer renders.** It formats precomputed values and performs no
arithmetic on an economic quantity. Two output paths that each computed
their own total would diverge silently, and agreement between them is
not evidence either is right.

**Legacy boundary.** Where an older computation is retained for output
stability it is labelled as such and structurally barred from
contributing to any total. A quarantine enforced only by a comment is
not enforced; reachability is established by import-graph trace.

Stage internals, registries and gate logic are not published.
