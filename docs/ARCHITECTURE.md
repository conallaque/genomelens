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
