# Data provenance

## `LifeTable_USA_Mx_2015.csv`

US 2015 age-specific all-cause mortality hazard rates (`Mx`), by sex and total.

* **Retrieved:** 2026-08-19
* **From:** https://raw.githubusercontent.com/DARTH-git/cohort-modeling-tutorial-timedep/main/data/LifeTable_USA_Mx_2015.csv
* **Used by:** `econ_engine.py` — background mortality in the cohort
  state-transition model, so that a prevented disease event is not credited
  with quality-adjusted life-years the person would not have lived to collect.
* **Modified:** no — vendored byte-for-byte so the model runs offline and so any
  future mismatch can be traced to the model rather than to the input.

The DARTH tutorial derives this file from the Human Mortality Database. The
model reads ages 25–99 at run time; the file itself is unfiltered.

Underlying method reference: Alarid-Escudero F, Krijkamp EM, Enns EA, et al.
*A Tutorial on Time-Dependent Cohort State-Transition Models in R Using a
Cost-Effectiveness Analysis Example.* Med Decis Making 2023;43(1):21-41.
PMID:35924564, doi:10.1177/0272989X221121747.

## `uniprot_gene_map.tsv.gz`

UniProt accession → primary gene symbol for the reviewed (Swiss-Prot) human
proteome. Two tab-separated columns, no header, 20,283 rows, 126 KB gzipped.

* **Retrieved:** 2026-09-06
* **From:** https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&fields=accession,gene_primary&format=tsv
* **Used by:** `risk/novel_variants.py` — `gene_for_uniprot()`, to name the gene
  behind a predicted-damaging variant. AlphaMissense reports a UniProt accession
  and an Ensembl transcript but no gene symbol, so without this join every
  whole-genome prediction outside the 18 curated ACMG accessions was reported
  with no gene name at all.
* **Modified:** yes — the header row is dropped, rows with an empty gene symbol
  are dropped, and the remainder is sorted and deduplicated. No symbol is
  altered. The retrieval and filter are one shell pipeline, reproducible from
  the URL above.
* **License:** CC BY 4.0 (UniProt). Attribution: The UniProt Consortium.
  *UniProt: the Universal Protein Knowledgebase in 2025.* Nucleic Acids Res
  2025;53(D1):D609-D617. PMID:39552041, doi:10.1093/nar/gkae1010.

Coverage measured against the accessions AlphaMissense actually emits: 822 of
826 distinct accessions in a 3M-row sample, 99.5%. The four misses are one
header artefact, one blank, and two TrEMBL-only accessions absent from
Swiss-Prot by construction. Curated entries in `UNIPROT_TO_GENE` take precedence
where they exist, because they were resolved against this repository's own
tables and pick canonical isoforms deliberately.
