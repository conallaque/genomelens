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
