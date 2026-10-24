# Reproduction kit: the arithmetic of one published table

`reproduce_table2.py` (Python standard library) and `reproduce_table2.R` (base R) recompute, from the printed totals of Kastrinos et al. 2021 (Gastroenterology, PMID 33839100), Table 2, the arithmetic the paper's conclusions rest on: which strategies are dominated, the efficiency frontier, each ratio recomputed beside the printed one, the optimal strategy at 100,000 USD per QALY and its net monetary benefit against two comparators, and a check of the paper's stated 3% discounting. `EXPECTED_OUTPUT.txt` is the output both scripts print on the table as printed.

## The table is not shipped

The table is the publisher's. The paper is an NIH author manuscript (PMC9330543): free to read, but not offered under an open license, so this repository does not redistribute it. To run the kit, transcribe five columns for the 20 printed strategies from Table 2 into `kastrinos2021_table2.tsv` in this folder, using `kastrinos2021_table2.TEMPLATE.tsv` for the layout: `gene`, `strategy_as_printed`, `QALYs`, `Cost` (with the `$` and thousands separators as printed) and `ICER` (`Dominated` or `--` where the table prints them). Then run either script; its output should equal `EXPECTED_OUTPUT.txt` exactly (`diff` prints nothing).

`pytest tests/public` runs both comparisons once the file is present and skips them, with that reason, until then.

## What the kit is and is not

It checks the arithmetic against the source; it does not check the source, and it does not rebuild the paper's Markov model, whose parameters are not in the table. It values nothing: the table compares colonoscopy schedules among carriers who are already known and has no "not identified" arm. One printed ratio does not reproduce from the printed totals, and the script says so.
