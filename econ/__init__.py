"""Health-economics subsystem.

Eight modules that used to sit loose at the repository root, grouped because
they are one thing: the evaluation layer that turns genomic findings into a
cost-effectiveness result. Nothing else in the codebase reaches into them
except through the two entry points, ``health_economics`` (the individual's
economic sheet) and ``value_of_information`` (the pooled payer analysis).

    params                parameter provenance registry — every number that
                          enters the model, with its tier, source and
                          distribution, validated at import
    engine                pooling and the cost-effectiveness calculation:
                          Finding, ConditionPool, PSA, CEAC, tornado
    decision              the decision-analytic layer — EVPI, EVPPI,
                          breakeven, efficiency frontier, equity weighting
    markov                cohort state-transition model and budget impact
    frontier              published-model comparison
    plain                 plain-language translation of the results
    health_economics      the individual's economic sheet
    value_of_information  the pooled payer analysis and report section

This package deliberately does no work at import time — it holds no state and
re-exports nothing, so importing one module never drags in the other seven.
"""
