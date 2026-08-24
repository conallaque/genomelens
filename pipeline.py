"""
Pipeline — DNA analysis orchestration.

V8 extracted this module out of analyze.py. The body of ``run_pipeline`` is
the exact orchestration code that previously lived inside ``analyze.main()``,
moved here verbatim. The split is:

    cli.py         — argparse parser
    pipeline.py    — orchestration (this file)
    analyze.py     — library functions (parse, tier1, tier2, AI) + back-compat shim
    renderers.py   — every build_*_html function
    snp_registry   — single source of truth for SNP metadata

``run_pipeline(args)`` is the entry point. ``cli.main()`` → ``analyze.main()``
→ ``run_pipeline()``. The body uses bare names (``log``, ``parse_dna_file``,
``tier1_lookup`` …) which are imported into this module's namespace via the
bulk re-export at the top — see the import block below.
"""

from __future__ import annotations

import argparse
import datetime
import html as _h
import json
import os
import sys
from pathlib import Path

from analyze import (
    PROFESSIONAL_MODULES_LOADED,
    REPORT_VERSION,
    SCRIPT_DIR,
    _build_module_context,
    ai_interpret_modules,
    analyze_addiction_genetics,
    analyze_ancestral_story,
    analyze_ancestry,
    analyze_blood_type,
    analyze_carriers,
    analyze_clinical_variants,
    analyze_deep_ancestry,
    analyze_detox,
    analyze_environmental_optimization,
    analyze_exercise,
    analyze_expanded_pgs,
    analyze_family_planning,
    analyze_genetic_age,
    analyze_gut_health,
    analyze_health_economics,
    analyze_hla,
    analyze_holistic_synthesis,
    analyze_immunogenetics,
    analyze_life_stage_playbook,
    analyze_local_ancestry,
    analyze_medications,
    analyze_metal_oxidative,
    analyze_mr,
    analyze_mt_haplogroup,
    analyze_multi_person,
    analyze_neurochemistry,
    analyze_novel_variants,
    analyze_nutrition,
    analyze_personal_economics,
    analyze_pgx,
    analyze_phewas,
    analyze_polygenic_scores,
    analyze_polygenic_traits,
    analyze_reproductive,
    analyze_urologic,
    analyze_value_of_information,
    analyze_wellness,
    # Optional module bindings (any of these can be None if the module
    # couldn't be imported — orchestration handles that via `is not None`).
    analyze_y_haplogroup,
    build_carrier_report,
    build_emergency_card,
    build_personalized_plan,
    build_supplement_stack,
    collect_references_used,
    compare_bloodwork,
    cross_category_synthesis,
    detect_interactions,
    detect_roh,
    diff_runs,
    evaluate_counseling_triggers,
    export_fhir,
    generate_narrative_report,
    html_to_pdf,
    imputation_available,
    impute_genotypes,
    load_snp_database,
    log,
    parse_dna_file,
    predict_traits,
    render_bloodwork_html,
    render_carrier_html,
    render_diff_text,
    render_economic_analysis_html,
    render_exercise_html,
    render_multi_person_html,
    render_nutrition_html,
    render_plan_html,
    render_supplements_html,
    retry_failed_categories,
    run_chat,
    run_qc,
    simulate_pgx,
    tier1_lookup,
    tier2_analysis,
    weasyprint_available,
    write_failed_categories,
)

# Bulk re-export: bring every symbol the orchestration body references into
# this module's namespace. The split between "library functions" (defined in
# analyze.py) and "orchestration" (this file) is intentional — analyze.py
# remains the canonical home for the small, well-tested helpers; pipeline.py
# owns the long sequencing logic.
from renderers import build_html_report

# Capture the module-load error for the PROFESSIONAL_MODULES_LOADED branch.
try:
    from analyze import _MODULE_LOAD_ERROR
except ImportError:
    _MODULE_LOAD_ERROR = ""

try:
    from pgx.pharmgkb import analyze_pharmgkb_clinical
except ImportError:
    analyze_pharmgkb_clinical = None

try:
    from pgx.top_drugs import analyze_top_drugs
except ImportError:
    analyze_top_drugs = None


def _render_longevity_html(integrated: dict, file_label: str = "") -> str:
    import html as _h
    def esc(s) -> str:
        return _h.escape(str(s) if s is not None else "")
    summary = integrated["executive_summary"]
    lon = integrated["longevity_composite"]
    plan = integrated["year_long_plan"]

    facts = "".join(f"<li>{esc(f)}</li>" for f in summary["key_facts"])
    actions = "".join(f"<li>{esc(a)}</li>" for a in summary["top_three_actions_this_week"])
    daily = "".join(f"<li>{esc(d)}</li>" for d in summary["non_negotiable_daily"])
    weekly = "".join(f"<li>{esc(d)}</li>" for d in summary["weekly_floors"])

    comp_rows = "".join(
        f"<tr><td>{esc(c['component'])}</td>"
        f"<td><strong>{c['score']}</strong></td>"
        f"<td>{int(c['weight']*100)}%</td></tr>"
        for c in lon["components"]
    )

    lever_rows = "".join(
        f"<tr><td><strong>{esc(lev['lever'])}</strong></td>"
        f"<td>{lev['current_score']}</td>"
        f"<td>{esc(lev['improvement_action'])}</td></tr>"
        for lev in lon["biggest_levers"]
    )

    plan_rows = "".join(
        f"<tr><td>{esc(b['mesocycle'])}</td><td>{esc(b['weeks'])}</td>"
        f"<td>{esc(b['exercise_focus'])}</td>"
        f"<td>{esc(b['nutrition_focus'])}</td>"
        f"<td>{esc('; '.join(b['labs_to_recheck']))}</td>"
        f"<td>{esc(b['goal_milestone'])}</td></tr>"
        for b in plan
    )

    return f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>Longevity Composite{(' — ' + esc(file_label)) if file_label else ''}</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
        color:#222; max-width:1100px; margin:24px auto; padding:0 16px; }}
h1 {{ font-size:1.6em; border-bottom:2px solid #333; padding-bottom:6px; }}
h2 {{ font-size:1.2em; margin-top:28px; padding-bottom:4px; border-bottom:1px solid #eee; }}
.card {{ background:#fcfcfd; border:1px solid #e2e2e6; border-radius:10px; padding:14px 16px; margin:10px 0; }}
.score-hero {{ font-size:3em; color:#1e6091; font-weight:600; }}
table {{ width:100%; border-collapse:collapse; }}
th,td {{ padding:8px 10px; border-bottom:1px solid #eee; text-align:left; vertical-align:top; }}
th {{ background:#f9f9f9; }}
</style></head><body>
<h1>Longevity Composite &amp; Integrated Plan</h1>

<h2>Executive Summary</h2>
<div class="card">
  <div class="score-hero">{lon['composite_score']}/100</div>
  <p><strong>{esc(lon['tier'])}</strong></p>
  <p>{esc(lon['interpretation'])}</p>
</div>

<h2>Key Facts</h2>
<div class="card"><ul>{facts}</ul></div>

<h2>Top 3 Actions This Week</h2>
<div class="card"><ol>{actions}</ol></div>

<h2>Non-Negotiable Daily / Weekly</h2>
<div class="card">
  <strong>Daily</strong><ul>{daily}</ul>
  <strong>Weekly floors</strong><ul>{weekly}</ul>
</div>

<h2>Component Breakdown</h2>
<div class="card">
  <table><tr><th>Component</th><th>Score</th><th>Weight</th></tr>{comp_rows}</table>
</div>

<h2>Biggest Improvable Levers</h2>
<div class="card">
  <table><tr><th>Lever</th><th>Current</th><th>Improvement Action</th></tr>{lever_rows}</table>
</div>

<h2>Year-Long Integrated Plan</h2>
<div class="card">
  <table style="font-size:0.88em">
    <tr><th>Mesocycle</th><th>Weeks</th><th>Exercise focus</th><th>Nutrition focus</th><th>Labs to recheck</th><th>Milestone</th></tr>
    {plan_rows}
  </table>
</div>

<p style="margin-top:30px;color:#888;font-size:0.85em">
Not medical advice. The composite score integrates genetic and capacity inputs only;
actual longevity outcomes are dominated by behaviour. Use this as a strategic map for
where to spend your effort over the next year.
</p>
</body></html>"""


def _stamp_html(html: str) -> str:
    """Append the build marker to a companion page before it is written.

    Every companion page lands in the same directory as the main report, and
    the pipeline does not clear that directory — deleting a user's files to
    guarantee freshness would be a worse cure than the disease. So instead each
    page says which commit produced it, and a leftover from an earlier run is
    identifiable by reading it rather than by trusting the timestamp. This is
    the write chokepoint: one place, so a page added later cannot miss it.
    """
    import build_stamp
    marker = build_stamp.build_stamp()["marker"]
    tag = (f'<div style="font-size:11px;color:#8a94a3;margin-top:18px;'
           f'font-family:ui-monospace,Menlo,monospace">{marker}</div>')
    return (html.replace("</body>", tag + "</body>", 1)
            if "</body>" in html else html + tag)


def _build_consolidated_econ_page(report_html: str, sheet_html: str,
                                  file_label: str) -> str:
    """One page carrying every economic result the run produced.

    WHY THIS EXISTS. The economics were split across two artefacts: the pooled
    payer analysis, the correction banners, the CEAC, the statistical-correction
    chain and the cohort projection all rendered inside ``report.html``, while
    the individual's sheet was its own two-page ``economic_analysis.html``.
    Anyone opening the sheet saw a correct but partial picture and had no way to
    know the rest existed — a reader cannot be expected to grep a 750 KB report
    for the other half of an answer.

    Sections are LIFTED from the generated report rather than re-rendered, which
    matters: re-rendering would create a second code path that could drift from
    the report, and drift between two views of one number is the exact defect
    this project spent a long time removing. Here there is one computation and
    two presentations of it, so they cannot disagree.
    """
    import re as _re

    def _section(html: str, sec_id: str) -> str:
        m = _re.search(rf'<section[^>]*id="{_re.escape(sec_id)}"', html)
        if not m:
            return ""
        i, depth = m.start(), 0
        for t in _re.finditer(r"<section\b|</section>", html[i:]):
            depth += 1 if t.group(0).startswith("<section") else -1
            if depth == 0:
                return html[i:i + t.end()]
        return ""

    styles = "\n".join(m.group(0) for m in
                       _re.finditer(r"<style[^>]*>.*?</style>", report_html, _re.S))
    sheet_styles = "\n".join(m.group(0) for m in
                             _re.finditer(r"<style[^>]*>.*?</style>", sheet_html, _re.S))
    voi = _section(report_html, "value-of-information")
    econ = _section(report_html, "health-economics")
    body = _re.search(r"<body[^>]*>(.*)</body>", sheet_html, _re.S)
    sheet_body = body.group(1) if body else ""

    nav = """
<div style="background:#f7f9fb;border:1px solid #e3e7ec;border-radius:10px;
     padding:12px 16px;margin:16px 0">
  <strong style="color:#12467a">What is on this page</strong>
  <ol style="margin:6px 0 0 18px;padding:0;color:#48545f;line-height:1.7">
    <li><a href="#econ-payer" style="color:#12467a">Pooled cost-effectiveness</a>
      &mdash; the payer view: findings combined once per condition, adherence
      charged, every correction shown beside the figure it replaced</li>
    <li><a href="#econ-personal" style="color:#12467a">Your economic sheet</a>
      &mdash; the same genome from your own perspective, health value kept
      separate from cash</li>
    <li><a href="#econ-detail" style="color:#12467a">Per-finding detail and
      cohort projection</a></li>
  </ol>
</div>"""

    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Health economics &mdash; {_h.escape(str(file_label))}</title>
{styles}{sheet_styles}
<style>
  @media print {{ details > summary {{ display:none }}
    .econ-break {{ break-before:page; page-break-before:always }}
    table, svg {{ break-inside:avoid }} thead {{ display:table-header-group }} }}
  .econ-h2 {{ font-size:19px; margin:20px 0 4px; padding-bottom:5px;
              border-bottom:2px solid #e3e7ec }}
</style></head>
<body class="wrap" style="max-width:1080px;margin:0 auto;padding:20px">
<h1 style="margin:0">Health economics</h1>
<div style="color:#57606a;font-size:13px;margin-top:4px">
  Every economic result from this run, in one place. Source:
  <code>{_h.escape(str(file_label))}</code></div>
{nav}
<h2 class="econ-h2" id="econ-payer">1 &middot; Pooled cost-effectiveness</h2>
{voi}
<div class="econ-break"></div>
<h2 class="econ-h2" id="econ-personal">2 &middot; Your economic sheet</h2>
{sheet_body}
<div class="econ-break"></div>
<h2 class="econ-h2" id="econ-detail">3 &middot; Per-finding detail and cohort projection</h2>
{econ}
</body></html>"""


def run_pipeline(args: argparse.Namespace) -> int:
    """Execute the analysis pipeline against parsed CLI args.

    Returns the exit code (0 = success, non-zero from sys.exit calls).
    Body below is the verbatim orchestration code that previously lived
    in analyze.main(); preserved unchanged to keep the golden snapshots
    byte-identical.
    """
    # ── --retry-failed: shortcut path that only re-runs the AI for previously
    # failed categories and patches them into the existing report.html.
    if args.retry_failed:
        log("=" * 60)
        log("DNA Analysis Tool — Retry Failed Categories")
        log("=" * 60)
        override = Path(args.output) if args.output else None
        sys.exit(retry_failed_categories(model=args.model, override_report=override))

    if not args.dna_file:
        log("ERROR: dna_file is required unless --retry-failed is used.")
        return 2

    output_path = (
        Path(args.output) if args.output else SCRIPT_DIR / "report.html"
    )

    log("=" * 60)
    log("DNA Analysis Tool — Local & Private")
    log("=" * 60)

    log("Loading SNP database ...")
    database = load_snp_database()
    log(f"  {len(database)} variants loaded")

    snps_df = parse_dna_file(args.dna_file)
    file_format = ""
    try:
        # snps library exposes source on the parser; recompute for the QC card
        from snps import SNPs as _SNPs
        file_format = _SNPs(args.dna_file, assign_par_snps=False).source
    except Exception:
        file_format = "unknown"

    # Capture the detected build NOW — pandas .attrs does not survive the concat
    # in VCF enrichment (or later transforms), so read it once and thread it as a
    # plain local rather than re-reading snps_df.attrs downstream.
    detected_build = snps_df.attrs.get("build", "unknown")

    # --assume-build overrides everything; otherwise, when probe-based detection
    # couldn't resolve the build (common for rsID-less whole-genome VCFs such as
    # GIAB benchmark callsets), fall back to reading it from the VCF header.
    _assume = getattr(args, "assume_build", None)
    if _assume:
        detected_build = _assume
        log(f"  Build forced by --assume-build: {detected_build}")
    elif detected_build not in ("grch37", "grch38"):
        try:
            import genome_input as _gi_hdr
            if _gi_hdr.looks_like_vcf(args.dna_file):
                _hb = _gi_hdr.detect_build_from_vcf_header(args.dna_file)
                if _hb in ("grch37", "grch38"):
                    log(f"  Build resolved from VCF header (rsID probes "
                        f"insufficient): {_hb}")
                    detected_build = _hb
        except Exception as _e:
            log(f"  WARNING: header-based build detection failed: {_e}")

    # ── Whole-genome VCF: back-fill curated registry positions the rsID reader
    #    dropped, and profile what the file contains that we can't yet interpret.
    vcf_profile = None
    try:
        import genome_input
        if genome_input.looks_like_vcf(args.dna_file):
            log(f"Input detected as VCF (build: {detected_build}). Enriching + profiling ...")
            snps_df, vcf_profile = genome_input.enrich_and_profile_vcf(
                snps_df, args.dna_file, detected_build, log=log)
    except Exception as e:
        log(f"  WARNING: VCF enrichment failed: {e}")

    # ── Optional two-genome comparison (--compare-genome) ──
    if getattr(args, "compare_genome", None) and analyze_multi_person is not None:
        try:
            log(f"Comparing against second genome: {args.compare_genome}")
            other_df = parse_dna_file(args.compare_genome)
            mp_result = analyze_multi_person(
                snps_df, other_df,
                label_a="You", label_b="Second genome")
            rel = mp_result["relationship"]
            log(f"  Relationship estimate: {rel['degree']} "
                f"(kinship φ={mp_result['king']['kinship']}, "
                f"{mp_result['king']['n_shared_snps']:,} shared SNPs, "
                f"{rel['confidence']} confidence)")
            # Write next to the main report output (gitignored filename).
            out_dir = Path(args.output).resolve().parent
            comp_path = out_dir / "genome_comparison.html"
            with open(comp_path, "w", encoding="utf-8") as f:
                f.write(render_multi_person_html(mp_result))
            log(f"  Wrote {comp_path}")
        except Exception as e:
            log(f"  WARNING: Genome comparison failed: {e}")

    # ── Optional imputation step ──
    imputation_info: dict | None = None
    if args.impute:
        if impute_genotypes is None:
            log("ERROR: imputation module not available")
        else:
            log("Imputation requested. Checking Beagle setup ...")
            ok, reason = imputation_available()
            log(f"  {reason}")
            if ok:
                snps_df, imputation_info = impute_genotypes(
                    snps_df, args.dna_file,
                    min_r2=args.impute_min_r2,
                )
                if imputation_info.get("from_cache"):
                    log(f"  Loaded cached imputation: "
                        f"{imputation_info['n_chip']:,} chip + "
                        f"{imputation_info['n_imputed']:,} imputed")
                else:
                    log(f"  Imputation complete: "
                        f"{imputation_info['n_chip']:,} chip + "
                        f"{imputation_info['n_imputed']:,} imputed")
            else:
                log("  Proceeding without imputation.")

    tier1_results, apoe_genotype = tier1_lookup(snps_df, database)

    # Y-DNA haplogroup analysis
    log("Analysing Y-DNA haplogroup ...")
    y_result = analyze_y_haplogroup(snps_df)
    log(f"  {y_result['message']}")

    # mtDNA haplogroup analysis (maternal lineage)
    mt_result: dict | None = None
    if analyze_mt_haplogroup is not None:
        log("Analysing mtDNA haplogroup ...")
        try:
            mt_result = analyze_mt_haplogroup(snps_df)
            log(f"  {mt_result.get('message', 'mtDNA analysis complete')}")
        except Exception as e:
            log(f"  WARNING: mtDNA analysis failed: {e}")
            mt_result = None
    else:
        log("mtDNA module not available — skipping maternal-lineage analysis")

    # ── Professional-grade analyses ──
    qc_result: dict | None = None
    prs_result: dict | None = None
    pgx_result: dict | None = None
    interactions_result: dict | None = None
    carrier_result: dict | None = None
    traits_result: dict | None = None
    counseling_result: dict | None = None
    references_used: list[dict] | None = None

    if PROFESSIONAL_MODULES_LOADED:
        log("Running professional-grade analyses ...")
        try:
            qc_result = run_qc(
                snps_df, filepath=args.dna_file,
                tier1_match_count=len(tier1_results),
                file_format=file_format, db_size=len(database),
            )
            log(f"  QC: {qc_result['grade']} (callability {qc_result['callability_pct']}%, "
                f"avg domain coverage {qc_result['average_domain_callability']}%)")
        except Exception as e:
            log(f"  WARNING: QC failed: {e}")
        try:
            prs_result = analyze_polygenic_scores(
                snps_df,
                sex=(qc_result["inferred_sex"] if qc_result else None),
            )
            n_high = sum(1 for p in prs_result["panels"].values()
                         if p["result"].get("tier") in ("Elevated", "High"))
            log(f"  PRS: {len(prs_result['panels'])} panels, {n_high} elevated/high tier")
        except Exception as e:
            log(f"  WARNING: PRS failed: {e}")
        try:
            pgx_result = analyze_pgx(snps_df)
            log(f"  PGx: {pgx_result['n_genes_tested']} genes, "
                f"{pgx_result['n_actionable_findings']} actionable drug findings, "
                f"{pgx_result.get('n_database_findings', 0)} drug-database matches")
        except Exception as e:
            log(f"  WARNING: PGx failed: {e}")
        try:
            interactions_result = detect_interactions(snps_df)
            log(f"  Compound interactions: {interactions_result['n_findings']} findings")
        except Exception as e:
            log(f"  WARNING: Interactions module failed: {e}")
        try:
            carrier_result = analyze_carriers(snps_df)
            log(f"  Carrier status: {carrier_result['n_affected']} affected, "
                f"{carrier_result['n_carriers']} carriers, "
                f"{carrier_result['n_untested']} untested")
        except Exception as e:
            log(f"  WARNING: Carrier module failed: {e}")
        try:
            traits_result = predict_traits(snps_df)
            log(f"  Traits: {traits_result['n_predictions']} predictions, "
                f"{traits_result['n_not_tested']} not tested")
        except Exception as e:
            log(f"  WARNING: Traits failed: {e}")
        try:
            references_used = collect_references_used(tier1_results)
            log(f"  References: {len(references_used)} catalogued variants matched")
        except Exception as e:
            log(f"  WARNING: References failed: {e}")
        try:
            counseling_result = evaluate_counseling_triggers(
                tier1_results=tier1_results,
                apoe_genotype=apoe_genotype or "",
                pgx_results=pgx_result or {"actionable_findings": []},
                prs_results=prs_result or {"panels": {}},
                carrier_results=carrier_result or {"affected": []},
                interactions_results=interactions_result or {"findings": []},
            )
            log(f"  Counseling triggers: {counseling_result['n_triggers']}")
        except Exception as e:
            log(f"  WARNING: Counseling module failed: {e}")
    else:
        log(f"Professional modules NOT loaded ({_MODULE_LOAD_ERROR}); using core analysis only.")

    # ── v3 modules: expanded PGS Catalog + ancestry PCA + medications ──
    expanded_pgs_result: dict | None = None
    ancestry_result: dict | None = None
    medications_result: dict | None = None

    if analyze_expanded_pgs is not None:
        try:
            log("Running expanded PGS Catalog panels ...")
            expanded_pgs_result = analyze_expanded_pgs(
                snps_df, sex=(qc_result or {}).get("inferred_sex"),
            )
            if expanded_pgs_result.get("available"):
                hf = expanded_pgs_result.get("headline_findings", [])
                log(f"  Expanded PGS: {expanded_pgs_result['n_panels']} panels, "
                    f"{len(hf)} elevated/high tier")
            else:
                log(f"  Expanded PGS not available: {expanded_pgs_result.get('reason','')}")
        except Exception as e:
            log(f"  WARNING: Expanded PGS failed: {e}")

    if analyze_ancestry is not None:
        try:
            log("Running ancestry estimation ...")
            # Pass the Y-DNA / mtDNA calls so the estimate can be cross-checked
            # against the deep paternal/maternal lineages (a small autosomal
            # panel that contradicts them is flagged, not trusted at face value).
            ancestry_result = analyze_ancestry(snps_df, y_result=y_result,
                                               mt_result=mt_result)
            if ancestry_result.get("available"):
                if ancestry_result.get("ancestry_call_suppressed"):
                    log("  Ancestry: no call — fallback marker panel is too small / "
                        "selection-confounded to classify ancestry (install the 1000G "
                        "PCA reference for a real call); trait markers shown instead.")
                else:
                    primary = ancestry_result.get("primary_population", "?")
                    conf = ancestry_result.get("confidence", "?")
                    amb = " [ambiguous]" if ancestry_result.get("ambiguous") else ""
                    log(f"  Ancestry: best match {primary} "
                        f"({conf} confidence{amb}; "
                        f"{ancestry_result.get('n_aims_independent', 0)}/"
                        f"{ancestry_result.get('n_aims_expected', '?')} markers)")
            else:
                log(f"  Ancestry: {ancestry_result.get('reason','unavailable')}")
        except Exception as e:
            log(f"  WARNING: Ancestry module failed: {e}")

    if args.medications and analyze_medications is not None:
        drugs = [d.strip() for d in args.medications.split(",") if d.strip()]
        if drugs:
            log(f"Running medication review for {len(drugs)} drugs ...")
            try:
                medications_result = analyze_medications(drugs, pgx_result or {})
                log(f"  Medications: {medications_result['n_known']} catalogued, "
                    f"{medications_result['n_unknown']} unknown")
            except Exception as e:
                log(f"  WARNING: Medications module failed: {e}")

    # ── V4: wellness predictions ──
    wellness_result: dict | None = None
    if analyze_wellness is not None:
        try:
            log("Generating wellness predictions ...")
            wellness_result = analyze_wellness(snps_df)
            log(f"  Wellness: {wellness_result['n_predictions']} predictions across "
                f"{len(wellness_result['categories'])} categories")
        except Exception as e:
            log(f"  WARNING: Wellness module failed: {e}")

    # ── V5: premium statistical genetics ──
    hla_result: dict | None = None
    roh_result: dict | None = None
    local_ancestry_result: dict | None = None
    phewas_result: dict | None = None
    mr_result: dict | None = None
    genetic_age_result: dict | None = None
    pgx_sim_result: dict | None = None
    reproductive_result: dict | None = None

    if analyze_hla is not None:
        try:
            log("HLA imputation (tag-SNP method) ...")
            hla_result = analyze_hla(snps_df)
            log(f"  HLA: {hla_result['n_alleles_called']}/{hla_result['n_alleles_tested']} alleles called, "
                f"{hla_result['n_carrier_alleles']} carrier")
        except Exception as e:
            log(f"  WARNING: HLA module failed: {e}")
    if detect_roh is not None:
        try:
            log("Runs of Homozygosity scan ...")
            roh_result = detect_roh(snps_df)
            log(f"  ROH: {roh_result['n_runs']} runs, total {roh_result['total_roh_mb']} Mb, "
                f"F_ROH = {roh_result['f_roh']} ({roh_result['context_tier']})")
        except Exception as e:
            log(f"  WARNING: ROH module failed: {e}")
    if analyze_local_ancestry is not None:
        try:
            log("Local-ancestry chromosome painting ...")
            global_props = (ancestry_result or {}).get("proportions", {})
            local_ancestry_result = analyze_local_ancestry(
                snps_df, global_proportions=global_props,
            )
            if local_ancestry_result.get("available"):
                log(f"  Local ancestry: {local_ancestry_result['n_windows_called']}/{local_ancestry_result['n_windows_total']} "
                    f"windows called, {local_ancestry_result['n_deviant']} deviant segments")
            else:
                log(f"  Local ancestry: {local_ancestry_result.get('reason','unavailable')}")
        except Exception as e:
            log(f"  WARNING: Local ancestry module failed: {e}")
    if analyze_phewas is not None:
        try:
            log("Phenome-wide trait/biomarker prediction ...")
            phewas_result = analyze_phewas(snps_df,
                                            sex=(qc_result or {}).get("inferred_sex"))
            log(f"  PheWAS: {phewas_result['n_scored']}/{phewas_result['n_traits']} traits scored, "
                f"{len(phewas_result['headline'])} extreme")
        except Exception as e:
            log(f"  WARNING: PheWAS module failed: {e}")
    if analyze_mr is not None:
        try:
            log("Mendelian randomization causal projections ...")
            mr_result = analyze_mr(snps_df, sex=(qc_result or {}).get("inferred_sex"))
            log(f"  MR: {mr_result['n_computed']}/{mr_result['n_total']} exposure-outcome pairs computed")
        except Exception as e:
            log(f"  WARNING: MR module failed: {e}")
    if analyze_genetic_age is not None:
        try:
            log("Genetic longevity / biological age proxy ...")
            genetic_age_result = analyze_genetic_age(snps_df)
            if genetic_age_result.get("available"):
                log(f"  Genetic longevity: {genetic_age_result['longevity']['percentile']}th percentile, "
                    f"{genetic_age_result['longevity_years_offset']:+.1f}y vs European mean")
        except Exception as e:
            log(f"  WARNING: Genetic age module failed: {e}")
    if simulate_pgx is not None:
        try:
            log("Quantitative PGx simulation ...")
            pgx_sim_result = simulate_pgx(pgx_result)
            n_drugs = len(pgx_sim_result.get("drugs", []))
            log(f"  PGx simulation: {n_drugs} drug models applied")
        except Exception as e:
            log(f"  WARNING: PGx simulation failed: {e}")
    if analyze_reproductive is not None:
        try:
            log("Reproductive inheritance simulator ...")
            reproductive_result = analyze_reproductive(carrier_result, roh_result)
            log(f"  Reproductive: {reproductive_result.get('n_scenarios', 0)} scenarios modeled")
        except Exception as e:
            log(f"  WARNING: Reproductive module failed: {e}")

    # ── V6: bloodwork comparison, supplements, exercise, nutrition ────────────
    bloodwork_result: dict | None = None
    if args.bloodwork:
        if compare_bloodwork is None:
            log("  Bloodwork comparison skipped: module not available")
        else:
            try:
                log(f"Comparing bloodwork JSON ({args.bloodwork}) to PheWAS predictions ...")
                _bw_meta = {"sex": (qc_result or {}).get("inferred_sex")}
                bloodwork_result = compare_bloodwork(
                    args.bloodwork, phewas_result, snps_df=snps_df, meta=_bw_meta)
                log(f"  Bloodwork: {bloodwork_result['n_matched']}/{bloodwork_result['n_labs_supplied']} "
                    f"labs compared, accuracy {bloodwork_result.get('accuracy_pct', 0)}% "
                    f"({bloodwork_result['n_confirmed']} confirmed, "
                    f"{bloodwork_result['n_partial']} partial, "
                    f"{bloodwork_result['n_diverged']} diverged)")
            except FileNotFoundError as e:
                log(f"  ERROR: {e}")
            except Exception as e:
                log(f"  WARNING: Bloodwork comparison failed: {e}")

    supplement_result: dict | None = None
    if build_supplement_stack is not None:
        try:
            log("Building personalised supplement stack ...")
            supplement_result = build_supplement_stack(
                snps_df=snps_df,
                pgx_result=pgx_result,
                wellness_result=wellness_result,
                carrier_result=carrier_result,
                phewas_result=phewas_result,
                tier1_results=tier1_results,
                bloodwork_result=bloodwork_result,
            )
            log(f"  Supplements: {supplement_result.get('n_supplements', 0)} items "
                f"(~${supplement_result.get('total_estimated_cost_usd_monthly', 0)}/mo)")
        except Exception as e:
            log(f"  WARNING: Supplements module failed: {e}")

    exercise_result: dict | None = None
    if analyze_exercise is not None:
        try:
            log("Generating personalised exercise programming ...")
            exercise_result = analyze_exercise(snps_df)
            if exercise_result.get("status") == "ok":
                log(f"  Exercise: {exercise_result['power_endurance']['bias']}; "
                    f"recovery {exercise_result['recovery']['speed'].lower()}; "
                    f"chronotype {exercise_result['chronotype']['chronotype'].lower()}")
        except Exception as e:
            log(f"  WARNING: Exercise module failed: {e}")

    nutrition_result: dict | None = None
    if analyze_nutrition is not None:
        try:
            log("Generating personalised nutrition plan ...")
            nutrition_result = analyze_nutrition(snps_df)
            if nutrition_result.get("status") == "ok":
                m = nutrition_result["macros"]
                log(f"  Nutrition: {m['pct_carbs']}C / {m['pct_fat']}F / {m['pct_protein']}P; "
                    f"caffeine {nutrition_result['caffeine']['metabolism'].lower()}")
        except Exception as e:
            log(f"  WARNING: Nutrition module failed: {e}")

    # Save Tier 1 JSON for reference / downstream use. The dict is hoisted
    # into a variable so the same structured findings feed both the on-disk
    # JSON and the health-economics module below (single source of truth).
    tier1_summary = {
        "report_version": REPORT_VERSION,
        "generated": datetime.datetime.now().isoformat(),
        "file_hash": (qc_result or {}).get("file_hash", "n/a"),
        "apoe_genotype": apoe_genotype,
        "y_haplogroup": y_result.get("haplogroup_path"),
        "y_haplogroup_status": y_result.get("status"),
        "mt_haplogroup": mt_result.get("haplogroup") if mt_result else None,
        "qc_grade": (qc_result or {}).get("grade"),
        "total_matched": len(tier1_results),
        "prs_summary": {
            name: {"tier": p["result"].get("tier"),
                   "percentile": p["result"].get("percentile")}
            for name, p in (prs_result or {}).get("panels", {}).items()
        } if prs_result else {},
        "pgx_summary": {
            g: {"phenotype": r["phenotype"], "activity_score": r.get("activity_score")}
            for g, r in (pgx_result or {}).get("per_gene", {}).items()
        } if pgx_result else {},
        "vo2max_tier": ((exercise_result or {}).get("vo2max") or {}).get("tier"),
        "longevity_percentile": (((genetic_age_result or {}).get("longevity") or {}).get("percentile")),
        "variants": tier1_results,
    }
    tier1_path = SCRIPT_DIR / "tier1_results.json"
    with open(tier1_path, "w") as f:
        json.dump(tier1_summary, f, indent=2)
    log(f"  Tier 1 results saved: {tier1_path}")

    # ── Health economics runs AFTER all analysis modules (wired below) ──
    economics_result: dict | None = None

    # ── Metal handling / oxidative-defense / neurodegeneration panel ──
    metal_oxidative_result: dict | None = None
    if analyze_metal_oxidative is not None:
        try:
            metal_oxidative_result = analyze_metal_oxidative(snps_df)
            log(f"  Metal/oxidative panel: "
                f"{metal_oxidative_result.get('n_predictions', 0)} findings")
        except Exception as e:
            log(f"  WARNING: Metal/oxidative module failed: {e}")

    # ── Detoxification & environmental resilience (smoke / PAH / metals) ──
    detox_result: dict | None = None
    if analyze_detox is not None:
        try:
            detox_result = analyze_detox(snps_df)
            if detox_result.get("available"):
                sr = detox_result.get("smoke_resilience", {})
                log(f"  Detoxification: {detox_result.get('n_findings', 0)} findings; "
                    f"smoke-resilience tier = {sr.get('tier', '?')}")
        except Exception as e:
            log(f"  WARNING: Detoxification module failed: {e}")

    # ── Urologic & genitourinary panel (OAB, BPH, prostate, stones, testis) ──
    urologic_result: dict | None = None
    if analyze_urologic is not None:
        try:
            urologic_result = analyze_urologic(snps_df)
            if urologic_result.get("available"):
                log(f"  Urologic panel: {urologic_result.get('n_findings', 0)} "
                    f"findings across {len(urologic_result.get('categories', []))} "
                    f"categories ({urologic_result.get('n_flagged', 0)} flagged)")
        except Exception as e:
            log(f"  WARNING: Urologic module failed: {e}")

    # ── Gut health (lactase, secretor, coeliac HLA tags, DAO, NOD2/IL23R) ──
    # THE MODULE THAT NEVER RAN. gut_health.py has existed with a registry
    # consistency test and no caller, so 327 lines of analysis produced nothing
    # in any report and could not reach the economic model no matter how the
    # economics were wired. Its coeliac haplotype call is the one trait with a
    # costable action behind it.
    gut_health_result: dict | None = None
    if analyze_gut_health is not None:
        try:
            gut_health_result = analyze_gut_health(snps_df)
            log(f"  Gut health: {gut_health_result.get('n_predictions', 0)} "
                f"predictions across "
                f"{len(gut_health_result.get('categories', []))} categories")
        except Exception as e:
            log(f"  WARNING: Gut-health module failed: {e}")

    # ── Deep ancestry (Neanderthal + ancient-pop + N/S Europe + timelines) ──
    deep_ancestry_result: dict | None = None
    if analyze_deep_ancestry is not None:
        try:
            deep_ancestry_result = analyze_deep_ancestry(
                snps_df, y_result=y_result, mt_result=mt_result)
            if deep_ancestry_result.get("available"):
                nea = (deep_ancestry_result.get("neanderthal") or {})
                log(f"  Deep ancestry: Neanderthal ~{nea.get('approx_pct','?')}% "
                    f"({nea.get('tier','?')}), ancient-pop top = "
                    f"{(deep_ancestry_result.get('ancient_populations') or {}).get('top','?')}")
        except Exception as e:
            log(f"  WARNING: Deep-ancestry module failed: {e}")

    # ── Blood type (ABO + Rh + secretor) ──
    blood_type_result: dict | None = None
    if analyze_blood_type is not None:
        try:
            blood_type_result = analyze_blood_type(snps_df)
            if blood_type_result.get("available"):
                log(f"  Blood type: {blood_type_result.get('combined') or '—'} "
                    f"({blood_type_result['abo']['genotype']}) · "
                    f"{blood_type_result['secretor'].get('secretor_status','?')}")
        except Exception as e:
            log(f"  WARNING: Blood-type module failed: {e}")

    # ── Neurochemistry (COMT axis + composite phenotype recommendations) ──
    neurochemistry_result: dict | None = None
    if analyze_neurochemistry is not None:
        try:
            neurochemistry_result = analyze_neurochemistry(snps_df)
            if neurochemistry_result.get("available"):
                c = neurochemistry_result["composite"]
                log(f"  Neurochemistry: COMT {c['comt_class']} · MAOA {c['maoa_class']} · "
                    f"BDNF {c['bdnf_class']}")
        except Exception as e:
            log(f"  WARNING: Neurochemistry module failed: {e}")

    # ── Holistic Synthesis is computed near the end (needs upstream outputs) ──

    # ── Addiction genetics (alcohol/opioid/nicotine/cannabis susceptibility) ──
    addiction_genetics_result: dict | None = None
    if analyze_addiction_genetics is not None:
        try:
            addiction_genetics_result = analyze_addiction_genetics(snps_df)
            if addiction_genetics_result.get("available"):
                c = addiction_genetics_result["composite"]
                log(f"  Addiction genetics: {addiction_genetics_result['n_findings']} findings · "
                    f"alcohol tier: {c['alcohol_tier']} · overall: {c['overall_tier']}")
        except Exception as e:
            log(f"  WARNING: Addiction-genetics module failed: {e}")

    # ── Trait genetics (genotype-level single-variant traits) ──
    polygenic_traits_result: dict | None = None
    if analyze_polygenic_traits is not None:
        try:
            polygenic_traits_result = analyze_polygenic_traits(snps_df)
            if polygenic_traits_result.get("available"):
                log(f"  Trait genetics: {polygenic_traits_result['n_findings']} "
                    f"single-variant trait calls")
        except Exception as e:
            log(f"  WARNING: Trait-genetics module failed: {e}")

    # ── TNRC18 rs117910193 novelty marker (single-variant genotype read-out) ──
    tnrc18_result: dict | None = None
    try:
        from risk.tnrc18_marker import analyze_tnrc18_marker
        tnrc18_result = analyze_tnrc18_marker(snps_df)
        if tnrc18_result.get("available"):
            log(f"  TNRC18 marker (rs117910193): {tnrc18_result['genotype_oriented']} "
                f"→ {tnrc18_result['marker']}")
        else:
            log("  TNRC18 marker (rs117910193): not typed in this file")
    except Exception as e:
        log(f"  WARNING: TNRC18 marker module failed: {e}")

    # ── Environmental optimization (light / exercise / vitamin-D × latitude) ──
    environmental_optimization_result: dict | None = None
    if analyze_environmental_optimization is not None:
        try:
            environmental_optimization_result = analyze_environmental_optimization(
                snps_df, latitude=getattr(args, "latitude", 40.0))
            if environmental_optimization_result.get("available"):
                ex = (environmental_optimization_result.get("exercise") or {})
                log(f"  Environmental optimization: exercise lean = "
                    f"{ex.get('lean','n/a')} · latitude "
                    f"{environmental_optimization_result['latitude_assumed']:.0f}°")
        except Exception as e:
            log(f"  WARNING: Environmental-optimization module failed: {e}")

    # ── Clinical variants (Phase 2 — ClinVar P/LP screen; VCF input only) ──
    clinical_variants_result: dict | None = None
    if analyze_clinical_variants is not None:
        try:
            import genome_input as _gi
            if _gi.looks_like_vcf(args.dna_file):
                clinical_variants_result = analyze_clinical_variants(
                    args.dna_file, detected_build,
                    inferred_sex=(qc_result or {}).get("inferred_sex"), log=log)
                if clinical_variants_result.get("available"):
                    log(f"  Clinical variants (ClinVar): {clinical_variants_result['n_plp']} "
                        f"P/LP · {clinical_variants_result['n_actionable']} actionable · "
                        f"{clinical_variants_result['n_carrier']} carrier · "
                        f"{clinical_variants_result['n_affected']} affected-consistent")
                else:
                    log(f"  Clinical variants: {clinical_variants_result.get('reason','n/a')}")
        except Exception as e:
            log(f"  WARNING: Clinical-variants module failed: {e}")

    # ── Novel/rare variants (Phase 3 — computational predictor screen; VCF only) ──
    novel_variants_result: dict | None = None
    if analyze_novel_variants is not None:
        try:
            import genome_input as _gi3
            if _gi3.looks_like_vcf(args.dna_file):
                novel_variants_result = analyze_novel_variants(
                    args.dna_file, detected_build,
                    clinvar_result=clinical_variants_result,
                    inferred_sex=(qc_result or {}).get("inferred_sex"),
                    commercial_safe=getattr(args, "commercial_safe", False),
                    log=log)
                if novel_variants_result.get("available"):
                    preds = [p["name"] for p in novel_variants_result.get("predictors_used", [])]
                    log(f"  Novel variants (predicted): "
                        f"{novel_variants_result['n_predicted_pathogenic']} predicted damaging · "
                        f"{novel_variants_result['n_rare_damaging']} rare/absent · predictors {preds}")
                else:
                    log(f"  Novel variants: {novel_variants_result.get('reason','n/a')}")
        except Exception as e:
            log(f"  WARNING: Novel-variants module failed: {e}")

    voi_result: dict | None = None
    # (VOI + health economics moved below, after all analysis modules)

    # ── Family planning (reproductive genetics — needs carrier + mt + sex) ──
    family_planning_result: dict | None = None
    if analyze_family_planning is not None:
        try:
            family_planning_result = analyze_family_planning(
                carrier_result=carrier_result,
                tier1_results=tier1_results,
                mt_result=mt_result,
                snps_df=snps_df,
                inferred_sex=(qc_result or {}).get("inferred_sex"),
                ancestry="European",
            )
            if family_planning_result.get("available"):
                log(f"  Family planning: {family_planning_result['n_recessive']} recessive · "
                    f"{family_planning_result['n_dominant']} dominant · mtDNA "
                    f"sex-gate = {family_planning_result['mtdna']['sex']}")
        except Exception as e:
            log(f"  WARNING: Family-planning module failed: {e}")

    # ── Immunogenetics (viral/bacterial/parasitic resistance + selection) ──
    immunogenetics_result: dict | None = None
    if analyze_immunogenetics is not None:
        try:
            immunogenetics_result = analyze_immunogenetics(snps_df)
            if immunogenetics_result.get("available"):
                log(f"  Immunogenetics: {immunogenetics_result['n_findings']} findings "
                    f"({immunogenetics_result['n_protective']} protective, "
                    f"{immunogenetics_result['n_susceptible']} susceptible)")
        except Exception as e:
            log(f"  WARNING: Immunogenetics module failed: {e}")

    # ── Ancestral Story (template narrative; AI-enhanced if AI enabled) ──
    ancestral_story_result: dict | None = None
    if analyze_ancestral_story is not None:
        try:
            ancestral_story_result = analyze_ancestral_story(
                y_result=y_result, mt_result=mt_result,
                deep_ancestry_result=deep_ancestry_result,
                immunogenetics_result=immunogenetics_result,
                ancestry_result=ancestry_result,
                model=args.model,
                use_ai=not args.no_ai,
            )
            if ancestral_story_result.get("available"):
                mode = "AI-enhanced" if ancestral_story_result.get("ai_used") else "template"
                log(f"  Ancestral story: {mode} narrative built with "
                    f"{len(ancestral_story_result['template']['chapters'])} chapters")
        except Exception as e:
            log(f"  WARNING: Ancestral-story module failed: {e}")

    # ── ClinPGx/PharmGKB clinical-variant annotations for typed rsIDs ──
    pharmgkb_result: dict | None = None
    if analyze_pharmgkb_clinical is not None:
        try:
            pharmgkb_result = analyze_pharmgkb_clinical(snps_df)
            if pharmgkb_result.get("available"):
                log(f"  PharmGKB annotations: {pharmgkb_result['n_typed_variants']} typed "
                    f"positions ({pharmgkb_result['n_high']} high-evidence), "
                    f"{pharmgkb_result['n_drugs']} drugs")
        except Exception as e:
            log(f"  WARNING: PharmGKB clinical module failed: {e}")

    # ── Top-prescribed-drugs pharmacogenomic screen ──
    top_drugs_result: dict | None = None
    if analyze_top_drugs is not None:
        try:
            top_drugs_result = analyze_top_drugs(snps_df, pgx_result)
            if top_drugs_result.get("available"):
                log(f"  Top-drugs screen: {top_drugs_result['n_screened']} drugs "
                    f"({top_drugs_result['n_actionable']} genotype-actionable, "
                    f"{top_drugs_result['n_with_pgx']} with PGx data)")
        except Exception as e:
            log(f"  WARNING: Top-drugs screen failed: {e}")

    # ── Dev affordance: snapshot the economics inputs for offline model work ──
    # Set GENOMELENS_ECON_SNAPSHOT=/abs/path.pkl to capture every module result
    # that feeds the economic model, so the econ layer can be iterated on (and
    # its numbers diffed) without re-running the whole genome pipeline. Never on
    # by default; the snapshot is genotype-derived, so it must be written
    # outside the repo.
    _econ_snap = os.environ.get("GENOMELENS_ECON_SNAPSHOT")
    if _econ_snap:
        try:
            import pickle as _pickle
            with open(_econ_snap, "wb") as _fh:
                _pickle.dump({
                    "tier1_summary": tier1_summary,
                    "expanded_pgs_result": expanded_pgs_result,
                    "hla_result": hla_result,
                    "carrier_result": carrier_result,
                    "interactions_result": interactions_result,
                    "addiction_result": addiction_genetics_result,
                    "metal_oxidative_result": metal_oxidative_result,
                    "mr_result": mr_result,
                    "neurochemistry_result": neurochemistry_result,
                    "urologic_result": urologic_result,
                    "clinical_variants_result": clinical_variants_result,
                    "phewas_result": phewas_result,
                    "immunogenetics_result": immunogenetics_result,
                    "wellness_result": wellness_result,
                    "detox_result": detox_result,
                    "family_planning_result": family_planning_result,
                    "top_drugs_result": top_drugs_result,
                    "novel_variants_result": novel_variants_result,
                    "genetic_age_result": genetic_age_result,
                    "roh_result": roh_result,
                    "bloodwork_result": bloodwork_result,
                    "qc_result": qc_result,
                }, _fh)
            log(f"  [dev] economics inputs snapshotted -> {_econ_snap}")
        except Exception as _se:
            log(f"  [dev] economics snapshot failed: {_se}")

    # ── Health economics (runs AFTER all analysis modules for maximum coverage) ──
    if analyze_health_economics is not None:
        try:
            economics_result = analyze_health_economics(
                tier1_summary, snps_df,
                expanded_pgs_result=expanded_pgs_result,
                hla_result=hla_result,
                carrier_result=carrier_result,
                interactions_result=interactions_result,
                addiction_result=addiction_genetics_result,
                metal_oxidative_result=metal_oxidative_result,
                mr_result=mr_result,
                neurochemistry_result=neurochemistry_result,
                urologic_result=urologic_result,
                clinical_variants_result=clinical_variants_result,
                phewas_result=phewas_result,
                immunogenetics_result=immunogenetics_result,
                wellness_result=wellness_result,
                detox_result=detox_result,
                family_planning_result=family_planning_result,
                top_drugs_result=top_drugs_result,
                    gut_health_result=gut_health_result)
            log(f"  Health economics: {economics_result.get('n_findings', 0)} findings "
                f"with modeled ROI")
        except Exception as e:
            log(f"  WARNING: Health economics module failed: {e}")

    # ── Value of Information (health economics — runs for chip AND WGS) ──
    if analyze_value_of_information is not None:
        try:
            import genome_input as _gi4
            _input_type = "wgs" if _gi4.looks_like_vcf(args.dna_file) else "chip"
            voi_result = analyze_value_of_information(
                economics_result=economics_result,
                clinical_variants_result=clinical_variants_result,
                novel_variants_result=novel_variants_result,
                genetic_age_result=genetic_age_result,
                roh_result=roh_result,
                input_type=_input_type, log=log)
            if voi_result.get("available"):
                log(f"  Value of Information ({_input_type}): expected genome value ≈ "
                    f"${voi_result.get('voi_expost_mean', 0):,.0f} · chip→WGS marginal ≈ "
                    f"${voi_result.get('marginal_chip_to_wgs', 0):,.0f}")
                try:
                    import cohort_simulator as _cs
                    from life_stage_playbook import resolve_age as _resolve_age
                    _pa = _resolve_age(getattr(args, "age", None), bloodwork_result)
                    voi_result["personalized"] = _cs.personalize_for_report(
                        voi_result, age=float(_pa or 40.0))
                    if voi_result["personalized"].get("available"):
                        _p = voi_result["personalized"]
                        log(f"    Personalized: efficient choice = "
                            f"{_p['frontier']['recommended_strategy']} · "
                            f"{_p['population_percentile']}th population percentile")
                except Exception as _pe:
                    log(f"    (personalized panel skipped: {_pe})")
                _cp = voi_result.get("carrier_panel_prior") or {}
                if _cp.get("available"):
                    log(f"    Carrier-panel prior (ROH, not monetised): "
                        f"{_cp['recommendation']}")
            else:
                log(f"  Value of Information: {voi_result.get('reason', 'n/a')}")
        except Exception as e:
            log(f"  WARNING: Value-of-Information module failed: {e}")

    ai_results: dict[str, str] = {}
    # ── Holistic Synthesis — computed BEFORE the AI block so its Genome
    #    Leverage Score can frame the executive summary + synthesis prompts. ──
    holistic_synthesis_result: dict | None = None
    if analyze_holistic_synthesis is not None:
        try:
            _sex = (qc_result or {}).get("inferred_sex")
            holistic_synthesis_result = analyze_holistic_synthesis(
                tier1_results={
                    "apoe_genotype": apoe_genotype,
                    "prs_summary": {p: r for p, r in
                                    ((prs_result or {}).get("panels") or {}).items()},
                    "variants": tier1_results,
                },
                bloodwork_result=bloodwork_result,
                immunogenetics_result=immunogenetics_result,
                neurochemistry_result=neurochemistry_result,
                deep_ancestry_result=deep_ancestry_result,
                ancestry_result=ancestry_result,
                prs_result=prs_result,
                pgx_result=pgx_result,
                meta={"sex": _sex, "smoker": False},
            )
            if holistic_synthesis_result.get("available"):
                gl = holistic_synthesis_result.get("genome_leverage") or {}
                log(f"  Holistic synthesis: {holistic_synthesis_result['n_insights']} "
                    f"cross-panel insights · Genome Leverage {gl.get('score','?')}/100 "
                    f"({gl.get('tier','?')})")
        except Exception as e:
            log(f"  WARNING: Holistic synthesis failed: {e}")

    # ── Shared module-context digest fed to every whole-report AI call ──
    module_context = ""
    try:
        module_context = _build_module_context(
            bloodwork_result=bloodwork_result,
            immunogenetics_result=immunogenetics_result,
            neurochemistry_result=neurochemistry_result,
            addiction_genetics_result=addiction_genetics_result,
            deep_ancestry_result=deep_ancestry_result,
            blood_type_result=blood_type_result,
            family_planning_result=family_planning_result,
            polygenic_traits_result=polygenic_traits_result,
            environmental_optimization_result=environmental_optimization_result,
            holistic_synthesis_result=holistic_synthesis_result,
            detox_result=detox_result,
            ancestry_result=ancestry_result,
            clinical_variants_result=clinical_variants_result,
            novel_variants_result=novel_variants_result,
            voi_result=voi_result,
            y_result=y_result,
            mt_result=mt_result,
        )
    except Exception as e:
        log(f"  WARNING: module-context digest failed: {e}")

    exec_summary: str | None = None
    cross_cat: str | None = None
    failed_categories: list[dict] = []
    module_ai: dict[str, str] = {}

    if not args.no_ai:
        log(f"Starting Tier 2 AI analysis (model: {args.model}) ...")
        try:
            ai_results, exec_summary, failed_categories = tier2_analysis(
                tier1_results, apoe_genotype, model=args.model,
                module_context=module_context,
            )
            log("  AI analysis complete.")
            # Cross-Category Interactions synthesis (now module-aware).
            try:
                cross_cat = cross_category_synthesis(
                    tier1_results, model=args.model, module_context=module_context)
            except Exception as e:
                log(f"  WARNING: cross-category synthesis failed: {e}")
                cross_cat = None
            # AI on ALL tiers — a per-module AI interpretation for each module
            # section (skippable with --no-module-ai for a faster run).
            if not getattr(args, "no_module_ai", False):
                try:
                    module_ai = ai_interpret_modules(
                        model=args.model, log=log,
                        bloodwork_result=bloodwork_result,
                        immunogenetics_result=immunogenetics_result,
                        neurochemistry_result=neurochemistry_result,
                        addiction_genetics_result=addiction_genetics_result,
                        deep_ancestry_result=deep_ancestry_result,
                        blood_type_result=blood_type_result,
                        family_planning_result=family_planning_result,
                        polygenic_traits_result=polygenic_traits_result,
                        environmental_optimization_result=environmental_optimization_result,
                        holistic_synthesis_result=holistic_synthesis_result,
                        clinical_variants_result=clinical_variants_result,
                        novel_variants_result=novel_variants_result,
                        voi_result=voi_result,
                    )
                    log(f"  Per-module AI: {len(module_ai)} module interpretations generated")
                except Exception as e:
                    log(f"  WARNING: per-module AI failed: {e}")
        except ConnectionError as e:
            log(f"  ERROR: {e}")
            log("  Generating Tier 1-only report.")
            args.no_ai = True
        except Exception as e:
            log(f"  ERROR during AI analysis: {e}")
            log("  Generating report with partial results.")
    else:
        log("Skipping AI analysis (--no-ai)")

    # ── Life-stage playbook (runs AFTER holistic_synthesis — needs leverage) ──
    life_stage_playbook_result: dict | None = None
    if analyze_life_stage_playbook is not None:
        try:
            # Age priority: --age flag, else bloodwork labs age; see resolve_age.
            from life_stage_playbook import resolve_age
            _age = resolve_age(getattr(args, "age", None), bloodwork_result)
            life_stage_playbook_result = analyze_life_stage_playbook(
                age=_age,
                holistic_synthesis_result=holistic_synthesis_result,
                immunogenetics_result=immunogenetics_result,
                addiction_genetics_result=addiction_genetics_result,
                neurochemistry_result=neurochemistry_result,
                family_planning_result=family_planning_result,
                tier1_results={"apoe_genotype": apoe_genotype},
            )
            log(f"  Life-stage playbook: current decade = "
                f"{life_stage_playbook_result.get('current_decade') or 'unknown (no age)'}")
        except Exception as e:
            log(f"  WARNING: Life-stage-playbook module failed: {e}")

    log("Generating HTML report ...")
    html = build_html_report(
        tier1_results=tier1_results,
        apoe_genotype=apoe_genotype,
        ai_results=ai_results,
        exec_summary=exec_summary,
        dna_filepath=args.dna_file,
        no_ai=args.no_ai,
        model=args.model,
        y_result=y_result,
        mt_result=mt_result,
        cross_cat_synthesis=cross_cat,
        qc_result=qc_result,
        prs_result=prs_result,
        pgx_result=pgx_result,
        interactions_result=interactions_result,
        carrier_result=carrier_result,
        traits_result=traits_result,
        counseling_result=counseling_result,
        references_used=references_used,
        imputation_info=imputation_info,
        expanded_pgs_result=expanded_pgs_result,
        ancestry_result=ancestry_result,
        medications_result=medications_result,
        wellness_result=wellness_result,
        hla_result=hla_result,
        roh_result=roh_result,
        local_ancestry_result=local_ancestry_result,
        phewas_result=phewas_result,
        mr_result=mr_result,
        genetic_age_result=genetic_age_result,
        pgx_sim_result=pgx_sim_result,
        reproductive_result=reproductive_result,
        economics_result=economics_result,
        pharmgkb_result=pharmgkb_result,
        top_drugs_result=top_drugs_result,
        metal_oxidative_result=metal_oxidative_result,
        detox_result=detox_result,
        urologic_result=urologic_result,
        deep_ancestry_result=deep_ancestry_result,
        blood_type_result=blood_type_result,
        immunogenetics_result=immunogenetics_result,
        ancestral_story_result=ancestral_story_result,
        neurochemistry_result=neurochemistry_result,
        holistic_synthesis_result=holistic_synthesis_result,
        addiction_genetics_result=addiction_genetics_result,
        family_planning_result=family_planning_result,
        polygenic_traits_result=polygenic_traits_result,
        tnrc18_result=tnrc18_result,
        environmental_optimization_result=environmental_optimization_result,
        life_stage_playbook_result=life_stage_playbook_result,
        clinical_variants_result=clinical_variants_result,
        novel_variants_result=novel_variants_result,
        voi_result=voi_result,
        module_ai=module_ai,
    )

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)

    log(f"  Report saved: {output_path}")

    # Persist the list of categories whose AI didn't fully succeed (always
    # overwrites — an empty list means the last run was clean) so the user
    # can re-run just those via `python analyze.py --retry-failed`.
    if not args.no_ai:
        write_failed_categories(
            failed_categories,
            model=args.model,
            report_path=output_path,
            dna_file=args.dna_file,
        )

    # ── v3: optional PDF export ──
    if args.pdf:
        if html_to_pdf is None or not weasyprint_available():
            log("  PDF skipped: weasyprint not installed. `pip install weasyprint`")
        else:
            pdf_path = output_path.with_suffix(".pdf")
            log(f"  Rendering PDF -> {pdf_path} ...")
            msg = html_to_pdf(
                html_path=output_path,
                pdf_path=pdf_path,
                file_label=Path(args.dna_file).name,
                file_hash=(qc_result or {}).get("file_hash", "n/a"),
                version=f"v{REPORT_VERSION}",
                qc_grade=(qc_result or {}).get("grade", ""),
            )
            log(f"  {msg}")

    # ── v3: optional standalone Carrier Report ──
    if args.carrier_report:
        if build_carrier_report is None or carrier_result is None:
            log("  Carrier report skipped: module or carrier data unavailable")
        else:
            log("Generating standalone carrier / family planning report ...")
            cr = build_carrier_report(carrier_result)
            cr_html = render_carrier_html(
                cr,
                file_label=Path(args.dna_file).name,
                version=f"v{REPORT_VERSION}",
            )
            cr_path = output_path.parent / "carrier_report.html"
            cr_path.write_text(_stamp_html(cr_html), encoding="utf-8")
            log(f"  Carrier report saved: {cr_path}")
            if args.pdf and html_to_pdf is not None and weasyprint_available():
                cr_pdf = cr_path.with_suffix(".pdf")
                msg = html_to_pdf(
                    html_path=cr_path, pdf_path=cr_pdf,
                    file_label=Path(args.dna_file).name,
                    file_hash=(qc_result or {}).get("file_hash", "n/a"),
                    version=f"v{REPORT_VERSION}",
                )
                log(f"  Carrier PDF: {msg}")

    # ── v5: optional Emergency Card ──
    if args.emergency_card:
        if build_emergency_card is None:
            log("  Emergency card skipped: module not loaded")
        else:
            log("Generating emergency medical genetics card ...")
            try:
                ec_html = build_emergency_card(
                    file_label=Path(args.dna_file).name,
                    file_hash=(qc_result or {}).get("file_hash", "n/a"),
                    version=f"v{REPORT_VERSION}",
                    apoe_genotype=apoe_genotype,
                    carrier_result=carrier_result,
                    pgx_result=pgx_result,
                    pgx_sim_result=pgx_sim_result,
                    hla_result=hla_result,
                    interactions_result=interactions_result,
                    report_link=str(output_path.name),
                )
                ec_path = output_path.parent / "emergency_card.html"
                ec_path.write_text(_stamp_html(ec_html), encoding="utf-8")
                log(f"  Emergency card saved: {ec_path}")
                if args.pdf and html_to_pdf is not None and weasyprint_available():
                    ec_pdf = ec_path.with_suffix(".pdf")
                    msg = html_to_pdf(html_path=ec_path, pdf_path=ec_pdf,
                                      file_label=Path(args.dna_file).name,
                                      file_hash=(qc_result or {}).get("file_hash", "n/a"),
                                      version=f"v{REPORT_VERSION}")
                    log(f"  Emergency PDF: {msg}")
            except Exception as e:
                log(f"  WARNING: Emergency card failed: {e}")

    # ── v5: optional Narrative report ──
    if args.narrative:
        if generate_narrative_report is None:
            log("  Narrative report skipped: module not loaded")
        elif args.no_ai:
            log("  Narrative report skipped: --no-ai also set")
        else:
            log("Generating natural-language narrative report ...")
            try:
                narrative_path = output_path.parent / "narrative_report.html"
                generate_narrative_report(
                    output_path=narrative_path,
                    file_label=Path(args.dna_file).name,
                    file_hash=(qc_result or {}).get("file_hash", "n/a"),
                    version=f"v{REPORT_VERSION}",
                    model=args.model,
                    tier1_results=tier1_results,
                    apoe_genotype=apoe_genotype,
                    y_result=y_result, mt_result=mt_result,
                    pgx_result=pgx_result, prs_result=prs_result,
                    interactions_result=interactions_result,
                    carrier_result=carrier_result, traits_result=traits_result,
                    wellness_result=wellness_result, ancestry_result=ancestry_result,
                    hla_result=hla_result, roh_result=roh_result,
                    phewas_result=phewas_result, mr_result=mr_result,
                    genetic_age_result=genetic_age_result,
                    counseling_result=counseling_result,
                )
                log(f"  Narrative report saved: {narrative_path}")
                if args.pdf and html_to_pdf is not None and weasyprint_available():
                    np_pdf = narrative_path.with_suffix(".pdf")
                    msg = html_to_pdf(html_path=narrative_path, pdf_path=np_pdf,
                                      file_label=Path(args.dna_file).name,
                                      file_hash=(qc_result or {}).get("file_hash", "n/a"),
                                      version=f"v{REPORT_VERSION}")
                    log(f"  Narrative PDF: {msg}")
            except Exception as e:
                log(f"  WARNING: Narrative report failed: {e}")

    # ── V6: standalone HTML outputs (supplements / exercise / nutrition / bloodwork) ──
    file_label = Path(args.dna_file).name if args.dna_file else ""

    if supplement_result and render_supplements_html is not None:
        try:
            supp_path = output_path.parent / "supplements.html"
            supp_path.write_text(
                _stamp_html(render_supplements_html(supplement_result, file_label=file_label)),
                encoding="utf-8",
            )
            log(f"  Supplements page saved: {supp_path}")
        except Exception as e:
            log(f"  WARNING: Supplements HTML failed: {e}")

    if exercise_result and render_exercise_html is not None:
        try:
            ex_path = output_path.parent / "exercise.html"
            ex_path.write_text(
                _stamp_html(render_exercise_html(exercise_result, file_label=file_label)),
                encoding="utf-8",
            )
            log(f"  Exercise page saved: {ex_path}")
        except Exception as e:
            log(f"  WARNING: Exercise HTML failed: {e}")

    if nutrition_result and render_nutrition_html is not None:
        try:
            nu_path = output_path.parent / "nutrition.html"
            nu_path.write_text(
                _stamp_html(render_nutrition_html(nutrition_result, file_label=file_label)),
                encoding="utf-8",
            )
            log(f"  Nutrition page saved: {nu_path}")
        except Exception as e:
            log(f"  WARNING: Nutrition HTML failed: {e}")

    if bloodwork_result and render_bloodwork_html is not None:
        try:
            bw_path = output_path.parent / "bloodwork.html"
            bw_path.write_text(
                _stamp_html(render_bloodwork_html(bloodwork_result, file_label=file_label)),
                encoding="utf-8",
            )
            log(f"  Bloodwork comparison saved: {bw_path}")
        except Exception as e:
            log(f"  WARNING: Bloodwork HTML failed: {e}")

    # Personal economic-impact sheet — standalone economic_analysis.html
    if analyze_personal_economics is not None and render_economic_analysis_html is not None:
        try:
            personal_econ = analyze_personal_economics(
                economics_result=economics_result,
                bloodwork_result=bloodwork_result,
                genetic_age_result=genetic_age_result,
                meta={"sex": (qc_result or {}).get("inferred_sex")},
                carrier_result=carrier_result,
                hla_result=hla_result,
                interactions_result=interactions_result,
                expanded_pgs_result=expanded_pgs_result,
                addiction_result=addiction_genetics_result,
                neurochemistry_result=neurochemistry_result,
                mr_result=mr_result,
                clinical_variants_result=clinical_variants_result,
                family_planning_result=family_planning_result,
                phewas_result=phewas_result,
                wellness_result=wellness_result,
            )
            if personal_econ.get("available"):
                econ_path = output_path.parent / "economic_analysis.html"
                econ_path.write_text(
                    _stamp_html(render_economic_analysis_html(personal_econ, file_label=file_label)),
                    encoding="utf-8",
                )
                # One consolidated page carrying every economic result, so
                # the pooled analysis and the individual sheet are no longer in
                # two separate artefacts a reader has to know to look for.
                try:
                    _all_econ = output_path.parent / "economics.html"
                    _all_econ.write_text(_stamp_html(
                        _build_consolidated_econ_page(
                            output_path.read_text(encoding="utf-8"),
                            econ_path.read_text(encoding="utf-8"),
                            file_label)), encoding="utf-8")
                    log(f"  Consolidated economics saved: {_all_econ}")
                except Exception as _ce:
                    log(f"  WARNING: consolidated economics page failed: {_ce}")
                log(f"  Economic-impact analysis saved: {econ_path} "
                    f"(modeled net benefit {personal_econ['total_net']:,} "
                    f"= {personal_econ['total_qaly_value']:,} health value + "
                    f"{personal_econ['net_cash']:,} net cash · "
                    f"{personal_econ['verdict']})")
        except Exception as e:
            log(f"  WARNING: Economic-analysis HTML failed: {e}")

    # Master dashboard — synthesises all V6 outputs into one page
    if build_personalized_plan is not None and render_plan_html is not None:
        try:
            plan = build_personalized_plan(
                supplement_result=supplement_result,
                exercise_result=exercise_result,
                nutrition_result=nutrition_result,
                bloodwork_result=bloodwork_result,
                phewas_result=phewas_result,
            )
            plan_path = output_path.parent / "personalized_plan.html"
            plan_path.write_text(
                _stamp_html(render_plan_html(plan, file_label=file_label,
                                 report_link=output_path.name)),
                encoding="utf-8",
            )
            log(f"  Master plan dashboard saved: {plan_path}")
        except Exception as e:
            log(f"  WARNING: Personalized plan failed: {e}")

    # ── Longevity composite & year-long integrated plan ──
    try:
        from longevity import integrated_longevity_plan
        if nutrition_result and nutrition_result.get("status") == "ok" \
                and exercise_result and exercise_result.get("status") == "ok":
            integrated = integrated_longevity_plan(nutrition_result, exercise_result)
            long_path = output_path.parent / "longevity.html"
            long_path.write_text(_stamp_html(_render_longevity_html(integrated, file_label)),
                                  encoding="utf-8")
            log(f"  Longevity composite + year-long plan saved: {long_path}")
    except Exception as e:
        log(f"  WARNING: Longevity composite failed: {e}")

    # ── V6: FHIR clinical export ──
    if args.fhir:
        if export_fhir is None:
            log("  FHIR export skipped: module not available")
        else:
            try:
                fhir_path = output_path.parent / "fhir_bundle.json"
                summary = export_fhir(
                    out_path=fhir_path,
                    pgx_result=pgx_result,
                    carrier_result=carrier_result,
                    hla_result=hla_result,
                    apoe_genotype=apoe_genotype,
                    inferred_sex=(qc_result or {}).get("inferred_sex"),
                    file_label=file_label,
                )
                log(f"  FHIR R4 bundle saved: {fhir_path} "
                    f"(PGx={summary['n_pgx_observations']}, "
                    f"Carrier={summary['n_carrier_observations']}, "
                    f"HLA={summary['n_hla_observations']}, "
                    f"APOE={'yes' if summary['apoe_included'] else 'no'})")
            except Exception as e:
                log(f"  WARNING: FHIR export failed: {e}")

    # ── v3: optional Compare vs previous run ──
    if args.compare:
        if diff_runs is None:
            log("  Compare skipped: module not loaded")
        else:
            prev = Path(args.compare)
            curr = SCRIPT_DIR / "tier1_results.json"
            if not prev.exists():
                log(f"  Compare: previous file not found: {prev}")
            else:
                log(f"Generating changelog ({prev} -> current) ...")
                d = diff_runs(prev, curr)
                print()
                print(render_diff_text(d))
                # Save JSON form alongside
                changelog_path = SCRIPT_DIR / "changelog.json"
                changelog_path.write_text(json.dumps(d, indent=2, default=str))
                log(f"  Changelog saved: {changelog_path}")

    log("=" * 60)
    log(f"Done! Open {output_path} in your browser.")
    if apoe_genotype:
        log(f"APOE Genotype:    {apoe_genotype}")
    log(f"Y-DNA Haplogroup: {y_result.get('haplogroup_path', 'Unknown')}")
    if mt_result:
        log(f"mtDNA Haplogroup: {mt_result.get('haplogroup', 'Unknown')}")
    risk_count = len([r for r in tier1_results if r["risk_copies"] > 0])
    log(f"Variants matched: {len(tier1_results)} | Risk alleles found: {risk_count}")
    log("=" * 60)

    # ── v3: optional Chat REPL ──
    if args.chat:
        if run_chat is None:
            log("Chat module not available")
        else:
            run_chat(
                model=args.model,
                tier1_path=SCRIPT_DIR / "tier1_results.json",
                pgx_result=pgx_result, prs_result=prs_result,
                traits_result=traits_result, carrier_result=carrier_result,
                interactions_result=interactions_result,
                ancestry_result=ancestry_result,
                bloodwork_result=bloodwork_result,
                detox_result=detox_result,
                urologic_result=urologic_result,
                deep_ancestry_result=deep_ancestry_result,
                blood_type_result=blood_type_result,
                immunogenetics_result=immunogenetics_result,
                neurochemistry_result=neurochemistry_result,
                addiction_genetics_result=addiction_genetics_result,
                holistic_synthesis_result=holistic_synthesis_result,
                clinical_variants_result=clinical_variants_result,
                novel_variants_result=novel_variants_result,
                voi_result=voi_result,
                polygenic_traits_result=polygenic_traits_result,
                environmental_optimization_result=environmental_optimization_result,
                family_planning_result=family_planning_result,
                economics_result=economics_result,
                y_result=y_result,
                mt_result=mt_result,
            )

    return 0
