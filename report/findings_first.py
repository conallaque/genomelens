"""The findings-first economics report — pages 1 and 2.

WHAT CHANGED AND WHY. The previous economics output opened with pooled
cost-effectiveness machinery and reached the genomic findings that produced it
about forty pages later. A reader had to work through the sophistication before
learning what was actually found in the genome. These pages invert that: the
findings come first, the economics explains their consequence, and the pooled
reference case arrives as the conclusion rather than the premise.

EVERY NUMBER COMES FROM THE PAYLOAD. Nothing here computes, derives or rounds a
quantity of its own — values arrive from ``EconomicsReportPayload`` and are
turned into strings by ``report.format``. That is the whole point of the
canonical payload: a figure on the page and the same figure in the JSON cannot
disagree, because there is only one of them.

TWO VISUAL RULES THE OLD OUTPUT BROKE.

*Confidence and economic direction are different dimensions.* Colour on the
evidence badge means evidence only, and nothing else on the page is tinted by
it. A finding is not shaded green for being high-confidence, because "well
established" and "good news" are not the same claim and a reader should never
have to guess which one a colour means.

*Cash and monetised health are never blended.* Expected net monetary benefit
includes health valued at a willingness-to-pay threshold and is not money
returned to anyone. It is labelled as such wherever it appears, and the words
"ROI", "return" and "savings" appear nowhere.
"""
from __future__ import annotations

import html

from report import format as fmt
from report.payload import EconomicsReportPayload, FindingEconomics

__all__ = ["CSS", "MAIN_REPORT_PAGES", "page_count",
           "render_findings_first"]

MAIN_REPORT_PAGES = 8
"""Sheets a typical genome renders to, and the floor for any genome.

Not a fixed total: the glance page continues onto further sheets when a genome
carries more findings than one holds. ``page_count`` reports the real number.
"""


def _e(s) -> str:
    return html.escape(str(s if s is not None else ""))


def _clip(s, n: int) -> str:
    """Truncate on a word boundary. Cutting mid-word ("published intervention
    ec") reads as a rendering fault rather than an editorial choice."""
    text = str(s or "").strip()
    if len(text) <= n:
        return text
    return text[:n].rsplit(" ", 1)[0].rstrip(" ,;:-—") + "…"


# ── design system ─────────────────────────────────────────────────────────────
#
# Deep teal for structure and headings, one mint tint for grouped blocks, and
# three badge colours that mean evidence and nothing else. Generous whitespace
# over density: an extra page costs less than an unreadable one.

CSS = """
:root{
  --teal-900:#0b3d47; --teal-700:#0f6b62; --teal-600:#12857a;
  --ink:#12303a; --ink-2:#5a6b74; --ink-3:#93a1aa;
  --mint:#f1f8f6; --mint-line:#d7e8e3;
  --amber:#a07818; --amber-bg:#fdf7ea; --amber-line:#ecdcb8;
  --rose:#a8475c; --rose-bg:#fdf0f2; --rose-line:#f0d3d9;
  --green:#127a55; --red:#a33a30;
  --line:#e4ecea; --paper:#ffffff;
}
*{box-sizing:border-box}
html,body{margin:0;padding:0}
body{color:var(--ink);background:var(--paper);
  font:14.5px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Inter,Roboto,
       "Helvetica Neue",Arial,sans-serif;-webkit-font-smoothing:antialiased}

.sheet{position:relative;width:960px;height:1240px;margin:0 auto;
  padding:0 62px 74px;display:flex;flex-direction:column;overflow:hidden}
.sheet + .sheet{break-before:page;page-break-before:always}
.topbar{position:absolute;left:0;right:0;top:0;height:9px;
  background:var(--teal-900)}

/* masthead */
.mast{padding-top:38px;display:flex;justify-content:space-between;
  align-items:flex-start}
.eyebrow{font-size:9.5px;letter-spacing:.17em;text-transform:uppercase;
  font-weight:700;color:var(--ink-3)}
.seclabel{font-size:16.5px;color:var(--teal-700);margin-top:3px;font-weight:500}
.flag{font-size:9px;font-weight:800;letter-spacing:.11em;text-transform:uppercase;
  color:var(--amber);background:var(--amber-bg);border:1px solid var(--amber-line);
  border-radius:20px;padding:4px 11px;white-space:nowrap;margin-top:4px}

h1.head{font-size:30px;line-height:1.16;letter-spacing:-.02em;font-weight:750;
  color:var(--teal-900);margin:19px 0 0;max-width:24em}
p.intro{margin:11px 0 0;color:var(--ink-2);font-size:13px;line-height:1.6;
  max-width:64em}

h2.sec{font-size:20px;font-weight:700;color:var(--teal-900);margin:0;
  letter-spacing:-.012em}
.sechead{display:flex;justify-content:space-between;align-items:baseline;
  gap:20px;margin-top:21px}
.sechead .aside{font-size:11.5px;color:var(--ink-2);text-align:right;
  max-width:27em;line-height:1.55}
.grouplabel{font-size:10px;letter-spacing:.13em;text-transform:uppercase;
  font-weight:800;color:var(--teal-700);margin:26px 0 0}
.grouplabel span{color:var(--ink-3);font-weight:600;letter-spacing:.02em;
  text-transform:none;font-size:11px;margin-left:7px}

/* metric cards */
.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:13px;
  margin-top:16px}
.mcard{border:1px solid var(--mint-line);background:var(--mint);
  border-radius:13px;padding:15px 16px 13px;break-inside:avoid;
  page-break-inside:avoid;display:flex;flex-direction:column}
/* 2.9em, not 2.6. At 9.5px/1.4 a two-line label is 26.6px while the old
   reserve was 24.7px, so "Healthcare cost change" overflowed it and sat
   its value 2px below the other three. The reserve exists to line the
   four figures up, so it has to clear the tallest label. */
.mcard .lbl{font-size:9.5px;letter-spacing:.11em;text-transform:uppercase;
  font-weight:800;color:var(--ink-3);line-height:1.4;min-height:2.9em}
.mcard .val{font-size:29px;font-weight:750;letter-spacing:-.03em;
  line-height:1.08;color:var(--teal-900);margin-top:2px}
.mcard .val.pos{color:var(--green)}
.mcard .val.neg{color:var(--red)}
.mcard .val.sm{font-size:19px;line-height:1.2;letter-spacing:-.02em}
.mcard .unit{font-size:10.5px;font-weight:600;color:var(--ink-2);
  line-height:1.35;margin-top:5px;letter-spacing:.005em}
.mcard .cap{font-size:11px;color:var(--ink-2);line-height:1.5;margin-top:8px}

/* breadth banner */
.banner{margin-top:17px;border:1px solid var(--mint-line);background:var(--mint);
  border-radius:14px;padding:15px 19px;display:flex;gap:17px;
  break-inside:avoid;page-break-inside:avoid}
.banner .count{flex:0 0 auto;width:58px;height:58px;border-radius:50%;
  background:var(--teal-700);color:#fff;display:flex;align-items:center;
  justify-content:center;font-size:24px;font-weight:750}
.banner .btxt{flex:1}
.banner .bh{font-size:14.5px;font-weight:700;color:var(--teal-700);
  line-height:1.4}
.banner .bs{font-size:11.5px;color:var(--ink-2);margin-top:3px;line-height:1.55}
.pills{display:flex;flex-wrap:wrap;gap:7px;margin-top:12px}
.pill{font-size:10.5px;background:#fff;border:1px solid var(--mint-line);
  border-radius:20px;padding:5px 12px;color:var(--ink-2);white-space:nowrap;
  line-height:1.3}

/* top-finding cards */
.tops{display:grid;grid-template-columns:repeat(4,1fr);gap:13px;margin-top:13px}
.tcard{border:1px solid var(--line);border-radius:13px;padding:13px 14px;
  background:var(--paper);break-inside:avoid;page-break-inside:avoid;
  display:flex;flex-direction:column}
.tcard .gene{font-size:9.5px;letter-spacing:.1em;text-transform:uppercase;
  font-weight:800;color:var(--ink-3);line-height:1.4;min-height:2.5em}
.tcard .amt{font-size:26px;font-weight:750;letter-spacing:-.03em;
  color:var(--teal-700);line-height:1.1;margin-top:3px}
.tcard .amt.neg{color:var(--red)}
.tcard .act{font-size:11px;color:var(--ink-2);line-height:1.5;margin-top:8px;
  flex:1}
.tcard .bwrap{margin-top:11px}

/* evidence badge — evidence only */
.badge{display:inline-block;font-size:9px;font-weight:800;letter-spacing:.1em;
  text-transform:uppercase;border-radius:20px;padding:4px 10px;line-height:1;
  white-space:nowrap}
.badge.high{color:var(--teal-600);background:#e8f5f1;border:1px solid #c9e6de}
.badge.moderate{color:var(--amber);background:var(--amber-bg);
  border:1px solid var(--amber-line)}
.badge.low{color:var(--rose);background:var(--rose-bg);
  border:1px solid var(--rose-line)}
.badge.none{color:var(--ink-3);background:#f3f6f6;border:1px solid #e4eaea}

/* waterfall */
.flow{margin-top:14px;border:1px solid var(--line);border-radius:14px;
  padding:16px 20px;break-inside:avoid;page-break-inside:avoid}
.flow h3{margin:0;font-size:15px;font-weight:700;color:var(--teal-900)}
.flow p.sub{margin:5px 0 0;font-size:11.5px;color:var(--ink-2);line-height:1.55}
.steps{display:flex;align-items:stretch;gap:0;margin-top:12px}
.step{flex:1;text-align:center;padding:0 6px}
.step .sl{font-size:9px;letter-spacing:.09em;text-transform:uppercase;
  font-weight:800;color:var(--ink-3);line-height:1.4;min-height:2.6em}
.step .sv{font-size:21px;font-weight:750;letter-spacing:-.025em;
  color:var(--teal-900);margin-top:3px}
.step .sn{font-size:10px;color:var(--ink-2);margin-top:5px;line-height:1.45}
.arrow{flex:0 0 auto;display:flex;flex-direction:column;
  align-items:center;justify-content:flex-start;padding:14px 2px 0;min-width:84px}
.arrow .amt{font-size:11.5px;font-weight:750;color:var(--red);white-space:nowrap}
.arrow .gl{font-size:9px;color:var(--ink-3);margin-top:3px;text-align:center;
  line-height:1.3}
.arrow .bar{width:100%;height:1.5px;background:var(--mint-line);margin-top:8px;
  position:relative}
.arrow .bar:after{content:"";position:absolute;right:-1px;top:-3px;
  border-left:6px solid var(--mint-line);border-top:4px solid transparent;
  border-bottom:4px solid transparent}
.flow .tail{margin-top:11px;padding-top:10px;border-top:1px solid var(--line);
  display:flex;justify-content:space-between;align-items:baseline;gap:24px}
.flow .tail .k{font-size:12px;color:var(--ink-2);line-height:1.55}
.flow .tail .v{font-size:21px;font-weight:750;color:var(--teal-900);
  letter-spacing:-.02em;white-space:nowrap}

/* info boxes */
.info{margin-top:14px;border:1px solid #cfe2ea;background:#f2f9fb;
  border-radius:12px;padding:15px 18px;font-size:12px;color:#1c4653;
  line-height:1.62;break-inside:avoid;page-break-inside:avoid}
.info b{color:var(--teal-900)}
.info + .info{margin-top:10px}

/* page-2 finding rows */
.rows{margin-top:6px}
.frow{display:grid;grid-template-columns:23% 1fr 132px 106px;gap:20px;
  padding:17px 4px;border-bottom:1px solid var(--line);
  break-inside:avoid;page-break-inside:avoid;align-items:start}
.frow:last-child{border-bottom:none}
.frow .fname{font-size:13.5px;font-weight:700;line-height:1.35;
  color:var(--teal-900)}
.frow .fsub{font-size:10.5px;color:var(--ink-3);margin-top:3px;line-height:1.4}
.frow .act{font-size:12.5px;color:var(--ink);line-height:1.55}
.frow .caveat{font-size:10.5px;color:var(--ink-3);margin-top:4px;
  line-height:1.45;font-style:italic}
.frow .basis{font-size:10.5px;color:var(--ink-3);margin-top:7px;line-height:1.5}
.frow .money{text-align:right}
.frow .money .v{font-size:20px;font-weight:750;letter-spacing:-.025em;
  color:var(--teal-700);line-height:1.15}
.frow .money .v.neg{color:var(--red)}
.frow .money .v.na{font-size:12.5px;font-weight:650;color:var(--ink-3);
  line-height:1.35}
.frow .money .vn{font-size:9px;letter-spacing:.09em;text-transform:uppercase;
  color:var(--ink-3);font-weight:800;margin-bottom:4px}
.frow .ev{text-align:left;padding-top:2px}

.foot{position:absolute;left:62px;right:62px;bottom:30px;padding-top:14px;
  border-top:1px solid var(--line);display:flex;justify-content:space-between;
  font-size:10px;color:var(--ink-3);background:var(--paper)}

/* The sheet is designed at 960x1240 CSS px, which is letter proportions at a
   comfortable reading size. Declaring the page box in those same units makes
   one .sheet map to exactly one printed page; leaving it as `size:letter`
   (816x1056 at 96dpi) silently split every sheet across two pages. */

/* glance rows (page 2) */
.glance{margin-top:6px}
.grow{display:grid;grid-template-columns:1.35fr 1fr 92px 112px;gap:14px;
  padding:10px 4px;border-bottom:1px solid var(--line);align-items:center;
  break-inside:avoid}
.grow:last-child{border-bottom:none}
.grow .n{font-size:12.5px;font-weight:700;color:var(--teal-900);line-height:1.35}
.grow .d{font-size:11.5px;color:var(--ink-2);line-height:1.4}
.grow .m{text-align:right;font-size:17px;font-weight:750;letter-spacing:-.025em;
  color:var(--teal-700);line-height:1.2}
.grow .m.neg{color:var(--red)}
.grow .m.na{font-size:11px;font-weight:650;color:var(--ink-3);line-height:1.35}

/* detail cards (pages 3-4) */
.dcard{border:1px solid var(--line);border-radius:13px;padding:11px 16px;
  margin-top:8px;break-inside:avoid}
.dcard .dh{display:flex;justify-content:space-between;align-items:baseline;gap:14px}
.dcard .dt{font-size:15px;font-weight:750;color:var(--teal-900);
  letter-spacing:-.015em;line-height:1.3}
.dcard .dsub{font-size:10px;color:var(--ink-3);margin-top:2px}
.dgrid{display:grid;grid-template-columns:1fr 226px;gap:20px;margin-top:10px}
.dfield{margin-top:5px}
.dfield:first-child{margin-top:0}
.dfield .fl{font-size:8.5px;letter-spacing:.11em;text-transform:uppercase;
  font-weight:800;color:var(--ink-3);margin-bottom:2px}
.dfield .fv{font-size:11px;color:var(--ink);line-height:1.5}
.econbox{border:1px solid var(--mint-line);background:var(--mint);
  border-radius:11px;padding:10px 13px}
.econbox .eh{font-size:8.5px;letter-spacing:.11em;text-transform:uppercase;
  font-weight:800;color:var(--ink-3)}
.econbox .ev{font-size:22px;font-weight:750;color:var(--teal-700);
  letter-spacing:-.03em;line-height:1.1;margin:2px 0 7px}
.econbox .ev.neg{color:var(--red)}
.econbox .ev.na{font-size:12px;line-height:1.3;color:var(--ink-3)}
.econbox dl{margin:0;display:grid;grid-template-columns:1fr auto;gap:3px 10px}
.econbox dt{font-size:10px;color:var(--ink-2)}
.econbox dd{margin:0;font-size:10px;font-weight:700;text-align:right;
  color:var(--ink);font-variant-numeric:tabular-nums}
.srcline{font-size:8.5px;color:var(--ink-3);margin-top:7px;line-height:1.45;
  padding-top:6px;border-top:1px solid var(--line)}

/* compact detail variant (page 4) — economics inline instead of boxed, so
   every finding fits one page without dropping any of them */
.dcard.c .dgrid{display:block}
.dcard.c .estrip{display:flex;gap:20px;align-items:baseline;margin-top:7px;
  padding-top:7px;border-top:1px solid var(--line);flex-wrap:wrap}
.dcard.c .estrip .ev{font-size:20px;font-weight:750;color:var(--teal-700);
  letter-spacing:-.03em;line-height:1}
.dcard.c .estrip .ev.neg{color:var(--red)}
.dcard.c .estrip .ev.na{font-size:12px;color:var(--ink-3);font-weight:650}
.dcard.c .estrip .es{font-size:10px;color:var(--ink-2)}
.cover{background:#fff}
.cv-wrap{padding:150px 84px 0 84px;display:flex;flex-direction:column}
.cv-brand{font-size:12px;letter-spacing:.16em;text-transform:uppercase;
  color:var(--accent,#0f6b5c);font-weight:700;margin-bottom:26px}
.cv-title{font-size:42px;line-height:1.12;letter-spacing:-.02em;margin:0 0 14px;
  font-weight:700;color:var(--ink,#12181f);max-width:15ch}
.cv-sub{font-size:15px;line-height:1.55;color:var(--ink-2,#4a5560);margin:0;
  max-width:52ch}
.cv-meta{margin-top:44px;border-top:1px solid var(--rule,#e3e6ea);
  border-bottom:1px solid var(--rule,#e3e6ea);padding:18px 0;
  display:grid;grid-template-columns:1fr 1fr;gap:11px 40px}
.cv-meta>div{display:flex;justify-content:space-between;gap:14px;font-size:12px}
.cv-meta span{color:var(--ink-3,#8a8f98);letter-spacing:.03em}
.cv-meta b{color:var(--ink,#12181f);text-align:right}
.cv-meta .mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
  font-size:10px;font-weight:500}
.cv-note{margin-top:30px;border-left:3px solid var(--accent,#0f6b5c);
  background:#f2f8f6;padding:13px 16px;font-size:12.5px;line-height:1.6;
  color:var(--ink-2,#4a5560);border-radius:0 4px 4px 0}
.cv-note.warn{border-left-color:#b8860b;background:#fdf6e7}
.cv-note b{color:var(--ink,#12181f);display:block;margin-bottom:3px}
.cv-disc{margin-top:16px;border:1px solid var(--rule,#e3e6ea);border-radius:4px;
  padding:14px 16px;font-size:11.5px;line-height:1.62;
  color:var(--ink-2,#4a5560)}
.cv-disc b{color:var(--ink,#12181f);display:block;margin-bottom:3px}
.dflt{font-size:9px;letter-spacing:.02em;color:var(--ink-3,#8a8f98);border:1px solid var(--rule,#e3e6ea);border-radius:3px;padding:0 3px;margin-left:3px;white-space:nowrap}
.dcard.c .estrip .es b{color:var(--ink);font-weight:700;
  font-variant-numeric:tabular-nums}
.dcard.c .estrip .el{font-size:8.5px;letter-spacing:.11em;text-transform:uppercase;
  font-weight:800;color:var(--ink-3);margin-right:-14px}

/* compact tables */
table.t{width:100%;border-collapse:collapse;margin-top:9px;font-size:10.5px}
table.t thead th{font-size:8px;letter-spacing:.09em;text-transform:uppercase;
  color:var(--ink-3);font-weight:800;text-align:left;padding:0 6px 5px;
  border-bottom:1.5px solid var(--mint-line);white-space:nowrap}
table.t td{padding:7px 6px;border-bottom:1px solid var(--line);
  vertical-align:top;line-height:1.4}
table.t tbody tr:last-child td{border-bottom:none}
table.t td.n{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
table.t td.b{font-weight:700;color:var(--teal-900)}

/* split panels */
.split{display:grid;grid-template-columns:1fr 1fr;gap:15px;margin-top:11px}
.split.wide{grid-template-columns:1.2fr 1fr}
.panel{border:1px solid var(--line);border-radius:13px;padding:14px 16px;
  break-inside:avoid}
.panel.mint{border-color:var(--mint-line);background:var(--mint)}
.panel h4{margin:0 0 3px;font-size:8.5px;letter-spacing:.12em;
  text-transform:uppercase;font-weight:800;color:var(--teal-700)}
.panel .ph{font-size:13.5px;font-weight:700;color:var(--teal-900);
  line-height:1.35;margin-bottom:5px}
.panel p{margin:0;font-size:10.5px;color:var(--ink-2);line-height:1.55}
.panel dl{margin:8px 0 0;display:grid;grid-template-columns:1fr auto;
  gap:5px 12px;align-items:baseline}
.panel dt{font-size:10.5px;color:var(--ink-2);line-height:1.4}
.panel dd{margin:0;font-size:13px;font-weight:750;text-align:right;
  color:var(--teal-900);white-space:nowrap}

/* interval bar */
.ibar{margin-top:10px;position:relative;height:54px}
.ibar .track{position:absolute;left:0;right:0;top:21px;height:9px;
  background:linear-gradient(90deg,#e6f2ee,#bfe0d6);border-radius:5px}
.ibar .zero{position:absolute;top:11px;width:2px;height:29px;background:var(--red)}
.ibar .tick{position:absolute;top:11px;width:1.5px;height:29px;
  background:var(--teal-700)}
.ibar .lab{position:absolute;font-size:9px;color:var(--ink-2);
  white-space:nowrap;transform:translateX(-50%)}
.ibar .lab.top{top:0} .ibar .lab.bot{top:42px}

/* tornado */
.tor{margin-top:8px}
.torow{display:grid;grid-template-columns:180px 1fr 72px;gap:11px;
  align-items:center;padding:3px 0;font-size:10px}
.torow .tn{color:var(--ink);line-height:1.3}
.torow .tn em{font-style:normal;color:var(--ink-3);font-size:8.5px;display:block;
  letter-spacing:.06em;text-transform:uppercase;font-weight:700}
.torow .tb{height:12px;background:var(--mint-line);border-radius:3px}
.torow .tb.assum{background:#f0dcc0}
.torow .tv{text-align:right;font-variant-numeric:tabular-nums;
  color:var(--ink-2);font-size:9.5px}

.checks{margin-top:8px;display:grid;grid-template-columns:1fr 1fr;gap:4px 18px}
.chk{font-size:10px;color:var(--ink-2);line-height:1.5}
.chk b{color:var(--green);font-weight:800;letter-spacing:.06em;font-size:8.5px;
  margin-right:6px;text-transform:uppercase}
ul.lim{margin:7px 0 0;padding-left:16px}
ul.lim li{font-size:10.5px;color:var(--ink-2);line-height:1.5;margin-bottom:3px}
.info.warn{border-color:var(--amber-line);background:var(--amber-bg);color:#6b5220}
h2.sec.lg{font-size:19px;letter-spacing:-.012em;text-transform:none;
  color:var(--teal-900);font-weight:700}
.metrics.three{grid-template-columns:repeat(3,1fr)}

@page{size:960px 1240px;margin:0}
@media print{ .sheet{height:1240px} }
"""

_SYNTH = ("Synthetic whole-genome input · illustrative model · "
          "not clinical or financial advice")


def _badge(conf: str) -> str:
    c = (conf or "").strip().lower()
    if c not in ("high", "moderate", "low"):
        return '<span class="badge none">no grade</span>'
    return f'<span class="badge {c}">{c}</span>'


def _mast(p: EconomicsReportPayload, section: str) -> str:
    flag = ('<div class="flag">synthetic input</div>'
            if p.metadata.is_synthetic else "")
    return f"""<div class="topbar"></div>
  <div class="mast">
    <div><div class="eyebrow">GenomeLens economics</div>
      <div class="seclabel">{_e(section)}</div></div>
    {flag}
  </div>"""


def _foot(p: EconomicsReportPayload, n: int,
          total: int = MAIN_REPORT_PAGES) -> str:
    note = (_SYNTH if p.metadata.is_synthetic
            else "Illustrative model · not clinical or financial advice")
    return (f'<div class="foot"><div>{_e(note)}</div>'
            f'<div>{n} / {total}</div></div>')


# ── page 1 ────────────────────────────────────────────────────────────────────

def _verdict(p: EconomicsReportPayload) -> str:
    """The headline, generated from the result rather than asserted.

    Never "your genome is worth $X": that turns a modelled expectation into a
    possession, and the quantity it names includes monetised health, which is
    not money anyone receives.
    """
    r = p.reference_case
    cheaper, healthier = r.incremental_cost < 0, r.incremental_qalys > 0
    if cheaper and healthier:
        return ("Acting on these findings is modelled to improve health and "
                "reduce healthcare spending.")
    if healthier:
        return ("Acting on these findings is modelled to improve health at "
                "additional healthcare cost.")
    if cheaper:
        return ("Acting on these findings is modelled to reduce healthcare "
                "spending without a measurable health gain.")
    return ("Acting on these findings is not modelled to improve health or "
            "reduce healthcare spending.")


def _confidence_mix(priced: list[FindingEconomics]) -> str:
    """Grade breakdown of the findings the caller counted in its headline.

    This used to take the whole payload and count every finding in it, while
    the headline above it counted only the priced ones — 6 high + 9 moderate +
    8 low against a headline of 22, plus 30 "ungraded" that were mostly the
    action-string phantoms. Two populations described by one sentence, with no
    way for a reader to reconcile them.

    It now takes THE LIST THE HEADLINE COUNTED. Sharing the list is what makes
    the totals agree; an assertion that they agree would only have restated the
    bug's absence, and this function cannot see the headline to check against.
    """
    counts: dict[str, int] = {}
    for f in priced:
        key = (f.evidence_confidence or "ungraded").lower()
        counts[key] = counts.get(key, 0) + 1
    order = {"high": 0, "moderate": 1, "low": 2, "ungraded": 3}
    parts = [f"{n} {c}" if c == "ungraded" else f"{n} {c}-confidence"
             for c, n in sorted(counts.items(), key=lambda kv: order.get(kv[0], 9))]
    return (", ".join(parts) + ".") if parts else ""



def render_page_one(p: EconomicsReportPayload) -> str:
    r = p.reference_case
    priced = [f for f in p.findings if f.canonical_expected_nmb is not None]
    top = sorted(priced, key=lambda f: -(f.canonical_expected_nmb or 0))[:4]

    # BOUNDED. This rendered one pill per finding with no limit, so a genome
    # with forty-three of them produced a 470px block that pushed "What rises to
    # the top" off the bottom of the cover page — the four headline figures the
    # page exists to show. The pills are a texture of what was found, not a
    # list to be read: past a couple of dozen they stop being scannable and
    # start being a wall. The remainder is counted, and every finding is listed
    # in full on the glance pages.
    _MAX_PILLS = 22
    _shown = p.findings[:_MAX_PILLS]
    pills = "".join(f'<span class="pill">{_e(f.short_name)}</span>'
                    for f in _shown)
    if len(p.findings) > _MAX_PILLS:
        pills += (f'<span class="pill">and {len(p.findings) - _MAX_PILLS} '
                  f'more</span>')

    tops = ""
    for f in top:
        v = f.canonical_expected_nmb or 0.0
        tops += f"""
      <div class="tcard">
        <div class="gene">{_e(f.short_name)}</div>
        <div class="amt{' neg' if v < 0 else ''}">{fmt.money(v)}</div>
        <div class="act">{_e(_clip(f.action_summary or f.category
                                   or 'no action recorded', 70))}</div>
        <div class="bwrap">{_badge(f.evidence_confidence)}</div>
      </div>"""

    u = p.uncertainty
    if u.psa_available and u.psa_iterations:
        unc_val = fmt.simulation_count(
            u.probability_cost_effective, u.psa_iterations
        ).replace(" simulations", "").replace("all ", "")
        # THE INTERVAL SITS BESIDE THE PROBABILITY, DELIBERATELY.
        #
        # A profile whose whole distribution lies above zero reports the same
        # 100% that a bug reports when the finding-level parameters are pinned
        # outside the sampling loop — and this report leads with that bug as a
        # self-caught error. The two are distinguishable only by whether the
        # spread is real, so the spread is printed in the same block rather
        # than elsewhere on the page: 100% next to a $9k-$43k interval reads as
        # a confident model, 100% next to no interval reads as a broken one.
        # `test_ceac_at_zero_threshold_is_below_certainty` is the machine
        # version of the same check.
        _ci = ""
        if u.nmb_ci_low or u.nmb_ci_high:
            _ci = (f" \u00b7 95% interval {fmt.money(u.nmb_ci_low)} to "
                   f"{fmt.money(u.nmb_ci_high)} across "
                   f"{fmt.count(u.n_parameters_varied)} varied parameters")
        # TWO SHARES OF ONE SET OF DRAWS, NOT A PARTITION. This used to read
        # "1,287 of 1,500 ... (86%); cost-saving in 1,241 of 1,500 (83%)" —
        # two same-denominator fractions separated by a semicolon, which a
        # reader naturally adds to 169% and treats as broken. They are the
        # same 1,500 simulations measured against two different bars:
        # cost-effective is INMB > 0 at the threshold, cost-saving is a
        # negative incremental cost, and the cost-saving draws sit inside the
        # cost-effective ones wherever the QALY gain is non-negative. Saying
        # "of the same draws" and naming the overlap is the whole fix; the
        # numbers were never wrong.
        _cs_n = round(u.probability_cost_saving * u.psa_iterations)
        unc_cap = (f"modelled simulations cost-effective "
                   f"({fmt.probability(u.probability_cost_effective)}); "
                   f"{fmt.count(_cs_n)} of the same draws are also cost-saving "
                   f"({fmt.probability(u.probability_cost_saving)}) "
                   f"\u2014 overlapping shares, not a split"
                   f"{_ci}")
    else:
        unc_val, unc_cap = fmt.MISSING, "no probabilistic analysis available"

    n_priced, n_total = len(priced), len(p.findings)
    unstd_line = (f" {n_total - n_priced} of the {n_total} findings shown are "
                  f"reported without a dollar figure rather than being assigned "
                  f"one the model cannot support." if n_total > n_priced else "")

    return f"""
<section class="sheet">
  {_mast(p, "Findings-first economics")}
  <h1 class="head">{_e(_verdict(p))}</h1>
  <p class="intro">This report starts with what was found in the genome. The
  economics explains the consequence of acting on each finding; the pooled
  reference case shows what remains after overlapping signals, imperfect
  adherence and the real cost of acting are all accounted for.</p>

  <div class="metrics">
    <div class="mcard"><div class="lbl">Health gain</div>
      <div class="val pos">~{fmt.healthy_days(r.incremental_qalys)}</div>
      <div class="unit">{fmt.qaly(r.incremental_qalys)} incremental QALYs</div>
      <div class="cap">of healthy life, on average, for a person with this
        pattern of findings</div></div>
    <div class="mcard"><div class="lbl">Healthcare cost change</div>
      <div class="val {'pos' if r.incremental_cost < 0 else 'neg'}"
        >{fmt.signed_money(r.incremental_cost)}</div>
      <div class="unit">healthcare-sector perspective</div>
      <div class="cap">negative means the modelled programme spends less than
        usual care</div></div>
    <div class="mcard"><div class="lbl">Modelled net benefit</div>
      <div class="val">{fmt.money(r.nmb)}</div>
      <div class="unit">reference-case NMB</div>
      <div class="cap">at {fmt.money(r.wtp)} per QALY. This values health at
        the threshold; it is not cash.</div></div>
    <div class="mcard"><div class="lbl">Decision uncertainty</div>
      <div class="val sm">{_e(unc_val)}</div>
      <div class="cap">{_e(unc_cap)}</div></div>
  </div>

  <div class="banner">
    <div class="count">{n_priced}</div>
    <div class="btxt">
      <div class="bh">findings carry a standardised economic estimate</div>
      <div class="bs">{_e(_confidence_mix(priced))}{_e(unstd_line)}</div>
      <div class="pills">{pills}</div>
    </div>
  </div>

  <div class="sechead">
    <h2 class="sec">What rises to the top</h2>
    <div class="aside">Ranked by expected net monetary benefit, each valued on
      its own.<br>Finding-level values are not additive; the pooled reference
      case is calculated separately.</div>
  </div>
  <div class="tops">{tops}</div>


  {_foot(p, 1)}
</section>"""


# ── page 2 ────────────────────────────────────────────────────────────────────

def _money_cell(f: FindingEconomics) -> str:
    if f.canonical_expected_nmb is None:
        label = ("Not routed to the<br>economic engine" if not f.is_monetized
                 else "No costed<br>pathway yet")
        return f'<div class="vn">Expected NMB</div><div class="v na">{label}</div>'
    v = f.canonical_expected_nmb
    return (f'<div class="vn">Expected NMB</div>'
            f'<div class="v{" neg" if v < 0 else ""}">{fmt.money(v)}</div>')


def _basis_line(f: FindingEconomics) -> str:
    """The secondary economic line. Deliberate absences read as choices."""
    if f.canonical_expected_nmb is None:
        if not f.is_monetized:
            if "reproduct" in (f.display_name + " " + f.category).lower():
                return ("Reproductive outcomes are deliberately not assigned a "
                        "dollar benefit.")
            return (_clip(f.reason_not_monetized, 96)
                    or "Deliberately not assigned a dollar benefit.")
        return "No registry-backed pathway yet."
    bits = []
    if f.medical_cost_averted:
        bits.append(f"{fmt.money(f.medical_cost_averted)} medical cost averted")
    if f.expected_qaly_gain:
        bits.append(f"{fmt.qaly(f.expected_qaly_gain)} QALY")
    if f.intervention_cost:
        bits.append(f"{fmt.money(f.intervention_cost)} to act")
    return " · ".join(bits) or "registry-backed expected NMB"


_GENERIC_CATEGORIES = {"pharmacogenomic / genomic", "genomic", "", "other"}


def _sub_label(f: FindingEconomics) -> str:
    """A truthful qualifier under the finding name, or nothing.

    `category` cannot be trusted here: analyze_personal_economics labels every
    record it reads from findings_with_economics "Pharmacogenomic / genomic",
    so APOE and PTPN22 both arrive tagged as pharmacogenomics. Printing that
    would put a wrong word under a correct number. The pathway id is
    structured and does not lie, so the drug comes from there; anything else
    gets no sub-label at all.
    """
    pid = f.economic_pathway_id or ""
    if pid.startswith("pgx:"):
        parts = pid.split(":")
        if len(parts) >= 3 and parts[2]:
            return parts[2].replace("_", " ")
    cat = (f.category or "").strip()
    return "" if cat.lower() in _GENERIC_CATEGORIES else _clip(cat, 42)


def _drug_of(f: FindingEconomics) -> str:
    pid = f.economic_pathway_id or ""
    if pid.startswith("pgx:"):
        parts = pid.split(":")
        if len(parts) >= 3 and parts[2]:
            return parts[2].replace("_", " ")
    return ""


def _decision_label(f: FindingEconomics) -> str:
    """The short decision label for the at-a-glance page.

    Page 2 scans; pages 3 and 4 explain. A full action sentence and an economic
    basis line on every row is what pushed this page onto a second sheet, so
    here it is the decision in a few words.
    """
    drug = _drug_of(f)
    if drug:
        return f"{drug.split()[0].capitalize()} prescribing"
    act = (f.action_summary or "").strip()
    if act:
        for sep in (" \u2014 ", " / ", " ("):
            if sep in act:
                act = act.split(sep, 1)[0].strip()
                break
        return _clip(act, 46)
    return _clip(f.category, 40) or "\u2014"


# How many finding rows fit on one glance sheet. The sheet is a fixed
# 960x1240px block with an absolutely-positioned footer, so content past this
# does not push the page taller — it runs underneath the footer and the page
# number. A genome with twenty findings printed its last rows with "Synthetic
# whole-genome input" struck through them, and the three closing notes ran 162px
# off the bottom entirely.
#
# Measured rather than guessed: rows render at ~50.6px, and the closing notes
# plus footer need ~250px, leaving room for twelve rows on the first sheet and
# a full sheet's worth on any continuation.
#
# OVERFLOW CONTINUES ONTO ANOTHER SHEET, it is not dropped. Truncating would
# have been easier and would have meant the summary page silently stopped
# listing findings once a genome had enough of them — on the page whose whole
# job is to show what was found. The report is 8 pages for a typical genome and
# longer for one carrying more; that is the correct direction for a document
# whose length should reflect its input.
_GLANCE_ROWS_FIRST_SHEET = 12
_GLANCE_ROWS_PER_SHEET = 16


def _glance_row(f) -> str:
    """One scannable row: finding, decision, evidence grade, standalone value."""
    if f.canonical_expected_nmb is None:
        cell = ('<div class="m na">Not routed to the engine</div>'
                if not f.is_monetized
                else '<div class="m na">Not yet standardised</div>')
    else:
        v = f.canonical_expected_nmb
        cell = (f'<div class="m{" neg" if v < 0 else ""}">{fmt.money(v)}</div>')
    return f"""
      <div class="grow">
        <div class="n">{_e(_clip(f.display_name, 72))}</div>
        <div class="d">{_e(_decision_label(f))}</div>
        <div>{_badge(f.evidence_confidence)}</div>
        {cell}
      </div>"""


def render_page_two(p: EconomicsReportPayload) -> str:
    """At a glance. Scanning, not explanation.

    Emits one sheet, or more when the genome carries more findings than a sheet
    holds. Every finding appears; none is dropped to make the page fit.
    """
    # Flatten to (group, finding) so a group can straddle a sheet boundary and
    # be re-labelled "continued" rather than forcing an early break.
    flat: list[tuple[str, object]] = []
    for group, items in p.findings_page_groups():
        for f in items:
            flat.append((group, f))

    # EVEN SPLIT, not fixed budgets. Filling each sheet to its cap and letting
    # the remainder fall onto the last one left a final sheet holding five rows
    # and 55% white space once the phantom action-rows were removed and the
    # count dropped from 53 findings to 33. The caps still decide HOW MANY
    # sheets are needed; the rows are then spread across those sheets so the
    # last one is never nearly empty.
    n = len(flat)
    if n <= _GLANCE_ROWS_FIRST_SHEET:
        n_sheets = 1
    else:
        rest = n - _GLANCE_ROWS_FIRST_SHEET
        n_sheets = 1 + -(-rest // _GLANCE_ROWS_PER_SHEET)   # ceil

    pages: list[list[tuple[str, object]]] = []
    i = 0
    for sheet in range(n_sheets):
        # Spread what remains over the sheets that remain, capped so the first
        # sheet (which carries the page heading) never exceeds its own budget.
        left = n_sheets - sheet
        take = -(-(n - i) // left)
        cap = _GLANCE_ROWS_FIRST_SHEET if sheet == 0 else _GLANCE_ROWS_PER_SHEET
        take = min(take, cap)
        pages.append(flat[i:i + take])
        i += take
    if i < n:                     # never drop a finding to make the page fit
        pages[-1].extend(flat[i:])
    if not pages:
        pages = [[]]

    n_unstd = sum(1 for f in p.findings if f.canonical_expected_nmb is None)
    unstd = ""
    if n_unstd:
        word = "One finding is" if n_unstd == 1 else f"{n_unstd} findings are"
        unstd = (f'<div class="info">{word} shown without a standardised '
                 f'expected NMB rather than being assigned a value the model '
                 f'cannot support.</div>')

    out = ""
    seen_groups: set[str] = set()
    for idx, chunk in enumerate(pages):
        rows, cur = "", None
        for group, f in chunk:
            if group != cur:
                if cur is not None:
                    rows += "</div>"
                n_in_group = sum(1 for g, _ in flat if g == group)
                n_here = sum(1 for g, _ in chunk if g == group)
                cont = " (continued)" if group in seen_groups else ""
                # THE REMAINDER, NOT THE TOTAL. A continuation header repeated
                # the group's full count — "Risk & prevention, 12 findings" on
                # the sheet showing 12 and again on the sheet showing the last
                # 4 — so the same twelve looked like twenty-four.
                if cont:
                    shown = sum(1 for pg in pages[:idx] for g, _ in pg
                                if g == group)
                    label = (f'{n_here} of {n_in_group} finding'
                             f'{"" if n_in_group == 1 else "s"}'
                             f' &middot; {shown} shown earlier')
                else:
                    label = (f'{n_in_group} finding'
                             f'{"" if n_in_group == 1 else "s"}')
                # WHY THESE FIGURES ARE IDENTICAL. On a real genome five
                # medication findings came back at exactly $210 each, because
                # PGX_CEA holds five gene-drug pairs and none of these was one
                # of them, so all five fell through to the generic
                # actionable-PGx flag. The figures were right; five identical
                # numbers in a column with nothing to explain them read as a
                # collapse bug. Said once per group rather than as a badge on
                # every row — the rows share one cause, not five.
                _gf = sum(1 for g, x in chunk if g == group
                          and getattr(x, "pgx_basis", "") == "generic_fallback")
                if _gf:
                    label += ('&nbsp;&middot; all priced off the generic PGx flag'
                              if _gf == n_here else
                              f'&nbsp;&middot; {_gf} priced off the generic '
                              f'PGx flag')
                rows += (f'<div class="grouplabel">{_e(group)}{cont}'
                         f'<span>{label}</span></div>'
                         f'<div class="glance">')
                seen_groups.add(group)
                cur = group
            rows += _glance_row(f)
        if cur is not None:
            rows += "</div>"

        last = idx == len(pages) - 1
        tail = ""
        if last:
            tail = f"""
  <div class="info" style="margin-top:16px">
    <b>Finding-level values are standalone expected values. They are not
    additive.</b> The reference-case total is calculated separately after
    overlapping signals are pooled by condition.</div>
  <div class="info">Expected NMB values health at
    {fmt.money(p.metadata.willingness_to_pay)} per QALY and is
    <b>not cash returned to anyone</b>. Ordering is by clinical group first and
    expected NMB second.</div>
  {unstd}"""
        head = ("Findings that could change a decision" if idx == 0
                else "Findings that could change a decision, continued")
        intro = ("""<p class="intro">Every finding the model examined, the
  decision it bears on, the strength of the evidence, and what it is worth on
  its own. No dollar figure appears without the genomic result that produced it
  &mdash; the detail behind each figure is on the pages that follow.</p>"""
                 if idx == 0 else "")
        out += f"""
<section class="sheet">
  {_mast(p, "Your findings at a glance")}
  <h1 class="head">{head}</h1>
  {intro}
  {rows}{tail}
  {_foot(p, 2 + idx)}
</section>"""
    return out


# ── pages 3 and 4: the detail ────────────────────────────────────────────────

def _econ_box(f: FindingEconomics) -> str:
    if f.canonical_expected_nmb is None:
        head = ("Not routed to the economic engine" if not f.is_monetized
                else "Economic pathway not yet standardised")
        why = _clip(f.reason_not_monetized, 62) or "no registry-backed pathway"
        return f"""<div class="econbox">
        <div class="eh">Expected NMB</div>
        <div class="ev na">{_e(head)}</div>
        <dl><dt style="font-size:9.5px">{_e(why)}</dt><dd></dd></dl></div>"""
    v = f.canonical_expected_nmb
    return f"""<div class="econbox">
        <div class="eh">Expected NMB</div>
        <div class="ev{' neg' if v < 0 else ''}">{fmt.money(v)}\
{_pgx_basis_note(f)}</div>
        <dl>
          <dt>Medical cost averted</dt><dd>{fmt.money(f.medical_cost_averted)}</dd>
          <dt>Health gain</dt><dd>{fmt.qaly(f.expected_qaly_gain)} QALY</dd>
          <dt>Cost to act</dt>
          <dd>{fmt.money(f.intervention_cost)}{_cost_basis_note(f)}</dd>
        </dl></div>"""


# The medication sheet bounds its CARD COUNT but never bounded its TEXT, so it
# sat 24px clear of the footer and adding ~60 characters to each of the three
# PGx action strings pushed it 9px past — a bleed onto the next page's
# masthead, which is the defect this document has already shipped once. The
# title was clipped at 62; the field values were not clipped at all.
_DETAIL_FIELD_CHARS = 230


def _detail_card(f: FindingEconomics, *, fields: list, source: str) -> str:
    body = "".join(
        f'<div class="dfield"><div class="fl">{_e(k)}</div>'
        f'<div class="fv">{_e(_clip(str(v), _DETAIL_FIELD_CHARS))}</div></div>'
        for k, v in fields if v)
    sub = _drug_of(f)
    return f"""
  <div class="dcard">
    <div class="dh">
      <div><div class="dt">{_e(_clip(f.display_name, 62))}</div>
        {f'<div class="dsub">{_e(sub)}</div>' if sub else ''}</div>
      {_badge(f.evidence_confidence)}
    </div>
    <div class="dgrid"><div>{body}</div>{_econ_box(f)}</div>
    <div class="srcline">{_e(source)}</div>
  </div>"""


_PGX_WHY = {
    "clopidogrel": ("Clopidogrel is a prodrug that CYP2C19 must activate. "
                    "Reduced activity means less active drug and a higher "
                    "chance the antiplatelet effect is inadequate."),
    "codeine": ("CYP2D6 activity governs conversion of codeine to morphine and "
                "clearance of several antidepressants, so one dose can be "
                "ineffective in some people and excessive in others."),
    "simvastatin": ("SLCO1B1 transports statins into the liver. Reduced "
                    "function raises systemic exposure and with it the risk of "
                    "muscle toxicity, most clearly for simvastatin."),
}


def render_page_three(p: EconomicsReportPayload) -> str:
    pgx = dict(p.findings_page_groups()).get("Medication & prescribing", [])[:3]
    cards = ""
    for f in pgx:
        drug = _drug_of(f)
        why = next((v for k, v in _PGX_WHY.items() if k in drug.lower()),
                   "Genotype affects how this medication is processed, so the "
                   "same dose does not have the same effect in everyone.")
        cards += _detail_card(f, fields=[
            ("What was found",
             f"Genotype predicts altered handling of "
             f"{drug or 'the relevant medication'}."),
            ("Why it matters", why),
            ("What could change",
             f.action_summary or "The prescribing decision this genotype "
                                 "informs, when it next arises."),
            ("Key dependency",
             f"Value accrues only if {drug or 'the relevant drug'} is actually "
             f"prescribed. The model prices that exposure explicitly rather "
             f"than assuming it."),
        ], source="Expected NMB = exposure probability x adverse-event "
                  "probability x relative risk reduction x cost of the avoided "
                  "event, less the cost of acting. Provenance on the methods page.")
    if not cards:
        cards = ('<div class="info">This run produced no medication-genotype '
                 'findings with a standardised economic pathway.</div>')
    return f"""
<section class="sheet">
  {_mast(p, "Medication-genotype findings")}
  <h1 class="head">Where a genotype could change a prescribing decision</h1>
  <p class="intro">These are the clearest cases of a genomic result bearing on
  a future medication choice. The clinical decision leads; the economics
  explains why preserving the genotype for that decision has value.</p>
  {cards}
  <div class="info warn"><b>Decision support, not a directive.</b> Nothing here
    is an instruction to start, stop or change a medication. Relevance depends
    on whether the drug is actually prescribed, and the prescribing decision
    belongs with a qualified clinician who can see the whole clinical
    picture.</div>
  {_foot(p, 3)}
</section>"""


def _estrip(f: FindingEconomics) -> str:
    if f.canonical_expected_nmb is None:
        label = ("Not routed to the economic engine" if not f.is_monetized
                 else "Economic pathway not yet standardised")
        why = _clip(f.reason_not_monetized, 70)
        return (f'<div class="estrip"><span class="el">Expected NMB</span>'
                f'<span class="ev na">{_e(label)}</span>'
                f'<span class="es">{_e(why)}</span></div>')
    v = f.canonical_expected_nmb
    return (f'<div class="estrip"><span class="el">Expected NMB</span>'
            f'<span class="ev{" neg" if v < 0 else ""}">{fmt.money(v)}</span>'
            f'<span class="es">medical cost averted '
            f'<b>{fmt.money(f.medical_cost_averted)}</b></span>'
            f'<span class="es">health gain '
            f'<b>{fmt.qaly(f.expected_qaly_gain)} QALY</b></span>'
            f'<span class="es">cost to act '
            f'<b>{fmt.money(f.intervention_cost)}</b>'
            f'{_cost_basis_note(f)}</span></div>')


def _cost_basis_note(f: FindingEconomics) -> str:
    """Mark a cost of acting that is the registered default, not this gene's.

    Eleven findings printed the same $500 with nothing to distinguish them, so
    risk-reducing surgery and passive symptom awareness read as separately
    costed at an identical price. They are one declared parameter used eleven
    times, and saying so is the difference between a default and a claim.
    """
    if f.intervention_cost_basis and f.intervention_cost_basis != "gene_specific":
        return ' <span class="dflt">registry default</span>'
    return ""


def _pgx_basis_note(f: FindingEconomics) -> str:
    """Mark a PGx figure that came from the generic flag, not a pair-specific CEA.

    Five findings on a real genome printed the same value because PGX_CEA holds
    five gene-drug pairs and none of them was one of these. Unlabelled that is
    indistinguishable from a bug; labelled it is a stated limit of the model.
    """
    if f.pgx_basis == "generic_fallback":
        return ' <span class="dflt">generic PGx flag</span>'
    return ""


def _compact_card(f: FindingEconomics, fields: list) -> str:
    body = "".join(
        f'<div class="dfield"><div class="fl">{_e(k)}</div>'
        f'<div class="fv">{_e(v)}</div></div>' for k, v in fields if v)
    return f"""
  <div class="dcard c">
    <div class="dh">
      <div class="dt">{_e(_clip(f.display_name, 66))}</div>
      {_badge(f.evidence_confidence)}
    </div>
    <div class="dgrid"><div>{body}</div></div>
    {_estrip(f)}
  </div>"""


def render_page_four(p: EconomicsReportPayload) -> str:
    groups = dict(p.findings_page_groups())
    items: list[FindingEconomics] = []
    for name in ("Risk & prevention", "Lower-confidence & exploratory",
                 "Reported, not costed"):
        items.extend(groups.get(name, []))
    # Detail cards are for findings the model can actually explain. The ones
    # with no standardised pathway have nothing to decompose, so a card would
    # be four empty fields; they are named in a closing line instead and appear
    # in full on page 2.
    priced = [f for f in items if f.canonical_expected_nmb is not None]
    unpriced = [f for f in items if f.canonical_expected_nmb is None]

    cards = ""
    generic = 0
    for f in priced[:5]:
        fields = [("Potential decision",
                   f.action_summary or "No standardised action pathway "
                                       "recorded.")]
        if f.action_caveat:
            fields.append(("Key assumption", f.action_caveat))
        else:
            generic += 1
        cards += _compact_card(f, fields)
    shared = ""
    if generic:
        shared = ('<p class="intro" style="font-size:10px;margin-top:9px">'
                  "Where no finding-specific caveat is shown: baseline event "
                  "probability and the effect of acting are registered "
                  "assumptions varied in the sensitivity analysis, and neither "
                  "is measured for this individual.</p>")

    tail = ""
    if unpriced:
        # BOUNDED. This listed every unpriced finding by name, so a genome with
        # twenty of them produced a paragraph that ran off the bottom of the
        # sheet and under the footer. The count is what carries the meaning —
        # "none of these was invented a value" — and the names are recoverable
        # in full on the glance pages, so the list is capped and the remainder
        # counted rather than the sentence being allowed to grow without limit.
        _SHOW = 6
        names = ", ".join(_clip(f.display_name, 40) for f in unpriced[:_SHOW])
        more = len(unpriced) - _SHOW
        if more > 0:
            names += f", and {more} more"
        tail = (f'<div class="info" style="margin-top:11px"><b>Reported without '
                f'a standardised value ({len(unpriced)}):</b> {_e(names)}. No '
                f'registry-backed economic pathway exists for these yet, so '
                f'none is invented. Each appears with its decision and evidence '
                f'grade on the glance pages.</div>')

    return f"""
<section class="sheet">
  {_mast(p, "Risk, prevention & exploratory findings")}
  <h1 class="head">Risk and prevention findings, in detail</h1>
  <p class="intro">The findings above carry an evidence grade rather than a
  single standard of proof: this page runs from well-established variants to
  exploratory ones, and the badge on each card says which is which. Lower
  confidence changes how a finding should be read, not whether it is shown. A
  negative expected NMB means the modelled cost of acting exceeds the modelled
  benefit &mdash; useful information, not an error. Evidence grade reflects the
  strength of the association, not the size of the value; provenance is on the
  methods page.</p>
  {cards}{shared}{tail}
  {_foot(p, 4)}
</section>"""


# ── page 5 ────────────────────────────────────────────────────────────────────

def render_page_five(p: EconomicsReportPayload) -> str:
    r, c = p.reference_case, p.corrections
    rows = ""
    for cond in sorted(p.condition_results, key=lambda x: -x.nmb):
        rows += f"""
      <tr><td class="b">{_e(cond.condition)}</td>
        <td class="n">{cond.n_contributing_findings}</td>
        <td class="n">{fmt.percentage(cond.baseline_risk, places=1)}</td>
        <td class="n">{fmt.percentage(cond.naive_additive_rrr, places=0)}</td>
        <td class="n">{fmt.percentage(cond.pooled_efficacy_rrr, places=0)}</td>
        <td class="n">{fmt.percentage(cond.adherence, places=0)}</td>
        <td class="n">{fmt.percentage(cond.adherence_adjusted_rrr, places=0)}</td>
        <td class="n">{fmt.money(cond.cost_averted)}</td>
        <td class="n">{fmt.qaly(cond.qaly_gain)}</td>
        <td class="n b">{fmt.money(cond.nmb)}</td></tr>"""
    # THE LABEL IS THE FIX, NOT SUPPRESSION. The Raw column adds the risk
    # reductions of findings that share a condition, which for a well-populated
    # pool exceeds 100% — 290% here, 475% before the duplicate paths were
    # collapsed. That is not a result the model reports; it is the artifact the
    # pooling correction exists to remove, and hiding it would leave the
    # correction unevidenced. Naming it is what turns the number into an
    # argument for the pooling logic instead of against it.
    raw_note = ("<p class=\"note\"><strong>Raw is an uncorrected sum, not a "
                "modelled result.</strong> It adds the risk reductions of "
                "findings that share a condition, which is not how "
                "probabilities combine. Any value above 100% is the artifact "
                "this correction removes. The pooled figure beside it is the "
                "model&rsquo;s answer.</p>")
    overlap = c.naive_cost_averted - c.efficacy_cost_averted
    adher = c.efficacy_cost_averted - c.effectiveness_cost_averted
    pct = (overlap / c.naive_cost_averted) if c.naive_cost_averted else 0.0
    apct = (adher / c.efficacy_cost_averted) if c.efficacy_cost_averted else 0.0
    # A second naive-to-realised walkthrough used to render here. It made the
    # same trip in one combined step; the "two
    # corrections" block below walks the same $60,291 to the same $28,824 and
    # additionally splits the reduction into its structural (pooling) and
    # behavioural (adherence) halves. Showing both put 1,293px on a 1,240px
    # sheet and printed the same journey twice, once less informatively.
    return f"""
<section class="sheet">
  {_mast(p, "How the findings combine")}
  <h1 class="head">Taken together, what do these findings mean?</h1>
  <p class="intro">This is the canonical reference case. Findings bearing on the
  same condition are combined on the risk scale, each condition's cost of
  illness is charged exactly once, and trial efficacy is discounted to
  real-world adherence before any value is credited.</p>
  <div class="metrics">
    <div class="mcard"><div class="lbl">Incremental healthcare cost</div>
      <div class="val {'pos' if r.incremental_cost < 0 else 'neg'}"
        >{fmt.signed_money(r.incremental_cost)}</div>
      <div class="cap">healthcare-sector perspective</div></div>
    <div class="mcard"><div class="lbl">Incremental QALYs</div>
      <div class="val pos">{fmt.qaly(r.incremental_qalys)}</div>
      <div class="cap">about {fmt.healthy_days(r.incremental_qalys)} of healthy
        life on average</div></div>
    <div class="mcard"><div class="lbl">Reference-case NMB</div>
      <div class="val">{fmt.money(r.nmb)}</div>
      <div class="cap">at {fmt.money(r.wtp)} per QALY</div></div>
    <div class="mcard"><div class="lbl">Classification</div>
      <div class="val sm">{_e(_clip(r.dominance_status, 38))}</div>
      <div class="cap">{_e(_clip(r.icer_note, 62))}</div></div>
  </div>
  <h2 class="sec">Condition pools &mdash; combined before value is credited</h2>
  <table class="t">
    <thead><tr><th>Condition</th><th style="text-align:right">Findings</th>
      <th style="text-align:right">Baseline</th>
      <th style="text-align:right">Raw<br><span class="thsub">uncorrected</span></th>
      <th style="text-align:right">Pooled</th><th style="text-align:right">Adherence</th>
      <th style="text-align:right">Realised</th><th style="text-align:right">Cost averted</th>
      <th style="text-align:right">QALYs</th><th style="text-align:right">NMB</th>
    </tr></thead><tbody>{rows}</tbody></table>
  {raw_note}
  <div class="flow">
    <h3>Two corrections, applied in order</h3>
    <p class="sub">The engine reports these as one combined reduction. They rest
      on different evidence and are separated here: overlap is a structural
      correction, adherence a behavioural one.</p>
    <div class="steps">
      <div class="step"><div class="sl">Naive avoided cost</div>
        <div class="sv">{fmt.money(c.naive_cost_averted)}</div>
        <div class="sn">every finding counted independently</div></div>
      <div class="arrow"><div class="amt">&minus;{fmt.money(overlap)}</div>
        <div class="gl">condition<br>pooling</div><div class="bar"></div></div>
      <div class="step"><div class="sl">Pooled, at trial efficacy</div>
        <div class="sv">{fmt.money(c.efficacy_cost_averted)}</div>
        <div class="sn">{fmt.percentage(pct, places=0)} removed &mdash; one
          liability, one cost of illness</div></div>
      <div class="arrow"><div class="amt">&minus;{fmt.money(adher)}</div>
        <div class="gl">real-world<br>adherence</div><div class="bar"></div></div>
      <div class="step"><div class="sl">Realised benefit</div>
        <div class="sv">{fmt.money(c.effectiveness_cost_averted)}</div>
        <div class="sn">{fmt.percentage(apct, places=0)} further removed
          &mdash; efficacy is not effectiveness</div></div>
    </div>
    <div class="tail"><div class="k">Adherence is charged to the benefit
      <em>and</em> to ongoing intervention cost, but not to the one-off test
      cost. {fmt.qaly(c.adherence_qaly_reduction)} QALYs are lost to imperfect
      adherence relative to trial efficacy.</div>
      <div class="v">{fmt.money(r.nmb)}</div></div>
  </div>
  {_foot(p, 5)}
</section>"""


# ── page 6 ────────────────────────────────────────────────────────────────────

def _ceac_at_zero_line(p: EconomicsReportPayload) -> str:
    """State the zero-threshold probability explicitly, whatever it is.

    At a willingness-to-pay of $0 a strategy is cost-effective only if it is
    cost-SAVING in that draw, so this is the number that distinguishes a
    genuinely dominant profile from one whose parameters stopped varying. It
    was computed and never printed, which meant the page showed the reassuring
    figure and withheld the diagnostic one.
    """
    pts = [c for c in (p.uncertainty.ceac or [])
           if float(c.get("wtp", c.get("lam", -1))) == 0]
    if not pts:
        return ""
    p0 = float(pts[0].get("p_cost_effective", pts[0].get("prob", 0.0)))
    if p0 >= 1.0:
        return ("At a $0 threshold this is still 100%, which on a dominant "
                "profile is expected rather than reassuring — read it with "
                "the interval above, not on its own.")
    return (f"At a $0 threshold it is {p0:.2%}, below certainty — the cash "
            f"arm genuinely varies.")


def _wgs_zero_note(p: EconomicsReportPayload) -> str:
    """Say why the observed contribution is zero, rather than printing a bare $0.

    A $0 beside "what did sequencing actually find that an array could not"
    reads as a failed calculation unless the page says otherwise. It is a real
    answer here: nothing in this genome is reportable ONLY by sequencing.
    """
    t = p.testing_decision
    if t.observed_wgs_only_findings:
        return ""
    return ('<p style="margin-top:7px"><b>Zero is an answer, not a gap.</b> '
            "Every finding in this genome sits at a position an array could "
            "have typed, so none of them is reportable only by sequencing. "
            "That makes the observed contribution genuinely nil for this "
            "input &mdash; it does not mean sequencing is worthless, which is "
            "what the prospective column beside it measures.</p>")


def _interval_bar(p: EconomicsReportPayload) -> str:
    u = p.uncertainty
    lo, hi, mean = u.nmb_ci_low, u.nmb_ci_high, u.psa_mean_nmb
    if hi <= lo:
        return ""
    span_lo, span_hi = min(0.0, lo), hi
    width = (span_hi - span_lo) or 1.0

    def x(v: float) -> float:
        return max(1.0, min(99.0, (v - span_lo) / width * 100.0))

    return f"""
  <div class="ibar">
    <div class="track"></div>
    <div class="zero" style="left:{x(0.0):.1f}%"></div>
    <div class="lab top" style="left:{x(0.0):.1f}%">$0 break-even</div>
    <div class="tick" style="left:{x(lo):.1f}%"></div>
    <div class="lab bot" style="left:{x(lo):.1f}%">5th pct {fmt.money(lo)}</div>
    <div class="tick" style="left:{x(mean):.1f}%"></div>
    <div class="lab bot" style="left:{x(mean):.1f}%">mean {fmt.money(mean)}</div>
    <div class="tick" style="left:{x(hi):.1f}%"></div>
    <div class="lab bot" style="left:{x(hi):.1f}%">95th pct {fmt.money(hi)}</div>
  </div>"""


def render_page_six(p: EconomicsReportPayload) -> str:
    u, pr = p.uncertainty, p.provenance
    tor = sorted(u.tornado, key=lambda t: -float(t.get("swing") or 0))[:6]
    mx = max((float(t.get("swing") or 0) for t in tor), default=1.0) or 1.0
    tor_rows = ""
    for t in tor:
        swing = float(t.get("swing") or 0)
        tier = str(t.get("tier") or "")
        tor_rows += f"""
      <div class="torow">
        <div class="tn">{_e(str(t.get('parameter', '')).replace('_', ' '))}
          <em>{_e(tier)}</em></div>
        <div class="tb{' assum' if tier == 'assumption' else ''}"
          style="width:{swing / mx * 100:.0f}%"></div>
        <div class="tv">{fmt.money(swing)}</div>
      </div>"""
    errs = sum(1 for f in p.report_validation if f.get("severity") == "ERROR")
    warns = sum(1 for f in p.report_validation if f.get("severity") == "WARNING")
    infos = sum(1 for f in p.report_validation if f.get("severity") == "INFO")
    plain = p.plain_language or {}
    share = plain.get("assumption_share")
    share_txt = (f"; they account for "
                 f"{fmt.percentage(float(share) / 100, places=1)} of total swing"
                 if share else "")
    return f"""
<section class="sheet">
  {_mast(p, "Uncertainty, evidence & model discipline")}
  <h1 class="head">How sure is this, and what rests on judgement?</h1>
  <p class="intro">Uncertainty has two layers: the spread of modelled outcomes,
  and the provenance of the inputs producing that spread. A robust decision does
  not make its inputs precise, so both are shown.</p>
  <h2 class="sec">Net monetary benefit across
    {fmt.count(u.psa_iterations)} simulations</h2>
  {_interval_bar(p)}
  <div class="split" style="margin-top:4px">
    <div class="panel mint"><h4>Cost-effective</h4>
      <div class="ph">{_e(fmt.simulation_count(u.probability_cost_effective,
                                               u.psa_iterations))}</div>
      <p>at {fmt.money(p.metadata.willingness_to_pay)} per QALY
        ({fmt.probability(u.probability_cost_effective)}). Reference-case NMB is
        {fmt.money(p.reference_case.nmb)}; the mean across simulations is
        {fmt.money(u.psa_mean_nmb)}. Both are correct and answer slightly
        different questions.</p>
      <p style="margin-top:6px"><b>95% interval
        {fmt.money(u.nmb_ci_low)} to {fmt.money(u.nmb_ci_high)}</b> across
        {fmt.count(u.n_parameters_varied)} varied parameters. A probability
        this high is only meaningful beside the spread that produced it:
        pinned parameters also report certainty, and the way to tell them
        apart is whether the interval moves.</p></div>
    <div class="panel mint"><h4>Cost-saving</h4>
      <div class="ph">{_e(fmt.simulation_count(u.probability_cost_saving,
                                               u.psa_iterations))}</div>
      <p>the modelled programme costs less than usual care
        ({fmt.probability(u.probability_cost_saving)}).
        {_e(_clip(u.note, 130))}</p>
      <p style="margin-top:6px">{_e(_ceac_at_zero_line(p))}</p></div>
  </div>
  <h2 class="sec">What drives the answer</h2>
  <div class="tor">{tor_rows}</div>
  <p class="intro" style="font-size:10px;margin-top:5px">Swing in net monetary
    benefit across each parameter's plausible range. Amber bars are declared
    assumptions rather than sourced values{share_txt}.</p>
  <div class="split wide" style="margin-top:12px">
    <div class="panel"><h4>Evidence foundation</h4>
      <dl>
        <dt>Registry sourced</dt>
        <dd>{pr.registry_n_parameters - pr.registry_n_assumption} / {pr.registry_n_parameters}
          &middot; {fmt.percentage(pr.registry_pct_sourced / 100, places=1)}</dd>
        <dt>Whole model &mdash; PMID/DOI</dt>
        <dd>{fmt.percentage(pr.model_pct_resolvable / 100, places=1)}</dd>
        <dt>Whole model &mdash; attributed</dt>
        <dd>{fmt.percentage(pr.model_pct_attributed_or_better / 100, places=1)}</dd>
        <dt>Whole model &mdash; judgement only</dt>
        <dd>{fmt.percentage(pr.model_pct_unsourced / 100, places=1)}</dd>
      </dl>
      <p style="margin-top:7px">Two different denominators &mdash;
        {pr.registry_n_parameters} registered parameters and
        {fmt.count(pr.model_n_total_known)} known model parameters. Reported
        separately, because collapsing them into one figure would create an
        inconsistency rather than remove one.</p></div>
    <div class="panel"><h4>What could change the answer</h4>
      <p>{_e(_clip(str(plain.get('what_would_change_it') or ''), 300))}</p>
      <div class="checks" style="grid-template-columns:1fr;margin-top:9px">
        <div class="chk"><b>{errs} errors</b> report-consistency validation
          &mdash; rendering permitted</div>
        <div class="chk">{warns} warnings &middot; {infos} information notes,
          listed in full in the technical appendix</div>
      </div></div>
  </div>
  {_foot(p, 6)}
</section>"""


# ── page 7 ────────────────────────────────────────────────────────────────────

def render_page_seven(p: EconomicsReportPayload) -> str:
    t = p.testing_decision
    rows = ""
    for s in t.strategies:
        icer = s.get("icer_vs_previous")
        rows += f"""
      <tr><td class="b">{_e(s.get('name', ''))}</td>
        <td class="n">{fmt.money(s.get('cost'))}</td>
        <td class="n">{float(s.get('qaly') or 0):.3f}</td>
        <td class="n">{fmt.money(icer) + '/QALY' if icer else '&mdash;'}</td>
        <td>{_e(s.get('status', ''))}</td></tr>"""
    return f"""
<section class="sheet">
  {_mast(p, "The testing decision")}
  <h1 class="head">Is whole-genome sequencing worth it?</h1>
  <p class="intro">Two different questions live on this page and they are never
  added together. One is a choice made <em>before</em> testing, averaged over a
  population. The other counts what sequencing actually found in this genome.</p>
  <div class="split">
    <div class="panel"><h4>Prospective &mdash; before testing</h4>
      <div class="ph">What is the expected value of choosing sequencing rather
        than an array?</div>
      <dl>
        <dt>Chance an array misses an actionable finding</dt>
        <dd>1 in {t.prospective_number_needed_to_sequence or '&mdash;'}</dd>
        <dt>Value if you are that person</dt>
        <dd>{fmt.money(t.prospective_value_per_finding)}</dd>
        <dt>Expected value, averaged over everyone</dt>
        <dd>{fmt.money(t.prospective_gross_expected_value)}</dd>
        <dt>Extra cost over an array</dt>
        <dd>{fmt.money(t.incremental_chip_to_wgs_cost)}</dd>
        <dt>Net expected value</dt>
        <dd>{fmt.money(t.prospective_net_expected_value)}</dd>
      </dl></div>
    <div class="panel mint"><h4>Observed &mdash; in this genome</h4>
      <div class="ph">What did sequencing actually find that an array could
        not?</div>
      <dl>
        <dt>Sequencing-only findings present</dt>
        <dd>{t.observed_wgs_only_findings}</dd>
        <dt>Observed standalone contribution</dt>
        <dd>{fmt.money(t.observed_wgs_only_value)}</dd>
      </dl>
      <p style="margin-top:9px">A count of findings this genome contains, not an
        expectation. For this input the population question is already
        settled.</p>
      {_wgs_zero_note(p)}
      <p style="margin-top:7px"><b>Not the same as the value of sequencing.</b>
        This line counts only findings an array <i>cannot</i> report. Sequencing
        also recovers findings an array could have carried and did not, which is
        the larger effect &mdash; so the full chip-to-whole-genome difference is
        bigger than the figure above.</p></div>
  </div>
  <div class="info"><b>Separation rule.</b> The prospective figures describe the
    population this genome was drawn from; the observed figure describes the
    genome itself. They answer different decision questions, come from different
    inputs, and are never combined into one headline total.</div>
  <h2 class="sec">Strategy frontier</h2>
  <table class="t">
    <thead><tr><th>Strategy</th><th style="text-align:right">Cost</th>
      <th style="text-align:right">QALYs</th>
      <th style="text-align:right">ICER vs previous</th><th>Status</th>
    </tr></thead><tbody>{rows}</tbody></table>
  <p class="intro" style="font-size:10px;margin-top:7px">Frontier values
    illustrate the testing choice for a population and are not this genome's
    reference case. {_e(_clip(t.caveat, 230))}</p>
  {_foot(p, 7)}
</section>"""


# ── page 8 ────────────────────────────────────────────────────────────────────

_LIMITATIONS = (
    "Population-average parameters are applied to one genome. The model values "
    "expected consequences for a person with this pattern of findings, not a "
    "prediction for any individual.",
    "Estimates are illustrative and decision-analytic. This is not a formal "
    "economic evaluation or a submission to any authority.",
    "Where a condition-specific baseline risk is unavailable a generic "
    "registered assumption is used; those assumptions dominate the sensitivity "
    "analysis and are labelled as such.",
    "Costs are not inflation-normalised to a single price year.",
    "Findings without a registry-backed pathway are excluded from net monetary "
    "benefit rather than assigned an invented value.",
    "Clinically relevant findings require confirmation by an appropriate "
    "professional before anything follows from them.",
)

_CHECKS = (
    ("NMB identity", "\u03bb \u00d7 \u0394QALY \u2212 \u0394Cost reconciles to the reported total"),
    ("Dominance classification", "sign of cost and effect agrees with the label"),
    ("Condition cost charged once", "no pool exceeds its cost of illness"),
    ("Pooled risk cap", "combined risk reduction stays within its bound"),
    ("WGS basis separation", "observed and prospective never combined"),
    ("Budget-impact consistency", "peak burden is the maximum net spend"),
    ("Report validation gate", "no errors; rendering permitted"),
)


def render_page_eight(p: EconomicsReportPayload) -> str:
    m, pr = p.metadata, p.provenance
    errs = sum(1 for f in p.report_validation if f.get("severity") == "ERROR")
    chk = "".join(f'<div class="chk"><b>pass</b>{_e(k)} &mdash; {_e(v)}</div>'
                  for k, v in _CHECKS) if errs == 0 else (
        f'<div class="chk">{errs} error(s) &mdash; see validation output</div>')
    lims = "".join(f"<li>{_e(x)}</li>" for x in _LIMITATIONS)
    return f"""
<section class="sheet">
  {_mast(p, "Methods & provenance")}
  <h1 class="head">What was assumed, where it came from, what was checked</h1>
  <p class="intro">The pages before this one tell the story; this page is their
  audit trail. Full parameter tables, references and the advanced analyses
  follow in the technical appendix.</p>
  <div class="split">
    <div class="panel"><h4>Model</h4>
      <dl>
        <dt>Perspective</dt><dd style="font-size:11px">{_e(m.perspective)}</dd>
        <dt>Comparator</dt><dd style="font-size:11px">usual care</dd>
        <dt>Outcomes</dt><dd style="font-size:11px">cost, QALYs, NMB</dd>
        <dt>Time horizon</dt><dd style="font-size:11px">{fmt.years(m.analysis_horizon_years)}</dd>
        <dt>Discounting</dt><dd style="font-size:11px">{fmt.percentage(m.discount_rate, places=0)}, costs and QALYs</dd>
        <dt>Willingness to pay</dt><dd style="font-size:11px">{fmt.money(m.willingness_to_pay)}/QALY</dd>
        <dt>Currency</dt><dd style="font-size:11px">{_e(m.currency)}, United States</dd>
        <dt>Pooling</dt><dd style="font-size:11px">complement of products</dd>
        <dt>Adherence</dt><dd style="font-size:11px">benefit and ongoing cost</dd>
        <dt>Uncertainty</dt><dd style="font-size:11px">{fmt.count(p.uncertainty.psa_iterations)}-run PSA</dd>
      </dl></div>
    <div class="panel"><h4>Provenance</h4>
      <dl>
        <dt>Registered parameters</dt><dd>{pr.registry_n_parameters}</dd>
        <dt>Published</dt><dd>{pr.registry_n_published}</dd>
        <dt>Derived</dt><dd>{pr.registry_n_derived}</dd>
        <dt>Declared assumptions</dt><dd>{pr.registry_n_assumption}</dd>
        <dt>Registry sourced</dt><dd>{fmt.percentage(pr.registry_pct_sourced / 100, places=1)}</dd>
        <dt>Model, resolvable</dt><dd>{fmt.percentage(pr.model_pct_resolvable / 100, places=1)}</dd>
        <dt>Model, attributed</dt><dd>{fmt.percentage(pr.model_pct_attributed_or_better / 100, places=1)}</dd>
        <dt>Model, judgement only</dt><dd>{fmt.percentage(pr.model_pct_unsourced / 100, places=1)}</dd>
      </dl>
      <p style="margin-top:8px">A parameter declared as an assumption may not
        cite a source &mdash; the registry refuses to load if one does.</p></div>
  </div>
  <h2 class="sec">Declared limitations</h2>
  <ul class="lim">{lims}</ul>
  <h2 class="sec">Validation</h2>
  <div class="checks">{chk}</div>
  <div class="info" style="margin-top:12px">
    <b>Audit trail.</b> Every figure in this report is serialised in
    <code>economics-payload.json</code>, the canonical object these pages render
    from. The parameter registry, references, cross-path reconciliation and the
    advanced analyses follow in the technical appendix.
    {f'Build {_e(m.build_id)}.' if m.build_id else ''}
    {'Synthetic input &mdash; no human genome was used.' if m.is_synthetic else ''}
  </div>
  {_foot(p, 8)}
</section>"""


# ── document ──────────────────────────────────────────────────────────────────

def page_count(p: EconomicsReportPayload) -> int:
    """Sheets this payload actually renders to.

    MAIN_REPORT_PAGES is the floor, not the answer: the glance page continues
    onto extra sheets when a genome carries more findings than one holds, and
    _renumber_footers already stamps the true total. Returning the constant
    here gave the document two different page counts.
    """
    return render_findings_first(p).count('<section class="sheet')


def _renumber_footers(html: str) -> str:
    """Number the footers by the sheets that actually exist.

    Each page function stamps its own number, which was correct while the report
    was exactly eight fixed sheets. Page 2 now continues onto another sheet when
    a genome carries more findings than one holds, so the hardcoded numbers
    collide — two sheets both claiming "3 / 8" — and the total is wrong besides.
    Counting after assembly is the only place that knows how many sheets there
    turned out to be.
    """
    import re as _re
    # PREFIX MATCH. The cover is `<section class="sheet cover">`, so an exact
    # match on `class="sheet"` counted 10 sheets while 11 footers were being
    # numbered — the last one would have read "11 / 10".
    total = html.count('<section class="sheet')
    n = 0

    def _sub(m):
        nonlocal n
        n += 1
        return f'<div>{n} / {total}</div></div>'

    return _re.sub(r"<div>\d+ / \d+</div></div>", _sub, html)


def _cover_date(iso: str) -> str:
    """ISO timestamp as a readable date, or the raw string if it will not parse."""
    from datetime import datetime
    t = str(iso or "").strip()
    if not t:
        return "\u2014"
    try:
        return datetime.fromisoformat(t).strftime("%d %B %Y")
    except ValueError:
        return t[:10] or t


def render_cover(p: EconomicsReportPayload) -> str:
    """The title sheet.

    NOT `core.pdf_export._build_cover_page`. That emits flowing-A4 `.pdf-cover`
    markup, and putting it in front of a self-paginated 960x1240 document is
    exactly what produced the stray cover and blank trailing page recorded in
    `core/pdf_export.py` and `pipeline.py` — it was tried and removed. This is
    a native sheet using the document's own geometry.

    WHAT IT EXISTS TO SAY. That the input is synthetic, and that this is a
    research and educational model rather than medical or financial advice.
    Those were carried in 8pt footer text and a corner chip; they are the most
    load-bearing sentences in the document and belong on the title page.
    """
    m = p.metadata
    synthetic = bool(m.is_synthetic)
    src = _e(m.input_label or "genome")
    kind = "whole-genome" if (m.input_type or "").lower() == "wgs" else "array"
    build = _e((m.build_id or "").replace("GENOMELENS-BUILD:", "").strip())
    return f"""
<section class="sheet cover">
  <div class="topbar"></div>
  <div class="cv-wrap">
    <div class="cv-brand">GenomeLens</div>
    <h1 class="cv-title">Economic analysis of a genomic report</h1>
    <p class="cv-sub">Cost&ndash;utility model of acting on the findings in one
      {kind} input, valued at
      {fmt.money(m.willingness_to_pay)} per QALY from a
      {_e(m.perspective or 'healthcare sector')} perspective.</p>

    <div class="cv-meta">
      <div><span>Source input</span><b>{src}</b></div>
      <div><span>Input type</span><b>{_e((m.input_type or '').upper() or '&mdash;')}</b></div>
      <div><span>Generated</span><b>{_e(_cover_date(m.generated_at))}</b></div>
      <div><span>Time horizon</span><b>{fmt.years(m.analysis_horizon_years)}</b></div>
      <div><span>Discounting</span><b>{fmt.percentage(m.discount_rate)}, costs and QALYs</b></div>
      <div><span>Build</span><b class="mono">{build or '&mdash;'}</b></div>
    </div>

    <div class="cv-note{'' if synthetic else ' warn'}">
      <b>{'Synthetic input &mdash; no human genome was used.' if synthetic
          else 'Personal genomic input.'}</b>
      {'Every finding in this report was produced from a purpose-built '
       'synthetic genome generated by this repository. No personal genomic or '
       'health data was read, stored or published to produce it.'
       if synthetic else
       'This report was produced from a personal genomic file. Treat it as '
       'health information and share it accordingly.'}
    </div>

    <div class="cv-disc">
      <b>Research and educational use &mdash; not medical or financial advice.</b>
      This is a decision-analytic model, not a clinical or economic evaluation
      submitted to any authority. Its parameters are population averages
      applied to one genome: it values the expected consequences of acting on a
      pattern of findings, and does not predict what will happen to any
      individual. Monetary figures include health valued at a
      willingness-to-pay threshold and are not cash returned to anyone.
      Clinically relevant findings require confirmation by an appropriate
      professional &mdash; a physician, clinical pharmacist or board-certified
      genetic counsellor &mdash; before anything follows from them.
    </div>
  </div>
  {_foot(p, 1)}
</section>"""


def render_findings_first(p: EconomicsReportPayload) -> str:
    """The findings-first report. Eight sheets for a typical genome, more when
    one carries enough findings that the glance page continues — the committed
    sample renders eleven."""
    label = p.metadata.input_label or "genome"
    pages = "\n".join(r(p) for r in (
        render_cover,
        render_page_one, render_page_two, render_page_three, render_page_four,
        render_page_five, render_page_six, render_page_seven, render_page_eight))
    pages = _renumber_footers(pages)
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GenomeLens economics &mdash; {_e(label)}</title>
<style>{CSS}</style></head><body>
{pages}
</body></html>"""
