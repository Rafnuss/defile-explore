#!/usr/bin/env python3
"""Species comparison pages of the QA viewer: every taxon on one figure, one row each, for the
abundance trend, the passage dates, and the shift of the passage dates. It reads `data/explore/`
only; the figures are drawn in the browser (Plotly), so the filters (group, tier, reliability
class, name) and the sort redraw them without a rebuild.

What each figure draws, and from which field:

- Abundance trend: the change of the smooth total (`trend.annual.smooth`) from the first to the
  last year, its class `key_numbers.trend`; or every year's smooth total against the taxon's own
  geometric mean (heat map). The export has no interval on the change itself: the one drawn
  combines the first and last years' 95% bands as if independent, so it is wider than the true
  interval (the level's uncertainty is shared by both years and cancels in the ratio).
- Passage: the 10/50/90% passage dates of `key_numbers.passage` (the smooth season in the last
  year with a trend, else the last `KEY_YEARS` counted seasons pooled), with the 80% band of the
  median date (`trend.passage`); or the share of the season passing on each day, the last
  `KEY_YEARS` seasons of `season.share` pooled (heat map).
- Passage shift: the smooth season's 10/50/90% dates (`trend.passage_q`) in the last year minus
  the first, class `trend.passage`; the 80% band combines the two years' bands (`trend.passage`)
  as if independent.

Usage:
    python scripts/species_compare.py          # -> logs/viewer/compare.html
"""

import argparse
import json
import math
import os

import numpy as np
from explore_viewer import ROOT, wrap

from defile_explore.pipeline import KEY_YEARS

Z95, Z80 = 1.959964, 1.281552  # normal quantiles of the 95% and 80% two-sided bands
GROUPS = {  # catalogue group: colour (validated all-pairs, light)
    "raptors": "#eb6834",
    "passerines": "#2a78d6",
    "waterbirds": "#1baf7a",
    "other": "#4a3aa7",
}
CLASS_OF = {  # figure: the element of `reliability.ELEMENTS` its rows rest on
    "trend": "key_numbers.trend",
    "passage": "key_numbers.passage",
    "shift": "trend.passage",
}


def num(x):
    return None if x is None or (isinstance(x, float) and math.isnan(x)) else x


def element_class(sp: dict, element: str) -> tuple[str, list[str]]:
    """The class defileViz gives `element`, and the reason codes of its claim; `counts` for a field
    from the counts alone (no trend)."""
    rel = sp.get("reliability")
    if not rel or element not in rel["elements"]:
        return "counts", []
    claim = {"key_numbers.trend": "trend", "trend.passage": "season"}.get(element, "season")
    return rel["elements"][element], [r["code"] for r in rel[claim]["reasons"]]


def trend_row(sp: dict) -> dict | None:
    t = sp.get("trend")
    if not t:
        return None
    a0, a1 = t["annual"][0], t["annual"][-1]
    if not (a0["smooth"] and a1["smooth"]):
        return None
    log_ratio = math.log(a1["smooth"] / a0["smooth"])
    if a0["smooth_q2.5"] and a1["smooth_q2.5"]:  # rounded to 0 for the rarest
        se = math.hypot(
            *(math.log(a["smooth_q97.5"] / a["smooth_q2.5"]) / (2 * Z95) for a in (a0, a1))
        )
    else:
        se = None
    smooth = [a["smooth"] for a in t["annual"]]
    gm = math.exp(np.mean(np.log(np.maximum(smooth, 1e-9))))
    cls, why = element_class(sp, CLASS_OF["trend"])
    return {
        "first": t["first_year"],
        "last": t["last_year"],
        "log": log_ratio,
        "lo": None if se is None else log_ratio - Z95 * se,
        "hi": None if se is None else log_ratio + Z95 * se,
        "s0": a0["smooth"],
        "s1": a1["smooth"],
        "years": [a["year"] for a in t["annual"]],
        "rel": [round(math.log2(max(s, 1e-9) / gm), 3) for s in smooth],
        "cls": cls,
        "why": why,
    }


def passage_row(sp: dict) -> dict | None:
    p = sp["key_numbers"].get("passage")
    if not p or num(p["q50"]) is None:
        return None
    cls, why = (
        element_class(sp, CLASS_OF["passage"]) if p["source"] == "smooth" else ("counts", [])
    )
    out = {k: round(p[k], 1) for k in ("q10", "q50", "q90")} | {"source": p["source"]}
    if sp.get("trend"):
        last = sp["trend"]["passage"][-1]
        out |= {"q50_lo": last["lo"], "q50_hi": last["hi"]}
    s = sp["season"]
    years = np.asarray(s["years"])
    share = np.asarray(s["share"], float)[years > years.max() - KEY_YEARS]
    pooled = np.nan_to_num(share).sum(axis=0)
    if pooled.max() > 0:
        out |= {"doy": s["doy"], "density": np.round(pooled / pooled.max(), 3).tolist()}
    return out | {"cls": cls, "why": why}


def shift_row(sp: dict) -> dict | None:
    t = sp.get("trend")
    if not t:
        return None
    q0, q1 = t["passage_q"][0], t["passage_q"][-1]
    b0, b1 = t["passage"][0], t["passage"][-1]
    bands = {"q10": ("q10_lo", "q10_hi"), "q50": ("lo", "hi"), "q90": ("q90_lo", "q90_hi")}
    out = {"first": t["first_year"], "last": t["last_year"]}
    for q, (lo, hi) in bands.items():
        d = q1[q] - q0[q]
        se = math.hypot(*((b[hi] - b[lo]) / (2 * Z80) for b in (b0, b1)))
        out[q] = {"d": round(d, 1), "lo": round(d - Z80 * se, 1), "hi": round(d + Z80 * se, 1),
                  "from": q0[q], "to": q1[q]}  # fmt: skip
    cls, why = element_class(sp, CLASS_OF["shift"])
    return out | {"cls": cls, "why": why}


def rows(taxa: list, export: str, pages: set) -> list:
    out = []
    for tx in taxa:
        path = os.path.join(export, "species", f"{tx['taxon_id']}.json")
        if not os.path.exists(path):
            continue
        sp = json.load(open(path))
        row = {
            "id": tx["taxon_id"],
            "en": tx["english_name"],
            "fr": tx["french_name"],
            "sci": tx["scientific_name"],
            "group": tx.get("group") or "other",
            "tier": tx["tier"],
            "order": tx.get("taxon_order") or 0,
            "birds": tx.get("season_birds") or 0,
            "page": tx["taxon_id"] in pages,
            "trend": trend_row(sp),
            "passage": passage_row(sp),
            "shift": shift_row(sp),
        }
        if row["trend"] or row["passage"]:
            out.append(row)
    return out


CSS = """
.bar{position:sticky;top:37px;z-index:8;background:#fff;border-bottom:1px solid #ddd;
 padding:8px 16px;display:flex;flex-wrap:wrap;gap:6px 18px;align-items:center;font-size:13px}
.bar fieldset{border:0;margin:0;padding:0;display:flex;gap:8px;align-items:center}
.bar legend{float:left;font-weight:600;color:#555;margin-right:4px;font-size:12px}
.bar label{white-space:nowrap;cursor:pointer}
.bar input[type=search]{padding:3px 6px;border:1px solid #ccc;border-radius:4px;width:150px}
.sw{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:3px;vertical-align:-1px}
.ctl{display:flex;flex-wrap:wrap;gap:6px 14px;align-items:center;font-size:13px;margin:4px 8px}
.ctl select{font-size:13px;padding:2px 4px}
.ctl .seg button{border:1px solid #ccc;background:#fff;padding:3px 10px;font-size:12px;cursor:pointer}
.ctl .seg button:first-child{border-radius:4px 0 0 4px}
.ctl .seg button:last-child{border-radius:0 4px 4px 0}
.ctl .seg button+button{border-left:0}
.ctl .seg button.on{background:#4C78A8;color:#fff;border-color:#4C78A8}
.count{color:#777;font-size:12px}
.key{font-size:12px;color:#555;margin:2px 8px}
.key span{margin-right:14px;white-space:nowrap}
"""

BODY = """
<div class="bar">
 <fieldset id="f-group"><legend>Group</legend>__GROUPS__</fieldset>
 <fieldset id="f-tier"><legend>Tier</legend>
  <label><input type="checkbox" value="full" checked> full</label>
  <label><input type="checkbox" value="short"> short</label>
  <label><input type="checkbox" value="rare"> rare</label></fieldset>
 <fieldset id="f-cls"><legend>On defileViz</legend>
  <label><input type="checkbox" value="show" checked> shown</label>
  <label><input type="checkbox" value="caveat" checked> caveated</label>
  <label><input type="checkbox" value="hide"> hidden</label>
  <label><input type="checkbox" value="counts" checked> from counts (no trend)</label></fieldset>
 <fieldset><legend>Name</legend><input id="f-name" type="search" placeholder="English, French, Latin"></fieldset>
</div>
<section>
<h2>Species compared</h2>
<p class="note">One row per taxon. The filters above apply to all three figures; the class
filter applies to the class of the figure's own claim (trend, season). Click a row's mark to open
its page (full tier). Hidden rows are greyed, caveated ones drawn hollow.</p>

<h3 id="trend">1. Abundance trend</h3>
<div class="panel">
 <div class="ctl">
  <span class="seg" id="t-view"><button data-v="dots" class="on">Change</button><button data-v="heat">Trajectory</button></span>
  <span class="seg" id="t-unit"><button data-v="total" class="on">first to last year</button><button data-v="year">per year</button></span>
  <label>Sort <select id="t-sort">
   <option value="value">by change</option><option value="order">taxonomic</option>
   <option value="birds">by birds per season</option><option value="name">by name</option>
   <option value="group">by group, then change</option></select></label>
  <span class="count" id="t-count"></span>
 </div>
 <div class="key" id="t-key"></div>
 <div id="fig-trend"></div>
 <p class="note">Change of the smooth total (trend without the year's level and weather episodes),
 from the taxon's first trend year (1993 or 2007) to 2025. Interval: 95%, conservative (the two
 years' bands combined as if independent; the export carries no interval on the change itself).
 Trajectory: each year's smooth total against the taxon's geometric mean over its years.</p>
</div>

<h3 id="passage">2. Passage dates</h3>
<div class="panel">
 <div class="ctl">
  <span class="seg" id="p-view"><button data-v="range" class="on">10-50-90%</button><button data-v="heat">Season</button></span>
  <label>Sort <select id="p-sort">
   <option value="q50">by median date</option><option value="q10">by 10% date</option>
   <option value="q90">by 90% date</option><option value="span">by span (10-90%)</option>
   <option value="order">taxonomic</option><option value="birds">by birds per season</option>
   <option value="name">by name</option><option value="group">by group, then median</option></select></label>
  <span class="count" id="p-count"></span>
 </div>
 <div class="key" id="p-key"></div>
 <div id="fig-passage"></div>
 <p class="note">Bar: 10% to 90% of the season passed; dot: half passed; black whisker: 80% band
 of the median date (taxa with a trend). With a trend, the smooth season of the last year
 (2025); without, the last __KEY_YEARS__ counted seasons pooled. Season: each day's share of the
 season, the last __KEY_YEARS__ seasons pooled, scaled to the taxon's peak day.</p>
</div>

<h3 id="shift">3. Shift of the passage dates</h3>
<div class="panel">
 <div class="ctl">
  <span class="seg" id="s-q"><button data-v="q10">10%</button><button data-v="q50" class="on">median</button><button data-v="q90">90%</button><button data-v="all">all three</button></span>
  <span class="seg" id="s-unit"><button data-v="total" class="on">first to last year</button><button data-v="decade">per decade</button></span>
  <label>Sort <select id="s-sort">
   <option value="value">by shift</option><option value="q50">by median date</option>
   <option value="order">taxonomic</option><option value="birds">by birds per season</option>
   <option value="name">by name</option><option value="group">by group, then shift</option></select></label>
  <span class="count" id="s-count"></span>
 </div>
 <div class="key" id="s-key"></div>
 <div id="fig-shift"></div>
 <p class="note">Days by which the smooth season's 10%, median or 90% passage date moved from the
 taxon's first trend year to 2025 (positive: later). Interval: 80%, conservative (the two years'
 bands combined as if independent).</p>
</div>
</section>
<script>
const DATA = __DATA__;
const GROUPS = __GROUP_COLORS__;
__JS__
</script>
"""

JS = r"""
const ROW = 15, TOP = 50, BOTTOM = 40, LEFT = 190;
const CLS_NAME = {show: "shown", caveat: "caveated", hide: "hidden", counts: "from counts"};
const MONTHS = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
const state = {t: {view: "dots", unit: "total"}, p: {view: "range"}, s: {q: "q50", unit: "total"}};

function date(doy) {  // day of year (non-leap) -> "12 Sep"
  const d = new Date(Date.UTC(2025, 0, 1) + (Math.round(doy) - 1) * 864e5);
  return d.getUTCDate() + " " + MONTHS[d.getUTCMonth()];
}
const DOY_TICKS = [182, 196, 213, 227, 244, 258, 274, 288, 305, 319, 335];
const pct = (lr) => { const v = Math.exp(lr) - 1; if (Math.abs(v) < 5e-4) return "0%"; return (v >= 0 ? "+" : "") + (Math.abs(v) < 0.1 ? (100 * v).toFixed(1) : Math.round(100 * v)) + "%"; };
const signed = (v, d = 0) => (v >= 0 ? "+" : "") + v.toFixed(d);
const esc = (s) => String(s).replace(/[&<>]/g, (c) => ({"&": "&amp;", "<": "&lt;", ">": "&gt;"}[c]));

function checked(id) { return new Set([...document.querySelectorAll(`#${id} input:checked`)].map((e) => e.value)); }
function filtered(key) {
  const g = checked("f-group"), t = checked("f-tier"), c = checked("f-cls");
  const q = document.getElementById("f-name").value.trim().toLowerCase();
  return DATA.filter((r) => r[key] && g.has(r.group) && t.has(r.tier) && c.has(r[key].cls)
    && (!q || [r.en, r.fr, r.sci].some((n) => String(n).toLowerCase().includes(q))));
}
function sorted(rows, how, value) {
  const by = {
    value: (a, b) => value(a) - value(b),
    order: (a, b) => a.order - b.order,
    birds: (a, b) => b.birds - a.birds,
    name: (a, b) => a.en.localeCompare(b.en),
    group: (a, b) => a.group.localeCompare(b.group) || value(a) - value(b),
  }[how] || ((a, b) => value(a) - value(b));
  return rows.slice().sort(by).reverse();  // Plotly draws the first category at the bottom
}
function layout(rows, extra) {
  const n = rows.length;
  return Object.assign({
    template: "plotly_white", height: TOP + BOTTOM + Math.max(n, 3) * ROW, showlegend: false,
    margin: {l: LEFT, r: 30, t: TOP, b: BOTTOM}, hovermode: "closest", hoverlabel: {align: "left"},
    yaxis: {type: "category", categoryorder: "array", categoryarray: rows.map((r) => r.label),
            automargin: true, tickfont: {size: 11}, fixedrange: true},
    xaxis: {side: "top", gridcolor: "#ececec", zeroline: false, tickangle: 0},
  }, extra);
}
// the class of a row: colour kept, hollow if caveated, grey if hidden
function markers(rows, color) {
  return {
    color: rows.map((r) => r.cls === "hide" ? "#bbb" : r.cls === "caveat" ? "#fff" : color(r)),
    line: {width: 1.5, color: rows.map((r) => r.cls === "hide" ? "#999" : color(r))},
    size: 8, symbol: "circle",
  };
}
function whiskers(rows, lo, hi, color, width) {  // one trace of segments, rows split by null
  const x = [], y = [];
  rows.forEach((r, i) => { if (lo(r) != null) { x.push(lo(r), hi(r), null); y.push(r.label, r.label, null); } });
  return {x, y, mode: "lines", line: {color, width}, hoverinfo: "skip", type: "scatter"};
}
function hoverhead(r) {
  return `<b>${esc(r.en)}</b> · ${esc(r.fr)} · <i>${esc(r.sci)}</i><br>${r.group}, tier ${r.tier}, ${Math.round(r.birds).toLocaleString()} birds per season`;
}
function why(c) { return CLS_NAME[c.cls] + (c.why.length ? ` (${c.why.join(", ")})` : ""); }
function draw(id, traces, lay, rows, countId, total) {
  const el = document.getElementById(id);
  Plotly.react(el, traces, lay, {responsive: true, displaylogo: false});
  document.getElementById(countId).textContent = `${rows.length} of ${total} taxa`;
  el.removeAllListeners && el.removeAllListeners("plotly_click");
  el.on("plotly_click", (e) => { const r = e.points[0].customdata; if (r && r.page) window.open(r.id + ".html", "_blank"); });
}
function key(id, items) {
  document.getElementById(id).innerHTML = items.map(([c, t, hollow]) =>
    `<span><span class="sw" style="background:${hollow ? "#fff" : c};border:1.5px solid ${c === "#bbb" ? "#999" : c}"></span>${t}</span>`).join("");
}
const CLASS_KEY = [["#555", "shown"], ["#555", "caveated", true], ["#bbb", "hidden"]];
const GROUP_KEY = Object.entries(GROUPS).map(([g, c]) => [c, g]);

// --- 1. abundance trend -------------------------------------------------------
const UP = "#2a78d6", DOWN = "#e34948";
function trendValue(r, v = r.trend.log) {
  return state.t.unit === "year" ? v / (r.trend.last - r.trend.first) : v;
}
function trendFig() {
  const total = DATA.filter((r) => r.trend).length;
  let rows = filtered("trend").map((r) => Object.assign({}, r, {cls: r.trend.cls, label: r.en}));
  rows = sorted(rows, document.getElementById("t-sort").value, (r) => trendValue(r));
  const color = (r) => r.trend.log >= 0 ? UP : DOWN;
  if (state.t.view === "heat") {
    const years = []; for (let y = 1993; y <= 2025; y++) years.push(y);
    const z = rows.map((r) => years.map((y) => { const i = r.trend.years.indexOf(y); return i < 0 ? null : r.trend.rel[i]; }));
    const text = rows.map((r) => years.map((y) => { const i = r.trend.years.indexOf(y); return i < 0 ? "" :
      `${hoverhead(r)}<br>${y}: ${r.trend.rel[i] >= 0 ? "×" + Math.pow(2, r.trend.rel[i]).toFixed(2) : "÷" + Math.pow(2, -r.trend.rel[i]).toFixed(2)} the taxon's mean<br>trend: ${why(r.trend)}`; }));
    const trace = {type: "heatmap", x: years, y: rows.map((r) => r.label), z, text, hovertemplate: "%{text}<extra></extra>",
      zmin: -3, zmax: 3, zmid: 0, xgap: 1, ygap: 1,
      colorscale: [[0, "#b8332f"], [0.25, "#ec9e9a"], [0.5, "#f0efec"], [0.75, "#86b6ef"], [1, "#1c5cab"]],
      colorbar: {title: {text: "vs mean", side: "right"}, tickvals: [-3, -2, -1, 0, 1, 2, 3], ticktext: ["÷8", "÷4", "÷2", "mean", "×2", "×4", "×8"], thickness: 10, len: Math.min(1, 400 / (rows.length * ROW + 1))}};
    key("t-key", [["#b8332f", "below the taxon's mean"], ["#1c5cab", "above"]]);
    draw("fig-trend", [trace], layout(rows, {xaxis: {side: "top", dtick: 4}}), rows, "t-count", total);
    return;
  }
  const v = (r, x) => trendValue(r, x);
  const fmt = (lr) => state.t.unit === "year" ? pct(lr) + "/yr" : pct(lr);
  const ticksTotal = [-0.9, -0.75, -0.5, 0, 1, 4, 9, 49].map((c) => Math.log(1 + c));
  const ticksYear = [-0.1, -0.05, 0, 0.05, 0.1, 0.2].map((c) => Math.log(1 + c));
  const ticks = state.t.unit === "year" ? ticksYear : ticksTotal;
  const traces = [
    whiskers(rows, (r) => r.trend.lo == null ? null : v(r, r.trend.lo), (r) => v(r, r.trend.hi), "rgba(80,80,80,0.45)", 1.5),
    {type: "scatter", mode: "markers", x: rows.map((r) => v(r, r.trend.log)), y: rows.map((r) => r.label),
     marker: markers(rows, color), customdata: rows,
     text: rows.map((r) => `${hoverhead(r)}<br><b>${pct(r.trend.log)}</b> ${r.trend.first}-${r.trend.last}${r.trend.lo == null ? "" : ` (95%, conservative: ${pct(r.trend.lo)} to ${pct(r.trend.hi)})`}<br>` +
       `${pct(r.trend.log / (r.trend.last - r.trend.first))} per year · smooth total ${Math.round(r.trend.s0).toLocaleString()} → ${Math.round(r.trend.s1).toLocaleString()}<br>trend: ${why(r.trend)}`),
     hovertemplate: "%{text}<extra></extra>"},
  ];
  const lay = layout(rows, {
    xaxis: {side: "top", tickangle: 0, gridcolor: "#ececec", tickvals: ticks, ticktext: ticks.map(fmt), zeroline: true, zerolinecolor: "#999",
            title: {text: state.t.unit === "year" ? "change per year" : "change from first to last year (log scale)", font: {size: 12}}},
  });
  key("t-key", [[UP, "increase"], [DOWN, "decline"], ...CLASS_KEY]);
  draw("fig-trend", traces, lay, rows, "t-count", total);
}

// --- 2. passage dates ---------------------------------------------------------
function passageFig() {
  const total = DATA.filter((r) => r.passage).length;
  let rows = filtered("passage").map((r) => Object.assign({}, r, {cls: r.passage.cls, label: r.en}));
  const how = document.getElementById("p-sort").value;
  const val = {q10: (r) => r.passage.q10, q90: (r) => r.passage.q90, span: (r) => r.passage.q90 - r.passage.q10}[how] || ((r) => r.passage.q50);
  rows = sorted(rows, ["q10", "q50", "q90", "span"].includes(how) ? "value" : how, val);
  const xaxis = {side: "top", tickangle: 0, gridcolor: "#ececec", tickvals: DOY_TICKS, ticktext: DOY_TICKS.map(date), range: [180, 336]};
  const head = (r) => `${hoverhead(r)}<br>10% ${date(r.passage.q10)} · <b>50% ${date(r.passage.q50)}</b> · 90% ${date(r.passage.q90)}` +
    (r.passage.q50_lo != null ? `<br>median 80% band ${date(r.passage.q50_lo)} to ${date(r.passage.q50_hi)}` : "") +
    `<br>${r.passage.source === "smooth" ? "smooth season, last year" : "counted seasons pooled"} · ${why(r.passage)}`;
  if (state.p.view === "heat") {
    const doy = []; for (let d = 180; d <= 336; d++) doy.push(d);
    const z = rows.map((r) => doy.map((d) => { if (!r.passage.density) return null; const i = d - r.passage.doy[0]; return i >= 0 && i < r.passage.density.length ? r.passage.density[i] : null; }));
    const trace = {type: "heatmap", x: doy, y: rows.map((r) => r.label), z, zmin: 0, zmax: 1, ygap: 1,
      colorscale: [[0, "#f7fafe"], [0.15, "#cde2fb"], [0.4, "#86b6ef"], [0.7, "#2a78d6"], [1, "#0d366b"]],
      customdata: rows.map((r) => doy.map(() => r)), text: rows.map((r) => doy.map((d) => head(r) + `<br>${date(d)}`)),
      hovertemplate: "%{text}<br>%{z:.0%} of the peak day<extra></extra>",
      colorbar: {title: {text: "of peak", side: "right"}, tickformat: ".0%", thickness: 10, len: Math.min(1, 400 / (rows.length * ROW + 1))}};
    const med = {type: "scatter", mode: "markers", x: rows.map((r) => r.passage.q50), y: rows.map((r) => r.label),
      marker: {symbol: "line-ns", size: 9, line: {width: 2, color: "#e34948"}}, hoverinfo: "skip"};
    key("p-key", [["#0d366b", "share of the season, against the peak day"], ["#e34948", "median date"]]);
    draw("fig-passage", [trace, med], layout(rows, {xaxis}), rows, "p-count", total);
    return;
  }
  const traces = [];
  for (const [g, c] of Object.entries(GROUPS)) {
    const gr = rows.filter((r) => r.group === g);
    const t = whiskers(gr, (r) => r.passage.q10, (r) => r.passage.q90, c, 6);
    t.opacity = 0.35; traces.push(t);
  }
  traces.push(whiskers(rows, (r) => r.passage.q50_lo, (r) => r.passage.q50_hi, "#222", 1.5));
  traces.push({type: "scatter", mode: "markers", x: rows.map((r) => r.passage.q50), y: rows.map((r) => r.label),
    marker: markers(rows, (r) => GROUPS[r.group]), customdata: rows, text: rows.map(head), hovertemplate: "%{text}<extra></extra>"});
  key("p-key", [...GROUP_KEY, ...CLASS_KEY, ["#999", "from counts (no trend): filled, no whisker"]]);
  draw("fig-passage", traces, layout(rows, {xaxis}), rows, "p-count", total);
}

// --- 3. passage shift ---------------------------------------------------------
const QCOL = {q10: "#86b6ef", q50: "#2a78d6", q90: "#104281"};
const QNAME = {q10: "10%", q50: "median", q90: "90%"};
function shiftFig() {
  const total = DATA.filter((r) => r.shift).length;
  let rows = filtered("shift").map((r) => Object.assign({}, r, {cls: r.shift.cls, label: r.en}));
  const qs = state.s.q === "all" ? ["q10", "q50", "q90"] : [state.s.q];
  const unit = (r, d) => state.s.unit === "decade" ? 10 * d / (r.shift.last - r.shift.first) : d;
  const how = document.getElementById("s-sort").value;
  const main = state.s.q === "all" ? "q50" : state.s.q;
  rows = sorted(rows, how === "q50" ? "value" : how, how === "q50" ? (r) => r.passage ? r.passage.q50 : r.shift.q50.to : (r) => unit(r, r.shift[main].d));
  const traces = [];
  const off = state.s.q === "all" ? {q10: -0.25, q50: 0, q90: 0.25} : {[main]: 0};
  for (const q of qs) {
    const color = (r) => state.s.q === "all" ? QCOL[q] : GROUPS[r.group];
    const y = rows.map((r, i) => i + off[q]);
    const wx = [], wy = [];
    rows.forEach((r, i) => { wx.push(unit(r, r.shift[q].lo), unit(r, r.shift[q].hi), null); wy.push(y[i], y[i], null); });
    traces.push({type: "scatter", mode: "lines", x: wx, y: wy, line: {color: "rgba(80,80,80,0.45)", width: 1.5}, hoverinfo: "skip"});
    traces.push({type: "scatter", mode: "markers", x: rows.map((r) => unit(r, r.shift[q].d)), y,
      marker: Object.assign(markers(rows, color), {size: state.s.q === "all" ? 7 : 8}), customdata: rows,
      text: rows.map((r) => { const s = r.shift[q]; return `${hoverhead(r)}<br>${QNAME[q]} date ${date(s.from)} (${r.shift.first}) → ${date(s.to)} (${r.shift.last})<br>` +
        `<b>${signed(s.d, 1)} days</b> (80%, conservative: ${signed(s.lo, 1)} to ${signed(s.hi, 1)}) · ${signed(10 * s.d / (r.shift.last - r.shift.first), 1)} days per decade<br>season: ${why(r.shift)}`; }),
      hovertemplate: "%{text}<extra></extra>"});
  }
  const lay = layout(rows, {
    yaxis: {tickvals: rows.map((r, i) => i), ticktext: rows.map((r) => r.en), tickfont: {size: 11}, fixedrange: true, range: [-0.7, rows.length - 0.3], automargin: true},
    xaxis: {side: "top", tickangle: 0, gridcolor: "#ececec", zeroline: true, zerolinecolor: "#999", ticksuffix: " d",
            title: {text: (state.s.unit === "decade" ? "days per decade" : "days, first to last year") + " (← earlier · later →)", font: {size: 12}}},
  });
  key("s-key", state.s.q === "all" ? [[QCOL.q10, "10% date"], [QCOL.q50, "median"], [QCOL.q90, "90% date"], ...CLASS_KEY] : [...GROUP_KEY, ...CLASS_KEY]);
  draw("fig-shift", traces, lay, rows, "s-count", total);
}

function all() { trendFig(); passageFig(); shiftFig(); }
function seg(id, obj, prop, redraw) {
  document.querySelectorAll(`#${id} button`).forEach((b) => b.addEventListener("click", () => {
    document.querySelectorAll(`#${id} button`).forEach((o) => o.classList.toggle("on", o === b));
    obj[prop] = b.dataset.v; redraw();
  }));
}
seg("t-view", state.t, "view", trendFig); seg("t-unit", state.t, "unit", trendFig);
seg("p-view", state.p, "view", passageFig);
seg("s-q", state.s, "q", shiftFig); seg("s-unit", state.s, "unit", shiftFig);
["t-sort"].forEach((i) => document.getElementById(i).addEventListener("change", trendFig));
document.getElementById("p-sort").addEventListener("change", passageFig);
document.getElementById("s-sort").addEventListener("change", shiftFig);
document.querySelectorAll(".bar input").forEach((e) => e.addEventListener("input", all));
all();
"""


def compare_page(data: list) -> str:
    groups = "".join(
        f'<label><input type="checkbox" value="{g}" checked><span class="sw" '
        f'style="background:{c}"></span>{g}</label>'
        for g, c in GROUPS.items()
    )
    body = (
        BODY.replace("__GROUPS__", groups)
        .replace("__KEY_YEARS__", str(KEY_YEARS))
        .replace("__DATA__", json.dumps(data, separators=(",", ":"), allow_nan=False))
        .replace("__GROUP_COLORS__", json.dumps(GROUPS))
        .replace("__JS__", JS)
    )
    return wrap("Species compared", f"<style>{CSS}</style>{body}")


def write(export: str, out: str, pages: set) -> str:
    taxa = json.load(open(os.path.join(export, "taxa.json")))
    path, page = os.path.join(out, "compare.html"), compare_page(rows(taxa, export, pages))
    with open(path, "w", encoding="utf-8") as f:
        f.write(page)
    return path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--export", default=os.path.join(ROOT, "data", "explore"))
    ap.add_argument("--out", default=os.path.join(ROOT, "logs", "viewer"))
    args = ap.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    pages = {f[:-5] for f in os.listdir(args.out) if f.endswith(".html")}
    print("->", write(args.export, args.out, pages))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
