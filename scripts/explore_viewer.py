#!/usr/bin/env python3
"""QA viewer of the Explore export: an index of every taxon with its diagnostics, and a page per
taxon with every block drawn (Plotly). It reads `data/explore/` only and computes nothing, so what
it shows is what defileViz gets; a panel that needs a calculation here belongs in the pipeline.

A developer's tool, focused on the model output: show everything that helps judge a fit, in full
detail, in the order of the pipeline (configuration, counts and effort, time of day and coverage,
season day by day, trend model, validation, derived blocks). The visitors' page is defileViz's,
trimmed and narrative, in the order of a visitor's questions; it does not have to look like this
one. Everything is drawn, including what defileViz will hide: a page opens with what defileViz
shows, caveats and hides (the `reliability` block's classes, applied to its `elements`), each
figure starts with a strip saying the same for its parts, and a hidden part is greyed and labelled
in the figure itself.

Usage:
    python scripts/explore_viewer.py                    # -> logs/viewer/index.html + one page per taxon
    python scripts/explore_viewer.py --taxa "Red Kite"  # only these pages (and the index)
"""

import argparse
import html
import json
import os

import numpy as np
import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from defile_explore.reliability import ELEMENTS

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COLORS = ["#4C78A8", "#F58518", "#54A24B", "#B279A2", "#E45756"]
FIG_H = 380
PROTOCOL = {1993: "daily counts start", 2007: "passerines counted", 2014: "hourly sheets",
            2021: "Trektellen"}  # fmt: skip
DOY_TICKS = [182, 196, 213, 227, 244, 258, 274, 288, 305, 319, 335]
LINK_NAMES = {
    "ebird": "eBird",
    "ebird_status": "eBird S&T",
    "birds_of_the_world": "Birds of the World",
    "ebba2": "EBBA2",
    "trektellen": "Trektellen",
    "vogelwarte": "Vogelwarte",
    "migration_atlas": "Migration Atlas",
}
CSS = """
body{font-family:-apple-system,Segoe UI,sans-serif;margin:0;background:#fafafa;color:#222}
nav{position:sticky;top:0;background:#fff;border-bottom:1px solid #ddd;padding:8px 16px;z-index:9}
nav a{margin-right:14px;color:#4C78A8;text-decoration:none;font-weight:600}
section{max-width:1200px;margin:0 auto;padding:16px}
h2{border-bottom:2px solid #4C78A8;padding-bottom:4px}
h3{margin-top:28px;color:#555;font-size:13px;text-transform:uppercase;letter-spacing:.06em}
.panel{background:#fff;border:1px solid #e3e3e3;border-radius:6px;padding:8px;margin:10px 0}
.note{color:#777;font-size:12px;margin:4px 8px}
.knrow{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:8px}
.kn{background:#fff;border:1px solid #e3e3e3;border-radius:6px;padding:8px}
.kl{font-size:11px;color:#777}.kv{font-size:20px;font-weight:700}.ks{font-size:11px;color:#999}
.lk{display:inline-block;background:#eef3f8;border-radius:12px;padding:2px 10px;margin:2px;
 font-size:12px;color:#335;text-decoration:none}
.flag{display:inline-block;background:#fdecea;color:#b71c1c;border-radius:3px;padding:1px 6px;
 margin:1px;font-size:11px}
ul.rep{max-height:320px;overflow:auto;font-size:13px}
table{border-collapse:collapse;font-size:13px;width:100%}
td,th{border-bottom:1px solid #eee;padding:3px 6px;text-align:left}
td.n{text-align:right}
p.rn{margin:2px 0;font-size:12px;max-height:120px;overflow:auto}
td{vertical-align:top}
.vis{font-size:12px;margin:2px 4px 6px;padding:4px 6px;background:#f6f6f6;border-radius:4px}
.vis b{color:#555}.why{color:#888;font-size:11px}
.toc{font-size:12px;margin:4px 0 8px}.toc a{color:#4C78A8}
"""


def doy_label(doy, fmt: str = "%d %b") -> str:
    """A (non-leap) day of year as a date; `"%d-%b"` for hovers."""
    import datetime as dt

    return (dt.date(2025, 1, 1) + dt.timedelta(days=int(round(doy)) - 1)).strftime(fmt)


def doy_of(mmdd: str) -> int:
    import datetime as dt

    return dt.date(2025, *map(int, mmdd.split("-"))).timetuple().tm_yday


def doy_axis(fig, row=None, col=None, axis="x"):
    upd = dict(tickvals=DOY_TICKS, ticktext=[doy_label(t) for t in DOY_TICKS])
    (fig.update_xaxes if axis == "x" else fig.update_yaxes)(**upd, row=row, col=col)


def arr(x):
    return np.array([np.nan if v is None else v for v in x], float)


def panel(fig, title, note="", top=40, strip=""):
    fig.update_layout(template="plotly_white", margin=dict(l=50, r=20, t=top, b=40),
                      title=dict(text=title, font=dict(size=14)))  # fmt: skip
    return (
        '<div class="panel">'
        + strip
        + fig.to_html(full_html=False, include_plotlyjs=False, config={"responsive": True})
        + (f'<p class="note">{note}</p>' if note else "")
        + "</div>"
    )


def method(block) -> str:
    return f" <span style='font-size:11px;color:#888'>[{block['method']}]</span>" if block else ""


# --- what defileViz shows ---------------------------------------------------------------------

CLASS_STYLE = {"show": "#e6f4ea;color:#1e6b34", "caveat": "#fff4e0;color:#8a5300",
               "hide": "#fdecea;color:#b71c1c"}  # fmt: skip
CLAIMS = ("totals", "trend", "season")
PARTS = {  # `reliability.ELEMENTS`: what this viewer draws of each
    "trend.annual.total": "season totals (points, 80%)",
    "trend.days.total": "filled days (mean, 80%)",
    "key_numbers.typical_season": "key number: typical season",
    "trend.annual.smooth": "smooth trend (line, 95% band)",
    "key_numbers.trend": "key number: trend",
    "trend.passage": "modelled passage dates per year",
    "trend.passage_q": "smooth passage dates (lines)",
    "key_numbers.passage": "key number: main passage",
    "season.share": "season day by day (gap-filled)",
    "season.passage": "passage dates per year (gap-filled, 80%)",
    "season.chances": "chances (gap-filled)",
}
HIDDEN_OPACITY = 0.3


def badge(cls: str, text: str | None = None) -> str:
    return (f"<span class='flag' style='background:{CLASS_STYLE[cls]}'>"
            f"{html.escape(text or cls)}</span>")  # fmt: skip


def fate(sp: dict, element: str) -> str:
    """The class defileViz gives a drawn element; `show` for counts, which rest on no claim."""
    if element == "key_numbers.passage" and sp["key_numbers"]["passage"]["source"] != "smooth":
        return "show"
    return ((sp.get("reliability") or {}).get("elements") or {}).get(element, "show")


def reasons(sp: dict, element: str) -> str:
    q = sp.get("reliability")
    if not q or element not in ELEMENTS or fate(sp, element) == "show":
        return ""
    claim = q[ELEMENTS[element]]
    codes = ", ".join(r["code"] for r in claim["reasons"])
    return f"<span class='why'>({ELEMENTS[element]}: {html.escape(codes)})</span>"


def strip(sp: dict, elements=(), counts: str = "") -> str:
    """One line over a figure: what defileViz does with each of its parts."""
    items = [f"{badge(fate(sp, e), PARTS[e])} {reasons(sp, e)}" for e in elements]
    if counts:
        items.insert(0, badge("show", counts))
    return f"<div class='vis'><b>defileViz:</b> {' '.join(items)}</div>"


def legend(sp: dict, element: str, name: str) -> dict:
    """A trace's name and opacity: tagged with its class unless shown, greyed when hidden."""
    cls = fate(sp, element)
    tag = "" if cls == "show" else f" [{cls.upper()}]"
    return {"name": name + tag, "opacity": HIDDEN_OPACITY if cls == "hide" else 1}


def hidden_note(fig, sp: dict, elements, row=None, col=None) -> None:
    """A watermark over a figure whose parts defileViz hides."""
    hidden = [PARTS[e] for e in elements if fate(sp, e) == "hide"]
    if hidden:
        fig.add_annotation(text="Hidden on defileViz: " + "; ".join(hidden), showarrow=False,
                           x=0.5, y=0.5, xref="x domain" if row else "paper",
                           yref="y domain" if row else "paper", bgcolor="rgba(253,236,234,0.85)",
                           font=dict(color="#b71c1c", size=13), row=row, col=col)  # fmt: skip


def display_panel(sp: dict, tx: dict) -> str:
    """The page's opening: every part of it, and what defileViz does with it."""
    q = sp.get("reliability")
    rows = []
    if q:
        for c in CLAIMS:
            parts = ", ".join(PARTS[e] for e, claim in ELEMENTS.items() if claim == c)
            why = ", ".join(f"{r['code']} ({r['level']})" for r in q[c]["reasons"]) or "–"
            rows.append((badge(q[c]["class"]), f"<b>{c}</b>: {parts}", why))
    else:
        rows.append((badge("hide"), "<b>trend</b>: no model (tier or settings)", "–"))
    d = sp["daytime"]
    rows += [
        (badge("show"), "birds counted per year, records, accounts, key numbers from counts"
                        + ("" if sp["season"].get("source") == "gam"
                           else ", phenology and chances (from the counts)"), "counts"),
        (badge("show") if d and d["hours"] is not None else badge("hide"),
         "passage during the day: " + ("not built" if not d or d["hours"] is None else
                                       "by date" if d["change"]["show"] else "one histogram"),
         "" if not d or d["hours"] is None or d["change"]["show"] else "no change over the season"),
        *((badge("show") if sp[f] else badge("hide"), f"{f}: " + (
            f"{len(sp[f]['years'])} years" if sp[f] else "not built"),
           "" if sp[f] else "too few birds classed") for f in ("age", "sex")),
    ]  # fmt: skip
    if q:
        est = ", ".join(map(str, q["estimated_years"])) or "none"
        rows.append((badge("caveat" if q["estimated_years"] else "show"),
                     f"years marked as mostly estimated: {est}", "under half counted"))  # fmt: skip
    body = "".join(
        f"<tr><td>{a}</td><td>{b}</td><td class='why'>{c}</td></tr>" for a, b, c in rows
    )
    return (f"<div class='panel'><table><tr><th>defileViz</th><th>part</th><th>why</th></tr>"
            f"{body}</table><p class='note'>Everything below is drawn; [CAVEAT] and [HIDE] mark "
            f"the parts defileViz caveats or leaves out.</p></div>")  # fmt: skip


# --- panels -----------------------------------------------------------------------------------


def key_numbers(sp: dict) -> str:
    k = sp["key_numbers"]
    items = []
    p = k["passage"]
    if p["q50"] is not None:
        items.append(("Main passage (80% of birds)", f"{doy_label(p['q10'])} – {doy_label(p['q90'])}",
                      f"median {doy_label(p['q50'])} ({p['source']})", "key_numbers.passage"))  # fmt: skip
    if "best_hours" in k:
        b = k["best_hours"]
        items.append(("Best hours (solar time)", f"{b['from']:02d}:00 – {b['to']:02d}:00",
                      f"{b['share']:.0%} of a peak day", None))  # fmt: skip
    if "typical_season" in k:
        t = k["typical_season"]
        items.append(("Typical season", f"{t['median']:,.0f}",
                      f"{t['min']:,.0f}–{t['max']:,.0f}, {t['years'][0]}–{t['years'][1]}",
                      "key_numbers.typical_season"))  # fmt: skip
    if "trend" in k:
        t = k["trend"]
        items.append(("Trend", f"{t['change']:+.0%}", f"smooth, {t['from']} → {t['to']}",
                      "key_numbers.trend"))  # fmt: skip
    if "chance" in k:
        c = k["chance"]
        items.append(("Chance of ≥1 / ≥10 in the main passage",
                      f"{c['at_least_1']:.0%} / {c['at_least_10']:.0%}", f"{c['days']} days", None))  # fmt: skip
    if "record" in k:
        items.append(("Record day", f"{k['record']['count']:,.0f}", k["record"]["date"], None))

    def tag(e):
        cls = "show" if e is None else fate(sp, e)
        return badge(cls, "shown" if cls == "show" else cls.upper()) + (
            reasons(sp, e) if e else ""
        )

    cells = "".join(
        f'<div class="kn"{" style=opacity:.5" if e and fate(sp, e) == "hide" else ""}>'
        f'<div class="kl">{html.escape(a)}</div><div class="kv">{b}</div>'
        f'<div class="ks">{html.escape(str(c))}</div><div>{tag(e)}</div></div>'
        for a, b, c, e in items
    )
    return f'<div class="knrow">{cells}</div>'


def trend_panel(sp: dict, effort_annual: list) -> str:
    t = sp["trend"]
    a = {k: arr([r[k] for r in t["annual"]]) for k in t["annual"][0]}
    start, last = t["first_year"], t["last_year"]
    raw = [r for r in sp["annual"] if r["year"] < start]
    texts = {r["year"]: r["en"] for r in (sp["accounts"] or {}).get("years", [])}
    hover = ["<br>".join(texts.get(int(y), "–")[i : i + 80] for i in range(0, 320, 80))
             for y in a["year"]]  # fmt: skip
    eff = [r for r in effort_annual if r["year"] <= last]
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.78, 0.22],
                        vertical_spacing=0.04)  # fmt: skip
    fig.add_trace(go.Bar(x=a["year"], y=a["observed"], name="birds counted", width=0.8,
                         marker_color="rgba(214,39,40,0.35)"), 1, 1)  # fmt: skip
    fig.add_trace(go.Bar(x=[r["year"] for r in raw], y=[r["window"] for r in raw], width=0.8,
                         name="historical counts (raw, not comparable)", visible=False,
                         marker_color="rgba(150,150,150,0.5)"), 1, 1)  # fmt: skip
    smooth, total = "trend.annual.smooth", "trend.annual.total"
    fig.add_trace(go.Scatter(x=a["year"], y=a["smooth_q97.5"], line=dict(width=0),
                             showlegend=False, hoverinfo="skip"), 1, 1)  # fmt: skip
    fig.add_trace(go.Scatter(x=a["year"], y=a["smooth_q2.5"], fill="tonexty", line=dict(width=0),
                             fillcolor="rgba(76,120,168,0.18)", hoverinfo="skip",
                             **legend(sp, smooth, "smooth trend 95%")), 1, 1)  # fmt: skip
    fig.add_trace(go.Scatter(x=a["year"], y=a["smooth"], line=dict(color=COLORS[0], width=2.5),
                             **legend(sp, smooth, "smooth trend (typical year)")), 1, 1)  # fmt: skip
    fig.add_trace(go.Scatter(
        x=a["year"], y=a["total"], mode="markers",
        **legend(sp, total, "season total (gap-filled, 80%)"),
        marker=dict(size=8, color="black"),
        error_y=dict(type="data", symmetric=False, array=a["q90"] - a["total"],
                     arrayminus=a["total"] - a["q10"], thickness=1, color="#555"),
        customdata=np.stack([a["observed"], a["observed_share"], hover], axis=1),
        hovertemplate="%{x}: %{y:,.0f} (counted %{customdata[0]:,.0f}, %{customdata[1]:.0%})"
                      "<br><i>%{customdata[2]}</i><extra></extra>"), 1, 1)  # fmt: skip
    fig.add_trace(go.Bar(x=[r["year"] for r in eff], y=[r["window_days"] for r in eff],
                         marker_color="#bbb", showlegend=False, name="days counted"), 2, 1)  # fmt: skip
    for k, (y, lab) in enumerate(PROTOCOL.items()):
        fig.add_vline(x=y - 0.5, line=dict(color="#aaa", dash="dot", width=1))
        fig.add_annotation(x=y - 0.5, y=1.0 - 0.05 * k, yref="paper", text=lab, showarrow=False,
                           font=dict(size=9, color="#777"), xanchor="left")  # fmt: skip
    modern = [start - 0.7, last + 0.7]
    full = [min([r["year"] for r in raw] + [start]) - 0.7, last + 0.7]
    vis = [True] * len(fig.data)
    vis_mod = vis.copy()
    vis_mod[1] = False
    fig.update_layout(
        height=FIG_H + 120, legend=dict(orientation="h", y=-0.12), barmode="overlay",
        xaxis=dict(range=modern), xaxis2=dict(range=modern),
        updatemenus=[dict(type="buttons", direction="right", x=0, y=1.13, xanchor="left", buttons=[
            dict(label=f"From {start}", method="update",
                 args=[{"visible": vis_mod}, {"xaxis.range": modern, "xaxis2.range": modern}]),
            dict(label="Include historical counts", method="update",
                 args=[{"visible": vis}, {"xaxis.range": full, "xaxis2.range": full}]),
        ])])  # fmt: skip
    fig.update_yaxes(title="birds / season", row=1, col=1)
    fig.update_yaxes(title="days", row=2, col=1)
    hidden_note(fig, sp, (total, smooth), 1, 1)
    est = set((sp.get("reliability") or {}).get("estimated_years", []))
    if est:  # the years defileViz marks as mostly estimated
        fig.add_trace(go.Scatter(x=[y for y in a["year"] if y in est],
                                 y=[v for y, v in zip(a["year"], a["total"]) if y in est],
                                 mode="markers", name="mostly estimated (marked)",
                                 marker=dict(size=13, color="rgba(0,0,0,0)",
                                             line=dict(color="#8a5300", width=2))), 1, 1)  # fmt: skip
    return panel(fig, "Long-term trend" + method(t),
                 "Bars: birds counted. Points: counted + posterior fill (80%); ringed: under half "
                 "counted. Line: smooth trend. Bottom: days counted in the window (all taxa).",
                 strip=strip(sp, (total, smooth), "birds counted"))  # fmt: skip


def reliability_panel(sp: dict) -> str:
    """The three claims' classes and reasons, and the numbers the rules read."""
    q = sp.get("reliability")
    if not q:
        return ""
    rows = "".join(
        f"<tr><td><b>{c}</b></td><td>{badge(q[c]['class'])}</td><td>"
        + (" ".join(badge(r["level"], r["code"]) for r in q[c]["reasons"]) or "")
        + "</td></tr>"
        for c in CLAIMS
    )
    est = ", ".join(map(str, q["estimated_years"])) or "none"
    x = {k: round(v, 3) if isinstance(v, float) else v for k, v in q["inputs"].items()}
    return (f"<div class='panel'><table><tr><th>claim</th><th>class</th><th>reasons</th></tr>{rows}"
            f"</table><p class='note'>Mostly estimated years: {est}.<br>{html.escape(json.dumps(x))}"
            f" [{q['method']}]</p></div>")  # fmt: skip


def benchmark_panel(sp: dict) -> str:
    """Each gap-transplant trial (a well-counted year given an older year's gaps): the fill against
    the birds it actually counted, with its 80% and 95% intervals; and the recent seasons predicted
    from the earlier ones."""
    b = sp.get("benchmark")
    if not b or not (b["gap"] or b["recent"]):
        return "<p class='note'>No benchmark.</p>"
    titles = ("gap transplant: filled vs. counted", "recent seasons: predicted vs. counted")
    fig = make_subplots(rows=1, cols=2, subplot_titles=titles)
    for col, key in ((1, "gap_trials"), (2, "recent_trials")):
        t = b[key]
        if not t:
            continue
        t = {k: arr([r[k] for r in t]) for k in t[0]}
        x, y = t["truth"], t["estimate"]
        label = t["target"] if "target" in t else t["year"]
        for q_lo, q_hi, w in (("q2.5", "q97.5", 1), ("q10", "q90", 3)):
            fig.add_trace(go.Scatter(x=x, y=y, mode="markers", marker=dict(size=1, color=COLORS[0]),
                                     error_y=dict(type="data", symmetric=False, array=t[q_hi] - y,
                                                  arrayminus=y - t[q_lo], width=0, thickness=w,
                                                  color=COLORS[0]),
                                     hoverinfo="skip", showlegend=False), 1, col)  # fmt: skip
        fig.add_trace(go.Scatter(x=x, y=y, mode="markers", marker=dict(color=COLORS[0], size=7),
                                 customdata=label, showlegend=False,
                                 hovertemplate="%{customdata}: counted %{x:,.0f}, estimated "
                                               "%{y:,.0f}<extra></extra>"), 1, col)  # fmt: skip
        lo, hi = np.nanmin([x.min(), y.min()]), np.nanmax([x.max(), y.max()])
        fig.add_trace(go.Scatter(x=[lo, hi], y=[lo, hi], mode="lines", showlegend=False,
                                 line=dict(color="#999", dash="dot")), 1, col)  # fmt: skip
        fig.update_xaxes(type="log", title="counted", row=1, col=col)
        fig.update_yaxes(type="log", title="estimated", row=1, col=col)
    fig.update_layout(height=FIG_H + 20)

    def line(name, s):
        if not s:
            return f"{name}: not run"
        return (f"{name}: error {s['abs_log_err']:.3f}, bias {s['bias']:+.3f}, 80% cover "
                f"{s['cover80']:.2f}, 95% cover {s['cover95']:.2f} ({s['n']:.0f} trials)")  # fmt: skip

    ratio = b["gap_ratio"]
    note = "<br>".join([
        line("Gap fill (GAM)", b["gap"]),
        f"Ratio estimator: error {ratio['abs_log_err']:.3f}" if ratio else "",
        line("Recent seasons", b["recent"]),
        "Error: median |log(estimate / counted)|. Intervals: 80% thick, 95% thin.",
    ])  # fmt: skip
    return panel(fig, "Benchmark" + method(b), note, top=75)


def season_panel(sp: dict) -> str:
    s = sp["season"]
    gam = s.get("source") == "gam"
    share = np.array([[np.nan if v is None else v for v in row] for row in s["share"]], float)
    zmax = float(np.nanquantile(share, 0.98)) if np.isfinite(share).any() else 1
    cov = s.get("c") or [[None] * len(s["doy"])] * len(s["years"])
    dates = [doy_label(x, "%d-%b") for x in s["doy"]]

    def cell(n, c):
        if n is None:
            return "not counted"
        return f"{n:,.0f} counted" + (f", c = {c:.2f}" if c else ", c < 0.1")

    custom = [[[dates[j], cell(n, c)] for j, (n, c) in enumerate(zip(nr, cr))]
              for nr, cr in zip(s["count"], cov)]  # fmt: skip
    fig = go.Figure(go.Heatmap(x=s["doy"], y=s["years"], z=share, zmin=0, zmax=zmax,
                               colorscale="YlOrRd", colorbar=dict(title="share<br>of year"),
                               customdata=custom,
                               hovertemplate="%{customdata[0]} %{y}: %{customdata[1]} "
                                             "(%{z:.1%} of the season)<extra></extra>"))  # fmt: skip
    if gam:  # days not counted, filled: hatched over
        filled = [[d for d, n in zip(s["doy"], nr) if n is None] for nr in s["count"]]
        fig.add_trace(go.Scatter(x=[d for f in filled for d in f],
                                 y=[y for y, f in zip(s["years"], filled) for _ in f],
                                 mode="markers", marker=dict(symbol="line-ew", size=4,
                                                             line=dict(color="rgba(0,0,0,0.35)",
                                                                       width=1)),
                                 name="not counted (filled)", hoverinfo="skip"))  # fmt: skip
    p = {k: arr([r[k] for r in s["passage"]]) for k in s["passage"][0]}
    smooth = sp["trend"]["passage_q"] if sp["trend"] else None
    for col, dash, lab in (("q10", "dot", "10%"), ("q50", "solid", "50%"), ("q90", "dot", "90%")):
        if smooth:
            fig.add_trace(go.Scatter(x=[r[col] for r in smooth], y=[r["year"] for r in smooth],
                                     mode="lines", legendgroup=lab,
                                     line=dict(color="black", dash=dash, width=1.5),
                                     **legend(sp, "trend.passage_q", f"{lab} smooth")))  # fmt: skip
        band = (dict(error_x=dict(type="data", symmetric=False, array=p["q50_hi"] - p["q50"],
                                  arrayminus=p["q50"] - p["q50_lo"], color="black",
                                  thickness=1.2, width=0))
                if gam and col == "q50" else {})  # fmt: skip
        opacity = 1 if gam else (0.3 + 0.7 * p["counted"]).tolist()
        fig.add_trace(go.Scatter(x=p[col], y=p["year"], mode="markers", legendgroup=lab,
                                 marker=dict(color="white" if col == "q50" else "#333", size=7,
                                             line=dict(color="black", width=1), opacity=opacity),
                                 **(legend(sp, "season.passage", f"{lab} that year") if gam
                                    else {"name": f"{lab} that year"}), **band))  # fmt: skip
    w = sp["window"]
    for i in (0, 1):  # window edges, where the model's differs from the default
        if w["model"][i] != w["default"][i]:
            for key, color in (("default", "#999"), ("model", "#1f77b4")):
                x = doy_of(w[key][i]) + (0.5 if i else -0.5)
                fig.add_vline(x=x, line=dict(color=color, dash="dash"))
    doy_axis(fig)
    fig.update_layout(height=FIG_H + 160, legend=dict(orientation="h", y=-0.1))
    shown = (("season.share", "season.passage") if gam else ()) + (
        ("trend.passage_q",) if smooth else ()
    )
    hidden_note(fig, sp, shown)
    title = "Phenology, every year: share of the season's birds per day"
    if smooth:
        title += (
            f" (median {smooth[-1]['q50'] - smooth[0]['q50']:+.0f} days since {smooth[0]['year']})"
        )
    source = (
        "Gap-filled by the trend GAM in the model window (mean of the draws; ticks: days not "
        "counted); beyond it, birds / c of days counted at c >= 0.5. Dots: that year's 10/50/90% "
        "passage dates over the draws, the median with its 80% band. "
        if gam
        else "From the counts; blank = not counted (c < 0.5). Dots: that year's 10/50/90% passage "
        "dates of the interpolated series (paler = more of the season uncounted). "
    )
    return panel(fig, title + method(s),
                 source + "Lines: the smooth. Dashed: the model "
                 f"window {' - '.join(w['model'])} (blue) where it differs from the default "
                 f"{' - '.join(w['default'])} (grey); shown: {' - '.join(w['view'])}. "
                 f"Passage beyond counting: {', '.join(w['beyond_counting']) or 'no'} (share of "
                 f"the passage in the envelope's first/last days: {w.get('edge_shares')}).",
                 strip=strip(sp, shown, "" if gam else "counted share per day"))  # fmt: skip


def chances_panel(sp: dict) -> str:
    ch = sp["season"]["chances"]
    gam = sp["season"].get("source") == "gam"
    fig = go.Figure()
    for i, k in enumerate((1, 10, 100, 1000)):
        if f"at_least_{k}" in ch:
            fig.add_trace(go.Scatter(x=ch["doy"], y=ch[f"at_least_{k}"], mode="lines+markers",
                                     name=f"≥ {k}", line=dict(color=COLORS[i], shape="spline")))  # fmt: skip
    doy_axis(fig)
    fig.update_yaxes(tickformat=".0%", range=[0, 1], title="share of days")
    fig.update_layout(height=FIG_H)
    if gam:
        hidden_note(fig, sp, ("season.chances",))
    note = (
        "Chance that a full day holds at least N birds, over the gap-filled days (every day of the "
        "model window, each weighing its probability over the draws; beyond it, well-counted days)."
        if gam
        else "Share of well-counted days with at least N birds."
    )
    return panel(fig, f"Your chances ({ch['years'][0]}–{ch['years'][1]})", note,
                 strip=strip(sp, ("season.chances",) if gam else (),
                             "" if gam else "counted days"))  # fmt: skip


def daytime_panel(sp: dict) -> str:
    d = sp["daytime"]
    if not d or d["hours"] is None:
        return "<p class='note'>No days timed to the hour.</p>"
    ch = d["change"]
    hist = [go.Bar(x=list(range(24)), y=d["hours"], marker_color="rgba(76,120,168,0.6)",
                   name="counted", hovertemplate="%{x}:00-%{x}:59: %{y:.1%}<extra></extra>")]  # fmt: skip
    if d["expected"] is not None:
        hist.append(go.Scatter(x=list(range(24)), y=d["expected"], line=dict(color="black"),
                               mode="lines+markers", line_shape="spline",
                               name="profile's prediction (same days, minutes counted)"))  # fmt: skip
    if ch["show"]:
        share = np.array([[np.nan if v is None else v for v in r] for r in d["share"]], float)
        fig = make_subplots(rows=1, cols=2, column_widths=[0.6, 0.4], horizontal_spacing=0.1,
                            subplot_titles=("by date (counted)", "whole season"))  # fmt: skip
        fig.add_trace(go.Heatmap(x=d["doy"], y=list(range(24)), z=share.T, colorscale="YlOrRd",
                                 showscale=False, customdata=[[doy_label(x, "%d-%b") for x in d["doy"]]] * 24,
                                 hovertemplate="%{customdata}, %{y}:00: %{z:.1%}<extra></extra>"),
                      1, 1)  # fmt: skip
        for t in hist:
            fig.add_trace(t, 1, 2)
        doy_axis(fig, 1, 1)
        fig.update_yaxes(title="solar hour (12 = solar noon)", range=[4.5, 19.5], row=1, col=1)
        fig.update_xaxes(range=[4.5, 19.5], title="solar hour (12 = solar noon)", row=1, col=2)
        fig.update_yaxes(tickformat=".0%", row=1, col=2)
        fig.update_annotations(yshift=6, font_size=12)
    else:
        fig = go.Figure(hist)
        fig.update_xaxes(range=[4.5, 19.5], title="solar hour (12 = solar noon)", dtick=1)
        fig.update_yaxes(tickformat=".0%", title="share of the day's birds")
    fig.update_layout(height=FIG_H, bargap=0.05)
    if "shift" in ch:
        change = (f"Late minus early season mean passage hour: {ch['shift']:+.1f} h "
                  f"(95% {ch['lo']:+.1f} to {ch['hi']:+.1f}, {ch['days'][0]} and {ch['days'][1]} "
                  f"days): {'shown by date' if ch['show'] else 'one histogram'}.")  # fmt: skip
    else:
        change = f"Change not tested (timed days early/late: {ch['days'][0]}/{ch['days'][1]})."
    return panel(fig, f"Passage during the day ({d['years'][0]}–{d['years'][1]}, {d['days']} "
                      f"days timed)" + method(d),
                 "Bars: birds per counted hour, as shares, pooled over the "
                 + (f"main passage ({doy_label(d['main_doy'][0], '%d %b')} to "
                    f"{doy_label(d['main_doy'][1], '%d %b')}, {d['main_days']} days)"
                    if d["main_doy"] else "season")
                 + ", hours counted at least half on each day; line: what the profile predicts "
                 "for the same days and minutes counted. Hours are solar time: add about 1.6 h "
                 "for summer time, 0.5 h for winter time. " + change,
                 top=75 if ch["show"] else 40,
                 strip=strip(sp, counts="by date" if ch["show"] else "one histogram"))  # fmt: skip


DEMOGRAPHY = {  # field: (what is plotted, its class labels, verb)
    "age": ("non-adult share", ("non-adult", "adult"), "aged"),
    "sex": ("male share", ("male", "female-type"), "sexed"),
}


def demography_panel(sp: dict, field: str) -> str:
    """Age or sex (same block): the first class's share per year, sized by the birds classed, and
    when each class passes."""
    a = sp[field]
    if not a:
        return ""
    what, labels, verb = DEMOGRAPHY[field]
    y, p, n = a["years"], arr(a["share"]), arr(a["n"])
    lo, hi, o = arr(a["lo"]), arr(a["hi"]), a["overall"]
    fig = make_subplots(rows=1, cols=2, subplot_titles=(f"{what} of birds {verb}",
                                                        f"timing: {labels[0]} vs. {labels[1]}"))  # fmt: skip
    x0, x1 = min(y) - 0.5, max(y) + 0.5
    fig.add_trace(go.Scatter(x=[x0, x1, x1, x0], y=[o["lo"], o["lo"], o["hi"], o["hi"]],
                             fill="toself", fillcolor="rgba(76,120,168,0.15)", line=dict(width=0),
                             hoverinfo="skip", name=f"all years: {o['share']:.0%}"), 1, 1)  # fmt: skip
    fig.add_trace(go.Scatter(x=[x0, x1], y=[o["share"]] * 2, mode="lines", showlegend=False,
                             line=dict(color=COLORS[0], dash="dash")), 1, 1)  # fmt: skip
    fig.add_trace(go.Scatter(x=y, y=p, mode="markers", name=f"per year (size = birds {verb})",
                             marker=dict(color=COLORS[0], size=6 + 22 * np.sqrt(n / n.max())),
                             error_y=dict(type="data", symmetric=False, array=hi - p,
                                          arrayminus=p - lo, width=0, color=COLORS[0]),
                             customdata=np.stack([n, arr(a["counted"])], axis=1),
                             hovertemplate=f"%{{x}}: %{{y:.0%}} of %{{customdata[0]:,.0f}} {verb} "
                                           "(of %{customdata[1]:,.0f} counted)<extra></extra>"),
                  1, 1)  # fmt: skip
    t = a["timing"]
    for i, (cls, lab) in enumerate(zip(a["classes"], labels)):
        if cls in t:
            fig.add_trace(go.Scatter(x=t["doy"], y=t[cls], line=dict(color=COLORS[i], shape="hv"),
                                     name=f"{lab} (n={t[cls + '_birds']:,.0f})"), 1, 2)  # fmt: skip
    fig.update_yaxes(tickformat=".0%", range=[0, 1])
    doy_axis(fig, 1, 2)
    fig.update_layout(height=FIG_H + 40, legend=dict(orientation="h", y=-0.15))
    note = "Female-type: female or juvenile, not separable in the field." if field == "sex" else ""
    return panel(fig, f"{field.capitalize()}: {o['share']:.0%} {labels[0]} ({o['lo']:.0%}–"
                      f"{o['hi']:.0%}), {o['n']:,.0f} {verb} of {o['counted']:,.0f} counted"
                      + method(a), note, top=75,
                 strip=strip(sp, counts=f"{len(y)} years"))  # fmt: skip


def records_panel(sp: dict) -> str:
    def note(n):
        ref = "" if n["ref"] is None else f"<i>{html.escape(n['ref'])}</i> — "
        return f"<p class='rn'>{ref}{html.escape(n['text'])}</p>"

    def day(r):
        if r.get("trektellen") is None:
            return r["date"]
        return f"<a href='{r['trektellen']}' target='_blank'>{r['date']}</a>"

    def notes(r):
        return "".join(note(n) for n in r["notes"])

    rows = "".join(
        f"<tr><td>{i + 1}</td><td>{day(r)}</td><td class='n'>{r['count']:,.0f}</td>"
        f"<td>{notes(r)}</td></tr>"
        for i, r in enumerate(sp["records"]["top_days"])
    )
    return (f"<div class='panel'><b>Top days</b><table><tr><th>#</th><th>date</th><th>birds</th>"
            f"<th>what was written that day</th></tr>{rows}</table></div>")  # fmt: skip


# --- inputs: counts, effort, time of day ------------------------------------------------------


def effort_panel(sp: dict, effort_annual: list) -> str:
    """Per year: the taxon's birds in the default window, the days counted (all taxa), the mean
    coverage of its counted days and the ratio index (birds per full-day equivalent)."""
    a = sp["annual"]
    year = [r["year"] for r in a]
    start = sp["settings"]["start_year"]["value"]
    eff = {r["year"]: r for r in effort_annual}
    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, row_heights=[0.45, 0.3, 0.25],
                        vertical_spacing=0.04)  # fmt: skip
    fig.add_trace(go.Bar(x=year, y=[r["window"] for r in a], name="birds counted (window)",
                         marker_color=["rgba(214,39,40,0.45)" if y >= start else "#bbb"
                                       for y in year]), 1, 1)  # fmt: skip
    fig.add_trace(go.Scatter(x=year, y=[r["index"] for r in a], mode="markers",
                             name="ratio index (birds / Σc)", marker=dict(color="black", size=5)),
                  2, 1)  # fmt: skip
    fig.add_trace(go.Scatter(x=year, y=[r["c_mean"] for r in a], mode="lines+markers",
                             name="mean coverage of counted days", line=dict(color=COLORS[0])),
                  3, 1)  # fmt: skip
    fig.add_trace(go.Bar(x=year, y=[eff.get(y, {}).get("window_days") for y in year],
                         name="days counted (window, all taxa)", marker_color="#ccc", opacity=0.5,
                         xaxis="x3", yaxis="y4"))  # fmt: skip
    fig.add_hline(y=0.5, line=dict(color="#b71c1c", dash="dot", width=1), row=3, col=1)
    for y, lab in PROTOCOL.items():
        fig.add_vline(x=y - 0.5, line=dict(color="#aaa", dash="dot", width=1))
    fig.add_vline(x=start - 0.5, line=dict(color="#b71c1c", width=1.5))
    fig.update_layout(height=FIG_H + 160, barmode="overlay",
                      legend=dict(orientation="h", y=-0.08),
                      yaxis4=dict(overlaying="y3", anchor="x3", side="right", showgrid=False,
                                  title="days"))  # fmt: skip
    fig.update_yaxes(title="birds", row=1, col=1)
    fig.update_yaxes(title="birds / day eq.", row=2, col=1)
    fig.update_yaxes(title="mean c", range=[0, 1.05], row=3, col=1)
    return panel(fig, "Counts and effort per year",
                 f"Raw, no model. Red line: start year {start} (grey bars before it are not "
                 "compared). Ratio index: Σ birds / Σ c over counted days (none below c = 0.5). "
                 "Dotted: protocol changes " + ", ".join(f"{y} {t}" for y, t in PROTOCOL.items())
                 + ".")  # fmt: skip


def profile_panel(sp: dict) -> str:
    """The time-of-day profile p(h | doy) the coverage is computed with."""
    pr = sp["profile"]
    z = np.array(pr["p"], float).T
    fig = go.Figure(go.Heatmap(x=pr["doy"], y=list(range(24)), z=z, colorscale="Blues",
                               colorbar=dict(title="p(h)", tickformat=".0%"),
                               customdata=[[doy_label(d, "%d-%b") for d in pr["doy"]]] * 24,
                               hovertemplate="%{customdata}, %{y}:00 solar: %{z:.1%}"
                                             "<extra></extra>"))  # fmt: skip
    doy_axis(fig)
    fig.update_yaxes(title="solar hour (12 = solar noon)", range=[3.5, 20.5])
    fig.update_layout(height=FIG_H)
    source = {"own": "its own", "uniform": "uniform over daylight"}.get(pr["source"],
                                                                         f"group: {pr['source']}")  # fmt: skip
    return panel(fig, f"Time-of-day profile ({source})",
                 "The share of a full day's passage in each solar hour, by date: a counted day's "
                 "coverage c is its sum over the hours counted. Fitted by the time-of-day GAM on "
                 "the days timed to the hour; a group's when the taxon has too few timed birds.")  # fmt: skip


# --- model components -------------------------------------------------------------------------


def components_panel(sp: dict) -> str:
    """The fitted season in the first and last year, and each year's expectation with and without
    its episodes over the counted birds per full day (one year at a time)."""
    t = sp["trend"]
    s, e = t["season"], t["episodes"]
    fig = make_subplots(rows=1, cols=2, column_widths=[0.35, 0.65], horizontal_spacing=0.08,
                        subplot_titles=("season: first and last year (share per day)",
                                        "one year: expected, with and without episodes"))  # fmt: skip
    for i, y in enumerate(k for k in s if k != "doy"):
        v = arr(s[y])
        fig.add_trace(go.Scatter(x=s["doy"], y=v / np.nansum(v), name=y, mode="lines",
                                 line=dict(color=COLORS[i], dash="dot" if i == 0 else "solid")),
                      1, 1)  # fmt: skip
    d = sp["days"]
    dates = np.array(d["date"], dtype="datetime64[D]")
    yr = dates.astype("datetime64[Y]").astype(int) + 1970
    # day of year in a non-leap year, as the trend's season days
    doy = (dates - dates.astype("datetime64[Y]")).astype(int) + 1
    leap = (yr % 4 == 0) & ((yr % 100 != 0) | (yr % 400 == 0))
    doy = np.where(leap & (doy > 59), doy - 1, doy)
    adj = arr(d["adjusted"])
    fd = t["days"]
    f_dates = np.array(fd["date"], dtype="datetime64[D]")
    f_yr = f_dates.astype("datetime64[Y]").astype(int) + 1970
    f_doy = (f_dates - f_dates.astype("datetime64[Y]")).astype(int) + 1
    f_leap = (f_yr % 4 == 0) & ((f_yr % 100 != 0) | (f_yr % 400 == 0))
    f_doy = np.where(f_leap & (f_doy > 59), f_doy - 1, f_doy)
    f_tot, f_lo, f_hi, f_c = arr(fd["total"]), arr(fd["q10"]), arr(fd["q90"]), arr(fd["c"])
    years = sorted(e["base"], key=int)
    first = len(fig.data)
    for y in years:
        m = (yr == int(y)) & np.isfinite(adj)
        vis = y == years[-1]
        k = f_yr == int(y)
        fig.add_trace(go.Scatter(x=f_doy[k], y=f_tot[k], mode="markers", visible=vis,
                                 **legend(sp, "trend.days.total", "filled day (mean, 80%)"),
                                 error_y=dict(type="data", symmetric=False,
                                              array=np.clip(f_hi[k] - f_tot[k], 0, None),
                                              arrayminus=np.clip(f_tot[k] - f_lo[k], 0, None),
                                              color="#4a90d9", thickness=1, width=0),
                                 marker=dict(color="#4a90d9", size=4, symbol="diamond"),
                                 customdata=f_c[k],
                                 hovertemplate="filled %{y:.0f} (c = %{customdata:.2f})"
                                               "<extra></extra>"), 1, 2)  # fmt: skip
        fig.add_trace(go.Scatter(x=doy[m], y=adj[m], mode="markers", visible=vis,
                                 name="counted / c", marker=dict(color="#555", size=5)), 1, 2)  # fmt: skip
        fig.add_trace(go.Scatter(x=e["doy"], y=e["base"][y], mode="lines", visible=vis,
                                 name="without episodes",
                                 line=dict(color="#555", dash="dash")), 1, 2)  # fmt: skip
        with_ep = arr(e["base"][y]) * np.exp(arr(e["episode"][y]))  # `episode` is on the log scale
        fig.add_trace(go.Scatter(x=e["doy"], y=with_ep, mode="lines", visible=vis,
                                 name="with episodes", line=dict(color="#e66c00")), 1, 2)  # fmt: skip
    n = len(fig.data)
    buttons = []
    for k, y in enumerate(years):
        vis = [True] * first + [first + 4 * k <= j < first + 4 * k + 4 for j in range(first, n)]
        buttons.append(dict(label=y, method="update", args=[{"visible": vis}]))
    fig.update_layout(height=FIG_H + 20, legend=dict(orientation="h", y=-0.15),
                      updatemenus=[dict(buttons=buttons, active=len(years) - 1, x=1, y=1.18,
                                        xanchor="right")])  # fmt: skip
    doy_axis(fig, 1, 1)
    doy_axis(fig, 1, 2)
    fig.update_yaxes(tickformat=".1%", row=1, col=1)
    fig.update_yaxes(title="birds per full day", row=1, col=2)
    return panel(fig, f"Model components (θ = {t['theta']:.2f}, κ = {t['kappa']:.1f})",
                 "Left: the smooth season (with its change of timing) of the first and last "
                 "year, each summing to 1. Right: a year's expected birds per full day without "
                 "(trend, year level, season, shift) and with its weather episodes, over the "
                 "counted birds / c (days with c ≥ 0.5) and each counted day's gap-filled total "
                 "(blue, any c: the count plus the hours not counted, drawn from the model). θ: "
                 "day-to-day over-dispersion; κ: hour-to-hour (flocks).", top=75,
                 strip=strip(sp, ["trend.days.total"], counts="counted / c"))  # fmt: skip


# --- pages ------------------------------------------------------------------------------------


# The page follows the pipeline: what defileViz does with it, then inputs, model, validation and
# the blocks derived from them.
PAGE_SECTIONS = {
    "display": "On defileViz",
    "config": "1 Configuration",
    "counts": "2 Counts and effort",
    "daytime": "3 Time of day",
    "season": "4 Season day by day",
    "trend": "5 Trend model",
    "validation": "6 Validation",
    "derived": "7 Derived blocks",
}
SECTION_NAMES = {
    "passage": "Passage at the Défilé",
    "evolution": "Changes in passage",
    "particularites": "Distinctive observations",
}


def accounts_panel(sp: dict) -> str:
    a = sp["accounts"]
    if not a:
        return "<p class='note'>No written account.</p>"
    general = (
        "".join(
            f"<p><b>{SECTION_NAMES[k]}.</b> {html.escape(v['en'])}"
            f"<br><span class='note'>{html.escape(v['fr'])}</span></p>"
            for k, v in a["general"].items()
        )
        or "<p class='note'>No general account.</p>"
    )
    years = "".join(
        f"<li><b>{r['year']}</b> — {html.escape(r['en'])}"
        f"<br><span class='note'>{html.escape(r['fr'])}</span></li>"
        for r in a["years"]
    )
    return (
        f"<div class='panel'>{general}<details><summary>{len(a['years'])} year accounts"
        f"</summary><ul class='rep'>{years}</ul></details>"
        f"<p class='note'>[{a['method']}] English, French below each.</p></div>"
    )


def page(tx: dict, sp: dict, effort_annual: list) -> str:
    settings = "".join(
        f"<tr><td>{k}</td><td>{html.escape(json.dumps(sp['window'][k.split('_')[0]] if k.endswith('_window') else v['value']))}</td><td>{v['source']}</td>"
        f"<td>{html.escape(v.get('reason') or '')}</td></tr>"
        for k, v in sp["settings"].items()
    )
    links = " ".join(
        f'<a class="lk" href="{u}" target="_blank">{LINK_NAMES.get(k, k)}</a>'
        for k, u in sp["links"].items()
    )
    diag = sp["diagnostics"]
    flags = "".join(f"<span class='flag'>{f}</span>" for f in diag["flags"]) or "none"
    no_trend = "<p class='note'>No trend.</p>"
    toc = " · ".join(f"<a href='#{k}'>{v}</a>" for k, v in PAGE_SECTIONS.items())
    body = f"""
<h2>{html.escape(tx['english_name'])} <small style="font-weight:normal;color:#777">
<i>{html.escape(str(tx['scientific_name']))}</i> · {html.escape(str(tx['french_name']))} ·
tier {tx['tier']} · from {tx['start_year']}</small></h2>
<div>{links}</div>
<div class="toc">{toc}</div>
<h3 id="display">On defileViz</h3>{display_panel(sp, tx)}
<h3 id="config">1. Configuration</h3><div class="panel">flags: {flags}<br><span class="note">
{html.escape(json.dumps({k: v for k, v in diag.items() if k != 'flags'}))}</span>
<table><tr><th>setting</th><th>value</th><th>source</th><th>reason</th></tr>{settings}</table></div>
<h3 id="counts">2. Counts and effort</h3>{effort_panel(sp, effort_annual)}
<h3 id="daytime">3. Time of day and coverage</h3>{profile_panel(sp)}{daytime_panel(sp)}
<h3 id="season">4. Season day by day</h3>{season_panel(sp)}{chances_panel(sp)}
<h3 id="trend">5. Trend model</h3>{trend_panel(sp, effort_annual) + components_panel(sp) if sp['trend'] else no_trend}
<h3 id="validation">6. Validation</h3>{reliability_panel(sp) + benchmark_panel(sp) if sp['trend'] else no_trend}
<h3 id="derived">7. Derived blocks</h3>{key_numbers(sp)}
{demography_panel(sp, "age") + demography_panel(sp, "sex") or "<p class='note'>Age and sex: not enough birds aged or sexed.</p>"}
{records_panel(sp)}{accounts_panel(sp)}"""
    return wrap(tx["english_name"], f"<section>{body}</section>")


def wrap(title: str, body: str) -> str:
    js = plotly.offline.get_plotlyjs_version()
    return f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(title)}</title>
<style>{CSS}</style>
<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@{js}/plotly.min.js"></script>
</head><body><nav><a href="index.html">Explore QA</a><a href="compare.html">Species compared</a>
{html.escape(title)}</nav>{body}</body></html>"""


def gap_error(sp: dict) -> str:
    g = (sp.get("benchmark") or {}).get("gap")
    return format(g["abs_log_err"], ".3f") if g else ""


def index(taxa: list, species: dict) -> str:
    rows = []
    for tx in taxa:
        sp = species.get(tx["taxon_id"])
        if sp is None:
            continue
        d, k = sp["diagnostics"], sp["key_numbers"]
        rows.append((len(d["flags"]), tx, d, k, sp))
    rows.sort(key=lambda r: (-(r[1]["tier"] == "full"), -r[0], r[1]["english_name"]))
    cells = "".join(
        f"<tr><td><a href='{tx['taxon_id']}.html'>{html.escape(tx['english_name'])}</a></td>"
        f"<td>{tx['tier']}</td><td>{tx.get('group', '')}{' ★' if tx.get('highlight') else ''}</td>"
        f"<td>{tx.get('story') or ''}</td>"
        f"<td>{tx['start_year']}</td><td class='n'>{tx['birds'] or 0:,.0f}</td>"
        f"<td>{d['profile']}</td>"
        f"<td class='n'>{'' if d['coverage_recent'] is None else format(d['coverage_recent'], '.2f')}</td>"
        f"<td class='n'>{d.get('interval_ratio', '') and format(d.get('interval_ratio'), '.2f')}</td>"
        f"<td class='n'>{d['timed_days']}</td><td class='n'>{d['age_years']}</td>"
        f"<td>{'' if 'trend' not in k else format(k['trend']['change'], '+.0%')}</td>"
        + "".join(
            f"<td>{badge(sp['reliability'][c]['class'])}</td>"
            if sp.get("reliability")
            else "<td></td>"
            for c in CLAIMS
        )
        + f"<td class='n'>{gap_error(sp)}</td>"
        f"<td>{''.join(f'<span class=flag>{f}</span>' for f in d['flags'])}</td></tr>"
        for _, tx, d, k, sp in rows
    )
    head = ("<tr><th>taxon</th><th>tier</th><th>group (★ highlight)</th><th>story</th><th>from</th><th>birds</th><th>profile</th>"
            "<th>coverage</th><th>q90/q10</th><th>timed days</th><th>age yrs</th><th>trend</th>"
            "<th>totals</th><th>trend</th><th>season</th><th>gap error</th><th>flags</th></tr>")  # fmt: skip
    return wrap("index", f"<section><h2>Explore QA ({len(rows)} taxa)</h2><div class='panel'>"
                         f"<table>{head}{cells}</table></div></section>")  # fmt: skip


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--export", default=os.path.join(ROOT, "data", "explore"))
    ap.add_argument("--out", default=os.path.join(ROOT, "logs", "viewer"))
    ap.add_argument("--taxa", nargs="+", help="pages for these taxa only (ids or English names)")
    ap.add_argument("--tiers", nargs="+", default=["full"], help="tiers to give a page")
    args = ap.parse_args(argv)
    taxa = json.load(open(os.path.join(args.export, "taxa.json")))
    effort_annual = json.load(open(os.path.join(args.export, "effort.json")))["annual"]
    os.makedirs(args.out, exist_ok=True)
    species = {}
    for tx in taxa:
        path = os.path.join(args.export, "species", f"{tx['taxon_id']}.json")
        if tx["tier"] in args.tiers and os.path.exists(path):
            species[tx["taxon_id"]] = json.load(open(path))
    wanted = set(args.taxa or [])
    pages = [tx for tx in taxa if tx["taxon_id"] in species
             and (not wanted or tx["taxon_id"] in wanted or tx["english_name"] in wanted)]  # fmt: skip
    for tx in pages:
        with open(os.path.join(args.out, f"{tx['taxon_id']}.html"), "w", encoding="utf-8") as f:
            f.write(page(tx, species[tx["taxon_id"]], effort_annual))
    with open(os.path.join(args.out, "index.html"), "w", encoding="utf-8") as f:
        f.write(index(taxa, species))
    from species_compare import write  # imports this module: here, not at the top

    write(args.export, args.out, set(species))
    print(f"{len(pages)} pages, the index and the comparison -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
