#!/usr/bin/env python3
"""QA viewer of the Explore export: an index of every taxon with its diagnostics, and a page per
taxon with every block drawn (Plotly). It reads `data/explore/` only and computes nothing, so what
it shows is what defileViz gets; a panel that needs a calculation here belongs in the pipeline.

A developer's tool, focused on the model output: show everything that helps judge a fit, in full
detail. The visitors' page is defileViz's, trimmed and narrative; it does not have to look like
this one.

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


def panel(fig, title, note="", top=40):
    fig.update_layout(template="plotly_white", margin=dict(l=50, r=20, t=top, b=40),
                      title=dict(text=title, font=dict(size=14)))  # fmt: skip
    return (
        '<div class="panel">'
        + fig.to_html(full_html=False, include_plotlyjs=False, config={"responsive": True})
        + (f'<p class="note">{note}</p>' if note else "")
        + "</div>"
    )


def method(block) -> str:
    return f" <span style='font-size:11px;color:#888'>[{block['method']}]</span>" if block else ""


# --- panels -----------------------------------------------------------------------------------


def key_numbers(k: dict) -> str:
    items = []
    p = k["passage"]
    if p["q50"] is not None:
        items.append(("Main passage (80% of birds)", f"{doy_label(p['q10'])} – {doy_label(p['q90'])}",
                      f"median {doy_label(p['q50'])} ({p['source']})"))  # fmt: skip
    if "best_hours" in k:
        b = k["best_hours"]
        items.append(("Best hours", f"{b['from']:02d}:00 – {b['to']:02d}:00",
                      f"{b['share']:.0%} of a peak day"))  # fmt: skip
    if "typical_season" in k:
        t = k["typical_season"]
        items.append(("Typical season", f"{t['median']:,.0f}",
                      f"{t['min']:,.0f}–{t['max']:,.0f}, {t['years'][0]}–{t['years'][1]}"))  # fmt: skip
    if "trend" in k:
        t = k["trend"]
        items.append(("Trend", f"{t['change']:+.0%}", f"smooth, {t['from']} → {t['to']}"))
    if "chance" in k:
        c = k["chance"]
        items.append(("Chance of ≥1 / ≥10 in the main passage",
                      f"{c['at_least_1']:.0%} / {c['at_least_10']:.0%}", f"{c['days']} days"))  # fmt: skip
    if "record" in k:
        items.append(("Record day", f"{k['record']['count']:,.0f}", k["record"]["date"]))
    cells = "".join(
        f'<div class="kn"><div class="kl">{html.escape(a)}</div><div class="kv">{b}</div>'
        f'<div class="ks">{html.escape(str(c))}</div></div>'
        for a, b, c in items
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
    fig.add_trace(go.Scatter(x=a["year"], y=a["smooth_q97.5"], line=dict(width=0),
                             showlegend=False, hoverinfo="skip"), 1, 1)  # fmt: skip
    fig.add_trace(go.Scatter(x=a["year"], y=a["smooth_q2.5"], fill="tonexty", line=dict(width=0),
                             fillcolor="rgba(76,120,168,0.18)", name="smooth trend 95%",
                             hoverinfo="skip"), 1, 1)  # fmt: skip
    fig.add_trace(go.Scatter(x=a["year"], y=a["smooth"], name="smooth trend (typical year)",
                             line=dict(color=COLORS[0], width=2.5)), 1, 1)  # fmt: skip
    fig.add_trace(go.Scatter(
        x=a["year"], y=a["total"], mode="markers", name="season total (gap-filled, 80%)",
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
    return panel(fig, "Long-term trend" + method(t),
                 "Bars: birds counted. Points: counted + posterior fill (80%). Line: smooth trend. "
                 "Bottom: days counted in the window (all taxa).")  # fmt: skip


CLASS_STYLE = {"show": "#e6f4ea;color:#1e6b34", "caveat": "#fff4e0;color:#8a5300",
               "hide": "#fdecea;color:#b71c1c"}  # fmt: skip
CLAIMS = ("totals", "trend", "season")


def badge(cls: str, text: str | None = None) -> str:
    return (f"<span class='flag' style='background:{CLASS_STYLE[cls]}'>"
            f"{html.escape(text or cls)}</span>")  # fmt: skip


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
    share = np.array([[np.nan if v is None else v for v in row] for row in s["share"]], float)
    zmax = float(np.nanquantile(share, 0.98)) if np.isfinite(share).any() else 1
    count = [[None if v is None else f"{v:,.0f}" for v in row] for row in s["count"]]
    dates = [doy_label(x, "%d-%b") for x in s["doy"]]
    custom = [[[dates[j], count[i][j]] for j in range(len(dates))] for i in range(len(count))]
    fig = go.Figure(go.Heatmap(x=s["doy"], y=s["years"], z=share, zmin=0, zmax=zmax,
                               colorscale="YlOrRd", colorbar=dict(title="share<br>of year"),
                               customdata=custom,
                               hovertemplate="%{customdata[0]} %{y}: %{customdata[1]} birds "
                                             "(%{z:.1%} of the season)<extra></extra>"))  # fmt: skip
    p = {k: arr([r[k] for r in s["passage"]]) for k in s["passage"][0]}
    smooth = sp["trend"]["passage_q"] if sp["trend"] else None
    for col, dash, lab in (("q10", "dot", "10%"), ("q50", "solid", "50%"), ("q90", "dot", "90%")):
        if smooth:
            fig.add_trace(go.Scatter(x=[r[col] for r in smooth], y=[r["year"] for r in smooth],
                                     mode="lines", name=f"{lab} smooth", legendgroup=lab,
                                     line=dict(color="black", dash=dash, width=1.5)))  # fmt: skip
        fig.add_trace(go.Scatter(x=p[col], y=p["year"], mode="markers", name=f"{lab} that year",
                                 legendgroup=lab,
                                 marker=dict(color="white" if col == "q50" else "#333", size=7,
                                             line=dict(color="black", width=1),
                                             opacity=(0.3 + 0.7 * p["counted"]).tolist())))  # fmt: skip
    w = sp["window"]
    for i in (0, 1):  # window edges, where the model's differs from the default
        if w["model"][i] != w["default"][i]:
            for key, color in (("default", "#999"), ("model", "#1f77b4")):
                x = doy_of(w[key][i]) + (0.5 if i else -0.5)
                fig.add_vline(x=x, line=dict(color=color, dash="dash"))
    doy_axis(fig)
    fig.update_layout(height=FIG_H + 160, legend=dict(orientation="h", y=-0.1))
    title = "Phenology, every year: share of the season's birds per day"
    if smooth:
        title += (
            f" (median {smooth[-1]['q50'] - smooth[0]['q50']:+.0f} days since {smooth[0]['year']})"
        )
    return panel(fig, title + method(s),
                 "Empirical; blank = not counted. Dots: that year's 10/50/90% passage dates "
                 "(paler = more of the season uncounted). Lines: the smooth. Dashed: the model "
                 f"window {' - '.join(w['model'])} (blue) where it differs from the default "
                 f"{' - '.join(w['default'])} (grey); shown: {' - '.join(w['view'])}. "
                 f"Passage beyond counting: {', '.join(w['beyond_counting']) or 'no'}.")  # fmt: skip


def chances_panel(sp: dict) -> str:
    ch = sp["season"]["chances"]
    fig = go.Figure()
    for i, k in enumerate((1, 10, 100, 1000)):
        if f"at_least_{k}" in ch:
            fig.add_trace(go.Scatter(x=ch["doy"], y=ch[f"at_least_{k}"], mode="lines+markers",
                                     name=f"≥ {k}", line=dict(color=COLORS[i], shape="spline")))  # fmt: skip
    doy_axis(fig)
    fig.update_yaxes(tickformat=".0%", range=[0, 1], title="share of counted days")
    fig.update_layout(height=FIG_H)
    return panel(fig, f"Your chances ({ch['years'][0]}–{ch['years'][1]})",
                 "Share of well-counted days with at least N birds.")  # fmt: skip


def daytime_panel(sp: dict) -> str:
    d = sp["daytime"]
    if not d or d["hours"] is None:
        return "<p class='note'>No days timed to the hour.</p>"
    ch = d["change"]
    hist = [go.Bar(x=list(range(24)), y=d["hours"], marker_color="rgba(76,120,168,0.6)",
                   name="counted", hovertemplate="%{x}:00-%{x}:59: %{y:.1%}<extra></extra>")]  # fmt: skip
    if d["profile"] is not None:
        hist.append(go.Scatter(x=list(range(24)), y=d["profile"], line=dict(color="black"),
                               name="smooth profile (main passage)"))  # fmt: skip
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
        fig.update_yaxes(title="local hour", range=[5.5, 20.5], row=1, col=1)
        fig.update_xaxes(range=[5.5, 20.5], title="local hour", row=1, col=2)
        fig.update_yaxes(tickformat=".0%", row=1, col=2)
        fig.update_annotations(yshift=6, font_size=12)
    else:
        fig = go.Figure(hist)
        fig.update_xaxes(range=[5.5, 20.5], title="local hour", dtick=1)
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
                 "Bars: birds per counted hour, as shares, pooled over the season. " + change,
                 top=75 if ch["show"] else 40)  # fmt: skip


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
                      + method(a), note, top=75)  # fmt: skip


def records_panel(sp: dict) -> str:
    def note(n):
        ref = "" if n["ref"] is None else f"<i>{html.escape(n['ref'])}</i> — "
        return f"<p class='rn'>{ref}{html.escape(n['text'])}</p>"

    def notes(r):
        return "".join(note(n) for n in r["notes"])

    rows = "".join(
        f"<tr><td>{i + 1}</td><td>{r['date']}</td><td class='n'>{r['count']:,.0f}</td>"
        f"<td>{notes(r)}</td></tr>"
        for i, r in enumerate(sp["records"]["top_days"])
    )
    return (f"<div class='panel'><b>Top days</b><table><tr><th>#</th><th>date</th><th>birds</th>"
            f"<th>what was written that day</th></tr>{rows}</table></div>")  # fmt: skip


# --- pages ------------------------------------------------------------------------------------


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
    body = f"""
<h2>{html.escape(tx['english_name'])} <small style="font-weight:normal;color:#777">
<i>{html.escape(str(tx['scientific_name']))}</i> · {html.escape(str(tx['french_name']))} ·
tier {tx['tier']} · from {tx['start_year']}</small></h2>
<div>{links}</div>
<h3>Key numbers</h3>{key_numbers(sp['key_numbers'])}
<h3>Diagnostics</h3><div class="panel">flags: {flags}<br><span class="note">
{html.escape(json.dumps({k: v for k, v in diag.items() if k != 'flags'}))}</span>
<table><tr><th>setting</th><th>value</th><th>source</th><th>reason</th></tr>{settings}</table></div>
<h3>Accounts</h3>{accounts_panel(sp)}
<h3>Trend</h3>{trend_panel(sp, effort_annual) if sp['trend'] else "<p class='note'>No trend.</p>"}
<h3>Reliability</h3>{reliability_panel(sp) or "<p class='note'>No trend.</p>"}{benchmark_panel(sp) if sp['trend'] else ""}
<h3>Phenology</h3>{season_panel(sp)}
<h3>Visiting</h3>{chances_panel(sp)}{daytime_panel(sp)}
<h3>Age and sex</h3>{demography_panel(sp, "age") + demography_panel(sp, "sex") or "<p class='note'>Not enough birds aged or sexed.</p>"}
<h3>Records</h3>{records_panel(sp)}"""
    return wrap(tx["english_name"], f"<section>{body}</section>")


def wrap(title: str, body: str) -> str:
    js = plotly.offline.get_plotlyjs_version()
    return f"""<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(title)}</title>
<style>{CSS}</style>
<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@{js}/plotly.min.js"></script>
</head><body><nav><a href="index.html">Explore QA</a> {html.escape(title)}</nav>{body}</body></html>"""


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
        f"<td>{tx['tier']}</td><td>{tx['start_year']}</td><td class='n'>{tx['birds'] or 0:,.0f}</td>"
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
    head = ("<tr><th>taxon</th><th>tier</th><th>from</th><th>birds</th><th>profile</th>"
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
    print(f"{len(pages)} pages and the index -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
