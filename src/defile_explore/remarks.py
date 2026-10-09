"""What the counters and the reports wrote about a day, for the taxon's record days.

Two sources in the release, both in French: a survey's `remark` (the counters' account of the day,
on Trektellen since 2013) and a count's `remark` (on one taxon's entry: plumage details, or the
report paragraphs about that day the dataset attached to it, each opening with its source in
brackets, `[Rapport annuel 2020; Contexte journalier 2020-10-14, Milan royal; ...]`).

`day_remarks` splits both into notes, one row per paragraph: `taxon_id` (None for a survey's),
`date`, `source` (`survey` or `count`), `ref` (the bracketed source's first part, e.g.
`Rapport annuel 2020`, or None) and `text`. A count remark's field label (`details:`, `remark:`,
...) is dropped, and so are empty paragraphs and the dataset's own placeholders (`PLACEHOLDERS`).
`notes_of` gives a day's notes, the day's survey notes first, each text once.
"""

import re

import pandas as pd

LABEL = re.compile(r"^(details?|remark|comment)\s*:\s*", re.IGNORECASE)
REF = re.compile(r"^\[([^\];]+)[^\]]*\]\s*")
PARAGRAPH = re.compile(r"\n\s*\n")
PLACEHOLDERS = ("An explicit 'no species' entry",)  # written by the dataset, not the counters
COLUMNS = ["taxon_id", "date", "source", "ref", "text"]


def split(remark: str) -> list[tuple[str | None, str]]:
    """A remark as `(ref, text)` paragraphs."""
    out = []
    for p in PARAGRAPH.split(remark):
        p = LABEL.sub("", p.strip()).strip()
        m = REF.match(p)
        ref, p = (m.group(1).strip(), p[m.end() :].strip()) if m else (None, p)
        if p and not p.startswith(PLACEHOLDERS):
            out.append((ref, p))
    return out


def _notes(frame: pd.DataFrame, source: str) -> pd.DataFrame:
    rows = [
        (t, d, source, ref, text)
        for t, d, r in frame[["taxon_id", "date", "remark"]].itertuples(index=False)
        for ref, text in split(r)
    ]
    return pd.DataFrame(rows, columns=COLUMNS)


def day_remarks(surveys: pd.DataFrame, counts: pd.DataFrame) -> pd.DataFrame:
    """The notes of every survey and count with a remark (see the module docstring)."""
    s = surveys.loc[surveys["remark"].notna(), ["date", "remark"]].assign(taxon_id=None)
    c = counts.loc[counts["remark"].notna(), ["taxon_id", "date", "remark"]]
    return pd.concat([_notes(s, "survey"), _notes(c, "count")], ignore_index=True)


def notes_of(remarks: pd.DataFrame, date) -> list[dict]:
    """The notes of one day (`remarks` already cut to the taxon), survey notes first, each text
    once: `[{source, ref, text}]`."""
    d = remarks[remarks["date"] == date].sort_values("source", ascending=False, kind="stable")
    d = d.drop_duplicates("text")
    return [
        {"source": s, "ref": r if isinstance(r, str) else None, "text": t}
        for s, r, t in d[["source", "ref", "text"]].values
    ]
