"""The written accounts of each taxon, in French and English, for its species file.

Authored in `content/accounts/` (its README has the editorial rules): `species-sections.tsv`, a
general account in up to three optional sections per taxon (`SECTIONS`), and `year-accounts.tsv`,
one short account of a taxon's season per year. Both are keyed by taxon id and carry `text_fr` and
`text_en`. The build reads only these two files: the report and paper extracts they were written
from are the editors' sources, not the page's.

`accounts_block` gives a taxon its own accounts only; a combined series has none of its own, and
its members' accounts stay on their own pages.
"""

import hashlib
import os

import pandas as pd

METHOD = "accounts@1"
ACCOUNTS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "content", "accounts")
SECTIONS_FILE = "species-sections.tsv"
YEARS_FILE = "year-accounts.tsv"
SECTIONS = ("passage", "evolution", "particularites")  # in page order
LANGUAGES = ("fr", "en")  # French is the editorial language
TEXT_COLUMNS = [f"text_{lang}" for lang in LANGUAGES]


def read(path: str, keys: list[str]) -> pd.DataFrame:
    t = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, encoding="utf-8-sig")
    missing = set(keys + TEXT_COLUMNS) - set(t.columns)
    if missing:
        raise ValueError(f"{path}: missing column(s) {sorted(missing)}")
    if t.duplicated(keys).any():
        raise ValueError(f"{path}: duplicated {keys}: {t.loc[t.duplicated(keys), keys].values}")
    empty = (t[TEXT_COLUMNS].apply(lambda c: c.str.strip()) == "").any(axis=1)
    if empty.any():
        raise ValueError(f"{path}: a text is missing in one language: {t.loc[empty, keys].values}")
    return t


def load_accounts(folder: str = ACCOUNTS_DIR, taxon_ids=None) -> dict:
    """`{"sections": DataFrame, "years": DataFrame}`, checked: unique keys, known sections, both
    languages present, and (given `taxon_ids`) every key a known taxon."""
    sections = read(os.path.join(folder, SECTIONS_FILE), ["key", "section"])
    years = read(os.path.join(folder, YEARS_FILE), ["year", "key"])
    bad = set(sections["section"]) - set(SECTIONS)
    if bad:
        raise ValueError(f"{SECTIONS_FILE}: unknown section(s) {sorted(bad)}")
    if taxon_ids is not None:
        unknown = (set(sections["key"]) | set(years["key"])) - set(taxon_ids)
        if unknown:
            raise ValueError(f"accounts for unknown taxa: {sorted(unknown)}")
    years["year"] = years["year"].astype(int)
    return {"sections": sections, "years": years}


def taxon_rows(accounts: dict, taxon_id: str) -> dict:
    """The rows of one taxon (what a job carries)."""
    return {k: t[t["key"] == taxon_id] for k, t in accounts.items()}


def accounts_block(rows: dict, last_year: int) -> dict | None:
    """The `accounts` block: `general` (`{section: {fr, en}}`, in `SECTIONS` order) and `years`
    (`[{year, fr, en}]`, newest first, none after `last_year`); None without any account."""
    s = rows["sections"].set_index("section")
    general = {
        sec: {lang: s.at[sec, f"text_{lang}"] for lang in LANGUAGES}
        for sec in SECTIONS
        if sec in s.index
    }
    y = rows["years"][rows["years"]["year"] <= last_year].sort_values("year", ascending=False)
    years = [
        {"year": int(r["year"]), **{lang: r[f"text_{lang}"] for lang in LANGUAGES}}
        for _, r in y.iterrows()
    ]
    if not general and not years:
        return None
    return {"method": METHOD, "general": general, "years": years}


def fingerprint(folder: str = ACCOUNTS_DIR) -> dict:
    """SHA-256 of the two authored files, for the export's manifest."""
    return {
        name: hashlib.sha256(open(os.path.join(folder, name), "rb").read()).hexdigest()
        for name in (SECTIONS_FILE, YEARS_FILE)
    }
