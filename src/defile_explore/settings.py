"""Per-taxon settings: a default for every taxon by rule, and exceptions by hand.

Every setting has a rule computed from the data (`resolve`), so a new taxon or a new release needs
no editing. Where a rule is wrong for a taxon, `overrides.yaml` sets the value, with a `reason`
(required): the exceptions stay few, visible and justified. Each value carries its `source`
(`rule` or `override`), written to the species file, so the page can say why a taxon looks the
way it does.

Settings:

- `start_year`: first year the taxon's counts are comparable (`export.start_year`).
- `trend`: whether the trend GAM is fitted (full tier, a species or a group with members).
- `age_years`, `sex_years`: years shown for age and sex (`demography.usable_years`). An override is
  a list of years, or `{from: Y}` / `{from: Y, to: Z}`: every year of the range with any bird aged
  (or sexed).
- `model_window`, `view_window`: the season days the trend is fitted on and the season panels
  show (`window.windows`). An override is `[MM-DD, MM-DD]`, inclusive.
- `links`: external pages. eBird, eBird Status and Trends, Birds of the World, EBBA2 and
  Trektellen by rule (from the taxonomy's codes); others (`vogelwarte`, `migration_atlas`) only
  from `overrides.yaml`, where they need no reason.
"""

import os
from dataclasses import dataclass

import pandas as pd
import yaml

from defile_explore import demography as G
from defile_explore import window as W
from defile_explore.export import has_members
from defile_explore.trend import TREND_RANKS

OVERRIDES_FILE = os.path.join(os.path.dirname(__file__), "overrides.yaml")
SETTINGS = ("start_year", "trend", "age_years", "sex_years", "model_window", "view_window")
WINDOW_SETTINGS = ("model_window", "view_window")  # set after the start year, from the data
LINK_KEYS = ("vogelwarte", "migration_atlas")  # set by hand only
TREKTELLEN_SITE = 2422  # Défilé de l'Écluse on trektellen.org
LINK_TEMPLATES = {
    "ebird": "https://ebird.org/species/{ebird}",
    "ebird_status": "https://science.ebird.org/status-and-trends/species/{ebird}/abundance-map",
    "birds_of_the_world": "https://birdsoftheworld.org/bow/species/{ebird}/cur/introduction",
    "ebba2": "https://ebba2.info/maps/species/{binomial}/ebba2/abundance/",
    "trektellen": "https://www.trektellen.org/species/graph/3/" + str(TREKTELLEN_SITE) + "/{tk}/0",
    "trektellen_day": "https://www.trektellen.org/count/view/" + str(TREKTELLEN_SITE) + "/{date}",
    "vogelwarte": "https://www.vogelwarte.ch/en/birds-of-switzerland/{value}/",
    "migration_atlas": "https://migrationatlas.org/node/{value}",
}


@dataclass
class Setting:
    value: object
    source: str = "rule"  # or "override"
    reason: str | None = None

    def as_dict(self) -> dict:
        d = {"value": self.value, "source": self.source}
        return d | ({"reason": self.reason} if self.reason else {})


def load_overrides(path: str = OVERRIDES_FILE) -> dict:
    """`{taxon_id: {setting: value, ..., "reason": str}}`, checked: known settings only, and a
    reason for any setting other than links."""
    with open(path, encoding="utf-8") as f:
        raw = yaml.safe_load(f) or {}
    for taxon_id, o in raw.items():
        unknown = set(o) - set(SETTINGS) - {"links", "reason"}
        if unknown:
            raise ValueError(f"{OVERRIDES_FILE}: {taxon_id}: unknown setting(s) {sorted(unknown)}")
        if set(o) & set(SETTINGS) and not o.get("reason"):
            raise ValueError(f"{OVERRIDES_FILE}: {taxon_id}: an override needs a `reason`")
        for name in set(o) & set(WINDOW_SETTINGS):
            W.from_dates(o[name])  # a malformed date fails here, not mid-build
        bad = set(o.get("links", {})) - set(LINK_KEYS)
        if bad:
            raise ValueError(f"{OVERRIDES_FILE}: {taxon_id}: unknown link(s) {sorted(bad)}")
    return raw


def year_range(rows: pd.DataFrame, field: str, spec) -> list[int]:
    """An override of `age_years`/`sex_years`: a list as given, or the years of `{from, to}` with
    any bird given a `field`."""
    if isinstance(spec, list):
        return [int(y) for y in spec]
    years = sorted(rows.loc[rows[field].notna(), "date"].dt.year.unique())
    return [int(y) for y in years if spec.get("from", 0) <= y <= spec.get("to", 9999)]


def rule_links(taxon: dict) -> dict:
    """Links derived from the taxonomy: by eBird code and Trektellen id, and EBBA2 by binomial."""
    out = {}
    if isinstance(taxon.get("ebird_code"), str) and taxon["taxon_rank"] == "species":
        for key in ("ebird", "ebird_status", "birds_of_the_world"):
            out[key] = LINK_TEMPLATES[key].format(ebird=taxon["ebird_code"])
    if taxon["taxon_rank"] == "species" and len(str(taxon["scientific_name"]).split()) == 2:
        out["ebba2"] = LINK_TEMPLATES["ebba2"].format(
            binomial=taxon["scientific_name"].replace(" ", "-")
        )
    tk = taxon.get("trektellen_species_id")
    if pd.notna(tk):  # unidentified groups ("harrier sp.") list several ids: the first
        out["trektellen"] = LINK_TEMPLATES["trektellen"].format(
            tk=int(float(str(tk).split(",")[0]))
        )
    return out


def resolve(
    taxon: dict, demography_rows: pd.DataFrame, counted: pd.Series, overrides: dict, last_year: int
):
    """The settings of one taxon (`taxon`: its `taxa.json` entry plus `ebird_code` and
    `trektellen_species_id`), as `{name: Setting}` plus `links`.

    Years after `last_year` (a season in progress) are never shown.
    """
    demography_rows = demography_rows[demography_rows["date"].dt.year <= last_year]
    s = {
        "start_year": Setting(int(taxon["start_year"])),
        "trend": Setting(
            taxon["tier"] == "full"
            and (taxon["taxon_rank"] in TREND_RANKS or has_members(taxon.get("members")))
        ),
        "age_years": Setting(G.usable_years(demography_rows, "age", counted)),
        "sex_years": Setting(G.usable_years(demography_rows, "sex", counted)),
    }
    o = overrides.get(taxon["taxon_id"], {})
    for name in SETTINGS:
        if name in o and name not in WINDOW_SETTINGS:
            value = o[name]
            if name in ("age_years", "sex_years"):
                value = year_range(demography_rows, name.split("_")[0], value)
            s[name] = Setting(value, "override", o["reason"])
    links = rule_links(taxon)
    for key, value in o.get("links", {}).items():
        links[key] = LINK_TEMPLATES[key].format(value=value)
    return s, links


def resolve_windows(taxon_id: str, rule: dict, overrides: dict) -> dict:
    """`model_window` and `view_window` (season days) from `rule` (`window.windows`), or from
    `overrides`; a view window never narrower than the model window."""
    o = overrides.get(taxon_id, {})
    s = {}
    for name in WINDOW_SETTINGS:
        if name in o:
            s[name] = Setting(W.from_dates(o[name]), "override", o["reason"])
        else:
            s[name] = Setting(rule[name.split("_")[0]])
    m, v = s["model_window"].value, s["view_window"].value
    s["view_window"].value = (min(m[0], v[0]), max(m[1], v[1]))
    return s
