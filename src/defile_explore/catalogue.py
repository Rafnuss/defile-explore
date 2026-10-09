"""Where each taxon sits in defileViz's species picker: the `taxa.json` fields that sort, group and
filter the list (`catalogue`).

- `taxon_order`: eBird's taxonomic sequence (`TAXON_ORDER`, by `ebird_code`); a combined series
  takes its first member's.
- `group`: a fixed group by order (`GROUPS`), for the picker's filters; a taxon without an order
  (a slash between orders) takes the first order of its scientific name, else `OTHER`.
- `season_birds`: the median birds counted per season in the default window over the last
  `RECENT_YEARS` complete seasons: how common the taxon is now, to sort by.
- `story`: how much of its page rests on the trend: `full` (a trend whose totals are shown),
  `caveat` (shown with a caveat), `counts` (no trend, or its totals hidden: the counts only); None
  below the full tier, which has no page.
- `highlight`: the default list, curated in `content/highlights.tsv` (`HIGHLIGHTS_FILE`, read
  by `load_highlights`): seeded with the species of the Défilé's 2019 (raptors) and 2020 (other
  monitored species) papers that have a page. A full-tier taxon only; the story says how good its
  data are. A rule (species with a `full` story and 100 birds a season or more) was tried first:
  it left out Common Crane and Grey Heron for their caveated totals and put in Brambling and
  Common Reed Bunting.
"""

import os

import numpy as np
import pandas as pd

GROUPS = {  # group: orders
    "raptors": ("Accipitriformes", "Falconiformes"),
    "waterbirds": (  # storks, herons, cranes, cormorants, ducks, geese, gulls, waders
        "Ciconiiformes",
        "Pelecaniformes",
        "Gruiformes",
        "Suliformes",
        "Anseriformes",
        "Charadriiformes",
        "Phoenicopteriformes",
        "Gaviiformes",
        "Podicipediformes",
    ),
    "passerines": ("Passeriformes",),
}
OTHER = "other"  # pigeons, swifts, bee-eaters, woodpeckers, ...
RECENT_YEARS = 10
HIGHLIGHTS_FILE = os.path.join(
    os.path.dirname(__file__), "..", "..", "content", "highlights.tsv"
)  # taxon_id, english_name (for the reader), source
STORY = {"show": "full", "caveat": "caveat", "hide": "counts"}  # totals class: story


def group_of(order, scientific_name: str) -> str:
    """The taxon's `GROUPS` key by its order, else by the first order named in its scientific name
    ("Accipitriformes/Falconiformes sp."), else `OTHER`."""
    if not isinstance(order, str):
        order = str(scientific_name).replace("/", " ").split(" ")[0]
    return next((g for g, orders in GROUPS.items() if order in orders), OTHER)


def taxon_order(taxa: pd.DataFrame, taxonomy: pd.DataFrame, ebird: pd.DataFrame) -> pd.Series:
    """EBird's `TAXON_ORDER` of each taxon in `taxa`; a combined series takes its first
    member's."""
    seq = taxonomy.set_index("taxon_id")["ebird_code"].map(
        ebird.set_index("SPECIES_CODE")["TAXON_ORDER"]
    )
    own = taxa["taxon_id"].map(seq)
    members = taxa["members"] if "members" in taxa else pd.Series(None, index=taxa.index)
    first = members.map(lambda m: seq.reindex(m).min() if isinstance(m, list) else np.nan)
    return own.fillna(first)


def season_birds(annual: list[dict], last_year: int) -> float:
    """Median birds counted in the default window over the last `RECENT_YEARS` seasons."""
    a = pd.DataFrame(annual)
    if a.empty:
        return 0.0
    recent = a.loc[a["year"].between(last_year - RECENT_YEARS + 1, last_year), "window"]
    return float(recent.median()) if len(recent) else 0.0


def load_highlights(taxon_ids, path: str = HIGHLIGHTS_FILE) -> set[str]:
    """The curated taxon ids; an id not in `taxon_ids` is refused."""
    ids = pd.read_csv(path, sep="\t")["taxon_id"]
    unknown = set(ids) - set(taxon_ids)
    if unknown:
        raise ValueError(f"{path}: unknown taxa {sorted(unknown)}")
    return set(ids)


def story(tier: str, reliability: dict | None) -> str | None:
    if tier != "full":
        return None
    return STORY[reliability["totals"]["class"]] if reliability else "counts"


def catalogue(
    taxa: pd.DataFrame,
    taxonomy: pd.DataFrame,
    ebird: pd.DataFrame,
    pages: dict,
    last_year: int,
    highlights: set[str],
) -> pd.DataFrame:
    """`taxa` with the picker's fields, from each taxon's species file in `pages` (taxon_id:

    `{annual, reliability}`) and the curated `highlights` (`load_highlights`), sorted
    taxonomically.
    """
    t = taxa.copy()
    t["taxon_order"] = taxon_order(t, taxonomy, ebird)
    t["group"] = [group_of(o, s) for o, s in zip(t["order"], t["scientific_name"])]
    page = t["taxon_id"].map(pages)
    t["season_birds"] = [season_birds(p["annual"], last_year) if p else 0.0 for p in page]
    t["story"] = [story(tr, p["reliability"] if p else None) for tr, p in zip(t["tier"], page)]
    t["highlight"] = t["taxon_id"].isin(highlights) & (t["tier"] == "full")
    return t.sort_values(["taxon_order", "scientific_name"], na_position="last")
