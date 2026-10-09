# Bilingual accounts for the Explore page

Prepared on 9 October 2026. These are editorial drafts for defileViz, based on the maintained annual and paper text extracts. French is the reference language; each English version follows the same observations, historical scope and level of uncertainty. Website prose states observations directly without citing reports. The original source extracts and release tables are copied unchanged from defile-dataset into `data/count/dataset/`. The accounts and their workflow are maintained here in defile-explore.

## Read the drafts

- [General accounts in French](reading/species-fr.md) and [English](reading/species-en.md).
- [Year accounts in French](reading/years-fr.md) and [English](reading/years-en.md).
- [General account CSV](review/species-accounts.csv): one row per taxon, paired optional `passage`, `evolution` and `particularites` fields in French and English.
- [Year account CSV](review/species-year-accounts.csv): one row per source taxon and year, with paired `text_fr` and `text_en` fields.

There are 63 populated general accounts and 563 populated year accounts in each language. All 615 keyed annual species entries have been reviewed, covering the available seasons from 2004 to 2025. There is no 2022 annual source in these extracts. The 52 annual entries containing only numerical details, dates or brief age/sex lists have empty narrative fields and status `numerical_or_list_only`; empty text does not mean no birds were observed. All 100 paper species/group entries were also reviewed.

The general CSV covers the 112 individual taxon keys present in the sources, including unidentified categories where applicable. Of these, 49 have no general prose because the source offers too little material for a useful synthesis. Three additional mixed-species paper groupings are retained in the coverage ledger, rather than exported as species. No annual account has been manufactured from a paper-only observation.

## Editorial approach

Annual accounts were drafted before general syntheses. They prioritise the character of the season, the sequence of meaningful movements and explanations explicitly supported by the source. Brief source material yields brief prose. Exact totals, peak counts and individual passage dates are generally omitted because the website will present those numerically. A rare encounter can remain an event-led account where it is the only available story.

The optional general sections describe passage at the Défilé, historical changes, and distinctive behaviour or observations. They remain absent where the source provides nothing useful. Historical trends retain their time frame: a conclusion through 2017 or 2019 is not silently extended to the present. Recent annual observations qualify older patterns where necessary. General statements about biology come from the supplied sources; no external natural-history material was added.

Weather explanations retain the original uncertainty. Temporal coincidence is not changed into a demonstrated cause, and a poor visible passage is not treated as proof of population decline. Improvements in dawn, dusk and distant-sector coverage are retained where relevant. Shared pigeon, swallow and corvid accounts are handled as shared observations, without allocating combined totals to individual species. Unidentified categories retain their meaning.

## Data checks and discrepancies

[Data checks](review/data-checks.csv) compare annual published totals with the current release; [coverage](review/coverage.csv) records every source species-year and general taxon/group reviewed. Coverage `source_rows` identifies the relevant original entries. General entries list the available same-taxon and explicitly shared inputs used during synthesis, rather than sentence-level citations. Annual `context_rows` lists available monitoring, weather and results entries for that season; it is supporting context, not a claim that each entry supplies a species-specific explanation.

The checks sum numerical values for `count_category=normal` over whole calendar years in Europe/Paris time. They include observations outside the official monitoring dates and hours when present in the release. Recorded numerical estimates or bounds are retained as recorded; this is not an estimate of unseen passage. Presence-only rows without a quantity are excluded, missing entries remain unknown, and counts are not adjusted for effort. Unidentified and identified taxa are not combined. Published totals can therefore legitimately differ from these sums. Differences are flags for review, not evidence that either source is wrong.

The descriptive rank compares available calendar-year totals from 1993 through the account year only. Later years are excluded. It is a cross-check, not a population trend or a substitute for an official-season record. Neither ranks nor totals are injected automatically into the prose.

Examples of decisions made during editing:

- Great Egret: the 2008 text claims a record, but the 2005 published total is higher. The record wording is omitted. The 2004 record claim is also omitted because the release has a higher earlier calendar-year total; the scope difference remains unresolved.
- Black Kite: the 2014 annual, paper and release totals differ. The account retains the exceptional passage and its seasonal pattern without choosing a precise figure. The 2025 annual total also differs from the calendar-year release sum.
- White Stork: 2025 is described as exceptional. Its published total and the calendar-year release total differ, so no precise figure is repeated. The older pilot’s record wording is superseded by the current draft.
- Honey Buzzard: the 2016/2017 annual texts describe decline since 2000, whereas the paper finds no clear trend over 1993–2017. The general account preserves that historical scope, alongside the lower levels noted in recent seasons.
- Grey Heron and Black Stork: some 2025 subtotals do not reconcile with their stated seasonal totals. The prose retains the sequence and observations without reproducing those figures or calculating age proportions.
- Rare birds: historical identification uncertainty and references to pending rarities acceptance are retained. Current committee decisions were not researched, so these drafts should not be read as a newly verified acceptance list.
- Apparent transcription problems or inconsistent dates are not turned into detailed narratives. The source extracts remain intact for review.

## Location and source snapshots

`content/accounts/` is versioned editorial content, outside the ignored `data/` directory. The two authored TSV files sit at its top; they are what the Explore build reads (`src/defile_explore/accounts.py`), into the `accounts` block of each species file. `reading/` holds the bilingual Markdown reading copies and `review/` the editorial records (coverage ledger, data checks, CSV exports, manifest, source snapshots, the pilot), all regenerated by `scripts/accounts/`. Those scripts read the original `report_text.csv`, `paper_text.csv`, count, survey and taxonomy tables from `data/count/dataset/`; the build does not.

`review/reference/annual-totals.csv` preserves the published totals used for editorial checks, and `review/reference/source_taxa.csv` preserves the original-name crosswalk used to resolve French names. These are unchanged supporting snapshots from defile-dataset, not additional observation tables. Update them deliberately when refreshing the source release. The paper-only Common Hoopoe name was resolved against the Clements 2025 reference checklist during the original drafting.

## Files and reproducibility

Edit [year-accounts.tsv](year-accounts.tsv) and [species-sections.tsv](species-sections.tsv). The CSV and Markdown files in `review/` and `reading/` are generated copies. The statuses are editorial metadata; the website shows non-empty text and omits empty sections, without the status or source-reference fields. The Explore build checks the TSVs again (unique keys, known sections and taxa, both languages filled), so a malformed edit fails the build rather than reaching the page.

From the repository root, run:

```sh
python3 scripts/accounts/check_website_account_totals.py
python3 scripts/accounts/build_website_accounts.py
```

The prose was selected, restructured and rephrased with AI assistance in a supervised editorial workflow: annual French drafts first, general French syntheses second, then paired English translations. The TSV files preserve that authored result. These scripts reproduce the exports and checks; they do not generate or translate prose automatically. To update an account, reread the relevant original text, check numerical claims where useful, edit both language fields in the TSV, and rebuild. Never infer a new biological explanation from a changed total alone.

The exporter checks unique taxon/year and taxon/section keys, source membership and complete language pairs. [manifest.json](review/manifest.json) records the counts, source years and SHA-256 fingerprints of the inputs used. Rebuilding after changing source data requires a fresh editorial review; it does not rewrite the authored narratives.

[pilot-fr.md](review/pilot-fr.md) and [pilot-editorial.md](review/pilot-editorial.md) preserve the earlier style pilot. The full drafts linked above supersede that pilot and incorporate the subsequent data checks. These exports do not change the released dataset schema.
