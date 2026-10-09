# Authored content

What the Explore build reads that is written by people, not computed:

- `accounts/`: the written accounts per species and season (see its README).
- `highlights.tsv`: the species the defileViz picker opens on (`catalogue.load_highlights`). One
  row per taxon: `taxon_id`, `english_name` (for the reader; not read) and `source`. Seeded with
  the species of the Défilé's 2019 and 2020 papers that have a page; add or remove rows by hand.
