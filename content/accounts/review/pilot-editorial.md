# Website species accounts — pilot editorial record

Prepared on 9 October 2026. Display copy is in [pilot-fr.md](pilot-fr.md). This pilot is for agreement on prose and structure, not publication or a complete species review. Maintained source CSVs are unchanged. French is the editorial reference version; translate into natural English after the French pilot is agreed, retaining the same claims, qualifiers and historical scope. Pair future language versions by taxon ID, year and section.

## Editorial rules

- Select, restructure and rephrase the accepted source text; do not add outside biology or infer a cause from an annual total.
- State observations directly, without mentioning reports in website prose. Preserve uncertainty in the source through words such as « probablement », « peut » and « semble ».
- Omit precise totals, peak counts and exact passage dates when the website already supplies them. Keep meaningful seasonal timing, historical periods and comparisons.
- Draft annual accounts before the general species synthesis. Use relevant annual context and papers; do not treat repeated passages as independent evidence.
- Keep an account short when the source is short. Leave optional species sections absent when they add no useful material.
- Avoid applying species-specific interpretations to mixed-species totals. Keep source disagreements visible in editorial notes rather than silently resolving them.

## Source rows

Sources are `raw/reports/report_text.csv` and `raw/reports/paper_text.csv`, relative to the repository root. Annual rows below all have `category=species`; paper rows also have `category=species`. Keys identify rows together with year or source_id.

| Species         | Key              | Annual accounts drafted | Paper rows used for general account    | Other annual material used for general account |
| --------------- | ---------------- | ----------------------- | -------------------------------------- | ---------------------------------------------- |
| Bondrée apivore | avibase-ED5A7E8F | 2014, 2016, 2025        | defile-paper-1996-I; defile-paper-2019 | 2025                                           |
| Milan noir      | avibase-06D9A2C8 | 2014, 2016, 2025        | defile-paper-1996-I; defile-paper-2019 | 2025                                           |
| Cigogne blanche | avibase-28825494 | 2017, 2025              | defile-paper-2020                      | 2025                                           |
| Pigeon ramier   | avibase-760F307A | 2019, 2025              | defile-paper-2020                      | 2014, 2025                                     |
| Geai des chênes | avibase-D8D10F2C | 2010, 2012              | defile-paper-2020                      | 2012                                           |

Monitoring context was read for 2014, 2016, 2017, 2019 and 2025, including the distinction between identified and unidentified pigeons in the 2025 methods. The 2025 weather entries and the 2012 weather/results context were also read. No general weather condition was assigned to a species as a cause unless its species entry made that connection.

## Decisions and limits

- Honey Buzzard: the 2016/2017 annual entries describe a decline since 2000, while the 2019 paper finds it difficult to identify a trend over 1993–2017. The general pilot uses the paper’s explicit time frame and interpretation. The short 2016 account retains its season-specific explanation without repeating that trend. The 2025 entry contains both an overall favourable assessment of counting conditions and discussion of weather-related weak visible passage; the pilot retains the concrete observations without claiming complete detection.
- Black Kite: annual and paper totals differ for the record 2014 season. No exact total is reproduced. The after-sunset group in 2025 is explicitly retained as outside the counting protocol. Its proposed motivation (« oiseaux pressés ») is omitted as speculative.
- White Stork: 2017 is described as a record at the time, not as the current record. Detection improvements are part of the general synthesis, not an invented explanation for the 2025 record. Predictions of still larger future seasons are omitted. The historical paper associates improved coverage particularly with staffing from 2017, but the 2014 methods already describe salaried coverage; the display text does not claim salaried observation began in 2017.
- Wood Pigeon: several annual passages combine Stock Doves, Wood Pigeons and unidentified pigeons. The 2019 pilot distinguishes identified Wood Pigeons from all-pigeon totals. The 2025 general account includes the recent exceptional passage so the historical decline is not presented as an uninterrupted current trend. Hunting, warmer winters and changes in migration routes remain proposed explanations, not quantified or independently established causes. No optional anecdote section is filled.
- Eurasian Jay: the annual CSV has species accounts for 2010 and 2012, but no keyed species account for 2019. No 2019 annual account was manufactured from the later paper. « Irruption » replaces « invasion » as an editorial term for the same movement. Food availability and breeding success are retained as general biological context from the paper, not asserted as demonstrated causes of the individual 2010 or 2012 events.

## Next steps

Review the French pilot for voice, length, section names and fidelity. Translate the agreed pilot into English, then extend drafting to the remaining source accounts. Before general accounts are finalised, read all available annual entries for each species and all relevant paper entries; this pilot uses selected years only. Decide the website CSV field names when the text structure is agreed rather than changing the dataset release schema during the prose pilot.
