# Archived score tables

The two archives contain 960,000 design-source records each for the shared structure generation × MPNN benchmark. They contain historical Rosetta and AlphaFold 2 scores and the full set of source-level AlphaFold 3 evaluations, including scores below the final thresholds. [IL7RA archive](IL7RA_scores.tsv.gz) · [PDGFRB archive](PDGFRB_scores.tsv.gz)

| Target | Design sources | Sources with AlphaFold 2 scores | Sources evaluated with AlphaFold 3 | Unique AlphaFold 3 tasks | Sources passing the final joint filter |
| --- | ---: | ---: | ---: | ---: | ---: |
| IL7RA | 960,000 | 326,908 | 1,572 | 1,551 | 105 |
| PDGFRB | 960,000 | 303,729 | 3,039 | 3,035 | 236 |

AlphaFold 3 scores come from the historical production records. IL7RA combines the initial 902 source records with the 670-record extension. PDGFRB uses the full production score table. These records match the manuscript's AlphaFold 3 figure-source IDs, task IDs and native top-model ipTM values. Counts by flow and MPNN match Figure 4b.

A blank cell indicates that the corresponding score record is absent. Some historical scores contain nonfinite values. `NaN` and `nan` indicate an undefined value, and `inf` indicates positive infinity. AlphaFold 3 scores are linked by design source, and a single prediction task can correspond to multiple sources.

## Design identifiers

| Column | Meaning |
| --- | --- |
| `design_id` | Archive ID of the design source, with one source per row. |
| `target` | Target code. `IL7RA` denotes IL‑7Rα and `PDGFRB` denotes PDGFRB. |
| `flow` | Backbone generator, one of `boltzgen`, `ppiflow`, `rfantibody`, or `iggm`. |
| `template` | Identifier of the starting VHH framework. |
| `mpnn` | Sequence design model, one of `proteinmpnn`, `solublempnn`, or `antibmpnn`. |
