# VHH Design Benchmark

This repository documents the code, configurations, historical run environments, and scores behind the PDGFRB and IL‑7Rα (IL7RA) VHH benchmark. PDGFRB is the reference workflow. It covers four backbone generators and three sequence design methods across four starting VHH frameworks. The BoltzGen native branch has its own analysis and selection.

1. [PDGFRB workflow](docs/workflow.md) — what ran at each stage and how its outputs were used downstream.
2. [Historical run environments](docs/dependencies.md) — the recorded software versions, commits, checkpoints, and database identifiers.
3. [Archived scores](data/README.md) — the two targets' evaluation records, field definitions, and file checksums.

The [IL7RA input package](configs/IL7RA/README.md) supplies the target, framework inputs, hotspot and CDR definitions. It uses the same stage entry points with its own AlphaFold 3 admission rule.

| Directory | Contents |
| --- | --- |
| [`configs/`](configs/) | Target and framework structures, hotspot and CDR definitions, and model-specific inputs |
| [`generate/`](generate/), [`mpnn/`](mpnn/) | Backbone generation, IgGM backbone repair, and sequence design |
| [`packing/`](packing/), [`score/`](score/) | Sequence backfill, side-chain packing, silent conversion, scoring, and selection |
| [`boltzgen_native/`](boltzgen_native/) | Native BoltzGen refolding, analysis, merging, and selection |
| [`data/`](data/) | Historical score archives and checksums |

Project licensing and review of third-party redistribution terms remain open for public release.

### Restore score files

The score archives are split into 40 MiB parts for GitHub. After cloning, run:

```bash
cd data
for target in IL7RA PDGFRB; do
  cat "${target}_scores.tsv.gz".part* > "${target}_scores.tsv.gz"
done
sha256sum -c scores.sha256
```

The reconstructed gzip files have the original names expected by the project.
