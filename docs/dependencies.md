# Historical run environments

The table records the environments used for the PDGFRB benchmark. Commits and version strings come from successful-run records or their launchers. CUDA toolkit versions loaded by launchers and CUDA build versions reported by PyTorch are identified separately.

Each software has one historical YAML record. For software used in multiple runtimes, the YAML reports one observed runtime and its applicable backbone groups. ProteinMPNN and SolubleMPNN share a record with separate checkpoint entries.

| Stage | Recorded software and inputs | Historical YAML |
| --- | --- | --- |
| BoltzGen | Commit `a3149cf18eeb58648d1abbb27539bd73f746cdda` for backbone generation and the native branch. The runs used the BoltzGen CLI and its configured model cache. [BoltzGen](https://github.com/HannesStark/boltzgen) | [boltzgen.yaml](../env/historical/boltzgen.yaml) |
| PPIFlow | Commit `920975e42fcaae80f6f10f7b9c21a4d57058d673`, the nanobody sampling entry point, and `nanobody.ckpt`. The launcher loaded CUDA Toolkit 13.2. | [ppiflow.yaml](../env/historical/ppiflow.yaml) |
| RFAntibody | Commit `8fe311415754e0276d1a39c87c57e69c88927a2d`. The production launcher used `rfdiffusion` from the RFAntibody environment. [RFAntibody](https://github.com/RosettaCommons/RFantibody) | [rfantibody.yaml](../env/historical/rfantibody.yaml) |
| IgGM and backbone repair | IgGM commit `06abc563b3fc8c7ea020543add16b69b6f8a1c8d` with the study's seed launcher. Backbone atom completion used PDBFixer and OpenMM on CPU with seed 42. [PDBFixer](https://github.com/openmm/pdbfixer) | [iggm.yaml](../env/historical/iggm.yaml) · [pdbfixer.yaml](../env/historical/pdbfixer.yaml) |
| ProteinMPNN and SolubleMPNN | ProteinMPNN commit `8907e6671bfbfc92303b5f79c4b5e6ce47cdef57`. The standard and soluble weight sets each contained a distinct `v_48_020.pt`. [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) | [proteinmpnn.yaml](../env/historical/proteinmpnn.yaml) |
| AntiB‑MPNN | Commit `7e6329ba3c317e7344c6b09add72dea2390af39a`, `Running_AntiBMPNN_run.py`, and `antibmpnn_000.pt`. The working directory was the tool's `example/` directory. | [antibmpnn.yaml](../env/historical/antibmpnn.yaml) |
| MPNN runtimes | For BoltzGen and PPIFlow backbones, the run logs report Python 3.11.16 and PyTorch 2.11.0+cu128. The inspected RFAntibody and IgGM run logs report Python 3.12.8 and PyTorch 2.6.0+cu124. | — |
| Side-chain packing and silent conversion | LigandMPNN-compatible packer with `ligandmpnn_sc_v_32_002_16.pt`. For BoltzGen and PPIFlow backbones, packing used PyTorch 2.8.0+cu128. For RFAntibody and IgGM backbones, it used PyTorch 2.2.1+cu121. Silent conversion used silent_tools with Rosetta or PyRosetta. [LigandMPNN](https://github.com/dauparas/LigandMPNN) | [ligandmpnn.yaml](../env/historical/ligandmpnn.yaml) · [silent_tools.yaml](../env/historical/silent_tools.yaml) |
| Rosetta | RosettaScripts from release-371, OpenMPI 3.1.6, DAlphaBall, and the standalone PSIPRED runner. The historical scoring protocol is in [score.xml](../score/rosetta/score.xml). [Rosetta downloads](https://downloads.rosettacommons.org/) | [rosetta.yaml](../env/historical/rosetta.yaml) |
| AlphaFold 2 initial guess | The historical initial-guess fork with `model_3_multimer_v2` weights, PyRosetta, and silent_tools. The batch launcher loaded CUDA 11.8 when available. | [af2.yaml](../env/historical/af2.yaml) |
| AlphaFold 3 | Version 3.0.0, the production model parameters and genetic databases, and the PDB database archive `pdb_2022_09_28_mmcif_files.tar`. [AlphaFold 3 installation](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md) | [af3.yaml](../env/historical/af3.yaml) |
