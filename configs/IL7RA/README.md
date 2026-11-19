# IL‑7Rα inputs

The released target is `common/IL7RA.pdb`. The starting VHH frameworks are 3EAK, 5jds, CAMV_U0058 and V7R-2E9. See [validation.json](validation.json) for hotspot positions and residue identities, the 1-based CDR ranges for each framework, allowed variable loop lengths, and fixed IgGM loop lengths. The method subdirectories provide the model-specific input formats.

Use the [shared workflow](../../docs/workflow.md) with `--target IL7RA`, `--config-dir configs/IL7RA` and an IL7RA framework. The AlphaFold 3 admission script takes the same target flag. Historical admission used AlphaFold 2 interaction PAE < 8 Å, plus designs with PAE < 10 Å from IgGM, AntiB‑MPNN, and PPIFlow × ProteinMPNN. The combined set contained 1,572 design sources and 1,551 unique AlphaFold 3 tasks. [Archived scores](../../data/README.md) Keep the model-specific chain IDs and residue numbering supplied in each input.
