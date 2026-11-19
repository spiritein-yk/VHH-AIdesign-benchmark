# PDGFRB reference workflow

PDGFRB is the reference for the 12 backbone generator and MPNN combinations. Each generator used four starting frameworks and produced 5,000 backbones per framework. The historical workflow below follows the files from those runs into the score archives. [Input package](../configs/PDGFRB/README.md) · [Historical run environments](dependencies.md)

## Historical processing order

| Stage | Historical entry point and settings | Output used downstream |
| --- | --- | --- |
| BoltzGen | `boltzgen run` with `nanobody-anything`, 500 backbone sampling steps, and diffusion batch size 1. The native inverse-folding branch generated four sequences per backbone. | Backbone CIFs and companion design masks for the shared workflow. Inverse-folded outputs for the separate native branch. |
| PPIFlow | `sample_antibody_nanobody.py` with the nanobody checkpoint, hotspot and framework inputs. The inference configuration used 100 sampling steps and batch size 1. | Generated PDBs and per-structure input metadata. |
| RFAntibody | `rfdiffusion` with 50 diffusion steps and final step 1. Each framework was run in five 1,000-backbone shards. | Generated PDBs with native CDR labels. |
| IgGM | `launch_iggm.py` called the upstream design code with 10 steps, temperature 1.0, chunk size 64 and setup seeds 42–61 across the twenty 1,000-backbone shards. PDBFixer completed missing backbone atoms before sequence design. | Repaired PDBs and the corresponding FASTA design mask. |
| MPNN | The normalized structures used A for VHH and B for target. ProteinMPNN, SolubleMPNN, and AntiB‑MPNN each sampled four sequences per backbone at temperature 0.1 and batch size 4, using zero backbone noise and `--omit_AAs CX` at designable sites. The recorded RFAntibody and IgGM jobs passed seed arguments 42, 43, and 44 for the three models, respectively. AntiB‑MPNN recorded its effective seed in the FASTA headers. The inspected BoltzGen and PPIFlow jobs used each model’s default seed handling. | Generated sequences and a source-to-backbone mapping. |
| Side-chain packing and silent conversion | The sampled VHH sequences were backfilled onto the corresponding backbones. The LigandMPNN packer used three denoising steps and 16 internal samples. The production batch size was 64 for BoltzGen and PPIFlow and 32 for RFAntibody and IgGM. The packed coordinates were converted to binary silent files and indexed with source tags. | Original packed structures for Rosetta and AlphaFold 2 initial guess, plus the tag mapping. |
| Rosetta and AlphaFold 2 | RosettaScripts ran the archived scoring XML in score-only mode. The AlphaFold 2 initial-guess predictor took structures that passed the four Rosetta gates and ran Multimer v2 model 3 with three recycles. | Rosetta scores and AlphaFold 2 predictions and scores, joined by design source. |
| AlphaFold 3 | Eligible sources had AlphaFold 2 interaction PAE < 10 Å. Each generator × MPNN combination contributed up to 500 candidates ranked by PAE across frameworks. Identical sequence pairs shared one AlphaFold 3 task. The run used seed 1, five samples and a shared target MSA. | The native top model and its ipTM were joined with original Rosetta ΔΔG and VHH Cα RMSD from target-aligned packed coordinates. |
| Final joint filter | Original-design ΔΔG < −30 REU, target-aligned VHH Cα RMSD < 5 Å, and native top-model ipTM > 0.8. | 236 passing sources out of 3,039 sources evaluated with AlphaFold 3. [Archived scores](../data/README.md) |

IL‑7Rα uses the same stages with its own inputs and AlphaFold 3 admission rule. [IL‑7Rα workflow inputs](../configs/IL7RA/README.md)

## Released command examples

These examples show how the released scripts pass files between stages. To run them locally, prepare the required software, weights, and databases, set the path environment variables for the stages you use, and run the commands from the repository root.

```bash
FLOW=ppiflow
SCAFFOLD=3EAK
MODEL=proteinmpnn
WORK=outputs/PDGFRB
```

## Backbone and sequence design

```bash
python generate/run.py --target PDGFRB --flow "$FLOW" --scaffold "$SCAFFOLD" \
  --output-dir "$WORK/raw/$FLOW/$SCAFFOLD"
python mpnn/prepare_inputs.py --target PDGFRB --flow "$FLOW" --scaffold "$SCAFFOLD" \
  --config-dir configs/PDGFRB --input-dir "$WORK/raw/$FLOW/$SCAFFOLD" \
  --output-dir "$WORK/mpnn-input/$FLOW/$SCAFFOLD"
python mpnn/run.py --model "$MODEL" --inputs "$WORK/mpnn-input/$FLOW/$SCAFFOLD" \
  --output-dir "$WORK/mpnn/$FLOW/$SCAFFOLD/$MODEL"
```

The generation and MPNN entry points support `--dry-run` to print the external command using the current environment variables. The historical backbone budget is 5,000 per flow and framework. The historical RFAntibody and IgGM runs used five independent 1,000-backbone shards. Use `--count 1000`, give each shard its own output directory, and repeat `--input-dir` in shard order when preparing MPNN inputs. The released IgGM driver calls the upstream entry point directly.

| Flow | Directory passed to MPNN input preparation |
| --- | --- |
| BoltzGen | `intermediate_designs/` within the generator output, with paired CIF and NPZ files |
| PPIFlow | The generator output directory, including the PDBs and `input/*_input.csv` metadata |
| RFAntibody | Generator output or independently generated shard directories, with native CDR PDB labels |
| IgGM | Repaired output or repaired shard directories, with the original IgGM structure names |

For IgGM, run `python generate/repair.py --input-dir RAW_SHARD --output-dir REPAIRED_SHARD` with the PDBFixer interpreter configured for that stage. The repair checks the backbone atom set and preserves supplied coordinates. The MPNN adapter normalizes structures to chain A (VHH) and chain B (target), writes the model JSONL files, and records actual CDR positions and source provenance.

## Packing and silent conversion

```bash
GROUP="$FLOW/PDGFRB/$SCAFFOLD/$MODEL"
python packing/prepare.py --inputs "$WORK/mpnn-input/$FLOW/$SCAFFOLD" \
  --sequences "$WORK/mpnn/$FLOW/$SCAFFOLD/$MODEL" --model "$MODEL" --output-root "$WORK/backfilled"
python packing/run.py --target PDGFRB --flow "$FLOW" --scaffold "$SCAFFOLD" --model "$MODEL" \
  --input-root "$WORK/backfilled" --output-root "$WORK/packed"
python packing/archive.py --manifest "$WORK/backfilled/$GROUP/manifest.tsv" \
  --packed-dir "$WORK/packed/$GROUP" --output-dir "$WORK/silent/$GROUP"
```

The backfill manifest includes sampled sequences and their scaffold, flow and model. The archive writes `packed.silent` and `tag_map.tsv`, with a canonical ID of the form `flow__target__scaffold__model__backbone__sample`.

## Rosetta and AlphaFold 2

```bash
SCORE="$WORK/scores/$GROUP"
python score/run.py rosetta "$WORK/silent/$GROUP/packed.silent" "$SCORE/rosetta.sc"
python score/filter_rosetta.py --scores "$SCORE/rosetta.sc" \
  --output "$SCORE/rosetta-passed.tsv" --annotated "$SCORE/rosetta-selection.tsv"
python score/handoff.py af2-input --passed "$SCORE/rosetta-passed.tsv" \
  --silent "$WORK/silent/$GROUP/packed.silent" --output "$SCORE/af2.silent"
python score/run.py af2 "$SCORE/af2.silent" "$SCORE/af2"
python score/filter_af2.py --scores "$SCORE/af2/out.sc" \
  --output "$SCORE/af2-passed.tsv" --annotated "$SCORE/af2-selection.tsv"
```

Rosetta uses the historical PDGFRB scoring XML. AlphaFold 2 receives the selected original packed structures. The AlphaFold 2 selection table is useful for inspection. The AlphaFold 3 handoff below joins the raw AlphaFold 2 scores to the original Rosetta scores and structures.

## AlphaFold 3 admission, inference and final scoring

```bash
ALPHAFOLD3_DIR="$WORK/af3"
python score/handoff.py af3-scores --target PDGFRB \
  --group "$WORK/silent/$GROUP/tag_map.tsv" "$SCORE/rosetta.sc" "$SCORE/af2/out.sc" \
  --output "$ALPHAFOLD3_DIR/source-scores.tsv" &&
python score/select_af3_candidates.py --target PDGFRB --scores "$ALPHAFOLD3_DIR/source-scores.tsv" \
  --output "$ALPHAFOLD3_DIR/candidates.tsv" --annotated "$ALPHAFOLD3_DIR/admission.tsv" &&
python score/handoff.py af3-input --target PDGFRB --candidates "$ALPHAFOLD3_DIR/candidates.tsv" \
  --output-dir "$ALPHAFOLD3_DIR/tasks" &&
(
  for task in "$ALPHAFOLD3_DIR/tasks/jobs/"*.json; do
    python score/run.py af3 "$task" "$ALPHAFOLD3_DIR/tasks/predictions" || exit $?
  done
) && python score/score_af3_metrics.py --sources "$ALPHAFOLD3_DIR/tasks/sources.tsv" --output "$ALPHAFOLD3_DIR/all-candidates.tsv" &&
python score/filter_af3.py --scores "$ALPHAFOLD3_DIR/all-candidates.tsv" \
  --output "$ALPHAFOLD3_DIR/passed.tsv" --annotated "$ALPHAFOLD3_DIR/final-selection.tsv"
```

For full PDGFRB admission, supply one `--group` for each completed route group before choosing candidates. The script applies the AlphaFold 2 PAE cut and groups by generator and MPNN across frameworks. Use `--target-msa` as documented by `af3-input --help` to reuse the historical target MSA when available. The command in the example lets AlphaFold 3 run its data pipeline. The AlphaFold 3 builder retains every admitted source ID and reuses predictions for identical sequences. Final RMSD is computed against original packed coordinates after target alignment.

## BoltzGen native branch

Use `boltzgen_native/run.py analyse` with one generator output for each scaffold and a separate output directory. Use `boltzgen_native/run.py select` with four repeated `--source` flags pointing to the four analysed directories. Those results use the native workflow's metrics and are reported separately.
