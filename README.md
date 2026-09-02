# bea2 — BioEnvAda v2

Nextflow pipeline for per-ortholog-group (OG) biophysical **and** evolutionary characterisation, version 2 of [BioEnvAda](https://github.com/Bio2Byte/bioenvada). Designed to stay simple: one process = one tool = one container, all logic in standalone `bin/*.py` scripts, per-OG fan-out with a single global statistics fan-in.

## What it does (per OG)

```
aa fasta ──mafft──> AA MSA ──cdskit backalign──> codon MSA
                       │
                 occupancy pruning (--occupancy 0.5)   <── THE canonical MSA:
                       │        pruned position == csubst codon_site == HyPhy site
        ┌──────────────┼───────────────────────┐
   biophysics      IQ-TREE (outgroup-      custom predictions table
   b2bTools        constrained) ──> rooted      (split per OG)
   AIUPred         tree ──> csubst (convergence/divergence)
   IPC                  └──> HyPhy  (FEL, FUBAR, MEME | Contrast-FEL, RELAX)
   DSSP
        └──────────────┬───────────────────────┘
              mapped onto the pruned MSA
                       │
        <og>_combined.tsv  — one row per residue per sequence per MSA position
              + per-OG plots + native tool outputs
```

Then one global stage collects all OGs: distributions of every predictor (residue-level and per-OG means, split by category), RELAX K distribution, csubst branch-pair summaries → `results/stats/`.

## Quickstart

```bash
conda activate nfconda    # or wherever nextflow lives

nextflow run pipeline.nf \
    -profile local \
    --aa_dir data/aa --nuc_dir data/nucs \
    --categories synechococcus_categories.tsv \
    --category_column temp_cat2 \
    --og_pattern 'CK_\d+' \
    --predictions_csv predictions_all.csv \
    -resume
```

Profiles: `local` (docker, laptop testing) · `hydra` (VUB HPC: Slurm + MAFFT/IQ-TREE **modules** + singularity — set the exact module strings in the `hydra` profile of `nextflow.config`) · `hpc` (generic Slurm + singularity for everything).

## Inputs

| Param | Meaning |
|---|---|
| `--aa_dir` / `--nuc_dir` | folders of per-OG amino-acid / CDS fastas, paired by filename (`<og>.fasta`) |
| `--categories` | TSV: `strain_name` first, then category columns (e.g. `temp_cat`, `temp_cat2`) |
| `--category_column` | which column defines the grouping — header name or 1-based index; swap to try different clusterings |
| `--fg_level` / `--outgroup_level` | values marking foreground (csubst/HyPhy) and rooting outgroup |
| `--og_pattern` | regex recognising an OG id **inside identifiers** (`'CK_\d+'`, `'OG\d+'`, …); used to split the custom predictions table |
| `--predictions_csv` | optional big table `block_id,residue_index,aa,pred_*` (`block_id` = `<og>_<strain>_<locus>`) |
| `--occupancy` | min MSA column occupancy (default 0.5) |

Sequence IDs must be `<og>_<strain>_<locus>` (strains may contain `_`). Strains absent from the categories table stay in trees/alignments as uncategorised background.

## Outputs

```
results/
  ogs/<OG>/
    alignment/            native mafft out, codon MSA, pruned MSAs, msa_columns.tsv
    biophysics/<tool>/    native/ (raw tool output) + <og>_<tool>_mapped.tsv
    evolution/
      iqtree/native/      full IQ-TREE output
      csubst/native/      search/ inspect/ csubst_iqtree/ (FULL IQ-TREE intermediates,
                          .state ancestral reconstructions — csubst has/keeps everything)
      csubst/             <og>.branch_pairs.tsv .sites.tsv .hotspots.tsv (convergent+divergent)
      hyphy/native/       all JSONs + logs        hyphy/  per-method persite + by_site tables
      <og>.rooted.nwk
    <og>_combined.tsv     ← THE table: row = residue per sequence per pruned-MSA position
    plots/                MSA heatmap + per-feature category line plots with hotspot overlay
  stats/                  global distributions & summaries over all OGs
  foreground/             foreground.tsv, strain lists, cleaned categories
  pipeline_info/          report / timeline / trace
```

## Containers

| Image | Used for | Recipe | Platforms |
|---|---|---|---|
| `slheidig/csubst:01` | cdskit, pruning, IQ-TREE, ete4 rooting, csubst | `docker/csubst/` | amd64 + arm64 |
| `slheidig/og_b2b_pca:latest` | b2bTools + all pandas/matplotlib steps | `docker/og_b2b_pca/` | amd64 + arm64 |
| `slheidig/hyphy:2.5.101` | HyPhy (FEL/FUBAR/MEME/contrast-fel/RELAX) | `docker/hyphy/` | amd64 + arm64 |
| `ghcr.io/doszilab/aiupred:cpu` | AIUPred disorder/binding | authors' official image, nothing to build | amd64 only |
| `slheidig/ipc:01` | IPC pI + per-residue charge | `docker/ipc/` | amd64 + arm64 |
| `slheidig/dssp3:2` | mkdssp secondary structure + SASA | shared with SIMSApiper, nothing to build | amd64 |
| `slheidig/simsapiper:06` | secstructartist secondary-structure plots | shared with SIMSApiper, nothing to build | amd64 |
| `quay.io/biocontainers/mafft` | mafft (`hpc`/`local`; module on hydra) | public | amd64 |

### Building them (multi-arch, from an Apple Silicon Mac)

```bash
docker/build.sh                  # build + push the images that need it
docker/build.sh hyphy            # a single target
docker/build.sh --load hyphy      # single-arch local image, for testing only
docker/build.sh --dry-run all     # show the plan
```

**Never use plain `docker build` for these.** On an arm Mac it produces an
arm64-only image that looks fine locally and cannot run on the clusters, and
Docker Desktop's default builder cannot write multi-arch manifests at all
(`docker exporter does not currently support exporting manifest lists`).
`docker/build.sh` creates a `docker-container` builder that can, takes the
per-image platform list from `docker/docker-bake.hcl`, pushes without keeping a
local copy (laptop disk), verifies the published manifest really contains
`linux/amd64`, and stops the builder again afterwards.

Disk: the builder's cache is capped at ~6 GB by `docker/buildkitd.toml`.
`docker buildx du` shows it, `docker buildx rm bea2` reclaims it.

Once pushed, the clusters pull from Docker Hub and singularity/apptainer
converts automatically, picking the amd64 entry from the manifest.

**AIUPred**: the `AIUPRED` process calls the `aiupred` CLI in the authors' image directly, with `--force-cpu` (the nodes have no GPU) and `-b`, which adds a Binding column next to Disorder in the same pass — so `aiupred_disorder` and `aiupred_binding` both come out of one run. Set `--aiupred_binding false` for disorder only. The CLI writes a `#` banner, then `#>id` per sequence followed by `position residue score...` rows; an awk step in the process pulls each id down onto its rows to produce the flat predictor table. Note the sequence headers are `#>id`, not `>id`. For the official IPC CLI instead of the built-in pKa computation, see `docker/ipc/Dockerfile`.

**DSSP** (`--dssp --structure_dir <dir>`, off by default): the only predictor that reads structures instead of sequences. It expects one predicted model per sequence at `<structure_dir>/<og>/<sequence_id>.pdb` — the layout ESMFold writes, single chain, residues numbered 1..N, pLDDT in the B-factor column. `DSSP_RUN` runs `mkdssp` over the OG's models and pulls the CA B-factors into a pLDDT table; `DSSP_PARSE` turns both into the flat predictor table. Because every sequence here has a model whose sequence matches the fasta exactly, the DSSP residue number is used directly as `residue_index` — none of SIMSApiper's sequence-reconciliation logic is needed. Columns: `dssp_ss8` (8-state — DSSP's H B E G I T S, with unassigned written as `X` as in SIMSApiper, so it reads as neither an alignment gap nor a loop), `dssp_ss3` (H/E/L — helix `HGI`, sheet `EB`, loop `TS` + unassigned; the same three buckets SIMSApiper uses), `dssp_acc` (SASA in Å²), `dssp_rsa` (ACC / max-ASA, Tien et al. 2013), `dssp_kd` (Kyte–Doolittle), `dssp_surface_hydrophobicity` (`rsa × kd` — hydrophobic *and* exposed), `dssp_plddt`. The two SS columns are strings, so `PLOT_OG` and `GLOBAL_STATS` skip them — `PLOT_DSSP` below plots them instead.

**Secondary-structure plots** (`PLOT_DSSP`, container `slheidig/simsapiper:06`, which carries `secstructartist`): `bin/plot_dssp_ss.py` adapts SIMSApiper's `2Dstructure_plot.py` / `DSSPcodesMSA_plot.py` / `dssp_seqview_plot.py` to the pruned-MSA coordinates and splits everything **by category**, which is what bea2 compares. Per OG:

- `<og>_ss_consensus.pdf` — a secstructartist cartoon of the consensus secondary structure, one ribbon track per category, above H/E/L frequency panels with the csubst hotspot sites overlaid.
- `<og>_ss_alignment.pdf` — per-sequence H/E/L, sequences grouped by category.
- `<og>_ss_consensus.tsv` — the numbers behind the cartoon. Published to `biophysics/dssp/` alongside the other DSSP tables, not to `plots/`.

A position gets a consensus H or E when at least `--ss_consensus` (default 0.5) of the *residues present* carry it, otherwise loop. Note the denominator: gaps are excluded rather than counted as loop, so the three frequencies sum to 1 — SIMSApiper instead divides by all sequences, which lets a gappy column drag every frequency down. A gap and a loop are different things here and are drawn differently (lightgrey vs orange in the per-sequence view).

## Extending with a new predictor

1. Write `bin/run_mytool.py` emitting the standard long format: `sequence_id  residue_index  residue  <features>` (1-based ungapped index).
2. Add a process in `modules/predictors.nf` emitting `tuple(og, 'mytool', tsv)`.
3. `preds = preds.mix(MYTOOL.out.pred)` in `pipeline.nf`. Mapping, combined table, plots and global stats pick it up automatically.

## Sizing (1500 OGs × 54 seqs, 40 cores)

Throughput-bound, not memory-bound: ~4 GB RAM/core is comfortable (peak concurrent use 60–100 GB). Aggregate ≈ 1–1.5 CPU-h per OG, dominated by MEME/RELAX/csubst → **≈ 2–3 days wall-time on 40 cores**. For a faster first pass: `--hyphy_methods FEL --hyphy_fg_methods relax`. Failed tasks retry twice with doubled resources, then are ignored (check `results/pipeline_info/trace_*.txt` for skipped OGs).

## Tips

- Rerun after a failure/parameter change: append `-resume`.
- Different clustering: `--category_column temp_cat --outdir results_strict` (with `-resume`, alignment/trees/predictors are reused; only foreground-dependent steps rerun).
- csubst arity is fixed at 2 (pairs); raise inside `modules/csubst.nf` if needed.
