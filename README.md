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

Profiles (all in `nextflow.config`; there is no `conf/`):

| Profile | What it does |
|---|---|
| `local` | docker, nothing submitted to a scheduler — laptop testing |
| `hydra` | VUB HPC: apptainer, Slurm for the nine heavy processes, MAFFT **module** for the alignment, and no container at all for the stdlib/module-covered steps. Set the module strings in the `hydra` profile |
| `hpc` | any cluster: Slurm + singularity, **no host software required** — every process declares a container and this profile overrides none of them |

By default only nine processes are submitted to Slurm; the eight seconds-long
ones run on the **local executor**, inside the driver's own allocation, because a
scheduler round-trip costs more than the work. So give the driver job real
cpus/memory — see [Running at scale](#running-at-scale).

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
                          tools: b2b, aiupred, ipc, dssp, custom
    evolution/
      iqtree/native/      .treefile .contree .iqtree .log
      csubst/native/      search/ (csubst_b, csubst_cb_2, csubst_cb_stats, log,
                          sites/<pair>/csubst.tsv)
      csubst/             <og>.branch_pairs.tsv .sites.tsv .hotspots.tsv (convergent+divergent)
      hyphy/native/       all JSONs + logs        hyphy/  per-method persite + by_site tables
      <og>.rooted.nwk
    <og>_combined.tsv     ← THE table: row = residue per sequence per pruned-MSA position
    plots/                only for OGs selected with --plot (see below)
  stats/                  global distributions & summaries over all OGs
  foreground/             foreground.tsv, strain lists, cleaned categories
```

Output volume is tuned for 1500 families — measured against the earlier layout,
**365 → 70 files per OG**. What changed:

- `--plot all | none | <og>,<og>` selects which families get figures; 27 PDFs per
  OG is ~40,000 files at scale. Plots live in their own processes, so changing
  `--plot` re-runs only the plot tasks and leaves every cached table intact under
  `-resume`. `./replot_og.sh <results_dir> <og>…` regenerates figures from
  published output if `-resume` will not cooperate.
- Only `csubst.tsv` is kept from each `search/sites/<pair>/` directory — the other
  seven files per pair are never read by the pipeline. `--publish_csubst_native`
  adds `inspect/` and csubst's own `csubst_iqtree/` back. Per-pair figures are off
  by default (`--csubst_site_plots`); `replot_og.sh` regenerates those too.
- DSSP publishes one concatenated `<og>.dssp.txt.gz` instead of a `.dssp` file per
  sequence, and IQ-TREE publishes four files instead of its whole prefix set.
- There is no `pipeline_info/`: `-with-trace`/`-with-report`/`-with-timeline`
  cannot be used here (see [Monitoring](#monitoring)).

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

**DSSP** (`--dssp --structure_dir <dir>`, off by default): the only predictor that reads structures instead of sequences. It expects one predicted model per sequence at `<structure_dir>/<og>/<sequence_id>.pdb` — the layout ESMFold writes, single chain, residues numbered 1..N, pLDDT in the B-factor column. The `DSSP` process runs `mkdssp` over the OG's models, pulls the CA B-factors into a pLDDT table, and turns both into the flat predictor table — all in one task, since the dssp3 image carries python3 as well as mkdssp. The per-sequence `.dssp` files therefore stay inside the task and only one concatenated `<og>.dssp.txt.gz` is published. Because every sequence here has a model whose sequence matches the fasta exactly, the DSSP residue number is used directly as `residue_index` — none of SIMSApiper's sequence-reconciliation logic is needed. Columns: `dssp_ss8` (8-state — DSSP's H B E G I T S, with unassigned written as `X` as in SIMSApiper, so it reads as neither an alignment gap nor a loop), `dssp_ss3` (H/E/L — helix `HGI`, sheet `EB`, loop `TS` + unassigned; the same three buckets SIMSApiper uses), `dssp_acc` (SASA in Å²), `dssp_rsa` (ACC / max-ASA, Tien et al. 2013), `dssp_kd` (Kyte–Doolittle), `dssp_surface_hydrophobicity` (`rsa × kd` — hydrophobic *and* exposed), `dssp_plddt`. The two SS columns are strings, so `PLOT_OG` and `GLOBAL_STATS` skip them — `PLOT_DSSP` below plots them instead.

**Secondary-structure plots** (`PLOT_DSSP`, container `slheidig/simsapiper:06`, which carries `secstructartist`): `bin/plot_dssp_ss.py` adapts SIMSApiper's `2Dstructure_plot.py` / `DSSPcodesMSA_plot.py` / `dssp_seqview_plot.py` to the pruned-MSA coordinates and splits everything **by category**, which is what bea2 compares. Per OG:

- `<og>_ss_consensus.pdf` — a secstructartist cartoon of the consensus secondary structure, one ribbon track per category, above H/E/L frequency panels with the csubst hotspot sites overlaid.
- `<og>_ss_alignment.pdf` — per-sequence H/E/L, sequences grouped by category.
- `<og>_ss_consensus.tsv` — the numbers behind the cartoon. Published to `biophysics/dssp/` alongside the other DSSP tables, not to `plots/`.

A position gets a consensus H or E when at least `--ss_consensus` (default 0.5) of the *residues present* carry it, otherwise loop. Note the denominator: gaps are excluded rather than counted as loop, so the three frequencies sum to 1 — SIMSApiper instead divides by all sequences, which lets a gappy column drag every frequency down. A gap and a loop are different things here and are drawn differently (lightgrey vs orange in the per-sequence view).

## Extending with a new predictor

1. Write `bin/run_mytool.py` emitting the standard long format: `sequence_id  residue_index  residue  <features>` (1-based ungapped index).
2. Add a process in `modules/predictors.nf` emitting `tuple(og, 'mytool', tsv)`.
3. `preds = preds.mix(MYTOOL.out.pred)` in `pipeline.nf`. Mapping, combined table, plots and global stats pick it up automatically.
4. Add a `withName: 'MYTOOL'` block in `nextflow.config`. **One block per process, complete** — container, and if it needs a Slurm job also executor, array, cpus, memory, time. Do not add a `label`: a `withLabel:`/`withName:` block does not merge with a same-named block elsewhere, it replaces it wholesale, which is how every resource setting once got silently dropped.
5. If the predictor is toggleable, add it to the `n_pred` count in `pipeline.nf` so `MAP_TO_MSA`'s `groupTuple(size:)` still knows how many predictors to expect.

## Running at scale

Measured on one complete 8-OG run (nothing cached), against the earlier pipeline:

| | 8 OGs | projected 1500 OGs |
|---|---|---|
| tasks | 208 → 127 | 39,000 → ~22,500 |
| **sbatch calls** | 208 → **10** | 39,000 → **~300** |
| files per OG | 365 → 70 | — |
| inodes (whole run) | 3,616 → 974 | — |

Three things compound: 26 → 15 tasks per OG, five of those fifteen no longer
touch Slurm, and the rest are batched into job arrays. The old
`submitRateLimit = '10 sec'` would have spent **108 hours purely submitting**
39,000 jobs.

Per-task resources are set per process in `nextflow.config`, each block complete
with executor, array size, cpus, memory, time and container. Peak RSS from
`sacct`: TREE 238 MB, CSUBST 472 MB, HyPhy ≤278 MB — the memory consumers are the
**predictors** (AIUPred 1.9 GB, the predictions splitter 2.9 GB, b2bTools 1.0 GB),
not the phylogenetics. Times are 2× the measured maximum; since `time` is
`X * task.attempt`, a retry doubles it again.

Knobs: `--array_size 0` disables array batching · `--heavy_executor local` runs
everything in one allocation · `--hyphy_methods FEL --hyphy_fg_methods relax` for
a faster first pass.

**The driver job** must outlive the whole run and holds task metadata for
~22,500 tasks. It also runs the local-executor tasks, so its cores cap how many
of those run concurrently. For a proteome run: days of walltime, ~8 cores, 16 GB,
`NXF_OPTS='-Xms1g -Xmx8g'`, and `NXF_ANSI_LOG=false` (otherwise the slurm `.out`
fills with ANSI redraw frames — that is what made them 137 KB). Expect Slurm's job
report to flag "suspiciously low CPU time": a workflow driver is wait-bound by
design, and the real compute is in the child jobs.

Failed tasks retry once with doubled resources, then are **ignored** — individual
families are expected to fail and that must never end a 1500-family run.

## Monitoring

`-with-trace`, `-with-report` and `-with-timeline` **cannot be used**: Nextflow
injects a `ps`-based metric collector *into each container*, and five of the seven
images ship no `procps`, so every containerised task exits 1 before running
anything. Two helpers replace them — both read files only, no JVM, no container,
safe on a login node and while the pipeline is running:

```bash
./progress.sh          # per-process [n of m] table that NXF_ANSI_LOG=false hides
./failed_tasks.sh      # per-process tally + failures + first failure's stderr
NXF_LOG=.nextflow.log.1 ./progress.sh    # an earlier run
```

`progress.sh` reads `.nextflow.log` plus `work/*/*/.exitcode`; `failed_tasks.sh`
queries `nextflow log` (which starts a JVM, so it is the heavier of the two).

## Tips

- Rerun after a failure/parameter change: append `-resume`.
- Different clustering: `--category_column temp_cat --outdir results_strict` (with `-resume`, alignment/trees/predictors are reused; only foreground-dependent steps rerun).
- csubst arity is fixed at 2 (pairs). Raising it would *increase* output rather than reduce it, and would break `select_significant_pairs.py` and `aggregate_csubst.py`, which assume `csubst_cb_2.tsv` and two-id pair names.
- Do not add `-asr`/`--rate` to the IQ-TREE call to save csubst its own run: csubst already uses `-te` on the rooted tree (fixed topology, no search) and relabels it before reconstructing states, so an external `.state` will not line up. Details in `modules/tree.nf`.
