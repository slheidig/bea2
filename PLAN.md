# BioEnvAda v2 — Nextflow pipeline plan

**Goal:** rework the csubst_trial analysis (csubst_pipeline.sh + selection_scan.sh + b2bplots/) into a simple, extensible Nextflow DSL2 pipeline for ~1500 ortholog groups (OGs) × 54 sequences, runnable on Hydra (VUB HPC, Slurm + modules + Singularity) and on any generic HPC (Slurm + Singularity only). Style guide: [Bio2Byte/bioenvada](https://github.com/Bio2Byte/bioenvada) v1 (`pipeline.nf` + `modules/*.nf` + `bin/*.py` + profile-based `nextflow.config`), with the resource-label / publishDir conventions of `proteinfam_pca_nf`.

---

## 1. Design principles

- **One process = one tool = one container.** No process assumes anything about the host beyond its container (or module on Hydra). All Python logic lives in `bin/` as standalone argparse scripts — nothing embedded in heredocs.
- **One coordinate system.** The occupancy-pruned MSA (default ≥ 50% occupancy, `--occupancy 0.5`) is *the* canonical alignment. Pruned-MSA position = csubst `codon_site` = HyPhy site = row key of every tabular output. `cdskit hammer` is dropped; a column map (original ↔ pruned position) is kept for traceability.
- **Per-OG fan-out, global fan-in.** Everything per-OG runs in parallel channels keyed by OG id; a final stage collects all OGs for global distributions.
- **Every tool step emits two things:** a `native/` folder (everything the tool wrote, untouched) and a tidy TSV keyed by `sequence_id` + `msa_position` (or per sequence where the value is sequence-level).

## 2. Inputs (params)

| Param | Default | Meaning |
|---|---|---|
| `--aa_dir` | required | unaligned amino-acid fastas, one per OG (`<og>.fasta`) |
| `--nuc_dir` | required | matching unaligned CDS fastas (`<og>.fasta`), same IDs |
| `--categories` | required | `synechococcus_categories.tsv` (strain_name + category columns) |
| `--category_column` | `temp_cat2` | column of interest (header name or 1-based index) — swap to try different clusterings |
| `--fg_level` | `warm` | value in that column marking foreground lineages |
| `--outgroup_level` | `Outgroup` | value marking rooting outgroup (never foreground) |
| `--predictions_csv` | optional | the big custom-tool table (1.9 GB; `block_id,residue_index,aa,pred_*`), split per OG once |
| `--occupancy` | `0.5` | min column occupancy kept in the pruned MSA |
| `--genetic_code` | `11` | bacterial |
| `--hyphy_methods` | `FEL,FUBAR,MEME` | whole-tree per-site methods |
| `--hyphy_fg_methods` | `contrast-fel,relax` | foreground methods |
| `--csubst` / `--hyphy` / `--aiupred` / `--b2btools` / `--ipc` | `true` | per-tool toggles (bioenvada style: easy to extend with new predictors) |
| `--dssp` / `--structure_dir` | `false` / `null` | DSSP on predicted models at `<structure_dir>/<og>/<sequence_id>.pdb` |
| `--ss_consensus` | `0.5` | min frequency of the residues present for a consensus H/E in the DSSP plots |
| `--ocn_cutoff`, `--omega_cutoff`, `--pvalue`, `--posterior` | as in current scripts | significance thresholds |
| `--outdir` | `results` | |

Strain parsing from sequence IDs (`<og>_<strain>_<locus>`, strains may contain `_`) reuses the robust rsplit logic of `01_find_seqsoi.py`. Strains absent from the categories table (e.g. WH8020) stay in tree/alignment as uncategorised background and are flagged in a warnings file, never silently dropped.

## 3. Workflow (modules and processes)

```
pipeline.nf
modules/
  align.nf        ALIGN_AA (mafft) · PRUNE (backalign + occupancy pruning)
  predictors.nf   B2BTOOLS · AIUPRED · IPC · DSSP · SPLIT_CUSTOM_PREDICTIONS
  tree.nf         MAKE_FOREGROUND · TREE (constraint + IQ-TREE + rooting)
  csubst.nf       CSUBST
  hyphy.nf        HYPHY_WHOLETREE · HYPHY_FOREGROUND
  collect.nf      MAP_TO_MSA · COLLECT (by-site + combined table) · PLOT_OG · PLOT_DSSP · GLOBAL_STATS
```

17 processes, 15 tasks per OG. Several steps were fused for proteome scale: the
seconds-long python steps are not worth a scheduler round-trip, and each task
costs ~30 work-dir inodes. See §7.

Per OG (channels joined on OG id):

1. **ALIGN_AA** — `mafft --auto` on the AA fasta → `<og>.aln.fa`. *(Hydra: `module MAFFT`; hpc: mafft container.)*
2. **PRUNE** (`bin/map_aa_to_nuc.py` + `bin/prune_msa.py`, fused) — back-translates the full AA alignment onto the CDS → in-frame codon alignment (before pruning, so codons stay consistent), then drops AA columns with occupancy < `--occupancy` and applies the *same* columns ×3 to the codon alignment. Outputs: pruned AA MSA, pruned codon MSA, `msa_columns.tsv` (pruned_pos ↔ original_pos), occupancy per column. Both steps are stdlib-only python, so on hydra this runs with no container at all.
3. **Predictors** (all on *ungapped* AA sequences except DSSP, mapped to the MSA later):
   - **B2BTOOLS** (container `slheidig/og_b2b_pca:latest`, all predictors) → long per-residue TSV (backbone, sidechain, helix, sheet, coil, ppII, earlyFolding, disoMine, agmata).
   - **AIUPRED** (new container, see §5) → per-residue disorder score (default mode; binding mode available via `ext.args`).
   - **IPC** (new container) → per-sequence pI (+ per-residue charge at pH 7, which is positional); pI repeated per residue row in the tables.
   - **DSSP** (container `slheidig/dssp3:2`, which also carries python3) — the one predictor that consumes **structures** rather than sequences: one predicted model per sequence at `<structure_dir>/<og>/<sequence_id>.pdb`. `mkdssp` gives 8-state secondary structure and absolute SASA; `bin/parse_dssp.py` adds RSA (Tien et al. 2013 max-ASA), Kyte–Doolittle hydropathy, surface hydrophobicity (`rsa × kd`) and the pLDDT read from the model's CA B-factors. Models cover exactly the fasta sequence, so the DSSP residue number is the ungapped `residue_index` directly. mkdssp and the parser are fused, which keeps the per-sequence `.dssp` files inside the task: only one concatenated `<og>.dssp.txt.gz` is published, instead of one file per sequence (57 → 1 on the test set).
   - **SPLIT_CUSTOM_PREDICTIONS** — runs **once**: streams the 1.9 GB CSV in chunks, splits rows by the `CK_########` prefix of `block_id` → one `<og>_custom.tsv` per OG, joined into the per-OG channel. (Note: `residue_index` there is 1-based ungapped.)
4. **TREE** (`bin/make_constraint.py` + IQ-TREE + `bin/root_tree.py`, fused; container `slheidig/csubst:01`, which carries iqtree3, python3 and ete4) — builds the outgroup-monophyly constraint from `--outgroup_level` (partial outgroup sets OK; a single taxon roots with `-o` alone; none → unconstrained + midpoint rooting), runs `-m ECMK07+F+R4` with 1000 UFBoot on the **pruned codon alignment**, then strips UFBoot labels and outgroup-roots (two-step `set_outgroup`). Publishes `.treefile`, `.contree`, `.iqtree`, `.log` and `<og>.rooted.nwk` — not the whole IQ-TREE prefix set.

   **Do not add `-asr`/`--rate` here** to try to save csubst its own IQ-TREE call. csubst already runs IQ-TREE with `-te` (fixed topology, no tree search, no bootstrap) on the rooted tree, and relabels that tree itself before reconstructing states, so an externally supplied `.state` will not match its node numbering. That pass costs ~13% of the tree search (197 s vs 1456 s measured), and `-asr` changes `.treefile`'s internal labels from `)98:` to `)Node8/99:`, which `root_tree.py` cannot parse.
5. **CSUBST** (container `slheidig/csubst:01`) — `doctor` → `search` (arity 2, foreground from `bin/make_foreground.py`) → `inspect` (ancestral seqs, branch maps) → `sites` on pairs significant for convergence OR divergence (`bin/select_significant_pairs.py`) → per-OG aggregation (`bin/aggregate_csubst.py`). `-safe` IQ-TREE wrapper kept. Only `csubst.tsv` is published from each `search/sites/<pair>/` directory, because that is the only file `aggregate_csubst.py` reads — the other seven per pair exist purely for eyeballing (160 → 22 files per OG). `inspect/` and csubst's own `csubst_iqtree/` are published only under `--publish_csubst_native`. Per-pair figures are off by default (`--csubst_site_plots false`); `./replot_og.sh` regenerates them from published output.
6. **HYPHY** (container from `docker/hyphy/` → pushed as `slheidig/hyphy:2.5.101`) — reprocessed from `selection_scan.sh`, split into `bin/mask_stops.py`, `bin/tag_foreground.py`, `bin/parse_hyphy_persite.py`, `bin/parse_relax.py`, `bin/hyphy_by_site.py`:
   - whole-tree: **FEL, FUBAR, MEME** on pruned codon MSA + rooted tree;
   - foreground: **Contrast-FEL** (warm tips tagged `{Foreground}`) and **RELAX** (`--models Minimal`, gene-level K).
   - one per-site TSV per method; the combined `by_site` table/plot is produced in COLLECT. Since the tree/alignment are the canonical pruned ones, HyPhy sites align 1:1 with csubst sites — no remapping.
7. **MAP_TO_MSA** (`bin/map_predictions_to_msa.py`, generalises `04a_filter_b2b.py`) — projects every per-residue prediction (b2b, AIUPred, IPC charge, DSSP, custom `pred_*`) onto pruned-MSA columns via the gap pattern; per-feature matrices (sequence × position). **One task per OG**, looping over that OG's predictors, rather than one task per (OG, predictor).
8. **COLLECT** (`bin/hyphy_by_site.py` + `bin/combine_og_table.py`, fused) — builds the per-OG HyPhy by-site table, then **the per-OG deliverable**: one row per (sequence, pruned-MSA position) with `og, sequence_id, strain, category, msa_position, original_msa_position, residue_index, aa`, all biophysical predictors, `pI, charge_pH7`, plus position-level evolutionary columns joined on `msa_position` (csubst: convergent/divergent hotspot flags + OCN/ω_C; HyPhy: α, β, dN/dS, p/q-values, FUBAR posteriors, contrast-fel q, selection class).
9. **PLOT_OG / PLOT_DSSP** — reworks `02_msa_cleanup.py` + `04b_plot_b2b.py`: MSA heatmap; per-feature line plots (mean/min–max band per category, colours assigned dynamically from the selected column's values, outgroup grey) with csubst hotspot overlays. Categories ordering/colouring no longer hardcoded to cool/warm. **Which OGs get figures is chosen with `--plot`** (`all` | `none` | a comma-separated list) — 27 PDFs per OG is ~40,000 files at proteome scale. Plotting is deliberately kept in its own processes so changing `--plot` does not invalidate the cached tables under `-resume`; `./replot_og.sh <results> <og>…` regenerates figures from published output.

Global (fan-in over all 1500 OGs):

10. **GLOBAL_STATS** (`bin/global_stats.py`, streaming accumulation like `plot_predictions.py`) — residue-level distributions per predictor (all + split by category), per-OG means, plus evolutionary distributions: ω_C across OGs, convergent/divergent pair counts, dN/dS and fraction of selected sites per method, RELAX K distribution and gene-level summary table; `stats.tsv` + PNGs. The per-OG tables are concatenated **driver-side** with `collectFile` before they reach it, so its command line is three fixed paths rather than one per OG — passing 1500 staged paths is what used to overflow it.

## 4. Output layout

```
results/
  ogs/<OG>/
    alignment/            native mafft out · pruned .aln.fa (AA+codon) · msa_columns.tsv
    biophysics/
      b2b/       native/ + <og>_b2b_mapped.tsv
      aiupred/   native/ + <og>_aiupred_mapped.tsv
      ipc/       native/ + <og>_ipc_mapped.tsv
      dssp/      native/ (<og>.dssp.txt.gz + plddt.tsv) + <og>_dssp_mapped.tsv
      custom/    <og>_custom_mapped.tsv
    evolution/
      iqtree/    native/            (.treefile .contree .iqtree .log)
      csubst/    native/search/     (csubst_b, csubst_cb_2, csubst_cb_stats, log,
                                     sites/<pair>/csubst.tsv)
                 + branch_pairs.tsv, sites.tsv, hotspots.tsv
      hyphy/     native/            (all JSONs + logs)
                 + <method>.persite.tsv, relax.tsv, by_site.tsv
      <OG>.rooted.nwk
    <OG>_combined.tsv     ← one row per residue per sequence per pruned-MSA position
    plots/                only for OGs selected with --plot: MSA heatmap, per-feature
                          line plots, and with --dssp the secstructartist consensus
                          cartoon + per-sequence SS view
  stats/                  global distributions, hotspot/RELAX/omega summaries
  foreground/             foreground.tsv, strain lists, cleaned categories
```

Measured on 8 test OGs, against the pre-scaling layout: **365 → 70 files per OG**
(5.2×) and 3,616 → 974 inodes over the whole run (3.7×). Directory count is the
remaining weak spot — 61 → 45 per OG, of which 22 are `csubst.branch_id*/`
directories holding a single `csubst.tsv` each. That is tracked separately (§9).

## 5. Containers

| Container | Status | Used for |
|---|---|---|
| `slheidig/csubst:01` | exists | cdskit backalign, prune, root_tree (ete4), csubst; iqtree inside |
| `slheidig/og_b2b_pca:latest` | exists | b2btools + all pandas/matplotlib steps (mapping, combine, plots, global stats) |
| `slheidig/hyphy:2.5.101` | built from `docker/hyphy/` | FEL/FUBAR/MEME/Contrast-FEL/RELAX (amd64 + arm64) |
| `ghcr.io/doszilab/aiupred:cpu` | authors' official CPU image (published with [AIUPred-NF](https://github.com/doszilab/AIUPred-NF)); the `aiupred` CLI is called directly with `--force-cpu -b`. Nothing to build; amd64-only | disorder + binding |
| `slheidig/ipc:01` | exists — built from `docker/ipc/` | pI + per-residue charge |
| `slheidig/dssp3:2`, `slheidig/simsapiper:06` | exist, shared with SIMSApiper, nothing to build | mkdssp; secstructartist SS plots |
| `quay.io/biocontainers/mafft` | public | alignment (**module on hydra**) |

IQ-TREE is no longer a module: TREE is fused with the constraint and rooting
steps, so it runs in `slheidig/csubst:01`, which carries iqtree3 alongside python3
and ete4. On hydra the five stdlib/module-covered local processes (ALIGN_AA,
MAKE_FOREGROUND, PRUNE, MAP_TO_MSA, PLOT_OG) run with **no container at all** --
that removes ~6 s of apptainer startup per task and a lot of shared-filesystem
metadata traffic. Every process still declares a container, which is what lets the
`hpc` profile run on any cluster with nothing but singularity installed.

All images must carry a **linux/amd64** entry (the clusters); arm64 is added
wherever every dependency has a native build, purely so the same tag runs on an
Apple Silicon laptop without emulation. No entrypoint tricks, tool on PATH — so
`singularity` conversion on the clusters is trivial. Dockerfiles live in
`docker/<tool>/`, platforms per image in `docker/docker-bake.hcl`, and
`docker/build.sh` is the only supported way to build them (see README).

## 6. Profiles

Everything lives in `nextflow.config`; there is no `conf/` directory. Three
profiles:

```
profiles {
  hydra { apptainer · slurm for the heavy processes · MAFFT module for ALIGN_AA ·
          no container at all for the stdlib/module-covered local steps }
  hpc   { singularity · slurm · containers for EVERYTHING, no host software
          beyond singularity: it overrides no process and declares no modules }
  local { docker · heavy_executor=local, array_size=0 — nothing is submitted }
}
```

**There are no resource labels.** Each process has one complete `withName`
block carrying executor, array batching, cpus, memory, time and container. This
is not a style choice: a `withLabel:`/`withName:` block does **not** merge with a
same-named block elsewhere, it *replaces* it wholesale. With resources on a label
in `nextflow.config` and `array`/`executor` on the same label in
`conf/hydra.config`, every cpus/memory/time setting was silently dropped — IQ-TREE
ran 1000 bootstraps on 2 cpus / 4 GB and `large`/`highmem_single` affected nothing.
Fixing it cut TREE from 2143 s to 1294 s on identical data.

The default executor is **local**, so the eight seconds-long processes run inside
the driver's allocation and are never submitted. Nine ask for a Slurm job:

| process | cpus | mem | time | arrayed |
|---|---|---|---|---|
| TREE | 4 | 4 GB | 1 h | yes |
| CSUBST | 4 | 4 GB | 30 min | yes |
| HYPHY_WHOLETREE / _FOREGROUND | 4 | 2 GB | 30 min | yes |
| B2BTOOLS | 2 | 4 GB | 30 min | yes |
| AIUPRED | 4 | 4 GB | 15 min | yes |
| COLLECT | 2 | 2 GB | 10 min | yes |
| SPLIT_CUSTOM_PREDICTIONS | 4 | 8 GB | 30 min | no (runs once) |
| GLOBAL_STATS | 4 | 16 GB | 1 h | no (runs once) |
| *the other eight* | 1 | 2 GB | 10 min | local executor |

`--array_size 0` disables batching; `--heavy_executor local` runs everything in
one allocation. `errorStrategy` is `retry` once then `ignore`, for **every**
process with no exceptions: individual families are expected to fail and that
must never end a 1500-family run.

## 7. Sizing (measured on 8 test OGs; nothing has run at proteome scale yet)

Throughput-bound, not memory-bound -- and memory was over-provisioned by up to
65x before it was measured. From `sacct` over one complete run (40m01, nothing
cached), peak RSS / max wall-time per task:

| process | peak RSS | max time |
|---|---|---|
| TREE (IQ-TREE + 1000 UFBoot, 4 cpus) | 238 MB | 1294 s |
| HYPHY_WHOLETREE (FEL/FUBAR/MEME) | 172 MB | 735 s |
| HYPHY_FOREGROUND (contrast-fel/RELAX) | 278 MB | 397 s |
| B2BTOOLS | 1000 MB | 377 s |
| SPLIT_CUSTOM_PREDICTIONS (once) | 2921 MB | 137 s |
| CSUBST | 472 MB | 97 s |
| AIUPRED | 1964 MB | 55 s |
| COLLECT | 232 MB | 15 s |
| GLOBAL_STATS (once) | 475 MB | 13 s |
| driver (`nf-bea`) | 1322 MB | 2408 s |

The real memory consumers are the **predictors** (AIUPred ~2 GB, the custom
predictions splitter ~3 GB, b2bTools ~1 GB), not the phylogenetics. Times in
`nextflow.config` are 2x the measured maximum, except B2BTOOLS and CSUBST (larger
OGs carry many more residues / branch pairs) and the two aggregate steps, whose
cost scales with the *number* of OGs rather than one family. Since `time` is
`X * task.attempt`, a retry doubles it again.

Scheduler load, same run vs the pre-scaling pipeline:

| | 8 OGs | projected 1500 OGs |
|---|---|---|
| tasks | 208 -> 127 | 39,000 -> ~22,500 |
| sbatch calls | 208 -> **10** | 39,000 -> **~300** |

Three effects compound: 26 -> 15 tasks per OG (fusion), 5 of those 15 no longer
touch Slurm (local executor), and arrays batch the remaining 10 per OG. The old
`submitRateLimit = '10 sec'` would have spent **108 hours purely submitting**
39,000 jobs.

The driver job must outlive the whole run and hold task metadata for ~22,500
tasks: give it days of walltime, `NXF_OPTS='-Xms1g -Xmx8g'`, 16 GB and ~8 cores
(the cores cap how many local tasks run concurrently, not driver speed). Set
`NXF_ANSI_LOG=false` or the slurm `.out` fills with redraw frames.

**Monitoring.** `-with-trace`, `-with-report` and `-with-timeline` cannot be used:
they inject a `ps`-based collector *into each container*, and five of the seven
images ship no `procps`, so every containerised task exits 1 before running.
Use the two helpers instead -- both read files only, no JVM, no container:

```
./progress.sh          # the per-process [n of m] table NXF_ANSI_LOG=false hides
./failed_tasks.sh      # per-process tally + failures + the first failure's stderr
```

## 8. Repo skeleton

```
bea2/
  pipeline.nf  nextflow.config          # no conf/ -- profiles live in nextflow.config
  modules/{align,predictors,tree,csubst,hyphy,collect}.nf
  bin/       ~22 argparse scripts split out of the two shell pipelines + b2bplots
  assets/    NO_RELAX, NO_PAIRS         # placeholders for optional GLOBAL_STATS inputs
  docker/{hyphy,ipc,csubst,og_b2b_pca}/ + build.sh + docker-bake.hcl
  testdata/  8 OGs (aa/, nuc/, structures/, categories, predictions csv)
  run_test.sh          launch on Hydra (sbatch)
  run_nf.sh            launch locally
  progress.sh          per-process [n of m] table, reads files only
  failed_tasks.sh      failure tally + first failure's stderr, via `nextflow log`
  replot_og.sh         regenerate one OG's figures from published output
```

## 9. Validation & open points

**Validated.** One complete 8-OG run (40m01, nothing cached) exercises every
process including csubst and all five HyPhy methods. Bootstrap confirmed running
(`Generating 1000 samples for ultrafast bootstrap`, `.contree` written).
`--plot <og>,<og>` confirmed to produce figures for exactly the named families.
The hydra module strings are verified by that run: `nxf_module_load` was called
with `MAFFT/7.526-GCC-14.2.0-with-extensions` and `SciPy-bundle`+`matplotlib`,
every such task exited 0 with no apptainer in its `.command.run`, and `PLOT_OG`
produced its 26 PDFs — which needs matplotlib, pandas and numpy to have resolved.

**Not yet validated.** Nothing has run beyond 8 small families. `GLOBAL_STATS`
and `SPLIT_CUSTOM_PREDICTIONS` have never seen more than 8 OGs, and TREE's
21-minute maximum will be exceeded by larger families. `replot_og.sh`'s apptainer
invocations and `failed_tasks.sh`'s no-failure / `all` paths are untested.

**Open, in priority order:**

1. **csubst `search/sites/` directory count** — 22 directories per OG holding a
   single `csubst.tsv` each, ~33,000 at proteome scale, now the dominant inode
   cost. Publish them flattened (`sites/<pair>.tsv`) or not at all, since
   `aggregate_csubst.py` already folds their content into `<og>.sites.tsv`.
   *Tracked in a separate branch.*
2. **Per-OG tables are not aggregated across OGs.** `collectFile` is applied only
   to what `GLOBAL_STATS` consumes; `<og>_combined.tsv`, `hotspots.tsv`,
   `by_site.tsv` and friends remain 1500 separate files (~9,000 files that could
   be ~6 partitioned tables).
3. **`publishDir` creates directories for patterns that match nothing** — 9 empty
   directories appeared in one run, ~13,500 at scale.
4. **OG-level batching** (`buffer(size: N)` so one task handles N OGs) would cut
   per-task overhead ~20x. Only worth it if the inode budget still does not fit
   after 1-3; costs resume granularity and per-OG retry.
5. `COLLECT` and `GLOBAL_STATS` were given containers back during the config
   merge. The green run shows both work container-free on the SciPy-bundle
   module, so `container = null` would save ~6 s of apptainer startup per
   `COLLECT` task (~2.5 CPU-hours at 1500 families).
6. csubst arity fixed at 2. Raising it would *increase* output, not reduce it
   (branch triplets are combinatorially more numerous), and would break both
   `select_significant_pairs.py` and `aggregate_csubst.py`, which hardcode
   `csubst_cb_2.tsv` and two-id pair names.
7. AIUPred default = disorder (`--aiupred_binding` adds binding in the same pass).
