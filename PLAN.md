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
| `--ocn_cutoff`, `--omega_cutoff`, `--pvalue`, `--posterior` | as in current scripts | significance thresholds |
| `--outdir` | `results` | |

Strain parsing from sequence IDs (`<og>_<strain>_<locus>`, strains may contain `_`) reuses the robust rsplit logic of `01_find_seqsoi.py`. Strains absent from the categories table (e.g. WH8020) stay in tree/alignment as uncategorised background and are flagged in a warnings file, never silently dropped.

## 3. Workflow (modules and processes)

```
pipeline.nf
modules/
  align.nf        ALIGN_AA (mafft) · BACKALIGN (cdskit) · PRUNE_MSA
  predictors.nf   B2BTOOLS · AIUPRED · IPC · SPLIT_CUSTOM_PREDICTIONS
  tree.nf         IQTREE · ROOT_TREE
  csubst.nf       CSUBST_SEARCH · CSUBST_SITES · CSUBST_AGGREGATE
  hyphy.nf        HYPHY_WHOLETREE · HYPHY_FOREGROUND · HYPHY_BYSITE
  collect.nf      MAP_TO_MSA · COMBINE_TABLE · PLOT_OG · GLOBAL_STATS
```

Per OG (channels joined on OG id):

1. **ALIGN_AA** — `mafft --auto` on the AA fasta → `<og>.aln.fa`. *(Hydra: `module MAFFT`; hpc: mafft container.)*
2. **BACKALIGN** — `cdskit backalign` maps the full AA alignment onto the CDS → in-frame codon alignment (before pruning, so codons stay consistent).
3. **PRUNE_MSA** (`bin/prune_msa.py`) — drops AA columns with occupancy < `--occupancy`; applies the *same* columns ×3 to the codon alignment. Outputs: pruned AA MSA, pruned codon MSA, `msa_columns.tsv` (pruned_pos ↔ original_pos), occupancy per column.
4. **Predictors** (all on *ungapped* AA sequences, mapped to the MSA later):
   - **B2BTOOLS** (container `slheidig/og_b2b_pca:latest`, all predictors) → long per-residue TSV (backbone, sidechain, helix, sheet, coil, ppII, earlyFolding, disoMine, agmata).
   - **AIUPRED** (new container, see §5) → per-residue disorder score (default mode; binding mode available via `ext.args`).
   - **IPC** (new container) → per-sequence pI (+ per-residue charge at pH 7, which is positional); pI repeated per residue row in the tables.
   - **SPLIT_CUSTOM_PREDICTIONS** — runs **once**: streams the 1.9 GB CSV in chunks, splits rows by the `CK_########` prefix of `block_id` → one `<og>_custom.tsv` per OG, joined into the per-OG channel. (Note: `residue_index` there is 1-based ungapped.)
5. **IQTREE** — GTR+G, 1000 UFBoot, on the **pruned codon alignment**, with outgroup-monophyly constraint built from `--outgroup_level` (partial outgroup sets OK; none → unconstrained). Publishes the **complete IQ-TREE native output**. *(Hydra: `module IQ-TREE`; hpc: container.)*
6. **ROOT_TREE** (`bin/root_tree.py`, ete4, in the csubst container) — strip UFBoot labels, outgroup-root (two-step `set_outgroup`), midpoint fallback — split out from csubst_pipeline step 3.
7. **CSUBST** (container `slheidig/csubst:01`) — `doctor` → `search` (arity 2, foreground from `bin/make_foreground.py`) → `inspect` (ancestral seqs, branch maps) → `sites` on pairs significant for convergence OR divergence (`bin/select_significant_pairs.py`) → per-OG aggregation (`bin/aggregate_csubst.py`). The process publishes its **entire work products including `csubst_iqtree/`** (the full IQ-TREE intermediates csubst generates — `.state` files, ancestral reconstruction — as you required). `-safe` IQ-TREE wrapper kept.
8. **HYPHY** (container from your `Dockerfile` → pushed as e.g. `slheidig/hyphy:2.5.1`) — reprocessed from `selection_scan.sh`, split into `bin/mask_stops.py`, `bin/tag_foreground.py`, `bin/parse_hyphy_persite.py`, `bin/parse_relax.py`, `bin/hyphy_by_site.py`:
   - whole-tree: **FEL, FUBAR, MEME** on pruned codon MSA + rooted tree;
   - foreground: **Contrast-FEL** (warm tips tagged `{Foreground}`) and **RELAX** (`--models Minimal`, gene-level K).
   - per-OG: one per-site TSV per method + the combined `by_site` table/plot. Since the tree/alignment are the canonical pruned ones, HyPhy sites align 1:1 with csubst sites — no remapping.
9. **MAP_TO_MSA** (`bin/map_predictions_to_msa.py`, generalises `04a_filter_b2b.py`) — projects every per-residue prediction (b2b, AIUPred, IPC charge, custom `pred_*`) onto pruned-MSA columns via the gap pattern; per-feature matrices (sequence × position).
10. **COMBINE_TABLE** (`bin/combine_og_table.py`) — **the per-OG deliverable**: one row per (sequence, pruned-MSA position) with: `og, sequence_id, strain, category, msa_position, original_msa_position, residue_index, aa`, all biophysical predictors, `pI, charge_pH7`, plus position-level evolutionary columns joined on `msa_position` (csubst: convergent/divergent hotspot flags + OCN/ω_C; HyPhy: α, β, dN/dS, p/q-values, FUBAR posteriors, contrast-fel q, selection class).
11. **PLOT_OG** — reworks `02_msa_cleanup.py` + `04b_plot_b2b.py`: MSA heatmap; per-feature line plots (mean/min–max band per category, colours assigned dynamically from the selected column's values, outgroup grey) with csubst hotspot overlays; HyPhy by-site plot. Categories ordering/colouring no longer hardcoded to cool/warm.

Global (fan-in over all 1500 OGs):

12. **GLOBAL_STATS** (`bin/global_stats.py`, streaming accumulation like `plot_predictions.py`) — residue-level distributions per predictor (all + split by category), per-OG means, plus evolutionary distributions: ω_C across OGs, convergent/divergent pair counts, dN/dS and fraction of selected sites per method, RELAX K distribution and gene-level summary table; `stats.tsv` + PNGs.

## 4. Output layout

```
results/
  ogs/<OG>/
    alignment/            native mafft out · pruned .aln.fa (AA+codon) · msa_columns.tsv
    biophysics/
      b2btools/  native/ + b2btools.tsv
      aiupred/   native/ + aiupred.tsv
      ipc/       native/ + ipc.tsv
      custom/    <og>_custom.tsv
    evolution/
      iqtree/    native/            (full intermediates)
      csubst/    native/            (search, inspect, sites, csubst_iqtree/ with .state)
                 + branch_pairs.tsv, sites.tsv, hotspots.tsv
      hyphy/     native/            (all JSONs + logs)
                 + <method>.persite.tsv, relax_summary.tsv, by_site.tsv
    <OG>_combined.tsv     ← one row per residue per sequence per pruned-MSA position
    plots/
  stats/                  global distributions, hotspot/RELAX/omega summaries, stats.tsv, warnings.tsv
```

## 5. Containers

| Container | Status | Used for |
|---|---|---|
| `slheidig/csubst:01` | exists | cdskit backalign, prune, root_tree (ete4), csubst; iqtree inside |
| `slheidig/og_b2b_pca:latest` | exists | b2btools + all pandas/matplotlib steps (mapping, combine, plots, global stats) |
| `slheidig/hyphy:2.5.1` | Dockerfile exists (`Dockerfile` in csubst_trial) → push | FEL/FUBAR/MEME/Contrast-FEL/RELAX (amd64-pinned) |
| `slheidig/aiupred` | **to build** — `python:3.11-slim` + torch-CPU + AIUPred from GitHub, models baked in at build time, plus a thin `bin/run_aiupred.py` CLI (fasta in → TSV out) | disorder |
| `slheidig/ipc` | **to build** — `python:3.11-slim` + IPC 2.x from isoelectric.org (free for academia) + a wrapper emitting per-sequence pI and per-residue charge | pI |
| mafft / iqtree | biocontainers images on `hpc`; **modules on hydra** | alignment, trees |

All images: linux/amd64, no entrypoint tricks, tool on PATH — so `singularity` conversion on the clusters is trivial. Dockerfiles live in `docker/<tool>/` in the repo.

## 6. Profiles

```
profiles {
  hydra { slurm · singularity (cacheDir $VSC_SCRATCH/.apptainer) ·
          withName MAFFT:  module 'MAFFT/...'
          withName IQTREE: module 'IQ-TREE/...'   // exact module strings from you
          all other processes: container }
  hpc   { slurm · singularity · containers for everything incl. mafft/iqtree }
  local { local executor · docker · reduced maxForks — for testing on your Mac }
  test  { local + 4 bundled OGs from csubst_trial as tests/data }
}
```

Resource labels (bioenvada/proteinfam style, `task.attempt` scaling, retry×2 then ignore so one bad OG never kills the run):

| Label | cpus | mem | time | processes |
|---|---|---|---|---|
| `small` | 1 | 2 GB | 30 min | mafft, backalign, prune, root, map, combine, plots, ipc |
| `medium` | 2 | 4 GB | 2 h | iqtree, b2btools, fubar |
| `large` | 4 | 8 GB | 8 h | csubst, fel/meme/contrast-fel/relax, aiupred |
| `highmem_single` | 4 | 16 GB | 4 h | split_custom_predictions, global_stats (run once) |

## 7. Sizing advice (40 cores, 1500 OGs × 54 seqs)

Per-OG tasks are small; the run is throughput-bound, not memory-bound. Plan on **~4 GB RAM per core ⇒ a 160 GB allocation comfortably saturates 40 cores** (peak concurrent usage realistically 60–100 GB; AIUPred and csubst are the fattest at ~4–8 GB each). Rough per-OG wall-times: mafft seconds; IQ-TREE 2–10 min; csubst search+sites 10–40 min (21 foreground lineages → 210 pairs); FEL/MEME 5–20 min each, FUBAR fast, RELAX-Minimal 10–30 min; b2btools 1–2 min; AIUPred <1 min. Aggregate ≈ 1–1.5 CPU-h per OG ⇒ **1500 OGs ≈ 1500–2200 CPU-h ≈ 2–3 days wall-time on 40 cores**. MEME/RELAX dominate — if that's too long, the toggles let you run a FEL-only first pass. Set `executor.queueSize ≈ 100` on Slurm and let the scheduler pack; no single job needs >4 cores.

## 8. Repo skeleton

```
bioenvada2/
  pipeline.nf  nextflow.config  conf/{hydra,hpc,local}.config
  modules/{align,predictors,tree,csubst,hyphy,collect}.nf
  bin/  (≈14 argparse scripts split out of the two shell pipelines + b2bplots)
  docker/{aiupred,ipc,hyphy}/Dockerfile
  tests/data/  (the 4 csubst_trial OGs)  run_nf.sh  README.md
```

## 9. Validation & open points

- **Validation:** run the `test` profile on the 4 csubst_trial families and diff csubst branch pairs / hotspot sites and HyPhy per-site classes against your existing `csubst_run` results (differences expected only where hammer- vs occupancy-pruning disagree on columns).
- Need from you: exact **Hydra module strings** for MAFFT/IQ-TREE (and cluster/partition names), and a decision on pushing the two existing Dockerfiles to Docker Hub under `slheidig/`.
- AIUPred default = disorder prediction (`binding` exposed as an option).
- csubst arity fixed at 2 for now (param-ready to raise later).
- IPC: confirm IPC 2.0 (isoelectric.org download) vs the pip `ipc2`-style reimplementations — I'll pin whichever you prefer in the Dockerfile.
