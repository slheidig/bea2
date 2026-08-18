#!/usr/bin/env bash
# Launch bea2. Adapt paths/profile and run:  ./run_nf.sh
# Local (Mac): conda activate nfconda   (nextflow lives in ~/miniconda3/envs/nfconda)
set -euo pipefail

nextflow run pipeline.nf \
    -profile local \
    --aa_dir      data/aa \
    --nuc_dir     data/nucs \
    --categories  synechococcus_categories.tsv \
    --category_column temp_cat2 \
    --og_pattern  'CK_\d+' \
    --predictions_csv results_haochen/synechococcus_predictions_all.csv \
    --outdir      results \
    -resume "$@"

# On the clusters:
#   -profile hydra     (Slurm + MAFFT/IQ-TREE modules + singularity)
#   -profile hpc       (Slurm + singularity for everything)
#     [--slurm_account ... --slurm_queue ...]
