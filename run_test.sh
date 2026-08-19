#!/usr/bin/env bash
# Launch bea2. Adapt paths/profile and run:  ./run_nf.sh
# Local (Mac): conda activate nfconda   (nextflow lives in ~/miniconda3/envs/nfconda)
set -euo pipefail

#!/bin/bash
#SBATCH --job-name=nf-bea
#SBATCH --time=1:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G


APPTAINERCACHE=$VSC_SCRATCH_VO_USER/.apptainer
export APPTAINER_CACHEDIR=$APPTAINERCACHE

module load Nextflow

house=$VSC_SCRATCH_VO_USER
nextflow run pipeline.nf \
    -profile hydra \
    --aa_dir      $house/testdata/aa \
    --nuc_dir     $house/testdata/nucs \
    --categories  $house/testdata/synechococcus_categories.tsv \
    --category_column temp_cat2 \
    --og_pattern  'CK_\d+' \
    --predictions_csv $house/results_haochen/synechococcus_predictions_all.csv \
    --outdir      $house/results \
    --apptainercache $APPTAINERCACHE


# On the clusters:
#   -profile hydra     (Slurm + MAFFT/IQ-TREE modules + singularity)
#   -profile hpc       (Slurm + singularity for everything)
#     [--slurm_account ... --slurm_queue ...]
