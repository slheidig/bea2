#!/bin/bash
#SBATCH --job-name=nf-bea
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G

# The 'small' processes now run on the local executor, inside THIS allocation
# (see conf/hydra.config); everything heavier goes to Slurm as job arrays.
# conf/hydra.config reads --cpus-per-task/--mem back from Slurm, so raising them
# here is all that is needed to run more short tasks concurrently.

APPTAINERCACHE=$VSC_SCRATCH_VO_USER/.apptainer
export APPTAINER_CACHEDIR=$APPTAINERCACHE

export NXF_OPTS='-Xms1g -Xmx8g'   # driver holds task metadata for the whole run
export NXF_ANSI_LOG=false         # otherwise the slurm .out fills with redraws

module load Nextflow/25.04.7

house=$VSC_SCRATCH_VO_USER/bea2
nextflow run pipeline.nf \
    -profile hydra \
    --aa_dir      $house/testdata/aa \
    --nuc_dir     $house/testdata/nuc \
    --categories  $house/testdata/synechococcus_categories.tsv \
    --category_column temp_cat2 \
    --og_pattern  'CK_\d+' \
    --predictions_csv $house/testdata/synechococcus_predictions_all.csv \
    --dssp \
    --structure_dir $house/testdata/structures \
    --outdir      $house/nnresults \
    --plot        CK_00000105,CK_00001561 \
    --apptainercache $APPTAINERCACHE \
    --csubst_reuse_iqtree -with-trace -resume

# ---------------------------------------------------------------------------
# New parameters
#
#   --plot all|none|<og>,<og>   which families get figures. 'all' is the old
#                               behaviour. Use a list for QC at scale:
#                                 --plot CK_00000105,CK_00001561
#                               Changing it re-runs only the plot tasks;
#                               every table stays cached under -resume.
#                               bin/replot_og.sh <results> <og>... regenerates
#                               plots from published output if resume is broken.
#
#   --csubst_site_plots         now defaults to false. bin/replot_og.sh
#                               regenerates the per-branch-pair figures, so
#                               toggling this does not re-run a csubst search.
#
#   --publish_csubst_native     also publish inspect/ and csubst's own
#                               csubst_iqtree/. Default false: only the
#                               csubst.tsv that aggregate_csubst.py reads is
#                               kept (160 -> 20 files per OG).
#
#   --csubst_reuse_iqtree       hand TREE's -asr/--rate output to csubst so it
#                               stops re-running IQ-TREE. DEFAULT OFF, still
#                               needs validating: run one OG with it on and
#                               confirm search/csubst_cb_2.tsv is unchanged
#                               against the same OG's current output.
#                               Requires --iqtree_model ECMK07+F+R4 (guarded).
#
# ---------------------------------------------------------------------------
# For the full ~1500-family run, change the header to roughly:
#   #SBATCH --time=3-00:00:00
#   #SBATCH --cpus-per-task=16
#   #SBATCH --mem=48G
# and add, to keep work/ off the project quota:
#   export NXF_WORK=$VSC_SCRATCH_VO_USER/bea2_work
# Consider --plot none (40k PDFs otherwise) and keep -with-trace: it is the
# only record of which OGs hit errorStrategy 'ignore'.
# ---------------------------------------------------------------------------

#nextflow clean -f -before $(nextflow log -q | tail -n 1)
