#!/usr/bin/env bash
# Regenerate an OG's plots from a published results folder, the same way the
# pipeline would. Use it for QC on families that ran with --plot none, or when
# -resume will not reuse the cached tables.
#
#   ./replot_og.sh <results_dir> <og> [<og> ...]
#
# Reads only published output. Writes into <results_dir>/ogs/<og>/plots/.
#
# matplotlib/pandas come from the cluster modules, so plot_og.py and
# hyphy_by_site.py run directly. secstructartist and csubst do not, so those two
# go through their apptainer images (override paths with SIMSAPIPER_IMG /
# CSUBST_IMG, or BEA2_APPTAINER_DIR to point at a different image directory).
set -euo pipefail

if [ $# -lt 2 ]; then
    echo "usage: $0 <results_dir> <og> [<og> ...]" >&2
    exit 1
fi

results=$1; shift
repo=$(cd "$(dirname "$0")" && pwd)
here=$repo/bin

# not APPTAINER_CACHEDIR: that can point at a staging dir rather than the images
cache=${BEA2_APPTAINER_DIR:-${VSC_SCRATCH_VO_USER:-$repo}/.apptainer}
SIMSAPIPER_IMG=${SIMSAPIPER_IMG:-$cache/slheidig-simsapiper-06.img}
CSUBST_IMG=${CSUBST_IMG:-$cache/slheidig-csubst-01.img}

cats=$results/foreground/categories_clean.tsv
column=${CATEGORY_COLUMN:-temp_cat2}
outgroup_level=${OUTGROUP_LEVEL:-Outgroup}
consensus=${SS_CONSENSUS:-0.5}
genetic_code=${GENETIC_CODE:-11}
plot_format=${CSUBST_PLOT_FORMAT:-pdf}

[ -f "$cats" ] || { echo "missing $cats" >&2; exit 1; }

# run a command inside an apptainer image, with the repo and results bound in
in_img() {
    local img=$1; shift
    apptainer exec --no-home --env PREPEND_PATH=/opt/conda/bin \
        -B "$repo" -B "$(cd "$results" && pwd)" "$img" "$@"
}

for og in "$@"; do
    d=$results/ogs/$og
    [ -d "$d" ] || { echo "!! no such OG: $d" >&2; continue; }
    echo "== $og"
    mkdir -p "$d/plots"

    mapped=$(find "$d/biophysics" -name "${og}_*_mapped.tsv" 2>/dev/null | sort | tr '\n' ' ')
    aln=$d/alignment/$og.aa.pruned.fa
    hotspots=$d/evolution/csubst/$og.hotspots.tsv
    hs=""
    [ -s "$hotspots" ] && hs="--hotspots $hotspots"

    if [ -n "$mapped" ] && [ -s "$aln" ]; then
        "$here/plot_og.py" --og "$og" --aln "$aln" --categories "$cats" \
            --column "$column" --outgroup-level "$outgroup_level" \
            --mapped $mapped $hs --outdir "$d/plots"
    else
        echo "  !! skipping plot_og.py (need mapped tables and $aln)"
    fi

    # secstructartist lives only in the simsapiper image
    dssp=$d/biophysics/dssp/${og}_dssp_mapped.tsv
    if [ -s "$dssp" ]; then
        if [ -f "$SIMSAPIPER_IMG" ]; then
            in_img "$SIMSAPIPER_IMG" "$here/plot_dssp_ss.py" --og "$og" --mapped "$dssp" \
                --categories "$cats" --column "$column" \
                --outgroup-level "$outgroup_level" --consensus "$consensus" \
                $hs --outdir "$d/plots"
        else
            echo "  !! skipping plot_dssp_ss.py (no image at $SIMSAPIPER_IMG)"
        fi
    fi

    persite=$(find "$d/evolution/hyphy" -name "$og.*persite.tsv" 2>/dev/null | sort | tr '\n' ' ')
    columns=$d/alignment/$og.msa_columns.tsv
    if [ -n "$persite" ] && [ -s "$columns" ]; then
        ( cd "$d/plots" && "$here/hyphy_by_site.py" --og "$og" \
            --columns "$columns" --persite $persite --plot yes )
    fi

    # csubst per-pair figures: needs the native search output and the codon MSA
    cb=$d/evolution/csubst/native/search/csubst_cb_2.tsv
    codon=$d/alignment/$og.codon.pruned.fa
    rooted=$d/evolution/$og.rooted.nwk
    if [ -s "$cb" ] && [ -s "$codon" ] && [ -s "$rooted" ] && [ -f "$CSUBST_IMG" ]; then
        "$here/select_significant_pairs.py" --cb "$cb" \
            --ocn "${OCN_CUTOFF:-0.5}" --omega "${OMEGA_CUTOFF:-1.0}" > "$d/plots/pairs.txt" || true
        while read -r p; do
            [ -n "$p" ] || continue
            in_img "$CSUBST_IMG" csubst sites --alignment_file "$codon" \
                --rooted_tree_file "$rooted" --genetic_code "$genetic_code" \
                --branch_id "$p" --cb_file "$cb" \
                --tree_site_plot yes --site_state_plot yes --site_summary_plot yes \
                --tree_site_plot_format "$plot_format" \
                --outdir "$d/plots/csubst_sites" >/dev/null 2>&1 \
                || echo "  !! csubst sites failed for pair $p"
        done < "$d/plots/pairs.txt"
    fi
done

echo "done -> $results/ogs/<og>/plots/"
