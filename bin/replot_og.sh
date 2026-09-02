#!/usr/bin/env bash
# Regenerate an OG's plots from a published results folder, the same way the
# pipeline would. Use it for QC on families that ran with --plot none, or when
# -resume will not reuse the cached tables.
#
#   bin/replot_og.sh <results_dir> <og> [<og> ...]
#
# Reads only published output. Writes into <results_dir>/ogs/<og>/plots/.
# csubst per-pair figures are regenerated too when the csubst native output is
# present (that is what --csubst_site_plots would have produced in-pipeline).
set -euo pipefail

if [ $# -lt 2 ]; then
    echo "usage: $0 <results_dir> <og> [<og> ...]" >&2
    exit 1
fi

results=$1; shift
here=$(cd "$(dirname "$0")" && pwd)
cats=$results/foreground/categories_clean.tsv
column=${CATEGORY_COLUMN:-temp_cat2}
outgroup_level=${OUTGROUP_LEVEL:-Outgroup}
consensus=${SS_CONSENSUS:-0.5}
genetic_code=${GENETIC_CODE:-11}
plot_format=${CSUBST_PLOT_FORMAT:-pdf}

[ -f "$cats" ] || { echo "missing $cats" >&2; exit 1; }

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

    dssp=$d/biophysics/dssp/${og}_dssp_mapped.tsv
    if [ -s "$dssp" ]; then
        "$here/plot_dssp_ss.py" --og "$og" --mapped "$dssp" --categories "$cats" \
            --column "$column" --outgroup-level "$outgroup_level" \
            --consensus "$consensus" $hs --outdir "$d/plots"
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
    if [ -s "$cb" ] && [ -s "$codon" ] && [ -s "$rooted" ] && command -v csubst >/dev/null; then
        "$here/select_significant_pairs.py" --cb "$cb" \
            --ocn "${OCN_CUTOFF:-0.5}" --omega "${OMEGA_CUTOFF:-1.0}" > "$d/plots/pairs.txt" || true
        while read -r p; do
            [ -n "$p" ] || continue
            csubst sites --alignment_file "$codon" --rooted_tree_file "$rooted" \
                --genetic_code "$genetic_code" --branch_id "$p" --cb_file "$cb" \
                --tree_site_plot yes --site_state_plot yes --site_summary_plot yes \
                --tree_site_plot_format "$plot_format" \
                --outdir "$d/plots/csubst_sites" >/dev/null 2>&1 \
                || echo "  !! csubst sites failed for pair $p"
        done < "$d/plots/pairs.txt"
    fi
done

echo "done -> $results/ogs/<og>/plots/"
