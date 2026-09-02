// csubst convergence/divergence analysis (one process per OG).
//
// Only csubst.tsv is read from each search/sites/<branch pair>/ directory
// (see aggregate_csubst.py); the other seven files per pair are unread, so
// they are neither declared nor published. The full native tree, including
// inspect/ and csubst's own csubst_iqtree/, is published only under
// --publish_csubst_native.
//
// With --csubst_reuse_iqtree, TREE's ancestral-state and rate files are handed
// to csubst so it skips its internal IQ-TREE run. The files are empty when the
// flag is off, which is what the [ -s ] test below detects.

process CSUBST {
    tag "$og"
    label 'large'
    publishDir "${params.outdir}/ogs/${og}/evolution/csubst", mode: 'copy',
        pattern: '*.tsv'
    publishDir "${params.outdir}/ogs/${og}/evolution/csubst/native", mode: 'copy',
        pattern: 'search/**'
    publishDir "${params.outdir}/ogs/${og}/evolution/csubst/native", mode: 'copy',
        pattern: '{inspect,csubst_iqtree}/**', enabled: params.publish_csubst_native

    input:
    tuple val(og), path(codon), path(rooted), path(treefile), path(state), path(rate), path(iqfile), path(iqlog)
    path foreground

    output:
    tuple val(og), path("${og}.branch_pairs.tsv"), path("${og}.sites.tsv"), path("${og}.hotspots.tsv"), emit: tables
    tuple val(og), path("${og}.hotspots.tsv"), emit: hotspots
    path 'search/csubst*'
    path 'search/sites/*/csubst.tsv',           optional: true
    path 'inspect/**',                          optional: true
    path 'csubst_iqtree/**',                    optional: true

    script:
    def plots = params.csubst_site_plots ? 'yes' : 'no'
    """
    # numerically-safe IQ-TREE wrapper for every call csubst makes internally
    IQ=\$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
    printf '#!/usr/bin/env bash\\nexec "%s" "\$@" -safe\\n' "\$IQ" > iqtree_safe.sh
    chmod +x iqtree_safe.sh
    IQX="\$PWD/iqtree_safe.sh"

    # reuse the pipeline's own IQ-TREE run instead of letting csubst redo it
    REUSE=""
    if [ -s ${state} ] && [ -s ${rate} ]; then
        REUSE="--iqtree_treefile ${treefile} --iqtree_state ${state} --iqtree_rate ${rate}"
        REUSE="\$REUSE --iqtree_iqtree ${iqfile} --iqtree_log ${iqlog}"
        echo "csubst: reusing IQ-TREE output from TREE (${treefile})"
    fi

    csubst doctor --alignment_file ${codon} --rooted_tree_file ${rooted} \\
        --iqtree_exe "\$IQX" --genetic_code ${params.genetic_code} \$REUSE | grep 'Doctor summary' || true

    csubst search --alignment_file ${codon} --rooted_tree_file ${rooted} \\
        --foreground ${foreground} --fg_format 1 --genetic_code ${params.genetic_code} \\
        --iqtree_exe "\$IQX" --threads ${task.cpus} --blas_threads 1 --outdir search \$REUSE

    csubst inspect --alignment_file ${codon} --rooted_tree_file ${rooted} \\
        --genetic_code ${params.genetic_code} --iqtree_exe "\$IQX" \\
        --plot_state_aa no --plot_state_codon no --outdir inspect \$REUSE || true

    # per-site analysis on every foreground pair significant for convergence OR divergence
    select_significant_pairs.py --cb search/csubst_cb_2.tsv \\
        --ocn ${params.ocn_cutoff} --omega ${params.omega_cutoff} > pairs.txt || true
    while read -r p; do
        [ -n "\$p" ] || continue
        csubst sites --alignment_file ${codon} --rooted_tree_file ${rooted} \\
            --genetic_code ${params.genetic_code} --iqtree_exe "\$IQX" --branch_id "\$p" \\
            --cb_file search/csubst_cb_2.tsv \\
            --tree_site_plot ${plots} --site_state_plot ${plots} --site_summary_plot ${plots} \\
            --tree_site_plot_format ${params.csubst_plot_format} \\
            --outdir search/sites \$REUSE > /dev/null 2>&1 || echo "  !! csubst sites failed for pair \$p"
    done < pairs.txt

    aggregate_csubst.py --og ${og} --search-dir search \\
        --ocn ${params.ocn_cutoff} --omega ${params.omega_cutoff}
    """
}
