// csubst convergence analysis (one process per OG).
//
// csubst runs on its own defaults: every flag below is a path, the genetic
// code, the thread count, or --output_stat (adds spe2spe/dif2spe). Only
// convergence and only is_fg pairs reach the tidy tables; divergence and the
// background stay in the native output.
//
// aggregate_csubst.py reads csubst.tsv from each search/sites/<pair>/; the
// rest of that directory is published only under --publish_csubst_native,
// together with inspect/ and csubst's own csubst_iqtree/.
//
// csubst runs its own IQ-TREE to reconstruct ancestral states, but it already
// does so with -te on the rooted tree from TREE: fixed topology, no tree
// search, no bootstrap (~13% of the tree search's wall time). It relabels that
// tree itself first, so an externally supplied .state would not match its node
// numbering -- do not try to feed it one.

process CSUBST {
    tag "$og"
    publishDir "${params.outdir}/ogs/${og}/evolution/csubst", mode: 'copy',
        pattern: '*.tsv'
    publishDir "${params.outdir}/ogs/${og}/evolution/csubst/native", mode: 'copy',
        pattern: 'search/csubst*', enabled: !params.publish_csubst_native
    publishDir "${params.outdir}/ogs/${og}/evolution/csubst/native", mode: 'copy',
        pattern: 'search/sites/*/csubst.tsv', enabled: !params.publish_csubst_native
    publishDir "${params.outdir}/ogs/${og}/evolution/csubst/native", mode: 'copy',
        pattern: 'search/**', enabled: params.publish_csubst_native
    publishDir "${params.outdir}/ogs/${og}/evolution/csubst/native", mode: 'copy',
        pattern: '{inspect,csubst_iqtree}/**', enabled: params.publish_csubst_native

    input:
    tuple val(og), path(codon), path(rooted)
    path foreground

    output:
    tuple val(og), path("${og}.branch_pairs.tsv"), path("${og}.sites.tsv"), path("${og}.hotspots.tsv"), emit: tables
    tuple val(og), path("${og}.hotspots.tsv"), emit: hotspots
    path 'search/csubst*'
    path 'search/sites/**',                     optional: true
    path 'inspect/**',                          optional: true
    path 'csubst_iqtree/**',                    optional: true

    script:
    """
    # numerically-safe IQ-TREE wrapper for every call csubst makes internally
    IQ=\$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
    printf '#!/usr/bin/env bash\\nexec "%s" "\$@" -safe\\n' "\$IQ" > iqtree_safe.sh
    chmod +x iqtree_safe.sh
    IQX="\$PWD/iqtree_safe.sh"

    csubst doctor --alignment_file ${codon} --rooted_tree_file ${rooted} \\
        --iqtree_exe "\$IQX" --genetic_code ${params.genetic_code} | grep 'Doctor summary' || true

    csubst search --alignment_file ${codon} --rooted_tree_file ${rooted} \\
        --foreground ${foreground} --genetic_code ${params.genetic_code} \\
        --iqtree_exe "\$IQX" --threads ${task.cpus} \\
        --output_stat ${params.csubst_output_stat} --outdir search

    csubst inspect --alignment_file ${codon} --rooted_tree_file ${rooted} \\
        --genetic_code ${params.genetic_code} --iqtree_exe "\$IQX" \\
        --outdir inspect || true

    # per-site analysis on every foreground pair convergent by --cutoff_stat
    select_significant_pairs.py --cb search/csubst_cb_2.tsv \\
        --stats search/csubst_cb_stats.tsv > pairs.txt || true
    while read -r p; do
        [ -n "\$p" ] || continue
        csubst sites --alignment_file ${codon} --rooted_tree_file ${rooted} \\
            --genetic_code ${params.genetic_code} --iqtree_exe "\$IQX" --branch_id "\$p" \\
            --cb_file search/csubst_cb_2.tsv \\
            --outdir search/sites > /dev/null 2>&1 || echo "  !! csubst sites failed for pair \$p"
    done < pairs.txt

    aggregate_csubst.py --og ${og} --search-dir search
    """
}
