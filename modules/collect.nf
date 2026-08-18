// Mapping predictions onto the pruned MSA, the per-OG combined table,
// per-OG plots, and the global statistics over all OGs.

process MAP_TO_MSA {
    tag "$og:$tool"
    label 'small'
    publishDir "${params.outdir}/ogs/${og}/biophysics/${tool}", mode: 'copy'

    input:
    tuple val(og), val(tool), path(pred), path(aln), path(columns)

    output:
    tuple val(og), val(tool), path("${og}_${tool}_mapped.tsv"), emit: mapped

    script:
    """
    map_predictions_to_msa.py --og ${og} --tool ${tool} --aln ${aln} \\
        --columns ${columns} --pred ${pred} --out ${og}_${tool}_mapped.tsv
    """
}

process COMBINE_TABLE {
    tag "$og"
    label 'medium'
    publishDir "${params.outdir}/ogs/${og}", mode: 'copy'

    input:
    tuple val(og), path(mapped), path(hotspots), path(bysite)
    path categories

    output:
    tuple val(og), path("${og}_combined.tsv"), emit: table

    script:
    def hs = hotspots ? "--hotspots ${hotspots}" : ''
    def bs = bysite ? "--bysite ${bysite}" : ''
    """
    combine_og_table.py --og ${og} --categories ${categories} --column '${params.category_column}' \\
        --mapped ${mapped} ${hs} ${bs} --out ${og}_combined.tsv
    """
}

process PLOT_OG {
    tag "$og"
    label 'small'
    publishDir "${params.outdir}/ogs/${og}", mode: 'copy'

    input:
    tuple val(og), path(mapped), path(pruned_aln), path(hotspots)
    path categories

    output:
    path 'plots/*'

    script:
    def hs = hotspots ? "--hotspots ${hotspots}" : ''
    """
    plot_og.py --og ${og} --aln ${pruned_aln} --categories ${categories} \\
        --column '${params.category_column}' --outgroup-level '${params.outgroup_level}' \\
        --mapped ${mapped} ${hs} --outdir plots
    """
}

process GLOBAL_STATS {
    label 'highmem_single'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path combined, stageAs: 'combined/*'
    path relax,    stageAs: 'relax/*'
    path pairs,    stageAs: 'pairs/*'

    output:
    path 'stats/*'

    script:
    def rel = relax ? "--relax relax/*" : ''
    def bp = pairs ? "--branch-pairs pairs/*" : ''
    """
    global_stats.py --combined combined/* ${rel} ${bp} --outdir stats
    """
}
