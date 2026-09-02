// Mapping predictions onto the pruned MSA, the per-OG combined table,
// per-OG plots, and the global statistics over all OGs.
//
// MAP_TO_MSA handles every predictor of one OG in a single task (it was one
// task per predictor), and COLLECT fuses the by-site table into the combined
// table. Plotting stays in its own processes so that --plot can be changed
// without invalidating the cached tables.

process MAP_TO_MSA {
    tag "$og"
    label 'small'
    publishDir "${params.outdir}/ogs/${og}/biophysics", mode: 'copy',
               saveAs: { f -> f.replaceFirst(/^mapped\//, '') }

    input:
    tuple val(og), val(tools), path(preds), path(aln), path(columns)

    output:
    tuple val(og), path('mapped/*/*'),    emit: mapped
    tuple val(og), path('mapped/dssp/*'), emit: dssp, optional: true

    script:
    def pairs = [tools, preds].transpose().collect { t, p -> "${t}:${p}" }.join(' ')
    """
    for tp in ${pairs} ; do
        tool=\${tp%%:*}
        pred=\${tp#*:}
        mkdir -p mapped/\$tool
        map_predictions_to_msa.py --og ${og} --tool \$tool --aln ${aln} \\
            --columns ${columns} --pred \$pred --out mapped/\$tool/${og}_\${tool}_mapped.tsv
    done
    """
}

process COLLECT {
    tag "$og"
    label 'medium'
    publishDir "${params.outdir}/ogs/${og}/evolution/hyphy", mode: 'copy', pattern: "${og}.by_site.*"
    publishDir "${params.outdir}/ogs/${og}", mode: 'copy', pattern: "${og}_combined.tsv"

    input:
    tuple val(og), path(mapped), path(columns), path(hotspots), path(persite)
    path categories

    output:
    tuple val(og), path("${og}_combined.tsv"), emit: table
    path "${og}.by_site.tsv", optional: true
    path "${og}.by_site.png", optional: true

    script:
    def hs = hotspots ? "--hotspots ${hotspots}" : ''
    def bysite = persite ? "hyphy_by_site.py --og ${og} --columns ${columns} --persite ${persite} --plot yes" : ''
    def bs = persite ? "--bysite ${og}.by_site.tsv" : ''
    """
    ${bysite}
    combine_og_table.py --og ${og} --categories ${categories} --column '${params.category_column}' \\
        --mapped ${mapped} --columns ${columns} ${hs} ${bs} --out ${og}_combined.tsv
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

process PLOT_DSSP {
    tag "$og"
    label 'small'
    publishDir "${params.outdir}/ogs/${og}", mode: 'copy', pattern: 'plots/*.pdf'
    // the consensus table is data, not a figure: it belongs with the other DSSP tables
    publishDir "${params.outdir}/ogs/${og}/biophysics/dssp", mode: 'copy',
               pattern: 'plots/*.tsv', saveAs: { f -> file(f).name }

    input:
    tuple val(og), path(mapped), path(hotspots)
    path categories

    output:
    path 'plots/*'

    script:
    def hs = hotspots ? "--hotspots ${hotspots}" : ''
    """
    plot_dssp_ss.py --og ${og} --mapped ${mapped} --categories ${categories} \\
        --column '${params.category_column}' --outgroup-level '${params.outgroup_level}' \\
        --consensus ${params.ss_consensus} ${hs} --outdir plots
    """
}

// The per-OG tables are concatenated on the driver side (collectFile) before
// they reach this process, so the command line stays a fixed three paths
// instead of one per OG.
process GLOBAL_STATS {
    label 'highmem_single'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path combined
    path relax
    path pairs

    output:
    path 'stats/*'

    script:
    def rel = relax.name != 'NO_RELAX'  ? "--relax ${relax}"        : ''
    def bp  = pairs.name != 'NO_PAIRS'  ? "--branch-pairs ${pairs}" : ''
    """
    global_stats.py --combined ${combined} ${rel} ${bp} --outdir stats
    """
}
