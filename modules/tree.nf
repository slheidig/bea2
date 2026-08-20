// Foreground/outgroup definitions from the categories table, per-OG constrained
// ML trees (IQ-TREE) and outgroup rooting (ete4).

process MAKE_FOREGROUND {
    label 'small'
    publishDir "${params.outdir}/foreground", mode: 'copy'

    input:
    path categories

    output:
    path 'foreground.tsv',        emit: foreground
    path 'fg_strains.txt',        emit: fg_strains
    path 'outgroup_strains.txt',  emit: outgroup
    path 'categories_clean.tsv',  emit: clean

    script:
    """
    make_foreground.py --categories ${categories} --column '${params.category_column}' \\
        --fg-level '${params.fg_level}' --outgroup-level '${params.outgroup_level}'
    """
}

process MAKE_CONSTRAINT {
    tag "$og"
    label 'small'

    input:
    tuple val(og), path(aln)
    path outgroup

    output:
    tuple val(og), path("${og}.constraint.nwk"), path("${og}.outgroup_labels.txt"), emit: files

    script:
    """
    make_constraint.py --og ${og} --aln ${aln} --outgroup ${outgroup}
    """
}

process IQTREE {
    tag "$og"
    label 'highmem_single'
    publishDir "${params.outdir}/ogs/${og}/evolution/iqtree/native", mode: 'copy'

    input:
    tuple val(og), path(codon), path(constraint), path(oglabels)

    output:
    tuple val(og), path("${og}.treefile"), emit: treefile
    path "${og}.*"

    script:
    def bb  = params.ufboot ? "-B ${params.ufboot}" : ''
    def opt = "-safe -s ${codon} --seqtype CODON${params.genetic_code} " +
              "-m ${params.iqtree_model} ${bb} -T ${task.cpus}"
    """
    IQ=\$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
    OGL=\$(cat ${oglabels})
    if [ -s ${constraint} ] && [ -n "\$OGL" ]; then
        "\$IQ" ${opt} -g ${constraint} -o "\$OGL" --prefix ${og} -redo -quiet
    elif [ -n "\$OGL" ]; then
        # single outgroup taxon: no monophyly constraint needed, root with -o
        "\$IQ" ${opt} -o "\$OGL" --prefix ${og} -redo -quiet
    else
        echo "NOTE: ${og}: no outgroup strains present -> unconstrained tree (midpoint root downstream)"
        "\$IQ" ${opt} --prefix ${og} -redo -quiet
    fi
    """
}

process ROOT_TREE {
    tag "$og"
    label 'small'
    publishDir "${params.outdir}/ogs/${og}/evolution", mode: 'copy'

    input:
    tuple val(og), path(treefile)
    path outgroup

    output:
    tuple val(og), path("${og}.rooted.nwk"), emit: rooted

    script:
    """
    root_tree.py --og ${og} --tree ${treefile} --outgroup ${outgroup} --out ${og}.rooted.nwk
    """
}
