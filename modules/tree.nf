// Foreground/outgroup definitions from the categories table, and the per-OG
// constrained ML tree. Constraint building, IQ-TREE and outgroup rooting are
// fused into TREE: the two python steps take ~4 s and 9 s around a job that
// takes minutes, and the csubst image carries iqtree3, python3 and ete4.
//
// With --csubst_reuse_iqtree, IQ-TREE additionally writes ancestral states
// (-asr) and per-site rates (--rate) so csubst can reuse them instead of
// running its own IQ-TREE. When the flag is off the two files are created
// empty, so CSUBST keeps a single input signature and decides in bash.

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

process TREE {
    tag "$og"
    label 'highmem_single'
    publishDir "${params.outdir}/ogs/${og}/evolution/iqtree/native", mode: 'copy',
        pattern: "${og}.{treefile,contree,iqtree,log}"
    publishDir "${params.outdir}/ogs/${og}/evolution", mode: 'copy',
        pattern: "${og}.rooted.nwk"

    input:
    tuple val(og), path(aln), path(codon)
    path outgroup

    output:
    tuple val(og), path("${og}.rooted.nwk"), emit: rooted
    tuple val(og), path("${og}.treefile"), path("${og}.state"), path("${og}.rate"),
          path("${og}.iqtree"), path("${og}.log"), emit: asr
    path "${og}.contree", optional: true

    script:
    def bb  = params.ufboot ? "-B ${params.ufboot}" : ''
    def asr = params.csubst_reuse_iqtree ? '-asr --rate' : ''
    def opt = "-safe -s ${codon} --seqtype CODON${params.genetic_code} " +
              "-m ${params.iqtree_model} ${bb} ${asr} -T ${task.cpus}"
    """
    make_constraint.py --og ${og} --aln ${aln} --outgroup ${outgroup}

    IQ=\$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
    OGL=\$(cat ${og}.outgroup_labels.txt)
    if [ -s ${og}.constraint.nwk ] && [ -n "\$OGL" ]; then
        "\$IQ" ${opt} -g ${og}.constraint.nwk -o "\$OGL" --prefix ${og} -redo -quiet
    elif [ -n "\$OGL" ]; then
        # single outgroup taxon: no monophyly constraint needed, root with -o
        "\$IQ" ${opt} -o "\$OGL" --prefix ${og} -redo -quiet
    else
        echo "NOTE: ${og}: no outgroup strains present -> unconstrained tree (midpoint root downstream)"
        "\$IQ" ${opt} --prefix ${og} -redo -quiet
    fi
    touch ${og}.state ${og}.rate

    root_tree.py --og ${og} --tree ${og}.treefile --outgroup ${outgroup} --out ${og}.rooted.nwk
    """
}
