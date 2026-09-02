// Alignment: mafft AA alignment -> back-translated codon alignment -> occupancy pruning.
// The pruned AA MSA is the canonical coordinate system of the whole pipeline.
// Back-translation and pruning are fused: both are stdlib-only python and run in
// seconds, so they are not worth two scheduler round-trips.

process ALIGN_AA {
    tag "$og"
    publishDir "${params.outdir}/ogs/${og}/alignment/native/mafft", mode: 'copy'

    input:
    tuple val(og), path(aa)

    output:
    tuple val(og), path("${og}.aln.fa"), emit: aln

    script:
    """
    mafft --auto --thread ${task.cpus} ${aa} > ${og}.aln.fa
    """
}

process PRUNE {
    tag "$og"
    publishDir "${params.outdir}/ogs/${og}/alignment", mode: 'copy'

    input:
    tuple val(og), path(aln), path(nuc)

    output:
    tuple val(og), path("${og}.aa.pruned.fa"),    emit: aa
    tuple val(og), path("${og}.codon.pruned.fa"), emit: codon
    tuple val(og), path("${og}.msa_columns.tsv"), emit: columns
    path "${og}.codon.fa"

    script:
    """
    map_aa_to_nuc.py ${aln} ${nuc} ${og}
    prune_msa.py --og ${og} --aa ${aln} --codon ${og}.codon.fa --occupancy ${params.occupancy}
    """
}
