// Biophysical predictors. Each runs on the OG's (de-gapped) amino-acid
// sequences and emits the standard long format:
//   sequence_id  residue_index(1-based, ungapped)  residue  <features...>
// New predictors: add a process here that emits tuple(og, 'toolname', tsv)
// and mix it into the `preds` channel in pipeline.nf — nothing else changes.

process B2BTOOLS {
    tag "$og"
    label 'medium'
    publishDir "${params.outdir}/ogs/${og}/biophysics/b2btools/native", mode: 'copy'

    input:
    tuple val(og), path(aa)

    output:
    tuple val(og), val('b2b'), path("${og}_b2b.tsv"), emit: pred

    script:
    """
    run_b2btools.py --fasta ${aa} --out ${og}_b2b.tsv --tools '${params.b2btools_tools}'
    """
}

process AIUPRED {
    tag "$og"
    label 'large'
    publishDir "${params.outdir}/ogs/${og}/biophysics/aiupred/native", mode: 'copy'

    input:
    tuple val(og), path(aa)

    output:
    tuple val(og), val('aiupred'), path("${og}_aiupred.tsv"), emit: pred

    script:
    def bind = params.aiupred_binding ? '-b' : ''
    def cols = params.aiupred_binding ? 'aiupred_disorder\\taiupred_binding' : 'aiupred_disorder'
    """
    aiupred -i ${aa} -o aiupred_raw.tsv --force-cpu ${bind}

    # AIUPred writes a '#' banner, then '#>id' per sequence followed by 'position residue score...' rows. 
    printf 'sequence_id\\tresidue_index\\tresidue\\t${cols}\\n' > ${og}_aiupred.tsv
    awk '/^#>/ { id = substr(\$0, 3); next }
         /^#/  { next }
         NF >= 3 { print id "\\t" \$0 }' aiupred_raw.tsv >> ${og}_aiupred.tsv
    """
}

process IPC {
    tag "$og"
    label 'small'
    publishDir "${params.outdir}/ogs/${og}/biophysics/ipc/native", mode: 'copy'

    input:
    tuple val(og), path(aa)

    output:
    tuple val(og), val('ipc'), path("${og}_ipc.tsv"), emit: pred

    script:
    """
    run_ipc.py --fasta ${aa} --out ${og}_ipc.tsv --ph ${params.ipc_ph} --pka-set ${params.ipc_pka_set}
    """
}

process SPLIT_CUSTOM_PREDICTIONS {
    label 'highmem_single'
    publishDir "${params.outdir}/custom_predictions", mode: 'copy', pattern: 'per_og/unmatched.tsv'

    input:
    path csv

    output:
    path 'per_og/*_custom.tsv', emit: tables
    path 'per_og/unmatched.tsv', optional: true

    script:
    """
    split_predictions.py --csv ${csv} --og-pattern '${params.og_pattern}' --outdir per_og
    """
}
