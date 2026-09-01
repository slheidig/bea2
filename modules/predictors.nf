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

// DSSP is the one predictor that reads structures instead of sequences:
// one predicted model per sequence, in ${params.structure_dir}/<og>/<sequence_id>.pdb.
process DSSP_RUN {
    tag "$og"
    label 'small'
    publishDir "${params.outdir}/ogs/${og}/biophysics/dssp/native", mode: 'copy'

    input:
    tuple val(og), path(structdir)

    output:
    tuple val(og), path('dssp'), path("${og}_plddt.tsv"), emit: raw

    script:
    """
    mkdir dssp
    printf 'sequence_id\\tresidue_index\\tplddt\\n' > ${og}_plddt.tsv
    for pdb in ${structdir}/*.pdb ; do
        id=\$(basename \$pdb .pdb)
        mkdssp -i \$pdb -o dssp/\$id.dssp
        # pLDDT = the CA B-factor of each residue (cols 61-66 of the ATOM record)
        awk -v id=\$id 'substr(\$0,1,4)=="ATOM" && substr(\$0,13,4)==" CA " { print id "\\t" substr(\$0,23,4)+0 "\\t" substr(\$0,61,6)+0 }' \$pdb >> ${og}_plddt.tsv
    done
    """
}

process DSSP_PARSE {
    tag "$og"
    label 'small'
    publishDir "${params.outdir}/ogs/${og}/biophysics/dssp/native", mode: 'copy'

    input:
    tuple val(og), path(dsspdir), path(plddt)

    output:
    tuple val(og), val('dssp'), path("${og}_dssp.tsv"), emit: pred

    script:
    """
    parse_dssp.py --dssp-dir ${dsspdir} --plddt ${plddt} --out ${og}_dssp.tsv
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
