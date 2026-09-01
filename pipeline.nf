#!/usr/bin/env nextflow
/*
================================================================================
  BioEnvAda v2 (bea2) — per-OG biophysical + evolutionary characterisation
================================================================================
  Per ortholog group (OG):
    mafft AA alignment -> back-translated codon alignment -> occupancy pruning
    (the pruned MSA is the canonical coordinate system for EVERYTHING below)
    biophysics : b2bTools, AIUPred, IPC, custom prediction table
    evolution  : IQ-TREE (outgroup-constrained) -> rooted tree ->
                 csubst (convergence/divergence) + HyPhy (FEL/FUBAR/MEME,
                 Contrast-FEL, RELAX)
    per-OG     : native tool outputs + mapped tables + ONE combined table
                 (row = residue per sequence per pruned-MSA position) + plots
    global     : distributions & summary statistics over all OGs
================================================================================
*/

nextflow.enable.dsl = 2

include { ALIGN_AA; BACKALIGN; PRUNE_MSA } from './modules/align'
include { MAKE_FOREGROUND; MAKE_CONSTRAINT; IQTREE; ROOT_TREE } from './modules/tree'
include { CSUBST } from './modules/csubst'
include { HYPHY_WHOLETREE; HYPHY_FOREGROUND; HYPHY_BYSITE } from './modules/hyphy'
include { B2BTOOLS; AIUPRED; IPC; DSSP_RUN; DSSP_PARSE; SPLIT_CUSTOM_PREDICTIONS } from './modules/predictors'
include { MAP_TO_MSA; COMBINE_TABLE; PLOT_OG; PLOT_DSSP; GLOBAL_STATS } from './modules/collect'

log.info """
================================================================================
 BioEnvAda v2  |  outdir: ${params.outdir}
================================================================================
 INPUT
   AA fastas (--aa_dir)              : ${params.aa_dir}
   CDS fastas (--nuc_dir)            : ${params.nuc_dir}
   categories (--categories)         : ${params.categories}
   category column (--category_column): ${params.category_column}
   foreground level (--fg_level)     : ${params.fg_level}
   outgroup level (--outgroup_level) : ${params.outgroup_level}
   OG id regex (--og_pattern)        : ${params.og_pattern}
   custom predictions (--predictions_csv): ${params.predictions_csv}
 MSA
   min occupancy (--occupancy)       : ${params.occupancy}
   genetic code (--genetic_code)     : ${params.genetic_code}
 TOOLS
   b2btools / aiupred / ipc          : ${params.b2btools} / ${params.aiupred} / ${params.ipc}
   dssp (--dssp)                     : ${params.dssp}
   structures (--structure_dir)      : ${params.structure_dir}
   csubst (--csubst)                 : ${params.csubst}
   hyphy (--hyphy)                   : ${params.hyphy}
   hyphy methods (--hyphy_methods)   : ${params.hyphy_methods}
   hyphy fg methods (--hyphy_fg_methods): ${params.hyphy_fg_methods}
================================================================================
"""

workflow {
    if (!params.aa_dir || !params.nuc_dir || !params.categories)
        error "Required: --aa_dir, --nuc_dir, --categories"
    if (params.dssp && !params.structure_dir)
        error "--dssp requires --structure_dir (<structure_dir>/<og>/<sequence_id>.pdb)"

    // ---- pair AA + CDS fastas by OG id (filename before the first dot) -----
    aa_ch  = Channel.fromPath("${params.aa_dir}/*.{fa,fasta,faa}", checkIfExists: true)
                    .map { f -> tuple(f.simpleName, f) }
    nuc_ch = Channel.fromPath("${params.nuc_dir}/*.{fa,fasta,fna}", checkIfExists: true)
                    .map { f -> tuple(f.simpleName, f) }
    ogs = aa_ch.join(nuc_ch)   // (og, aa, nuc)

    // ---- categories -> foreground / outgroup definitions (once) ------------
    MAKE_FOREGROUND(Channel.fromPath(params.categories, checkIfExists: true).first())
    fg_file    = MAKE_FOREGROUND.out.foreground.first()
    fg_strains = MAKE_FOREGROUND.out.fg_strains.first()
    og_strains = MAKE_FOREGROUND.out.outgroup.first()
    cats_clean = MAKE_FOREGROUND.out.clean.first()

    // ---- alignment + canonical pruned MSA -----------------------------------
    ALIGN_AA(ogs.map { og, aa, nuc -> tuple(og, aa) })
    aln = ALIGN_AA.out.aln
    BACKALIGN(aln.join(ogs.map { og, aa, nuc -> tuple(og, nuc) }))
    PRUNE_MSA(aln.join(BACKALIGN.out.codon))
    pruned_aa    = PRUNE_MSA.out.aa
    pruned_codon = PRUNE_MSA.out.codon
    columns      = PRUNE_MSA.out.columns

    // ---- evolutionary analyses ----------------------------------------------
    hotspots_ch   = Channel.empty()   // (og, hotspots.tsv)
    csubst_tables = Channel.empty()   // (og, branch_pairs, sites, hotspots)
    bysite_ch     = Channel.empty()   // (og, by_site.tsv)
    relax_ch      = Channel.empty()   // (og, relax.tsv)

    if (params.csubst || params.hyphy) {
        MAKE_CONSTRAINT(aln, og_strains)
        IQTREE(pruned_codon.join(MAKE_CONSTRAINT.out.files))
        ROOT_TREE(IQTREE.out.treefile, og_strains)
        rooted = ROOT_TREE.out.rooted
        evo_in = pruned_codon.join(rooted)   // (og, codon, rooted)

        if (params.csubst) {
            CSUBST(evo_in, fg_file)
            csubst_tables = CSUBST.out.tables
            hotspots_ch   = CSUBST.out.hotspots
        }
        if (params.hyphy) {
            methods    = Channel.fromList(params.hyphy_methods.tokenize(',').collect { it.trim() }.findAll { it })
            fg_methods = Channel.fromList(params.hyphy_fg_methods.tokenize(',').collect { it.trim() }.findAll { it })
            HYPHY_WHOLETREE(evo_in.combine(methods))
            HYPHY_FOREGROUND(evo_in.combine(fg_methods), fg_strains)
            relax_ch = HYPHY_FOREGROUND.out.relax
            persite  = HYPHY_WHOLETREE.out.persite.mix(HYPHY_FOREGROUND.out.persite)
            HYPHY_BYSITE(persite.groupTuple().join(columns))
            bysite_ch = HYPHY_BYSITE.out.table
        }
    }

    // ---- biophysical predictors ---------------------------------------------
    // each emits (og, toolname, long-format tsv); add new predictors here
    aa_in = ogs.map { og, aa, nuc -> tuple(og, aa) }
    preds = Channel.empty()
    if (params.b2btools) { B2BTOOLS(aa_in); preds = preds.mix(B2BTOOLS.out.pred) }
    if (params.aiupred)  { AIUPRED(aa_in);  preds = preds.mix(AIUPRED.out.pred) }
    if (params.ipc)      { IPC(aa_in);      preds = preds.mix(IPC.out.pred) }
    if (params.dssp) {
        DSSP_RUN(ogs.map { og, aa, nuc -> tuple(og, file("${params.structure_dir}/${og}", checkIfExists: true)) })
        DSSP_PARSE(DSSP_RUN.out.raw)
        preds = preds.mix(DSSP_PARSE.out.pred)
    }
    if (params.predictions_csv) {
        SPLIT_CUSTOM_PREDICTIONS(Channel.fromPath(params.predictions_csv, checkIfExists: true))
        custom = SPLIT_CUSTOM_PREDICTIONS.out.tables.flatten()
                    .map { f -> tuple(f.name.replaceAll(/_custom\.tsv$/, ''), 'custom', f) }
        preds = preds.mix(custom)
    }

    // ---- map every predictor onto the pruned MSA ----------------------------
    MAP_TO_MSA(preds.combine(aln.join(columns), by: 0))
    mapped_grouped = MAP_TO_MSA.out.mapped
        .map { og, tool, f -> tuple(og, f) }
        .groupTuple()   // (og, [mapped tsvs])

    // ---- per-OG combined table + plots --------------------------------------
    combine_in = mapped_grouped
        .join(hotspots_ch, remainder: true)
        .join(bysite_ch, remainder: true)
        .join(columns)
        .filter { og, m, hs, bs, col -> m }
        .map { og, m, hs, bs, col -> tuple(og, m, hs ?: [], bs ?: [], col) }
    COMBINE_TABLE(combine_in, cats_clean)

    plot_in = mapped_grouped
        .join(pruned_aa)
        .join(hotspots_ch, remainder: true)
        .filter { og, m, aln_p, hs -> m && aln_p }
        .map { og, m, aln_p, hs -> tuple(og, m, aln_p, hs ?: []) }
    PLOT_OG(plot_in, cats_clean)

    // secondary-structure cartoon + per-sequence view (secstructartist)
    if (params.dssp) {
        ss_in = MAP_TO_MSA.out.mapped
            .filter { og, tool, f -> tool == 'dssp' }
            .map { og, tool, f -> tuple(og, f) }
            .join(hotspots_ch, remainder: true)
            .filter { og, m, hs -> m }
            .map { og, m, hs -> tuple(og, m, hs ?: []) }
        PLOT_DSSP(ss_in, cats_clean)
    }

    // ---- global statistics over all OGs --------------------------------------
    GLOBAL_STATS(
        COMBINE_TABLE.out.table.map { it[1] }.collect(),
        relax_ch.map { it[1] }.collect().ifEmpty([]),
        csubst_tables.map { it[1] }.collect().ifEmpty([])
    )
}

workflow.onComplete {
    log.info "bea2 finished: ${workflow.success ? 'OK' : 'FAILED'} | duration ${workflow.duration} | outdir: ${params.outdir}"
}
