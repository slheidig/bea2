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

include { ALIGN_AA; PRUNE } from './modules/align'
include { MAKE_FOREGROUND; TREE } from './modules/tree'
include { CSUBST } from './modules/csubst'
include { HYPHY_WHOLETREE; HYPHY_FOREGROUND } from './modules/hyphy'
include { B2BTOOLS; AIUPRED; IPC; DSSP; SPLIT_CUSTOM_PREDICTIONS } from './modules/predictors'
include { MAP_TO_MSA; COLLECT; PLOT_OG; PLOT_DSSP; GLOBAL_STATS } from './modules/collect'

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
   csubst reuses IQ-TREE (--csubst_reuse_iqtree): ${params.csubst_reuse_iqtree}
   hyphy (--hyphy)                   : ${params.hyphy}
   hyphy methods (--hyphy_methods)   : ${params.hyphy_methods}
   hyphy fg methods (--hyphy_fg_methods): ${params.hyphy_fg_methods}
 OUTPUT
   plots (--plot)                    : ${params.plot}
   publish csubst native (--publish_csubst_native): ${params.publish_csubst_native}
================================================================================
"""

workflow {
    if (!params.aa_dir || !params.nuc_dir || !params.categories)
        error "Required: --aa_dir, --nuc_dir, --categories"
    if (params.dssp && !params.structure_dir)
        error "--dssp requires --structure_dir (<structure_dir>/<og>/<sequence_id>.pdb)"
    // csubst reconstructs ancestral states under its own default model; reusing
    // our IQ-TREE run is only valid when we ran the same model.
    if (params.csubst_reuse_iqtree && params.iqtree_model != 'ECMK07+F+R4')
        error "--csubst_reuse_iqtree requires --iqtree_model 'ECMK07+F+R4' (csubst's default), got '${params.iqtree_model}'"

    // ---- --plot all | none | <og>[,<og>...] --------------------------------
    def plot_mode = params.plot.toString().trim()
    def plot_all  = plot_mode.equalsIgnoreCase('all')
    def plot_none = plot_mode.equalsIgnoreCase('none')
    def plot_set  = (plot_all || plot_none) ? [] : plot_mode.tokenize(',').collect { it.trim() }.findAll { it }
    if (!plot_all && !plot_none && !plot_set)
        error "--plot must be 'all', 'none' or a comma-separated list of OG ids"

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
    PRUNE(aln.join(ogs.map { og, aa, nuc -> tuple(og, nuc) }))
    pruned_aa    = PRUNE.out.aa
    pruned_codon = PRUNE.out.codon
    columns      = PRUNE.out.columns

    // ---- evolutionary analyses ----------------------------------------------
    hotspots_ch   = Channel.empty()   // (og, hotspots.tsv)
    csubst_tables = Channel.empty()   // (og, branch_pairs, sites, hotspots)
    persite_ch    = Channel.empty()   // (og, persite.tsv)
    relax_ch      = Channel.empty()   // (og, relax.tsv)

    if (params.csubst || params.hyphy) {
        TREE(aln.join(pruned_codon), og_strains)
        rooted = TREE.out.rooted
        evo_in = pruned_codon.join(rooted)   // (og, codon, rooted)

        if (params.csubst) {
            CSUBST(evo_in.join(TREE.out.asr), fg_file)
            csubst_tables = CSUBST.out.tables
            hotspots_ch   = CSUBST.out.hotspots
        }
        if (params.hyphy) {
            methods    = Channel.fromList(params.hyphy_methods.tokenize(',').collect { it.trim() }.findAll { it })
            fg_methods = Channel.fromList(params.hyphy_fg_methods.tokenize(',').collect { it.trim() }.findAll { it })
            HYPHY_WHOLETREE(evo_in.combine(methods))
            HYPHY_FOREGROUND(evo_in.combine(fg_methods), fg_strains)
            relax_ch   = HYPHY_FOREGROUND.out.relax
            persite_ch = HYPHY_WHOLETREE.out.persite.mix(HYPHY_FOREGROUND.out.persite)
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
        DSSP(ogs.map { og, aa, nuc -> tuple(og, file("${params.structure_dir}/${og}", checkIfExists: true)) })
        preds = preds.mix(DSSP.out.pred)
    }
    if (params.predictions_csv) {
        SPLIT_CUSTOM_PREDICTIONS(Channel.fromPath(params.predictions_csv, checkIfExists: true))
        custom = SPLIT_CUSTOM_PREDICTIONS.out.tables.flatten()
                    .map { f -> tuple(f.name.replaceAll(/_custom\.tsv$/, ''), 'custom', f) }
        preds = preds.mix(custom)
    }

    // ---- map every predictor of an OG onto the pruned MSA in one task -------
    // size: lets each OG start as soon as its own predictors are done instead of
    // waiting for all of them; remainder: keeps OGs whose predictors were ignored
    n_pred = [params.b2btools, params.aiupred, params.ipc, params.dssp,
              params.predictions_csv as boolean].count { it }
    MAP_TO_MSA(preds.groupTuple(size: n_pred, remainder: true).join(aln.join(columns)))
    mapped_grouped = MAP_TO_MSA.out.mapped

    // ---- per-OG combined table (with the hyphy by-site table folded in) ----
    collect_in = mapped_grouped
        .join(columns)
        .join(hotspots_ch,            remainder: true)
        .join(persite_ch.groupTuple(), remainder: true)
        .filter { it[1] }
        .map { og, m, col, hs, ps -> tuple(og, m, col, hs ?: [], ps ?: []) }
    COLLECT(collect_in, cats_clean)

    // ---- per-OG plots (subset with --plot) ----------------------------------
    plot_in = mapped_grouped
        .join(pruned_aa)
        .join(hotspots_ch, remainder: true)
        .filter { it[1] && it[2] && (plot_all || (!plot_none && it[0] in plot_set)) }
        .map { og, m, aln_p, hs -> tuple(og, m, aln_p, hs ?: []) }
    PLOT_OG(plot_in, cats_clean)

    // secondary-structure cartoon + per-sequence view (secstructartist)
    if (params.dssp) {
        ss_in = MAP_TO_MSA.out.dssp
            .join(hotspots_ch, remainder: true)
            .filter { it[1] && (plot_all || (!plot_none && it[0] in plot_set)) }
            .map { og, m, hs -> tuple(og, m, hs ?: []) }
        PLOT_DSSP(ss_in, cats_clean)
    }

    // ---- global statistics over all OGs --------------------------------------
    // concatenated driver-side: GLOBAL_STATS gets three paths regardless of OG count
    combined_all = COLLECT.out.table.map { it[1] }
        .collectFile(name: 'combined_all.tsv', keepHeader: true, skip: 1, sort: true)
    relax_all = relax_ch.map { it[1] }
        .collectFile(name: 'relax_all.tsv', keepHeader: true, skip: 1, sort: true)
        .ifEmpty(file("${projectDir}/assets/NO_RELAX"))
    pairs_all = csubst_tables.map { it[1] }
        .collectFile(name: 'branch_pairs_all.tsv', keepHeader: true, skip: 1, sort: true)
        .ifEmpty(file("${projectDir}/assets/NO_PAIRS"))

    GLOBAL_STATS(combined_all, relax_all, pairs_all)
}

workflow.onComplete {
    log.info "bea2 finished: ${workflow.success ? 'OK' : 'FAILED'} | duration ${workflow.duration} | outdir: ${params.outdir}"
}
