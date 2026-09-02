// HyPhy selection scans, reprocessed from selection_scan.sh.
// Whole-tree per-site methods (FEL/FUBAR/MEME/SLAC) and foreground methods
// (Contrast-FEL, RELAX) run as separate parallel tasks per (OG, method).
// Sites == pruned-MSA positions == csubst codon_site (canonical coordinates).

process HYPHY_WHOLETREE {
    tag "$og:$method"
    publishDir "${params.outdir}/ogs/${og}/evolution/hyphy/native", mode: 'copy',
        pattern: '*.{json,log}'
    publishDir "${params.outdir}/ogs/${og}/evolution/hyphy", mode: 'copy',
        pattern: '*.persite.tsv'

    input:
    tuple val(og), path(codon), path(rooted), val(method)

    output:
    tuple val(og), path("${og}.${method}.persite.tsv"), emit: persite
    path "${og}.${method}.{json,log}"

    script:
    def m = method.toLowerCase()
    def branches = (m == 'fubar') ? '' : "--branches ${params.hyphy_branches}"
    """
    mask_stops.py --in ${codon} --out nostop.fa
    hyphy CPU=${task.cpus} ${m} --alignment nostop.fa --tree ${rooted} \\
        --code ${params.hyphy_code} ${branches} --output ${og}.${method}.json \\
        < /dev/null > ${og}.${method}.log 2>&1
    parse_hyphy_persite.py --json ${og}.${method}.json --out ${og}.${method}.persite.tsv \\
        --og ${og} --method ${method} --pvalue ${params.pvalue} --posterior ${params.posterior}
    """
}

process HYPHY_FOREGROUND {
    tag "$og:$method"
    publishDir "${params.outdir}/ogs/${og}/evolution/hyphy/native", mode: 'copy',
        pattern: '*.{json,log,nwk}'
    publishDir "${params.outdir}/ogs/${og}/evolution/hyphy", mode: 'copy',
        pattern: '*.{persite.tsv,relax.tsv}'

    input:
    tuple val(og), path(codon), path(rooted), val(method)
    path fg_strains

    output:
    tuple val(og), path("${og}.${method}.fg.persite.tsv"), emit: persite, optional: true
    tuple val(og), path("${og}.relax.tsv"),                emit: relax,   optional: true
    path "${og}.${method}.fg.{json,log}"
    path "${og}.foreground.nwk"

    script:
    def m = method.toLowerCase()
    if (m == 'relax')
        """
        mask_stops.py --in ${codon} --out nostop.fa
        tag_foreground.py --og ${og} --tree ${rooted} --fg-strains ${fg_strains} --out ${og}.foreground.nwk
        # --models Minimal keeps the K=1 vs K-free test, skips the slow descriptive model
        hyphy CPU=${task.cpus} relax --alignment nostop.fa --tree ${og}.foreground.nwk \\
            --code ${params.hyphy_code} --test Foreground --models Minimal \\
            --output ${og}.${method}.fg.json < /dev/null > ${og}.${method}.fg.log 2>&1
        parse_relax.py --json ${og}.${method}.fg.json --og ${og} --out ${og}.relax.tsv
        """
    else if (m == 'contrast-fel')
        """
        mask_stops.py --in ${codon} --out nostop.fa
        tag_foreground.py --og ${og} --tree ${rooted} --fg-strains ${fg_strains} --out ${og}.foreground.nwk
        # contrast-fel's branch menu: with one custom label, 5 = Foreground, d = done
        printf '5\\nd\\n' | hyphy CPU=${task.cpus} contrast-fel --alignment nostop.fa \\
            --tree ${og}.foreground.nwk --code ${params.hyphy_code} \\
            --output ${og}.${method}.fg.json > ${og}.${method}.fg.log 2>&1
        parse_hyphy_persite.py --json ${og}.${method}.fg.json --out ${og}.${method}.fg.persite.tsv \\
            --og ${og} --method ${method} --pvalue ${params.pvalue} --posterior ${params.posterior}
        """
    else
        """
        mask_stops.py --in ${codon} --out nostop.fa
        tag_foreground.py --og ${og} --tree ${rooted} --fg-strains ${fg_strains} --out ${og}.foreground.nwk
        hyphy CPU=${task.cpus} ${m} --alignment nostop.fa --tree ${og}.foreground.nwk \\
            --code ${params.hyphy_code} --branches Foreground \\
            --output ${og}.${method}.fg.json < /dev/null > ${og}.${method}.fg.log 2>&1
        parse_hyphy_persite.py --json ${og}.${method}.fg.json --out ${og}.${method}.fg.persite.tsv \\
            --og ${og} --method ${method} --pvalue ${params.pvalue} --posterior ${params.posterior}
        """
}

// The per-site tables are folded into COLLECT (modules/collect.nf) so that the
// by-site step does not cost its own scheduler round-trip.
