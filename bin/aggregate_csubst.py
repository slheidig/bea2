#!/usr/bin/env python3
"""Aggregate one OG's foreground csubst convergence into three tidy tables.

Convergence only; divergence (any2dif) stays in the native output. Reported in
three flavours: any2spe (main statistic), spe2spe (parallelism) and dif2spe
(convergence sensu stricto), the last derived per site as any2spe - spe2spe.

Cutoffs are csubst's own: --cutoff_stat for pairs, --min-combinat-prob for sites.

Outputs (tab-separated), codon_site == pruned-MSA position:
  <og>.branch_pairs.tsv  convergent pairs, ranked by OCNany2spe
  <og>.sites.tsv         per-site convergent events
  <og>.hotspots.tsv      sites ranked by number of recurring strains
"""
import argparse
import csv
import glob
import os
from collections import defaultdict

PAIR_STATS = ["OCNany2spe", "OCNspe2spe", "OCNdif2spe", "ECNany2spe", "omegaCany2spe"]


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def z(x):
    return 0.0 if x != x else x


def fmt(x, nd=3):
    return "" if x != x else round(x, nd)


def read_cutoffs(stats_tsv, arity):
    if not os.path.exists(stats_tsv):
        return []
    for r in csv.DictReader(open(stats_tsv), delimiter="\t"):
        if int(float(r["arity"])) != arity:
            continue
        spec = r.get("cutoff_stat", "").strip()
        if spec:
            return [(s.split(",")[0], float(s.split(",")[1])) for s in spec.split("|")]
    return []


def passes(row, cutoffs):
    vals = [num(row.get(stat)) for stat, _ in cutoffs]
    if any(v != v for v in vals):
        return False
    return all(v >= cut for v, (_, cut) in zip(vals, cutoffs))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--og", required=True)
    ap.add_argument("--search-dir", required=True, help="csubst search output dir (with csubst_cb_2.tsv, csubst_b.tsv, sites/)")
    ap.add_argument("--min-combinat-prob", type=float, default=0.5)
    args = ap.parse_args()
    og = args.og

    def strain_of(bn):
        s = bn[len(og) + 1:] if bn.startswith(og + "_") else bn
        return s.rsplit("_", 1)[0] if "_" in s else s

    # branch_id -> strain
    bmap = {}
    b_tsv = os.path.join(args.search_dir, "csubst_b.tsv")
    if os.path.exists(b_tsv):
        for row in csv.DictReader(open(b_tsv), delimiter="\t"):
            bmap[row["branch_id"]] = strain_of(row["branch_name"])

    cutoffs = read_cutoffs(os.path.join(args.search_dir, "csubst_cb_stats.tsv"), 2)
    cb = os.path.join(args.search_dir, "csubst_cb_2.tsv")
    cb_rows = list(csv.DictReader(open(cb), delimiter="\t")) if os.path.exists(cb) else []

    with open(f"{og}.branch_pairs.tsv", "w") as bp:
        bp.write("og\tstrain_1\tstrain_2\t" + "\t".join(PAIR_STATS) + "\n")
        out = []
        for r in cb_rows:
            if r.get("is_fg") != "Y" or not passes(r, cutoffs):
                continue
            b1, b2 = r["branch_id_1"], r["branch_id_2"]
            out.append([bmap.get(b1, b1), bmap.get(b2, b2)]
                       + [fmt(num(r.get(s)), 2 if s == "omegaCany2spe" else 3) for s in PAIR_STATS])
        out.sort(key=lambda x: -z(num(x[2])))
        for x in out:
            bp.write(og + "\t" + "\t".join(map(str, x)) + "\n")

    with open(f"{og}.sites.tsv", "w") as sf, open(f"{og}.hotspots.tsv", "w") as hf:
        sf.write("og\tstrain_1\tstrain_2\tcodon_site\tsub_strain1\tsub_strain2"
                 "\tOCNany2spe\tOCNspe2spe\tOCNdif2spe\n")
        hf.write("og\tcodon_site\tn_strains\tstrains\n")
        hot = defaultdict(set)
        for dd in sorted(glob.glob(os.path.join(args.search_dir, "sites", "csubst.branch_id*"))):
            if not os.path.isdir(dd):
                continue
            b1, b2 = os.path.basename(dd).split("branch_id")[1].split(",")
            tsv = os.path.join(dd, "csubst.tsv")
            if not os.path.exists(tsv):
                continue
            for r in csv.DictReader(open(tsv), delimiter="\t"):
                any2spe = z(num(r.get("OCNany2spe")))
                if any2spe < args.min_combinat_prob:
                    continue
                spe2spe = z(num(r.get("OCNspe2spe")))
                s1, s2 = bmap.get(b1, b1), bmap.get(b2, b2)
                site = r["codon_site_alignment"]
                sf.write(f"{og}\t{s1}\t{s2}\t{site}\t"
                         f"{r.get('aa_' + b1 + '_anc', '?')}->{r.get('aa_' + b1, '?')}\t"
                         f"{r.get('aa_' + b2 + '_anc', '?')}->{r.get('aa_' + b2, '?')}\t"
                         f"{round(any2spe, 3)}\t{round(spe2spe, 3)}\t{round(any2spe - spe2spe, 3)}\n")
                hot[site] |= {s1, s2}
        for site, strains in sorted(hot.items(), key=lambda kv: -len(kv[1])):
            hf.write(f"{og}\t{site}\t{len(strains)}\t{','.join(sorted(strains))}\n")
    print(f"{og}: wrote branch_pairs / sites / hotspots tables")


if __name__ == "__main__":
    main()
