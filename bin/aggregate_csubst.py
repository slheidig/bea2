#!/usr/bin/env python3
"""Aggregate one OG's csubst outputs into three tidy tables.

Both signals are reported in the same style:
  spe = convergence (both branches substitute to the SAME amino acid)
  dif = divergence  (same site, substitutions to DIFFERENT amino acids)

Outputs (tab-separated, 'signal' column = convergent|divergent):
  <og>.branch_pairs.tsv  significant foreground pairs, ranked by OCN
  <og>.sites.tsv         per-site events from csubst sites
  <og>.hotspots.tsv      sites ranked by number of recurring foreground strains

Site numbering: codon_site == pruned-MSA position (the pipeline's canonical
coordinate system).
"""
import argparse
import csv
import glob
import os
from collections import defaultdict

MODES = {"spe": "convergent", "dif": "divergent"}


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def z(x):
    return 0.0 if x != x else x


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--og", required=True)
    ap.add_argument("--search-dir", required=True, help="csubst search output dir (with csubst_cb_2.tsv, csubst_b.tsv, sites/)")
    ap.add_argument("--ocn", type=float, default=0.5)
    ap.add_argument("--omega", type=float, default=1.0)
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

    cb = os.path.join(args.search_dir, "csubst_cb_2.tsv")
    cb_rows = list(csv.DictReader(open(cb), delimiter="\t")) if os.path.exists(cb) else []

    with open(f"{og}.branch_pairs.tsv", "w") as bp:
        bp.write("og\tsignal\tstrain_1\tstrain_2\tOCN\tECN\tomegaC\n")
        out = []
        for suf, label in MODES.items():
            for r in cb_rows:
                if r.get("is_fg", "") != "Y":
                    continue
                oc, wc = num(r.get(f"OCNany2{suf}")), num(r.get(f"omegaCany2{suf}"))
                if oc >= args.ocn and (wc > args.omega or wc != wc):
                    out.append((label, bmap.get(r["branch_id_1"], r["branch_id_1"]),
                                bmap.get(r["branch_id_2"], r["branch_id_2"]),
                                round(oc, 3), round(z(num(r.get(f"ECNany2{suf}"))), 3),
                                "inf" if wc != wc else round(wc, 2)))
        out.sort(key=lambda x: (x[0], -x[3]))
        for x in out:
            bp.write(og + "\t" + "\t".join(map(str, x)) + "\n")

    with open(f"{og}.sites.tsv", "w") as sf, open(f"{og}.hotspots.tsv", "w") as hf:
        sf.write("og\tsignal\tstrain_1\tstrain_2\tcodon_site\tsub_strain1\tsub_strain2\tOCN\n")
        hf.write("og\tsignal\tcodon_site\tn_strains\tstrains\n")
        for suf, label in MODES.items():
            hot = defaultdict(set)
            for dd in sorted(glob.glob(os.path.join(args.search_dir, "sites", "csubst.branch_id*"))):
                if not os.path.isdir(dd):
                    continue
                b1, b2 = os.path.basename(dd).split("branch_id")[1].split(",")
                tsv = os.path.join(dd, "csubst.tsv")
                if not os.path.exists(tsv):
                    continue
                for r in csv.DictReader(open(tsv), delimiter="\t"):
                    if z(num(r.get(f"OCNany2{suf}", 0))) >= args.ocn:
                        s1, s2 = bmap.get(b1, b1), bmap.get(b2, b2)
                        site = r["codon_site_alignment"]
                        sf.write(f"{og}\t{label}\t{s1}\t{s2}\t{site}\t"
                                 f"{r.get('aa_' + b1 + '_anc', '?')}->{r.get('aa_' + b1, '?')}\t"
                                 f"{r.get('aa_' + b2 + '_anc', '?')}->{r.get('aa_' + b2, '?')}\t"
                                 f"{round(z(num(r.get(f'OCNany2{suf}'))), 3)}\n")
                        hot[site] |= {s1, s2}
            for site, strains in sorted(hot.items(), key=lambda kv: -len(kv[1])):
                hf.write(f"{og}\t{label}\t{site}\t{len(strains)}\t{','.join(sorted(strains))}\n")
    print(f"{og}: wrote branch_pairs / sites / hotspots tables")


if __name__ == "__main__":
    main()
