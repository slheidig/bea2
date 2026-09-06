#!/usr/bin/env python3
"""Select convergent foreground branch pairs; prints 'branch_id_1,branch_id_2'.

Cutoffs are csubst's own, read back from csubst_cb_stats.tsv. Non-finite
statistics (omegaC = inf on a near-zero denominator) do not pass.
"""
import argparse
import csv


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def read_cutoffs(stats_tsv, arity):
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
    ap.add_argument("--cb", required=True, help="csubst_cb_2.tsv")
    ap.add_argument("--stats", required=True, help="csubst_cb_stats.tsv")
    ap.add_argument("--arity", type=int, default=2)
    args = ap.parse_args()

    cutoffs = read_cutoffs(args.stats, args.arity)
    if not cutoffs:
        return

    for r in csv.DictReader(open(args.cb), delimiter="\t"):
        if r.get("is_fg") == "Y" and passes(r, cutoffs):
            print(f"{r['branch_id_1']},{r['branch_id_2']}")


if __name__ == "__main__":
    main()
