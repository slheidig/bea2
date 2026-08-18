#!/usr/bin/env python3
"""Select foreground branch pairs significant for convergence OR divergence.

Reads csubst_cb_2.tsv; prints 'branch_id_1,branch_id_2' lines (stdout) for
is_fg == Y pairs with OCNany2{spe,dif} >= --ocn and omegaCany2{spe,dif} > --omega
(omega NaN counts as significant: near-zero denominator == inf).
"""
import argparse
import csv


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cb", required=True, help="csubst_cb_2.tsv")
    ap.add_argument("--ocn", type=float, default=0.5)
    ap.add_argument("--omega", type=float, default=1.0)
    args = ap.parse_args()

    def sig(oc, wc):
        return oc >= args.ocn and (wc > args.omega or wc != wc)

    for r in csv.DictReader(open(args.cb), delimiter="\t"):
        if r.get("is_fg", "") != "Y":
            continue
        conv = sig(num(r.get("OCNany2spe")), num(r.get("omegaCany2spe")))
        div = sig(num(r.get("OCNany2dif")), num(r.get("omegaCany2dif")))
        if conv or div:
            print(f"{r['branch_id_1']},{r['branch_id_2']}")


if __name__ == "__main__":
    main()
