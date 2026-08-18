#!/usr/bin/env python3
"""Parse a HyPhy RELAX JSON into a one-row gene-level summary TSV."""
import argparse
import json


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", required=True)
    ap.add_argument("--og", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fmt = lambda x: f"{x:.4g}" if isinstance(x, (int, float)) else "NA"
    K = p = LRT = None
    verdict = "parse_failed"
    try:
        j = json.load(open(args.json))
        tr = j.get("test results", {})

        def g1(*keys):
            for k in keys:
                for kk, v in tr.items():
                    if kk.lower() == k.lower():
                        return v
            return None

        K = g1("relaxation or intensification parameter",
               "relaxation or intensification test statistic")
        p = g1("p-value")
        LRT = g1("LRT")
        verdict = "NA"
        if isinstance(K, (int, float)):
            verdict = "intensified(K>1)" if K > 1 else ("relaxed(K<1)" if K < 1 else "none")
    except Exception:
        pass
    with open(args.out, "w") as fh:
        fh.write("og\tK\tLRT\tp_value\tverdict\n")
        fh.write(f"{args.og}\t{fmt(K)}\t{fmt(LRT)}\t{fmt(p)}\t{verdict}\n")
    print(f"{args.og}: RELAX K={fmt(K)} p={fmt(p)} -> {verdict}")


if __name__ == "__main__":
    main()
