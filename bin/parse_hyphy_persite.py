#!/usr/bin/env python3
"""Parse a HyPhy per-site MLE JSON (FEL / FUBAR / MEME / Contrast-FEL) to TSV.

Adds a dN_dS column (where alpha & beta exist) and a 'selection' call:
  p-value methods : diversifying / purifying at p <= --pvalue
  FUBAR           : posterior Prob[alpha<beta] / Prob[alpha>beta] >= --posterior
"""
import argparse
import json
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--json", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--og", required=True)
    ap.add_argument("--method", required=True)
    ap.add_argument("--pvalue", type=float, default=0.1)
    ap.add_argument("--posterior", type=float, default=0.9)
    args = ap.parse_args()

    try:
        j = json.load(open(args.json))
    except Exception as e:
        open(args.out, "w").write(f"# no MLE table parsed: {e}\n")
        print(f"  !! parse failed: {e}")
        sys.exit(0)

    mle = j.get("MLE", {})
    headers = [h[0] for h in mle.get("headers", [])]
    rows = mle.get("content", {}).get("0", [])
    # some MLE tables (e.g. FUBAR) list fewer header names than row columns
    width = max((len(r) for r in rows), default=len(headers))
    while len(headers) < width:
        headers.append(f"col{len(headers) + 1}")

    def idx(name):
        for i, h in enumerate(headers):
            if h.lower() == name.lower():
                return i
        return None

    ai, bi, pi = idx("alpha"), idx("beta"), idx("p-value")
    agi, abi = idx("prob[alpha>beta]"), idx("prob[alpha<beta]")
    have_ab = ai is not None and bi is not None
    mode = ("pval" if (have_ab and pi is not None)
            else ("post" if (agi is not None and abi is not None) else None))
    extra = (["dN_dS"] if have_ab else []) + (["selection"] if mode else [])

    with open(args.out, "w") as fh:
        fh.write("og\tcodon_site\t" + "\t".join(headers + extra) + "\n")
        for i, row in enumerate(rows, 1):
            cells = [args.og, str(i)] + [f"{v:.6g}" if isinstance(v, (int, float)) else str(v) for v in row]
            if have_ab:
                a, b = row[ai], row[bi]
                dnds = (b / a) if a > 0 else (float("inf") if b > 0 else 0.0)
                cells.append("inf" if dnds == float("inf") else f"{dnds:.4f}")
            if mode == "pval":
                a, b, p = row[ai], row[bi], row[pi]
                cells.append("diversifying" if (b > a and p <= args.pvalue)
                             else ("purifying" if (b < a and p <= args.pvalue) else "neutral"))
            elif mode == "post":
                pag, pab = row[agi], row[abi]
                cells.append("diversifying" if pab >= args.posterior
                             else ("purifying" if pag >= args.posterior else "neutral"))
            fh.write("\t".join(cells) + "\n")
    print(f"  wrote {args.out}: {len(rows)} sites")


if __name__ == "__main__":
    main()
