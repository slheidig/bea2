#!/usr/bin/env python3
"""Tag foreground tips with {Foreground} in a Newick tree (for HyPhy)."""
import argparse
import re


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--og", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--fg-strains", required=True, help="one strain name per line")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fg = {s.strip() for s in open(args.fg_strains) if s.strip()}
    txt = open(args.tree).read()

    def strain_of(lbl):
        rest = lbl[len(args.og) + 1:] if lbl.startswith(args.og + "_") else lbl
        return rest.rsplit("_", 1)[0] if "_" in rest else rest

    # leaf labels = tokens immediately followed by ':' (branch length)
    labels = set(re.findall(r"[(),]([^(),:]+):", txt))
    n = 0
    for lbl in labels:
        if strain_of(lbl) in fg:
            txt = txt.replace(lbl + ":", lbl + "{Foreground}:")
            n += 1
    open(args.out, "w").write(txt)
    print(f"{args.og}: tagged {n} foreground tips")


if __name__ == "__main__":
    main()
