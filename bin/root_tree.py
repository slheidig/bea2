#!/usr/bin/env python3
"""Root a gene tree on the outgroup clade (ete4 two-step), with midpoint fallback.

- strips UFBoot support labels (csubst requires unique/absent internal labels)
- >=2 outgroup leaves: root inside the ingroup first, then on the outgroup MRCA
  (gives a clean bifurcating root)
- 1 outgroup leaf: root on it directly
- 0 outgroup leaves: midpoint root (fallback: longest terminal branch)
"""
import argparse
import re

from ete4 import Tree


def strain_of(label, og):
    rest = label[len(og) + 1:] if label.startswith(og + "_") else label
    return rest.rsplit("_", 1)[0] if "_" in rest else rest


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--og", required=True)
    ap.add_argument("--tree", required=True, help="IQ-TREE .treefile")
    ap.add_argument("--outgroup", required=True, help="outgroup_strains.txt")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ogset = {s.strip() for s in open(args.outgroup) if s.strip()}
    nwk = re.sub(r"\)[0-9.]+:", "):", open(args.tree).read())  # strip support labels
    t = Tree(nwk)
    for n in t.traverse():
        if "support" in n.props:
            n.del_prop("support")

    og_leaves = [l.name for l in t.leaves() if strain_of(l.name, args.og) in ogset]
    in_leaves = [l.name for l in t.leaves() if strain_of(l.name, args.og) not in ogset]

    if not og_leaves:
        print(f"{args.og}: no outgroup leaves -> midpoint root")
        try:
            t.set_outgroup(t.get_midpoint_outgroup())
        except Exception:
            far = max(t.leaves(), key=lambda l: l.dist or 0)
            t.set_outgroup(t[far.name])
    elif len(og_leaves) == 1:
        t.set_outgroup(t[og_leaves[0]])
    else:
        t.set_outgroup(t[in_leaves[0]])            # step 1: root inside the ingroup
        t.set_outgroup(t.common_ancestor(og_leaves))  # step 2: root on outgroup clade

    t.write(outfile=args.out)
    # belt-and-braces: remove any residual support-style internal labels
    txt = re.sub(r"\)[0-9.]+:", "):", open(args.out).read())
    open(args.out, "w").write(txt)
    print(f"{args.og}: rooted on {len(og_leaves)} outgroup leaves -> {args.out}")


if __name__ == "__main__":
    main()
