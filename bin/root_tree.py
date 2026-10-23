#!/usr/bin/env python3
"""Root a gene tree on the outgroup clade (ete4 two-step), with midpoint fallback.

- strips UFBoot support labels (csubst requires unique/absent internal labels)
- >=2 outgroup leaves: root inside the ingroup first, then on the outgroup MRCA
  (gives a clean bifurcating root)
- 1 outgroup leaf: root on it directly
- 0 outgroup leaves: midpoint root (fallback: longest terminal branch)
"""
import argparse
import io
import re

from ete4 import Tree


def strain_of(label, og):
    rest = label[len(og) + 1:] if label.startswith(og + "_") else label
    return rest.rsplit("_", 1)[0] if "_" in rest else rest


def plot(t, og_leaves, og, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from Bio import Phylo

    tree = Phylo.read(io.StringIO(t.write(parser=0)), "newick")
    fig, ax = plt.subplots(figsize=(10, max(4, 0.22 * len(tree.get_terminals()))))
    Phylo.draw(tree, axes=ax, do_show=False,
               branch_labels=lambda c: f"{c.confidence:g}" if c.confidence is not None else "",
               label_colors=lambda name: "grey" if name in og_leaves else "black")
    ax.set_title(og)
    fig.tight_layout()
    fig.savefig(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--og", required=True)
    ap.add_argument("--tree", required=True, help="IQ-TREE .treefile")
    ap.add_argument("--outgroup", required=True, help="outgroup_strains.txt")
    ap.add_argument("--out", required=True)
    ap.add_argument("--plot", default=None, help="also draw the rooted tree with supports")
    args = ap.parse_args()

    ogset = {s.strip() for s in open(args.outgroup) if s.strip()}
    t = Tree(open(args.tree).read())
    t.del_prop("support")

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

    if args.plot:
        plot(t, set(og_leaves), args.og, args.plot)
    t.write(outfile=args.out, parser=5)
    # belt-and-braces: remove any residual support-style internal labels
    txt = re.sub(r"\)[0-9.]+:", "):", open(args.out).read())
    open(args.out, "w").write(txt)
    print(f"{args.og}: rooted on {len(og_leaves)} outgroup leaves -> {args.out}")


if __name__ == "__main__":
    main()
