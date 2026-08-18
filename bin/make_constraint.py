#!/usr/bin/env python3
"""Per-OG IQ-TREE outgroup files: monophyly constraint tree + outgroup leaf labels.

Not every OG contains every outgroup strain; a partial set is used as-is. If an
OG contains none of the outgroup strains, both outputs are empty files and the
tree is built unconstrained (midpoint-rooted downstream).

Outputs:
  <og>.outgroup_labels.txt  comma-joined leaf labels of present outgroup strains (may be empty)
  <og>.constraint.nwk       ((outgroup...),ingroup...); or empty
"""
import argparse


def read_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0].strip())
    return ids


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--og", required=True)
    ap.add_argument("--aln", required=True, help="alignment (leaf label source)")
    ap.add_argument("--outgroup", required=True, help="outgroup_strains.txt (one strain per line)")
    args = ap.parse_args()

    strains = [s.strip() for s in open(args.outgroup) if s.strip()]
    ids = read_ids(args.aln)

    labels = []
    for s in strains:
        hit = next((i for i in ids if f"_{s}_" in i), None)
        if hit:
            labels.append(hit)
    ingroup = [i for i in ids if i not in labels]

    with open(f"{args.og}.outgroup_labels.txt", "w") as fh:
        fh.write(",".join(labels))
    with open(f"{args.og}.constraint.nwk", "w") as fh:
        # a monophyly constraint only makes sense for >=2 outgroup taxa;
        # with 1 taxon, rooting via -o alone is sufficient (empty constraint)
        if len(labels) >= 2 and ingroup:
            fh.write(f"(({','.join(labels)}),{','.join(ingroup)});\n")
    print(f"{args.og}: {len(labels)}/{len(strains)} outgroup strains present")


if __name__ == "__main__":
    main()
