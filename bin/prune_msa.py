#!/usr/bin/env python3
"""Prune MSA columns below a minimum occupancy; apply the same columns to the codon MSA.

The pruned amino-acid MSA is the CANONICAL coordinate system of the pipeline:
pruned position i (1-based) == csubst codon_site i == HyPhy site i.

Outputs:
  <og>.aa.pruned.fa       pruned amino-acid alignment
  <og>.codon.pruned.fa    pruned codon alignment (same columns x3)
  <og>.msa_columns.tsv    pruned_position <-> original_position + occupancy
"""
import argparse


def read_fasta(path):
    seqs, name, order = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
                order.append(name)
            elif line and name:
                seqs[name].append(line.strip())
    return order, {k: "".join(v) for k, v in seqs.items()}


def write_fasta(path, order, seqs, width=80):
    with open(path, "w") as fh:
        for n in order:
            fh.write(f">{n}\n")
            s = seqs[n]
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--og", required=True)
    ap.add_argument("--aa", required=True, help="full amino-acid alignment")
    ap.add_argument("--codon", required=True, help="full in-frame codon alignment (backaligned)")
    ap.add_argument("--occupancy", type=float, default=0.5,
                    help="min fraction of non-gap residues per column (default 0.5)")
    args = ap.parse_args()

    order, aa = read_fasta(args.aa)
    corder, codon = read_fasta(args.codon)

    L = len(next(iter(aa.values())))
    assert all(len(s) == L for s in aa.values()), "AA alignment is not flush"
    n = len(order)

    # occupancy per column (gap characters: - and .)
    occ = []
    for j in range(L):
        filled = sum(1 for name in order if aa[name][j] not in "-.")
        occ.append(filled / n)
    kept = [j for j in range(L) if occ[j] >= args.occupancy]
    if not kept:
        raise SystemExit(f"ERROR: {args.og}: no columns pass occupancy >= {args.occupancy}")

    aa_pruned = {name: "".join(aa[name][j] for j in kept) for name in order}

    # codon alignment: same columns, x3. Codon MSA may lack sequences dropped by
    # backalign; keep the intersection but warn.
    missing = [x for x in order if x not in codon]
    if missing:
        print(f"WARNING: {args.og}: {len(missing)} sequences absent from codon alignment: {missing[:5]}")
    corder2 = [x for x in order if x in codon]
    Lc = len(next(iter(codon.values())))
    assert Lc == 3 * L, f"codon alignment length {Lc} != 3 x AA alignment length {L}"
    codon_pruned = {name: "".join(codon[name][3 * j:3 * j + 3] for j in kept) for name in corder2}

    write_fasta(f"{args.og}.aa.pruned.fa", order, aa_pruned)
    write_fasta(f"{args.og}.codon.pruned.fa", corder2, codon_pruned)
    with open(f"{args.og}.msa_columns.tsv", "w") as fh:
        fh.write("pruned_position\toriginal_position\toccupancy\n")
        for i, j in enumerate(kept, 1):
            fh.write(f"{i}\t{j + 1}\t{occ[j]:.4f}\n")
    print(f"{args.og}: kept {len(kept)}/{L} columns (occupancy >= {args.occupancy})")


if __name__ == "__main__":
    main()
