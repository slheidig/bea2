#!/usr/bin/env python3
"""Mask in-frame stop codons (TAA/TAG/TGA) to gaps — HyPhy rejects them."""
import argparse

STOPS = {"TAA", "TAG", "TGA"}


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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    order, seqs = read_fasta(args.inp)
    with open(args.out, "w") as fh:
        for n in order:
            s = seqs[n]
            cod = [s[i:i + 3] for i in range(0, len(s) - len(s) % 3, 3)]
            fh.write(f">{n}\n{''.join('---' if c.upper() in STOPS else c for c in cod)}\n")


if __name__ == "__main__":
    main()
