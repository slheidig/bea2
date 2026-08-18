#!/usr/bin/env python3
"""Isoelectric point (IPC) + per-residue charge for one OG's sequences.

Whole-protein pI is sequence-level; it is repeated on every residue row so the
table joins cleanly with the other per-residue predictors. The per-residue
side-chain charge at --ph (default 7.0) IS positional.

pI is computed by bisection of the Henderson-Hasselbalch net charge using the
IPC pKa sets of Kozlowski 2016 (http://isoelectric.org). If the official IPC
CLI is present in the container ($IPC_CMD), it can be preferred with --use-cli.

Output: sequence_id  residue_index  residue  ipc_pI  ipc_charge
"""
import argparse
import os
import subprocess
import sys

# pKa sets: (Cterm, Nterm, {residue: (pKa, sign)})
PKA_SETS = {
    # Kozlowski L.P. 2016, IPC — Isoelectric Point Calculator, Biol Direct 11:55
    "IPC_protein": (2.869, 9.094, {
        "C": (7.555, -1), "D": (3.872, -1), "E": (4.412, -1), "Y": (10.85, -1),
        "H": (5.637, +1), "K": (9.052, +1), "R": (11.84, +1)}),
    "IPC_peptide": (2.383, 9.564, {
        "C": (8.297, -1), "D": (3.887, -1), "E": (4.317, -1), "Y": (10.071, -1),
        "H": (6.018, +1), "K": (10.517, +1), "R": (12.503, +1)}),
    "EMBOSS": (3.6, 8.6, {
        "C": (8.5, -1), "D": (3.9, -1), "E": (4.1, -1), "Y": (10.1, -1),
        "H": (6.5, +1), "K": (10.8, +1), "R": (12.5, +1)}),
}


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
    return order, {k: "".join(v).replace("-", "").replace(".", "").upper() for k, v in seqs.items()}


def net_charge(seq, ph, cterm, nterm, pka):
    q = 1.0 / (1.0 + 10 ** (ph - nterm))          # N-terminus (+)
    q -= 1.0 / (1.0 + 10 ** (cterm - ph))         # C-terminus (-)
    for aa in seq:
        if aa in pka:
            pk, sign = pka[aa]
            q += sign / (1.0 + 10 ** (sign * (ph - pk)))
    return q


def isoelectric_point(seq, cterm, nterm, pka):
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if net_charge(seq, mid, cterm, nterm, pka) > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-4:
            break
    return (lo + hi) / 2


def residue_charge(aa, ph, pka):
    if aa not in pka:
        return 0.0
    pk, sign = pka[aa]
    return sign / (1.0 + 10 ** (sign * (ph - pk)))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ph", type=float, default=7.0)
    ap.add_argument("--pka-set", default="IPC_protein", choices=sorted(PKA_SETS))
    ap.add_argument("--use-cli", action="store_true",
                    help="prefer the official IPC CLI ($IPC_CMD) for pI values")
    args = ap.parse_args()

    order, seqs = read_fasta(args.fasta)
    cterm, nterm, pka = PKA_SETS[args.pka_set]

    cli_pi = {}
    if args.use_cli and os.environ.get("IPC_CMD"):
        try:
            res = subprocess.run([*os.environ["IPC_CMD"].split(), args.fasta],
                                 capture_output=True, text=True, timeout=600)
            # IPC prints '>header ... pI: <value>' style lines; parse leniently
            cur = None
            for line in (res.stdout + res.stderr).splitlines():
                if line.startswith(">"):
                    cur = line[1:].split()[0]
                elif cur and "pI" in line:
                    for tok in line.replace(":", " ").split():
                        try:
                            cli_pi[cur] = float(tok)
                            break
                        except ValueError:
                            continue
        except Exception as e:
            print(f"WARNING: IPC CLI failed ({e}); using built-in {args.pka_set}")

    with open(args.out, "w") as fh:
        fh.write("sequence_id\tresidue_index\tresidue\tipc_pI\tipc_charge\n")
        for sid in order:
            seq = seqs[sid]
            pi = cli_pi.get(sid, isoelectric_point(seq, cterm, nterm, pka))
            for i, aa in enumerate(seq):
                fh.write(f"{sid}\t{i + 1}\t{aa}\t{pi:.3f}\t{residue_charge(aa, args.ph, pka):+.4f}\n")
    print(f"wrote {args.out}: {len(order)} sequences (pKa set {args.pka_set}, pH {args.ph})")


if __name__ == "__main__":
    main()
