#!/usr/bin/env python3
"""Run AIUPred disorder (or binding) prediction on one OG's sequences.

Drives the `aiupred` command-line tool shipped in the authors' official CPU
image, ghcr.io/doszilab/aiupred:cpu (AIUPred 3.1.2, package `aiupred`). Earlier
AIUPred releases exposed a flat `aiupred_lib` module with a Python API; that
module no longer exists, so this wrapper goes through the CLI instead — the
interface the authors version and test.

Sequences are de-gapped before prediction, so residue_index is 1-based over the
ungapped sequence.

Output (standard predictor format):
  sequence_id  residue_index  residue  aiupred_<mode>
"""
import argparse
import os
import subprocess
import sys
import tempfile


def read_fasta(path):
    seqs, order, name = {}, [], None
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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", default="disorder", choices=["disorder", "binding"])
    ap.add_argument("--aiupred-bin", default=os.environ.get("AIUPRED_BIN", "aiupred"),
                    help="AIUPred executable (default: aiupred, from the container's PATH)")
    args = ap.parse_args()

    order, seqs = read_fasta(args.fasta)
    if not order:
        sys.exit(f"ERROR: no sequences read from {args.fasta}")

    with tempfile.TemporaryDirectory(dir=".") as tmp:
        # Feed AIUPred the de-gapped sequences under normalised ids, so the ids
        # it echoes back line up with what the rest of the pipeline expects.
        in_fa = os.path.join(tmp, "degapped.fasta")
        with open(in_fa, "w") as fh:
            for sid in order:
                fh.write(f">{sid}\n{seqs[sid]}\n")

        raw = os.path.join(tmp, "aiupred.tsv")
        # --force-cpu: the cluster nodes have no GPU, and without it AIUPred
        # defaults to GPU index 0 and dies. -b switches to AIUPred-binding.
        cmd = [args.aiupred_bin, "-i", in_fa, "-o", raw, "--force-cpu"]
        if args.mode == "binding":
            cmd.append("-b")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            sys.stderr.write(proc.stdout + proc.stderr)
            sys.exit(f"ERROR: {' '.join(cmd)} failed with exit {proc.returncode}")

        scores = parse_aiupred(raw, args.mode)

    col = f"aiupred_{args.mode}"
    with open(args.out, "w") as fh:
        fh.write(f"sequence_id\tresidue_index\tresidue\t{col}\n")
        for sid in order:
            seq = seqs[sid]
            vals = scores.get(sid)
            if vals is None:
                print(f"WARNING: {sid}: no AIUPred output, writing NaN")
                vals = []
            if len(vals) != len(seq):
                print(f"WARNING: {sid}: {len(vals)} scores for {len(seq)} residues")
            for i, aa in enumerate(seq):
                v = vals[i] if i < len(vals) else float("nan")
                fh.write(f"{sid}\t{i + 1}\t{aa}\t{float(v):.6g}\n")
    print(f"wrote {args.out}: {len(order)} sequences ({col})")


def parse_aiupred(path, mode):
    """Read AIUPred's TSV into {sequence_id: [score, ...]}.

    The file opens with a '#' banner, then one '>id' line per sequence followed
    by per-residue rows: position, residue, then one or more score columns.
    With -b AIUPred appends the binding score to the disorder one, so take the
    last numeric column in binding mode and the first otherwise.
    """
    out, current = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                out[current] = []
                continue
            if current is None:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            nums = []
            for p in parts[2:]:
                try:
                    nums.append(float(p))
                except ValueError:
                    pass
            if not nums:
                continue
            out[current].append(nums[-1] if mode == "binding" and len(nums) > 1 else nums[0])
    return out


if __name__ == "__main__":
    main()
