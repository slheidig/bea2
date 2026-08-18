#!/usr/bin/env python3
"""Run b2bTools (DynaMine/EFoldMine/DisoMine/AgMata) on one OG's sequences.

Sequences are de-gapped before prediction. Output is a long TSV in the
pipeline's standard predictor format:
  sequence_id  residue_index(1-based, ungapped)  residue  <feature columns>
"""
import argparse
import os
import sys


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
    return order, {k: "".join(v).replace("-", "").replace(".", "") for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tools", default="dynamine,efoldmine,disomine,agmata")
    args = ap.parse_args()

    order, seqs = read_fasta(args.fasta)
    degapped = "degapped.fasta"
    with open(degapped, "w") as fh:
        for n in order:
            fh.write(f">{n}\n{seqs[n]}\n")

    from b2bTools import SingleSeq
    try:
        from b2bTools import constants
        toolmap = {"dynamine": getattr(constants, "TOOL_DYNAMINE", "dynamine"),
                   "efoldmine": getattr(constants, "TOOL_EFOLDMINE", "efoldmine"),
                   "disomine": getattr(constants, "TOOL_DISOMINE", "disomine"),
                   "agmata": getattr(constants, "TOOL_AGMATA", "agmata")}
    except ImportError:
        toolmap = {t: t for t in ("dynamine", "efoldmine", "disomine", "agmata")}

    tools = [toolmap[t.strip().lower()] for t in args.tools.split(",") if t.strip()]
    ss = SingleSeq(os.path.abspath(degapped))
    ss.predict(tools=tools)
    preds = ss.get_all_predictions()
    if isinstance(preds, dict) and "proteins" in preds:
        preds = preds["proteins"]

    # collect feature names (skip non-list entries like 'seq')
    def clean(vals):
        out = []
        for v in vals:
            if isinstance(v, (list, tuple)):   # (residue, value) pairs in some versions
                v = v[-1]
            try:
                out.append(float(v))
            except (TypeError, ValueError):
                out.append(float("nan"))
        return out

    feats = set()
    table = {}
    for sid, d in preds.items():
        if not isinstance(d, dict):
            continue
        table[sid] = {}
        for k, v in d.items():
            if k in ("seq", "sequence") or not isinstance(v, (list, tuple)):
                continue
            table[sid][k] = clean(v)
            feats.add(k)
    feats = sorted(feats)
    if not feats:
        sys.exit("ERROR: b2bTools returned no per-residue features")

    with open(args.out, "w") as fh:
        fh.write("sequence_id\tresidue_index\tresidue\t" + "\t".join(feats) + "\n")
        for sid in order:
            d = table.get(sid, {})
            seq = seqs[sid]
            for i, aa in enumerate(seq):
                row = [sid, str(i + 1), aa]
                for f in feats:
                    vals = d.get(f)
                    row.append(f"{vals[i]:.6g}" if vals and i < len(vals) and vals[i] == vals[i] else "")
                fh.write("\t".join(row) + "\n")
    print(f"wrote {args.out}: {len(order)} sequences, features: {', '.join(feats)}")


if __name__ == "__main__":
    main()
