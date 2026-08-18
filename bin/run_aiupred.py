#!/usr/bin/env python3
"""Run AIUPred disorder (or binding) prediction on one OG's sequences.

Uses the aiupred_lib API from https://github.com/doszilab/AIUPred (expected at
$AIUPRED_PATH, default /opt/aiupred — see docker/aiupred/Dockerfile). Runs on
CPU. Sequences are de-gapped first.

Output (standard predictor format):
  sequence_id  residue_index  residue  aiupred_<mode>
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
    return order, {k: "".join(v).replace("-", "").replace(".", "").upper() for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", default="disorder", choices=["disorder", "binding"])
    args = ap.parse_args()

    sys.path.insert(0, os.environ.get("AIUPRED_PATH", "/opt/aiupred"))
    import aiupred_lib

    # AIUPred API: init the transformer + regression networks once, then predict
    # per sequence. Fall back through the known entry-point names.
    predict_one = None
    if hasattr(aiupred_lib, "init_models") and hasattr(aiupred_lib, "predict"):
        models = aiupred_lib.init_models(args.mode) if _accepts_arg(aiupred_lib.init_models) \
            else aiupred_lib.init_models()
        models = models if isinstance(models, tuple) else (models,)
        predict_one = lambda seq: aiupred_lib.predict(seq, *models)
    elif hasattr(aiupred_lib, "aiupred_disorder"):
        predict_one = aiupred_lib.aiupred_disorder
    if predict_one is None:
        sys.exit("ERROR: unrecognised aiupred_lib API — check the AIUPred version in the container")

    order, seqs = read_fasta(args.fasta)
    col = f"aiupred_{args.mode}"
    with open(args.out, "w") as fh:
        fh.write(f"sequence_id\tresidue_index\tresidue\t{col}\n")
        for sid in order:
            seq = seqs[sid]
            scores = predict_one(seq)
            scores = list(getattr(scores, "tolist", lambda: scores)())
            if len(scores) != len(seq):
                print(f"WARNING: {sid}: {len(scores)} scores for {len(seq)} residues")
            for i, aa in enumerate(seq):
                v = scores[i] if i < len(scores) else float("nan")
                fh.write(f"{sid}\t{i + 1}\t{aa}\t{float(v):.6g}\n")
    print(f"wrote {args.out}: {len(order)} sequences ({col})")


def _accepts_arg(fn):
    try:
        import inspect
        return len(inspect.signature(fn).parameters) >= 1
    except (TypeError, ValueError):
        return False


if __name__ == "__main__":
    main()
