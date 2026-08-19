#!/usr/bin/env python3
"""Run AIUPred disorder (or binding) prediction on one OG's sequences.

Uses the aiupred_lib API from https://github.com/doszilab/AIUPred. Runs inside
the authors' official CPU image, ghcr.io/doszilab/aiupred:cpu, where the code
lives under /opt/aiupred; $AIUPRED_PATH overrides that root. Runs on CPU.
Sequences are de-gapped first.

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

    _add_aiupred_to_syspath(os.environ.get("AIUPRED_PATH", "/opt/aiupred"))
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


def _add_aiupred_to_syspath(root):
    """Put the directory that actually contains aiupred_lib.py on sys.path.

    The AIUPred distribution does not keep aiupred_lib.py at its top level in
    every layout (a plain clone root does not import), so rather than assume,
    walk the tree once and use whatever directory holds it. Falls back to the
    root itself, which is correct when the module is already importable.
    """
    sys.path.insert(0, root)
    if os.path.isfile(os.path.join(root, "aiupred_lib.py")):
        return
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in {".git", "__pycache__"}]
        if "aiupred_lib.py" in filenames:
            sys.path.insert(0, dirpath)
            return


def _accepts_arg(fn):
    try:
        import inspect
        return len(inspect.signature(fn).parameters) >= 1
    except (TypeError, ValueError):
        return False


if __name__ == "__main__":
    main()
