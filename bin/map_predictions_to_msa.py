#!/usr/bin/env python3
"""Project a per-residue prediction table onto the pruned-MSA coordinate system.

Uses the FULL alignment's gap pattern to map each sequence's ungapped
residue_index to an original MSA column, keeps only the pruned columns
(msa_columns.tsv), and renumbers to pruned positions 1..N (== csubst
codon_site == HyPhy site).

Input predictor TSV (standard format): sequence_id, residue_index (1-based
ungapped), residue, <feature columns>.

Output: sequence_id, msa_position, original_msa_position, residue_index,
residue, <tool>_<feature>...  (features of tool 'custom' keep their names).
"""
import argparse
import csv


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif line and name:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--og", required=True)
    ap.add_argument("--tool", required=True)
    ap.add_argument("--aln", required=True, help="FULL amino-acid alignment")
    ap.add_argument("--columns", required=True, help="msa_columns.tsv")
    ap.add_argument("--pred", required=True, help="predictor long TSV")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    aln = read_fasta(args.aln)
    # original position (1-based) -> pruned position (1-based)
    kept = {}
    for r in csv.DictReader(open(args.columns), delimiter="\t"):
        kept[int(r["original_position"])] = int(r["pruned_position"])

    # per sequence: ungapped residue_index (1-based) -> original column (1-based)
    res2col = {}
    for sid, s in aln.items():
        res2col[sid] = [i + 1 for i, ch in enumerate(s) if ch not in "-."]

    with open(args.pred) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        H = {h: i for i, h in enumerate(header)}
        for req in ("sequence_id", "residue_index"):
            if req not in H:
                raise SystemExit(f"ERROR: {args.pred} lacks required column '{req}'")
        feat_cols = [h for h in header if h not in ("sequence_id", "residue_index", "residue")]
        rename = {h: (h if args.tool == "custom" or h.startswith(args.tool + "_")
                      else f"{args.tool}_{h}") for h in feat_cols}

        n_out = n_skip = 0
        mismatches = 0
        with open(args.out, "w") as out:
            out.write("sequence_id\tmsa_position\toriginal_msa_position\tresidue_index\tresidue\t"
                      + "\t".join(rename[h] for h in feat_cols) + "\n")
            for row in reader:
                sid = row[H["sequence_id"]]
                if sid not in res2col:
                    n_skip += 1
                    continue
                ridx = int(row[H["residue_index"]])
                cols = res2col[sid]
                if not (1 <= ridx <= len(cols)):
                    n_skip += 1
                    continue
                orig = cols[ridx - 1]
                pruned = kept.get(orig)
                if pruned is None:      # column dropped by occupancy pruning
                    continue
                res = row[H["residue"]] if "residue" in H else aln[sid][orig - 1]
                if "residue" in H and aln[sid][orig - 1].upper() != res.upper():
                    mismatches += 1
                out.write(f"{sid}\t{pruned}\t{orig}\t{ridx}\t{res}\t"
                          + "\t".join(row[H[h]] for h in feat_cols) + "\n")
                n_out += 1
        if mismatches:
            print(f"WARNING: {args.og}/{args.tool}: {mismatches} residue mismatches vs alignment")
        if n_skip:
            print(f"NOTE: {args.og}/{args.tool}: {n_skip} rows skipped (sequence/index not in alignment)")
        print(f"{args.og}/{args.tool}: {n_out} mapped rows -> {args.out}")


if __name__ == "__main__":
    main()
