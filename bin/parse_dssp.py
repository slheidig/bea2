#!/usr/bin/env python3
"""Turn one OG's DSSP output into the pipeline's standard predictor format:
  sequence_id  residue_index(1-based, ungapped)  residue  <feature columns>

Every sequence in this pipeline has a predicted structure covering exactly its
fasta sequence, so the DSSP residue number is used directly as residue_index.

Inputs: a directory of <sequence_id>.dssp files and the per-residue pLDDT table
extracted from the model B-factors (sequence_id, residue_index, plddt).
"""
import argparse
import csv
import os

# theoretical maximum solvent accessibility, Tien et al. 2013 (PLoS ONE 8:e80635)
MAX_ASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
           "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
           "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174}

# Kyte & Doolittle 1982 hydropathy
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "E": -3.5, "Q": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

# 8-state -> 3-state: helix / sheet / loop. T, S and unassigned all become loop,
# the same three buckets SIMSApiper uses.
SS3 = {"H": "H", "G": "H", "I": "H", "E": "E", "B": "E"}

COLS = ["dssp_ss8", "dssp_ss3", "dssp_acc", "dssp_rsa", "dssp_kd",
        "dssp_surface_hydrophobicity", "dssp_plddt"]


def parse_dssp(path):
    """(residue_index, residue, ss8, acc) per residue of one .dssp file."""
    rows = []
    body = False
    with open(path) as fh:
        for line in fh:
            if not body:
                body = line.startswith("  #  RESIDUE")
                continue
            if len(line) < 38 or line[13] == "!":   # chain break
                continue
            aa = line[13]
            aa = "C" if aa.islower() else aa        # lowercase = disulfide cysteine
            ss = line[16]
            # unassigned is blank in the file; 'X' as in SIMSApiper, so it can be
            # confused with an alignment gap ('-'); it becomes ss3 loop below
            rows.append((int(line[5:10]), aa, "X" if ss == " " else ss,
                         int(line[34:38])))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dssp-dir", required=True)
    ap.add_argument("--plddt", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    plddt = {}
    for r in csv.DictReader(open(args.plddt), delimiter="\t"):
        plddt[(r["sequence_id"], int(r["residue_index"]))] = r["plddt"]

    files = sorted(f for f in os.listdir(args.dssp_dir) if f.endswith(".dssp"))
    n_res = 0
    with open(args.out, "w") as out:
        out.write("sequence_id\tresidue_index\tresidue\t" + "\t".join(COLS) + "\n")
        for f in files:
            sid = f[: -len(".dssp")]
            for idx, aa, ss8, acc in parse_dssp(os.path.join(args.dssp_dir, f)):
                rsa = acc / MAX_ASA[aa] if aa in MAX_ASA else None
                kd = KD.get(aa)
                surf = rsa * kd if rsa is not None and kd is not None else None
                out.write("\t".join([
                    sid, str(idx), aa, ss8, SS3.get(ss8, "L"), str(acc),
                    "" if rsa is None else f"{rsa:.4f}",
                    "" if kd is None else f"{kd:g}",
                    "" if surf is None else f"{surf:.4f}",
                    plddt.get((sid, idx), ""),
                ]) + "\n")
                n_res += 1
    print(f"wrote {args.out}: {len(files)} structures, {n_res} residues")


if __name__ == "__main__":
    main()
