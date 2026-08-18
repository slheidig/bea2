#!/usr/bin/env python3
"""Split the big custom-predictions CSV into per-OG TSVs (streaming).

Input columns: block_id,residue_index,aa,pred_*  where block_id is
<OG>_<strain>_<locus>. The OG id is recognised with --og-pattern, a regex
anchored at the START of block_id (e.g. 'CK_\\d+' or 'OG\\d+').
residue_index is 1-based over the ungapped sequence.

Output: <outdir>/<og>_custom.tsv in the standard predictor format
(sequence_id, residue_index, residue, pred_*). Rows whose block_id does not
match the pattern are collected in <outdir>/unmatched.tsv.
"""
import argparse
import os
import re

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--og-pattern", default=r"CK_\d+",
                    help="regex matching the OG id at the start of block_id")
    ap.add_argument("--outdir", default="per_og")
    ap.add_argument("--chunksize", type=int, default=2_000_000)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    pat = re.compile(f"^({args.og_pattern})")
    seen = set()
    n_rows = n_unmatched = 0
    for chunk in pd.read_csv(args.csv, chunksize=args.chunksize):
        chunk = chunk.rename(columns={"block_id": "sequence_id", "aa": "residue"})
        chunk["_og"] = chunk["sequence_id"].str.extract(pat, expand=False)
        bad = chunk[chunk["_og"].isna()]
        if len(bad):
            path = os.path.join(args.outdir, "unmatched.tsv")
            bad.drop(columns=["_og"]).to_csv(path, sep="\t", index=False,
                                             mode="a", header=not os.path.exists(path))
            n_unmatched += len(bad)
            chunk = chunk.dropna(subset=["_og"])
        for og, sub in chunk.groupby("_og", sort=False):
            path = os.path.join(args.outdir, f"{og}_custom.tsv")
            new = og not in seen
            sub.drop(columns=["_og"]).to_csv(path, sep="\t", index=False,
                                             mode="w" if new else "a", header=new)
            seen.add(og)
        n_rows += len(chunk)
        print(f"  ... {n_rows:,} rows, {len(seen)} OGs, {n_unmatched:,} unmatched")
    print(f"done: {n_rows:,} rows split into {len(seen)} per-OG tables in {args.outdir}/ "
          f"({n_unmatched:,} unmatched rows)")


if __name__ == "__main__":
    main()
