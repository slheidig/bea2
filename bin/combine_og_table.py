#!/usr/bin/env python3
"""Combine all mapped predictor tables + evolutionary results into ONE per-OG table.

One row per (sequence, pruned-MSA position): identity columns, category, every
biophysical predictor, and the position-level evolutionary annotations
(identical for all sequences at a position):
  csubst : convergent_n_strains / divergent_n_strains (hotspot recurrence)
  hyphy  : every column of the by-site table (per-method stats)

Output: <og>_combined.tsv
"""
import argparse
import sys

import pandas as pd


def strain_of(sid, og):
    rest = sid[len(og) + 1:] if sid.startswith(og + "_") else sid
    return rest.rsplit("_", 1)[0] if "_" in rest else rest


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--og", required=True)
    ap.add_argument("--categories", required=True, help="categories_clean.tsv")
    ap.add_argument("--column", default="temp_cat2")
    ap.add_argument("--mapped", nargs="+", required=True, help="*_mapped.tsv predictor tables")
    ap.add_argument("--hotspots", default=None, help="<og>.hotspots.tsv (csubst)")
    ap.add_argument("--bysite", default=None, help="<og>.by_site.tsv (hyphy)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    og = args.og

    ID = ["sequence_id", "msa_position", "original_msa_position", "residue_index", "residue"]
    KEY = ["sequence_id", "msa_position"]

    base = None
    for f in args.mapped:
        df = pd.read_csv(f, sep="\t")
        if base is None:
            base = df
            continue
        aux = [c for c in ID if c not in KEY and c in df.columns]
        base = base.merge(df.drop(columns=aux), on=KEY, how="outer")
    if base is None or base.empty:
        sys.exit(f"ERROR: {og}: no mapped predictor rows")
    base = base.sort_values(KEY).reset_index(drop=True)

    # strain + category
    cats = pd.read_csv(args.categories, sep="\t", dtype=str)
    col = args.column
    if col.isdigit():
        col = cats.columns[int(col) - 1]
    cat_map = dict(zip(cats.iloc[:, 0], cats[col]))
    base.insert(1, "strain", base["sequence_id"].map(lambda s: strain_of(s, og)))
    base.insert(2, "category", base["strain"].map(cat_map))
    base.insert(0, "og", og)

    # csubst hotspot recurrence per position
    if args.hotspots:
        try:
            hs = pd.read_csv(args.hotspots, sep="\t")
            for signal in ("convergent", "divergent"):
                sub = hs[hs["signal"] == signal][["codon_site", "n_strains"]]
                m = dict(zip(sub["codon_site"].astype(int), sub["n_strains"]))
                base[f"csubst_{signal}_n_strains"] = (
                    base["msa_position"].map(m).fillna(0).astype(int))
        except Exception as e:
            print(f"WARNING: {og}: hotspots not joined ({e})")

    # hyphy by-site columns per position
    if args.bysite:
        try:
            bs = pd.read_csv(args.bysite, sep="\t")
            bs = bs.drop(columns=[c for c in ("og", "original_msa_position") if c in bs.columns])
            bs = bs.rename(columns={"codon_site": "msa_position"})
            bs.columns = ["msa_position"] + [f"hyphy_{c}" for c in bs.columns[1:]]
            base = base.merge(bs, on="msa_position", how="left")
        except Exception as e:
            print(f"WARNING: {og}: hyphy by-site not joined ({e})")

    base.to_csv(args.out, sep="\t", index=False)
    print(f"{og}: combined table {base.shape[0]} rows x {base.shape[1]} cols -> {args.out}")


if __name__ == "__main__":
    main()
