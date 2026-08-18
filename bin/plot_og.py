#!/usr/bin/env python3
"""Per-OG plots: pruned-MSA heatmap + per-feature line plots by category.

Line plots (one PDF per feature): per-category mean with min/max band across
the pruned-MSA positions, overlaid with csubst convergent/divergent hotspot
sites. Category colours are assigned dynamically from the values present in
the chosen categories column (outgroup level always grey, drawn first).

Inputs: pruned AA alignment, categories_clean.tsv, mapped predictor TSVs,
optional <og>.hotspots.tsv.
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

STANDARD_AA = list("ACDEFGHIKLMNPQRSTVWY")
PALETTE = ["blue", "darkred", "darkorange", "green", "purple", "teal", "olive", "magenta"]
HOTSPOT_STYLE = {
    "convergent": dict(color="darkorange", linestyle="-", linewidth=1.0, alpha=0.8),
    "divergent": dict(color="purple", linestyle=":", linewidth=1.0, alpha=0.8),
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
    return order, {k: "".join(v) for k, v in seqs.items()}


def strain_of(sid, og):
    rest = sid[len(og) + 1:] if sid.startswith(og + "_") else sid
    return rest.rsplit("_", 1)[0] if "_" in rest else rest


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--og", required=True)
    ap.add_argument("--aln", required=True, help="pruned AA alignment")
    ap.add_argument("--categories", required=True)
    ap.add_argument("--column", default="temp_cat2")
    ap.add_argument("--outgroup-level", default="Outgroup")
    ap.add_argument("--mapped", nargs="*", default=[])
    ap.add_argument("--hotspots", default=None)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    og = args.og
    os.makedirs(args.outdir, exist_ok=True)

    cats = pd.read_csv(args.categories, sep="\t", dtype=str)
    col = cats.columns[int(args.column) - 1] if args.column.isdigit() else args.column
    cat_map = dict(zip(cats.iloc[:, 0], cats[col]))

    order, seqs = read_fasta(args.aln)
    seq_cat = {sid: cat_map.get(strain_of(sid, og), "unknown") for sid in order}

    levels = [l for l in pd.unique(cats[col].dropna()) if l]
    # outgroup grey & drawn first; everything else from the palette
    colours, i = {}, 0
    for l in levels:
        if l == args.outgroup_level:
            colours[l] = "grey"
        else:
            colours[l] = PALETTE[i % len(PALETTE)]
            i += 1
    colours["unknown"] = "black"
    draw_order = ([args.outgroup_level] if args.outgroup_level in levels else []) + \
                 [l for l in levels if l != args.outgroup_level] + ["unknown"]

    hotspots = {}
    if args.hotspots and os.path.exists(args.hotspots):
        hs = pd.read_csv(args.hotspots, sep="\t")
        for signal in ("convergent", "divergent"):
            hotspots[signal] = sorted(hs.loc[hs["signal"] == signal, "codon_site"].astype(int))

    # ---- MSA heatmap -------------------------------------------------------
    L = len(next(iter(seqs.values())))
    tab20 = plt.cm.tab20(np.linspace(0, 1, 20))
    aa2c = {aa: tab20[i] for i, aa in enumerate(STANDARD_AA)}
    img = np.ones((len(order), L, 4))
    for r, sid in enumerate(order):
        for c, ch in enumerate(seqs[sid]):
            img[r, c] = aa2c.get(ch.upper(), (0.7, 0.7, 0.7, 1) if ch not in "-." else (1, 1, 1, 1))
    fig, ax = plt.subplots(figsize=(max(6, L / 25), max(3, len(order) / 8)))
    ax.imshow(img, aspect="auto", interpolation="nearest", extent=(0.5, L + 0.5, len(order) - 0.5, -0.5))
    for r, sid in enumerate(order):
        ax.plot([0.2], [r], marker="s", ms=3, color=colours.get(seq_cat[sid], "black"), clip_on=False)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=4)
    ax.set_xlabel("pruned-MSA position")
    ax.set_title(f"{og}: pruned MSA ({len(order)} seqs x {L} positions)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, f"{og}_msa_heatmap.pdf"))
    plt.close(fig)

    # ---- per-feature line plots -------------------------------------------
    for f in args.mapped:
        df = pd.read_csv(f, sep="\t")
        df["category"] = df["sequence_id"].map(seq_cat)
        feat_cols = [c for c in df.columns if c not in
                     ("sequence_id", "msa_position", "original_msa_position",
                      "residue_index", "residue", "category")]
        for feat in feat_cols:
            vals = pd.to_numeric(df[feat], errors="coerce")
            if vals.notna().sum() == 0:
                continue
            sub = df[["msa_position", "category"]].copy()
            sub["v"] = vals
            fig, ax = plt.subplots(figsize=(max(6, L / 25), 3.2))
            for catlev in draw_order:
                g = sub[sub["category"] == catlev]
                if g.empty:
                    continue
                agg = g.groupby("msa_position")["v"].agg(["mean", "min", "max"])
                c = colours.get(catlev, "black")
                ax.fill_between(agg.index, agg["min"], agg["max"], color=c, alpha=0.15, lw=0)
                ax.plot(agg.index, agg["mean"], color=c, lw=1.0, label=catlev)
            for signal, sites in hotspots.items():
                for s in sites:
                    ax.axvline(s, **HOTSPOT_STYLE[signal])
            ax.set_xlabel("pruned-MSA position")
            ax.set_ylabel(feat)
            ax.set_title(f"{og}: {feat} by {col}", fontsize=10)
            handles = [Line2D([0], [0], color=colours[l], lw=1.5, label=l)
                       for l in draw_order if l in sub["category"].values]
            handles += [Line2D([0], [0], label=s, **HOTSPOT_STYLE[s]) for s in hotspots if hotspots[s]]
            ax.legend(handles=handles, fontsize=6, ncol=3, framealpha=0.8)
            fig.tight_layout()
            fig.savefig(os.path.join(args.outdir, f"{og}_{feat}_plot.pdf"))
            plt.close(fig)
    print(f"{og}: plots written to {args.outdir}/")


if __name__ == "__main__":
    main()
