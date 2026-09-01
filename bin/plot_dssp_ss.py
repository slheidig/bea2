#!/usr/bin/env python3
"""Secondary-structure plots for one OG, from the mapped DSSP table.

Adapted from SIMSApiper's 2Dstructure_plot.py / DSSPcodesMSA_plot.py /
dssp_seqview_plot.py, but on the pruned-MSA coordinate system and split by
category, which is what bea2 compares.

  <og>_ss_consensus.pdf   secstructartist cartoon of the consensus secondary
                          structure per category, over H/E/C frequency panels,
                          with csubst hotspots overlaid
  <og>_ss_alignment.pdf   per-sequence H/E/C, sequences grouped by category
  <og>_ss_consensus.tsv   the numbers behind the cartoon
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import secstructartist as ssa

PALETTE = ["blue", "darkred", "darkorange", "green", "purple", "teal", "olive", "magenta"]
HOTSPOT_STYLE = {
    "convergent": dict(color="darkorange", linestyle="-", linewidth=1.0, alpha=0.8),
    "divergent": dict(color="purple", linestyle=":", linewidth=1.0, alpha=0.8),
}
SS3 = ["H", "E", "L"]
SS_NAME = {"H": "Helix", "E": "Sheet", "L": "Loop"}   # loop = DSSP T, S, unassigned
SS_COLOUR = {"H": "green", "E": "blue", "L": "darkorange"}
SS_ARTIST = {"H": "H", "E": "S", "L": "L"}     # dssp_ss3 -> secstructartist code


def strain_of(sid, og):
    rest = sid[len(og) + 1:] if sid.startswith(og + "_") else sid
    return rest.rsplit("_", 1)[0] if "_" in rest else rest


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--og", required=True)
    ap.add_argument("--mapped", required=True, help="<og>_dssp_mapped.tsv")
    ap.add_argument("--categories", required=True)
    ap.add_argument("--column", default="temp_cat2")
    ap.add_argument("--outgroup-level", default="Outgroup")
    ap.add_argument("--hotspots", default=None)
    ap.add_argument("--consensus", type=float, default=0.5,
                    help="min frequency for a position to get a consensus H/E")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    og = args.og
    os.makedirs(args.outdir, exist_ok=True)
    out = lambda n: os.path.join(args.outdir, n)

    df = pd.read_csv(args.mapped, sep="\t")
    if "dssp_ss3" not in df.columns:
        raise SystemExit(f"ERROR: {args.mapped} has no dssp_ss3 column")

    cats = pd.read_csv(args.categories, sep="\t", dtype=str)
    col = cats.columns[int(args.column) - 1] if args.column.isdigit() else args.column
    cat_map = dict(zip(cats.iloc[:, 0], cats[col]))
    df["category"] = df["sequence_id"].map(lambda s: cat_map.get(strain_of(s, og), "unknown"))

    # outgroup grey & drawn first; everything else from the palette (as in plot_og.py)
    levels = [l for l in pd.unique(cats[col].dropna()) if l]
    colours, i = {}, 0
    for l in levels:
        if l == args.outgroup_level:
            colours[l] = "grey"
        else:
            colours[l] = PALETTE[i % len(PALETTE)]
            i += 1
    colours["unknown"] = "black"
    order = ([args.outgroup_level] if args.outgroup_level in levels else []) + \
            [l for l in levels if l != args.outgroup_level] + ["unknown"]
    present = [l for l in order if (df["category"] == l).any()]
    if not present:
        raise SystemExit(f"ERROR: {og}: no sequence matched a category")

    L = int(df["msa_position"].max())
    x = np.arange(1, L + 1)
    mat = (df.pivot_table(index="sequence_id", columns="msa_position",
                          values="dssp_ss3", aggfunc="first")
             .reindex(columns=x))

    hotspots = {}
    if args.hotspots and os.path.exists(args.hotspots):
        hs = pd.read_csv(args.hotspots, sep="\t")
        for signal in ("convergent", "divergent"):
            hotspots[signal] = sorted(hs.loc[hs["signal"] == signal, "codon_site"].astype(int))

    # ---- per-category frequency + consensus ---------------------------------
    seq_cat = df.drop_duplicates("sequence_id").set_index("sequence_id")["category"]
    freq, consensus, records = {}, {}, []
    for lev in present:
        sub = mat.loc[seq_cat[seq_cat == lev].index]
        n = sub.notna().sum(axis=0)                       # residues present per position
        f = {s: (sub == s).sum(axis=0) / n.replace(0, np.nan) for s in SS3}
        freq[lev] = f
        best = pd.DataFrame(f).idxmax(axis=1)
        top = pd.DataFrame(f).max(axis=1)
        consensus[lev] = best.where(top >= args.consensus, "L").fillna("L")
        for p in x:
            records.append([og, lev, p, int(n[p]), f["H"][p], f["E"][p], f["L"][p],
                            consensus[lev][p]])
    pd.DataFrame(records, columns=["og", "category", "msa_position", "n_residues",
                                   "freq_helix", "freq_sheet", "freq_loop",
                                   "consensus"]).to_csv(out(f"{og}_ss_consensus.tsv"),
                                                        sep="\t", index=False, na_rep="")

    fig_width = min(max(L / 25, 12), 80)

    # ---- figure 1: consensus cartoon + frequency panels ---------------------
    fig, axes = plt.subplots(4, 1, figsize=(fig_width, 3 + 1.2 * len(present)),
                             sharex=True,
                             gridspec_kw={"height_ratios": [0.6 * len(present), 1, 1, 1]})
    ax = axes[0]
    artist = ssa.SecStructArtist()
    artist.height = 0.7
    artist["H"].fillcolor = SS_COLOUR["H"]
    artist["S"].fillcolor = SS_COLOUR["E"]
    for row, lev in enumerate(present):
        artist.draw("".join(SS_ARTIST[s] for s in consensus[lev]), x, ypos=row, ax=ax)
    ax.set_ylim(len(present) - 0.5, -0.5)
    ax.set_yticks(range(len(present)))
    ax.set_yticklabels(present, fontsize=7)
    for tick, lev in zip(ax.get_yticklabels(), present):
        tick.set_color(colours.get(lev, "black"))
    ax.set_title(f"{og}: consensus secondary structure "
                 f"(>= {args.consensus:g} of residues present), by category", fontsize=10)

    for panel, s in zip(axes[1:], SS3):
        for lev in present:
            panel.plot(x, freq[lev][s].values, linewidth=1,
                       color=colours.get(lev, "black"), label=lev)
        panel.set_ylim(-0.05, 1.05)
        panel.set_ylabel(f"{SS_NAME[s]}\nfrequency", fontsize=8)
        panel.tick_params(labelsize=7)
    for signal, sites in hotspots.items():
        for panel in axes[1:]:
            for site in sites:
                panel.axvline(site, **HOTSPOT_STYLE[signal])
    handles = [plt.Line2D([], [], color=colours.get(l, "black"), label=l) for l in present] + \
              [plt.Line2D([], [], label=f"{s} site", **HOTSPOT_STYLE[s]) for s in hotspots]
    axes[1].legend(handles=handles, fontsize=6, ncol=len(handles), loc="upper right")
    axes[-1].set_xlabel("pruned-MSA position")
    axes[-1].set_xlim(0.5, L + 0.5)
    fig.tight_layout()
    fig.savefig(out(f"{og}_ss_consensus.pdf"))
    plt.close(fig)

    # ---- figure 2: per-sequence secondary structure, grouped by category ----
    rows = [sid for lev in present for sid in sorted(seq_cat[seq_cat == lev].index)]
    code = {s: i for i, s in enumerate(SS3)}
    img = mat.loc[rows].map(lambda v: code.get(v, len(SS3))).to_numpy()
    cmap = ListedColormap([SS_COLOUR[s] for s in SS3] + ["lightgrey"])

    fig, ax = plt.subplots(figsize=(fig_width, min(max(len(rows) / 8, 3), 80)))
    ax.imshow(img, aspect="auto", interpolation="nearest", cmap=cmap,
              vmin=0, vmax=len(SS3), extent=(0.5, L + 0.5, len(rows) - 0.5, -0.5))
    for r, sid in enumerate(rows):
        ax.plot([0.2], [r], marker="s", ms=3, clip_on=False,
                color=colours.get(seq_cat[sid], "black"))
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=4)
    ax.set_xlabel("pruned-MSA position")
    ax.set_title(f"{og}: secondary structure per sequence ({len(rows)} seqs x {L} positions)",
                 fontsize=10)
    ax.legend(handles=[Patch(color=SS_COLOUR[s], label=SS_NAME[s]) for s in SS3]
                      + [Patch(color="lightgrey", label="gap")],
              fontsize=6, bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(out(f"{og}_ss_alignment.pdf"))
    plt.close(fig)

    print(f"{og}: secondary-structure plots for {len(rows)} sequences x {L} positions "
          f"-> {args.outdir}")


if __name__ == "__main__":
    main()
