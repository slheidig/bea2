#!/usr/bin/env python3
"""Global statistics & distributions across ALL OGs (streaming, one OG at a time).

Inputs: the per-OG combined tables (plus optional RELAX summaries and csubst
branch-pair tables). Produces, in --outdir:
  stats.tsv                     per feature x category: n, mean, sd, quantiles
  feature_distributions.png     residue-level histogram per feature, by category
  og_means.png                  distribution of per-OG means per feature, by category
  og_means.tsv                  per-OG per-category mean of every feature
  relax_summary.tsv (+ .png)    all RELAX rows + K distribution
  branch_pairs_all.tsv (+ .png) all csubst pairs + OCN/omega distributions
  csubst_hotspot_counts.tsv     hotspot sites per OG
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ID_COLS = {"og", "sequence_id", "strain", "category", "msa_position",
           "original_msa_position", "residue_index", "residue"}
MAX_SAMPLE = 2_000_000  # residue-level values kept per feature x category


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--combined", nargs="+", required=True)
    ap.add_argument("--relax", nargs="*", default=[])
    ap.add_argument("--branch-pairs", nargs="*", default=[])
    ap.add_argument("--outdir", default="stats")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    out = lambda n: os.path.join(args.outdir, n)

    samples = {}      # (feature, category) -> list of arrays
    counts = {}
    og_means = []

    for f in sorted(args.combined):
        try:
            df = pd.read_csv(f, sep="\t", low_memory=False)
        except Exception as e:
            print(f"WARNING: skipping {f}: {e}")
            continue
        og = str(df["og"].iloc[0]) if "og" in df.columns and len(df) else os.path.basename(f)
        feats = [c for c in df.columns if c not in ID_COLS
                 and not c.startswith("hyphy_") and pd.api.types.is_numeric_dtype(
                     pd.to_numeric(df[c], errors="coerce"))]
        cat = df["category"].fillna("unknown") if "category" in df.columns else "all"
        for feat in feats:
            v = pd.to_numeric(df[feat], errors="coerce")
            for lev, g in v.groupby(cat):
                g = g.dropna()
                if g.empty:
                    continue
                key = (feat, lev)
                counts[key] = counts.get(key, 0) + len(g)
                kept = samples.setdefault(key, [])
                if sum(len(a) for a in kept) < MAX_SAMPLE:
                    kept.append(g.to_numpy())
                og_means.append({"og": og, "feature": feat, "category": lev,
                                 "mean": g.mean(), "n": len(g)})

    om = pd.DataFrame(og_means)
    if om.empty:
        raise SystemExit("ERROR: no data accumulated from combined tables")
    om.to_csv(out("og_means.tsv"), sep="\t", index=False)

    # stats.tsv
    rows = []
    for (feat, lev), arrs in sorted(samples.items()):
        a = np.concatenate(arrs)
        rows.append({"feature": feat, "category": lev, "n_residues": counts[(feat, lev)],
                     "mean": a.mean(), "sd": a.std(),
                     "q05": np.quantile(a, .05), "q25": np.quantile(a, .25),
                     "median": np.quantile(a, .5), "q75": np.quantile(a, .75),
                     "q95": np.quantile(a, .95)})
    pd.DataFrame(rows).to_csv(out("stats.tsv"), sep="\t", index=False, float_format="%.5g")

    feats = sorted({k[0] for k in samples})
    levels = sorted({k[1] for k in samples})
    cmap = plt.cm.tab10(np.linspace(0, 1, max(len(levels), 1)))
    lev_c = {l: cmap[i] for i, l in enumerate(levels)}

    def grid(n):
        nc = int(np.ceil(np.sqrt(n)))
        return int(np.ceil(n / nc)), nc

    # residue-level distributions
    nr, nc = grid(len(feats))
    fig, axes = plt.subplots(nr, nc, figsize=(4 * nc, 2.8 * nr), squeeze=False)
    for ax, feat in zip(axes.flat, feats):
        for lev in levels:
            arrs = samples.get((feat, lev))
            if not arrs:
                continue
            a = np.concatenate(arrs)
            ax.hist(a, bins=100, density=True, histtype="step", color=lev_c[lev], label=lev)
        ax.set_title(feat, fontsize=8)
        ax.tick_params(labelsize=6)
    for ax in axes.flat[len(feats):]:
        ax.axis("off")
    axes.flat[0].legend(fontsize=6)
    fig.suptitle("Residue-level distributions by category")
    fig.tight_layout()
    fig.savefig(out("feature_distributions.png"), dpi=140)
    plt.close(fig)

    # per-OG means distributions
    fig, axes = plt.subplots(nr, nc, figsize=(4 * nc, 2.8 * nr), squeeze=False)
    for ax, feat in zip(axes.flat, feats):
        sub = om[om["feature"] == feat]
        for lev in levels:
            v = sub.loc[sub["category"] == lev, "mean"].dropna()
            if v.empty:
                continue
            ax.hist(v, bins=50, density=True, histtype="step", color=lev_c[lev], label=lev)
        ax.set_title(feat, fontsize=8)
        ax.tick_params(labelsize=6)
    for ax in axes.flat[len(feats):]:
        ax.axis("off")
    axes.flat[0].legend(fontsize=6)
    fig.suptitle("Per-OG mean distributions by category")
    fig.tight_layout()
    fig.savefig(out("og_means.png"), dpi=140)
    plt.close(fig)

    # RELAX
    if args.relax:
        rel = pd.concat([pd.read_csv(f, sep="\t") for f in args.relax], ignore_index=True)
        rel.to_csv(out("relax_summary.tsv"), sep="\t", index=False)
        K = pd.to_numeric(rel["K"], errors="coerce").dropna()
        if len(K):
            fig, ax = plt.subplots(figsize=(5, 3))
            ax.hist(np.clip(K, 0, 10), bins=60, color="#4477aa")
            ax.axvline(1, color="k", ls="--", lw=0.8)
            ax.set_xlabel("RELAX K (clipped at 10)")
            ax.set_ylabel("OGs")
            ax.set_title("Selection intensity on foreground (K<1 relaxed, K>1 intensified)")
            fig.tight_layout()
            fig.savefig(out("relax_K_hist.png"), dpi=140)
            plt.close(fig)

    # csubst branch pairs + hotspot counts
    if args.branch_pairs:
        bp = pd.concat([pd.read_csv(f, sep="\t") for f in args.branch_pairs], ignore_index=True)
        bp.to_csv(out("branch_pairs_all.tsv"), sep="\t", index=False)
        fig, axes = plt.subplots(1, 2, figsize=(9, 3))
        for ax, colname in zip(axes, ("OCN", "omegaC")):
            v = pd.to_numeric(bp[colname], errors="coerce").dropna()
            if len(v):
                ax.hist(np.clip(v, 0, 20), bins=60, color="#cc6677")
            ax.set_xlabel(f"{colname} (clipped at 20)")
            ax.set_ylabel("pairs")
        fig.suptitle("csubst significant foreground pairs, all OGs")
        fig.tight_layout()
        fig.savefig(out("branch_pairs_hist.png"), dpi=140)
        plt.close(fig)
        (bp.groupby(["og", "signal"]).size().rename("n_pairs").reset_index()
           .to_csv(out("csubst_hotspot_counts.tsv"), sep="\t", index=False))

    print(f"global stats written to {args.outdir}/")


if __name__ == "__main__":
    main()
