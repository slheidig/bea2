#!/usr/bin/env python3
"""Join all HyPhy per-site tables of one OG into a wide by-site table (+ plot).

Input persite TSVs are named <og>.<METHOD>[.fg].persite.tsv; the label becomes
'<METHOD>' or '<METHOD>_fg'. codon_site == pruned-MSA position (canonical).
The original MSA position is added from msa_columns.tsv.

Outputs: <og>.by_site.tsv and (optionally) <og>.by_site.png — one subplot per
method; bar height = dN/dS (capped) or -log10(p); colour by significance.
"""
import argparse
import csv
import math
import os


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def read_tsv(p):
    with open(p) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return (rows[0], rows[1:]) if rows else ([], [])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--og", required=True)
    ap.add_argument("--columns", required=True, help="msa_columns.tsv")
    ap.add_argument("--persite", nargs="+", required=True)
    ap.add_argument("--plot", default="yes")
    ap.add_argument("--dnds-cap", type=float, default=10.0)
    args = ap.parse_args()
    og = args.og

    # pruned -> original MSA position
    orig = {}
    for r in csv.DictReader(open(args.columns), delimiter="\t"):
        orig[int(r["pruned_position"])] = int(r["original_position"])

    per = {}
    for f in args.persite:
        base = os.path.basename(f)
        lab = base[len(og) + 1:-len(".persite.tsv")].replace(".fg", "_fg")
        hdr, rows = read_tsv(f)
        if not hdr or hdr[0].startswith("#"):
            continue
        H = {h.lower(): i for i, h in enumerate(hdr)}

        def first(*names):
            for n in names:
                if n in H:
                    return H[n]
            return None

        want = {"dNdS": first("dn_ds"), "class": first("selection"),
                "p": first("p-value", "p-value (overall)"), "q": first("q-value (overall)"),
                "Ppos": first("prob[alpha<beta]")}
        cols = {k: v for k, v in want.items() if v is not None}
        ci = H.get("codon_site")
        if ci is None:
            continue
        d = {}
        for row in rows:
            d[int(row[ci])] = {k: row[i] for k, i in cols.items()}
        per[lab] = (list(cols.keys()), d)

    sites = sorted({s for _, (ks, d) in per.items() for s in d})
    labels = sorted(per)
    out = f"{og}.by_site.tsv"
    with open(out, "w") as fh:
        header = ["og", "codon_site", "original_msa_position"]
        for lab in labels:
            header += [f"{lab}.{k}" for k in per[lab][0]]
        fh.write("\t".join(header) + "\n")
        for s in sites:
            cells = [og, str(s), str(orig.get(s, ""))]
            for lab in labels:
                ks, d = per[lab]
                rec = d.get(s, {})
                cells += [rec.get(k, "") for k in ks]
            fh.write("\t".join(cells) + "\n")
    print(f"{og}: wrote {out} ({len(sites)} sites x {len(labels)} analyses)")

    if args.plot.lower() not in ("yes", "true", "1"):
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
    except ImportError:
        print("matplotlib unavailable — table written, plot skipped")
        return

    SIG1, SIG2, CAP = 0.05, 0.1, args.dnds_cap
    GREEN, YELLOW, GREY = "#4daf4a", "#ffcc33", "#bdbdbd"

    def sig_of(rec):
        for k in ("q", "p"):
            v = num(rec.get(k))
            if v is not None:
                return v
        v = num(rec.get("Ppos"))
        return (1.0 - v) if v is not None else None

    def colour(s):
        return GREEN if (s is not None and s <= SIG1) else (YELLOW if (s is not None and s <= SIG2) else GREY)

    usable = [lab for lab in labels if ("dNdS" in per[lab][0] or "p" in per[lab][0])]
    if not usable or not sites:
        return
    fig, axes = plt.subplots(len(usable), 1,
                             figsize=(max(7, len(sites) / 22), 1.7 * len(usable) + 1),
                             sharex=True, squeeze=False)
    axes = axes[:, 0]
    for ax, lab in zip(axes, usable):
        ks, d = per[lab]
        use_dnds = "dNdS" in ks
        xs, hs, cs = [], [], []
        for s in sites:
            rec = d.get(s, {})
            if use_dnds:
                v = num(rec.get("dNdS"))
                if v is None:
                    continue
                h = min(v, CAP)
            else:
                pv = num(rec.get("p"))
                if pv is None or pv <= 0:
                    continue
                h = -math.log10(pv)
            xs.append(s)
            hs.append(h)
            cs.append(colour(sig_of(rec)))
        ax.bar(xs, hs, width=1.0, align="center", color=cs, edgecolor="none")
        if use_dnds:
            ax.axhline(1.0, color="k", lw=0.6, ls="--")
            ax.set_ylabel(f"{lab}\ndN/dS (cap {CAP:g})", fontsize=8)
        else:
            ax.set_ylabel(f"{lab}\n-log10(p)", fontsize=8)
        ax.tick_params(labelsize=7)
    axes[0].set_title(f"{og}: per-site selection (pruned-MSA position)", fontsize=10)
    axes[-1].set_xlabel("pruned-MSA position (codon)", fontsize=9)
    fig.legend(handles=[Patch(color=GREEN, label="sig<=0.05"),
                        Patch(color=YELLOW, label="sig<=0.1"),
                        Patch(color=GREY, label="ns")],
               loc="upper right", fontsize=7, framealpha=0.9,
               title="significance (q if avail, else p)", title_fontsize=6)
    fig.tight_layout()
    fig.savefig(f"{og}.by_site.png", dpi=140)
    print(f"{og}: wrote {og}.by_site.png")


if __name__ == "__main__":
    main()
