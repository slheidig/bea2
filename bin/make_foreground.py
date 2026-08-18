#!/usr/bin/env python3
"""Build csubst foreground file + strain lists from the categories table.

Reads the categories TSV (CRLF tolerated), resolves the requested column by
header name or 1-based index, and writes:
  foreground.tsv        csubst foreground (fg_format 1): lineage_id <TAB> .*_STRAIN_.*
  fg_strains.txt        foreground strain names, one per line
  outgroup_strains.txt  outgroup strain names, one per line
  categories_clean.tsv  CRLF-stripped copy of the categories table
"""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--categories", required=True)
    ap.add_argument("--column", default="temp_cat2", help="header name or 1-based index")
    ap.add_argument("--fg-level", default="warm")
    ap.add_argument("--outgroup-level", default="Outgroup")
    args = ap.parse_args()

    rows = []
    with open(args.categories) as fh:
        for line in fh:
            rows.append(line.rstrip("\n").replace("\r", "").split("\t"))
    header = rows[0]

    if args.column.isdigit():
        idx = int(args.column) - 1
        if not (0 <= idx < len(header)):
            sys.exit(f"ERROR: column index {args.column} out of range 1..{len(header)}")
    else:
        try:
            idx = header.index(args.column)
        except ValueError:
            sys.exit(f"ERROR: no column '{args.column}'. Available: {', '.join(header)}")
    print(f"Category column: {header[idx]} (index {idx + 1})")

    fg, og = [], []
    with open("categories_clean.tsv", "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows[1:]:
            if not r or not r[0]:
                continue
            fh.write("\t".join(r) + "\n")
            val = r[idx] if idx < len(r) else ""
            if val == args.fg_level:
                fg.append(r[0])
            elif val == args.outgroup_level:
                og.append(r[0])

    if len(fg) < 2:
        sys.exit(f"ERROR: need >=2 foreground strains for arity-2 convergence; "
                 f"got {len(fg)} with {header[idx]} == '{args.fg_level}'")

    with open("foreground.tsv", "w") as fh:
        for i, s in enumerate(fg, 1):
            fh.write(f"{i}\t.*_{s}_.*\n")
    open("fg_strains.txt", "w").write("\n".join(fg) + "\n")
    open("outgroup_strains.txt", "w").write(("\n".join(og) + "\n") if og else "")
    print(f"Foreground: {len(fg)} strains | Outgroup: {len(og)} strains")


if __name__ == "__main__":
    main()
