# -*- coding: utf-8 -*-

"""Phase 2c: combine the two pseudobulk files, roll up to cluster/supercluster level, and
compute cell-type specificity (a gene's share of its total CPM across cell types).

Run from the repo root:  python scripts/06_specificity.py

Outputs (data/atlas/):
    specificity_{subcluster,cluster,supercluster}.parquet   genes (rows) x cell types (columns)
    celltypes_{level}.csv                                   label, name, parent supercluster, cells
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ap = argparse.ArgumentParser()
ap.add_argument("--npz", nargs="+",
                default=["data/atlas/WHB-10Xv3-Neurons-raw.pseudobulk.npz",
                         "data/atlas/WHB-10Xv3-Nonneurons-raw.pseudobulk.npz"])
ap.add_argument("--taxonomy", default="data/atlas/taxonomy_map.csv")
ap.add_argument("--genes", default="data/atlas/metadata/gene.csv")
ap.add_argument("--outdir", default="data/atlas")
args = ap.parse_args()
out = Path(args.outdir)

LEVELS = ["subcluster", "cluster", "supercluster"]
CHECK_TYPES = ["Hippocampal dentate gyrus", "Astrocyte", "Oligodendrocyte", "Medium spiny neuron"]

tax = pd.read_csv(args.taxonomy, index_col="cluster_alias").sort_index()
n_sub = int(tax.index.max()) + 1
assert len(tax) == n_sub and (tax.index.to_numpy() == np.arange(n_sub)).all(), \
    "taxonomy table must cover cluster_alias 0..N-1 without gaps"

# combine the pseudobulk files
counts = n_cells = total = None
seen = np.zeros(n_sub, dtype=bool)
gene_ids = None
for f in args.npz:
    d = np.load(f, allow_pickle=False)
    if gene_ids is None:
        gene_ids = d["gene_ids"]
        counts = np.zeros((n_sub, len(gene_ids)))
        n_cells = np.zeros(n_sub, dtype=np.int64)
        total = np.zeros(n_sub)
    elif not np.array_equal(gene_ids, d["gene_ids"]):
        raise SystemExit(f"Gene order in {f} differs from the first file")
    a = d["alias"]
    if seen[a].any():
        raise SystemExit(f"{f}: some subclusters were already present in an earlier file")
    seen[a] = True
    counts[a] = d["counts"]
    n_cells[a] = d["n_cells"]
    total[a] = d["total_counts"]
    print(f"{Path(f).name}: {len(a):,} subclusters, {int(d['n_cells'].sum()):,} cells")

print(f"\nSubclusters covered: {int(seen.sum()):,}/{n_sub:,}; cells: {int(n_cells.sum()):,}; "
      f"genes: {len(gene_ids):,}")
assert seen.all(), "some subclusters are missing from the inputs"
assert (n_cells == tax["n_cells"].to_numpy()).all(), "cell counts disagree with the taxonomy"

symbols = pd.read_csv(args.genes).set_index("gene_identifier")["gene_symbol"]
sym = symbols.reindex(gene_ids).fillna("").to_numpy()

# roll up and compute specificity 
for level in LEVELS:
    labels = tax[f"{level}_label"].to_numpy()
    uniq, codes = np.unique(labels, return_inverse=True)
    G = sp.csr_matrix((np.ones(n_sub), (codes, np.arange(n_sub))), shape=(len(uniq), n_sub))
    agg_counts = np.asarray(G @ counts)
    agg_total = np.asarray(G @ total)
    agg_cells = np.asarray(G @ n_cells).astype(np.int64)

    cpm = (agg_counts / agg_total[:, None] * 1e6).astype(np.float32)
    del agg_counts
    colsum = cpm.sum(axis=0, keepdims=True)
    spec = np.divide(cpm, colsum, out=np.zeros_like(cpm), where=colsum > 0)

    names = tax.groupby(f"{level}_label")[f"{level}_name"].first().reindex(uniq).to_numpy()
    parent = tax.groupby(f"{level}_label")["supercluster_name"].first().reindex(uniq).to_numpy()
    pd.DataFrame(spec.T, index=gene_ids, columns=uniq).to_parquet(out / f"specificity_{level}.parquet")
    pd.DataFrame({"label": uniq, "name": names, "supercluster_name": parent,
                  "n_cells": agg_cells, "total_counts": agg_total}) \
        .to_csv(out / f"celltypes_{level}.csv", index=False)
    print(f"{level}: {len(uniq):,} cell types, cells per type min/median "
          f"{agg_cells.min():,}/{int(np.median(agg_cells)):,}, "
          f"genes with no expression anywhere: {int((colsum == 0).sum()):,}")

    if level == "supercluster":
        print("\nTop 5 specific genes per supercluster (restricted to CPM >= 20 in that type):")
        for t in CHECK_TYPES:
            i = np.flatnonzero(names == t)
            if len(i) == 0:
                print(f"  {t}: not found")
                continue
            i = i[0]
            s = np.where(cpm[i] >= 20, spec[i], -1)
            top = np.argsort(s)[::-1][:5]
            print(f"  {t}: " + ", ".join(f"{sym[j] or gene_ids[j]} ({spec[i, j]:.2f})" for j in top))
    del cpm, spec

print(f"\nWrote specificity and cell-type tables to {out}/")