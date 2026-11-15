# -*- coding: utf-8 -*-

"""Phase 2b: aggregate single-nucleus counts into one pseudobulk profile per subcluster.

Reads the raw-count .h5ad in backed mode, one chunk of cells at a time, so only a
(subclusters x protein-coding genes) accumulator is ever held in memory.

Usage (from the repo root):
    python scripts/04_pseudobulk.py data/atlas/WHB-10Xv3-Nonneurons-raw.h5ad
    python scripts/04_pseudobulk.py data/atlas/WHB-10Xv3-Nonneurons-raw.h5ad --limit 100000   # quick test

Output: data/atlas/<file stem>.pseudobulk.npz with
    counts       (n_present_subclusters x n_genes) summed raw counts, protein-coding genes only
    alias        cluster_alias of each row
    n_cells      cells contributing to each row
    total_counts all-gene UMI total of each row (for CPM later)
    gene_ids     Ensembl IDs of the columns
"""
import argparse
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ap = argparse.ArgumentParser()
ap.add_argument("h5ad")
ap.add_argument("--metadata", default="data/atlas/metadata/cell_metadata.csv")
ap.add_argument("--genes", default="data/atlas/metadata/gene.csv")
ap.add_argument("--taxonomy", default="data/atlas/taxonomy_map.csv")
ap.add_argument("--out", default=None)
ap.add_argument("--chunk", type=int, default=20000, help="cells per chunk")
ap.add_argument("--limit", type=int, default=None, help="only process the first N cells (testing)")
args = ap.parse_args()

h5 = Path(args.h5ad)
out = Path(args.out) if args.out else Path("data/atlas") / f"{h5.stem}.pseudobulk.npz"

# cell -> subcluster mapping 
need = ["cell_label", "cluster_alias"]
cols = pd.read_csv(args.metadata, nrows=0).columns.tolist()
missing = [c for c in need if c not in cols]
if missing:
    sys.exit(f"{args.metadata} lacks columns {missing}. Columns found: {cols}")
meta = pd.read_csv(args.metadata, usecols=need,
                   dtype={"cell_label": "string", "cluster_alias": "Int64"})
alias = meta.set_index("cell_label")["cluster_alias"]
if not alias.index.is_unique:
    sys.exit("cell_label is not unique in the metadata file")

tax = pd.read_csv(args.taxonomy, index_col="cluster_alias")
n_groups = int(tax.index.max()) + 1

# open the matrix, choose genes 
adata = ad.read_h5ad(h5, backed="r")
n = adata.n_obs if args.limit is None else min(args.limit, adata.n_obs)
print(f"{h5.name}: {adata.n_obs:,} cells x {adata.n_vars:,} genes; processing {n:,} cells")

genes = pd.read_csv(args.genes)
coding_ids = set(genes.loc[genes["biotype"] == "protein_coding", "gene_identifier"])
keep_idx = np.flatnonzero(adata.var_names.isin(coding_ids))
gene_ids = np.array(adata.var_names[keep_idx], dtype=str)
print(f"Keeping {len(keep_idx):,} protein-coding genes")

cell_alias = alias.reindex(adata.obs_names[:n]).to_numpy(dtype=float, na_value=np.nan)
n_unmapped = int(np.isnan(cell_alias).sum())
print(f"Cells without a subcluster in the metadata: {n_unmapped:,}")
if np.nanmax(cell_alias) >= n_groups:
    sys.exit("cluster_alias larger than the taxonomy table; check the files match")

# chunked aggregation 
acc = np.zeros((n_groups, len(keep_idx)), dtype=np.float64)
n_cells = np.zeros(n_groups, dtype=np.int64)
total_counts = np.zeros(n_groups, dtype=np.float64)
t0 = time.time()
for start in range(0, n, args.chunk):
    stop = min(start + args.chunk, n)
    g = cell_alias[start:stop]
    valid = ~np.isnan(g)
    if not valid.any():
        continue
    x = sp.csr_matrix(adata.X[start:stop])
    if not valid.all():
        x = x[np.flatnonzero(valid)]
        g = g[valid]
    g = g.astype(np.int64)

    lib = np.asarray(x.sum(axis=1)).ravel()
    total_counts += np.bincount(g, weights=lib, minlength=n_groups)
    n_cells += np.bincount(g, minlength=n_groups)

    G = sp.csr_matrix((np.ones(len(g)), (g, np.arange(len(g)))), shape=(n_groups, len(g)))
    res = (G @ x[:, keep_idx].astype(np.float64)).tocoo()
    acc[res.row, res.col] += res.data

    done = stop / n
    print(f"  {stop:>9,}/{n:,} cells ({done:5.1%}), {time.time() - t0:6.0f}s elapsed", flush=True)

# integrity checks 
present = np.flatnonzero(n_cells > 0)
print(f"\nSubclusters with cells: {len(present):,}")
print(f"Cells aggregated: {int(n_cells.sum()):,}")
if args.limit is None:
    expected = tax.loc[present, "n_cells"].to_numpy()
    mism = int((n_cells[present] != expected).sum())
    print(f"Subclusters whose cell count differs from the taxonomy table: {mism}")

np.savez(out, counts=acc[present], alias=present, n_cells=n_cells[present],
         total_counts=total_counts[present], gene_ids=gene_ids)
print(f"Wrote {out}")