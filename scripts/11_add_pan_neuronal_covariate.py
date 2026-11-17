# -*- coding: utf-8 -*-

"""Phase 6a: add a 'pan_neuronal' gene property to the MAGMA covariate files, so that
cell-type enrichment can be conditioned on general neuronal expression.

    pan_neuronal = CPM in all neuronal nuclei / (CPM in neuronal + CPM in non-neuronal nuclei)

Run from the repo root:  python scripts/11_add_pan_neuronal_covariate.py
Writes data/atlas/magma/covar_{level}_cond.txt (the original columns plus pan_neuronal).
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--neurons", default="data/atlas/WHB-10Xv3-Neurons-raw.pseudobulk.npz")
ap.add_argument("--nonneurons", default="data/atlas/WHB-10Xv3-Nonneurons-raw.pseudobulk.npz")
ap.add_argument("--levels", nargs="+", default=["supercluster", "cluster"])
ap.add_argument("--magma-dir", default="data/atlas/magma")
args = ap.parse_args()


def group_cpm(path):
    d = np.load(path, allow_pickle=False)
    return d["gene_ids"], d["counts"].sum(axis=0) / d["total_counts"].sum() * 1e6


genes_n, cpm_n = group_cpm(args.neurons)
genes_nn, cpm_nn = group_cpm(args.nonneurons)
assert np.array_equal(genes_n, genes_nn), "gene order differs between the two files"
den = cpm_n + cpm_nn
pan = pd.Series(np.where(den > 0, cpm_n / np.where(den > 0, den, 1), np.nan), index=genes_n)
print(f"pan_neuronal: {int(pan.notna().sum()):,} genes, median {pan.median():.2f}, "
      f"share of genes above 0.5: {(pan > 0.5).mean():.0%}")

for level in args.levels:
    f = Path(args.magma_dir) / f"covar_{level}.txt"
    covar = pd.read_csv(f, sep="\t", index_col="GENE")
    covar["pan_neuronal"] = pan.reindex(covar.index)
    covar = covar.dropna(subset=["pan_neuronal"])
    rho = covar.drop(columns="pan_neuronal").corrwith(covar["pan_neuronal"], method="spearman")
    covar.to_csv(Path(args.magma_dir) / f"covar_{level}_cond.txt", sep="\t",
                 index_label="GENE", float_format="%.5g")
    print(f"{level}: {covar.shape[1] - 1} cell types, {len(covar):,} genes; Spearman correlation with "
          f"pan_neuronal: min {rho.min():.2f}, median {rho.median():.2f}, max {rho.max():.2f}")