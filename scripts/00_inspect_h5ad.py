# -*- coding: utf-8 -*-

"""Inspect an .h5ad without loading the expression matrix into memory.

Usage: python scripts/00_inspect_h5ad.py data/atlas/WHB-10Xv3-Nonneurons-raw.h5ad
"""
import sys
import anndata as ad

path = sys.argv[1]
adata = ad.read_h5ad(path, backed="r")

print(f"Shape (cells x genes): {adata.shape}")
print(f"Matrix storage: {adata.file['X'].attrs.get('encoding-type', 'dense or unknown')}")
print(f"X dtype: {adata.X.dtype if hasattr(adata.X, 'dtype') else 'n/a'}")
print("\nobs columns:", list(adata.obs.columns))
print(adata.obs.head(3).to_string())
print("\nvar columns:", list(adata.var.columns))
print(adata.var.head(3).to_string())
print("\nobs index example:", list(adata.obs_names[:3]))
print("var index example:", list(adata.var_names[:3]))
print("uns keys:", list(adata.uns.keys()))