# -*- coding: utf-8 -*-

"""Phase 3b: turn specificity matrices into MAGMA inputs.

For each level it writes, into data/atlas/magma/:
    covar_{level}.txt         gene-property file for --gene-covar (header; first column = gene ID)
    sets_{level}.txt          row-based gene sets for --set-annot (top fraction of genes per cell type)
    celltype_map_{level}.csv  variable name <-> label, name, parent supercluster, cell count

Cell types with fewer than --min-cells cells are excluded, and specificity is re-normalised
across the cell types that remain (so each gene's values again sum to 1).

Run from the repo root:  python scripts/07_make_magma_celltype_inputs.py
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--levels", nargs="+", default=["supercluster", "cluster"])
ap.add_argument("--min-cells", type=int, default=100)
ap.add_argument("--top-frac", type=float, default=0.10)
ap.add_argument("--atlas-dir", default="data/atlas")
ap.add_argument("--outdir", default="data/atlas/magma")
args = ap.parse_args()

atlas = Path(args.atlas_dir)
out = Path(args.outdir)
out.mkdir(parents=True, exist_ok=True)

for level in args.levels:
    spec = pd.read_parquet(atlas / f"specificity_{level}.parquet")
    ct = pd.read_csv(atlas / f"celltypes_{level}.csv").set_index("label")

    keep = ct.index[ct["n_cells"] >= args.min_cells]
    dropped = len(ct) - len(keep)
    spec = spec[keep].astype(np.float64)
    total = spec.sum(axis=1)
    spec = spec[total > 0].div(total[total > 0], axis=0)

    # Variable names: readable names for the 31 superclusters, labels otherwise
    if level == "supercluster":
        var = ct.loc[keep, "name"].map(lambda s: re.sub(r"\W+", "_", s).strip("_"))
    else:
        var = pd.Series(keep, index=keep)
    assert var.is_unique, "variable names are not unique"
    spec.columns = var.loc[spec.columns].to_numpy()

    spec.to_csv(out / f"covar_{level}.txt", sep="\t", index_label="GENE", float_format="%.5g")

    n_top = int(round(args.top_frac * len(spec)))
    with open(out / f"sets_{level}.txt", "w") as fh:
        for name in spec.columns:
            top = spec[name].nlargest(n_top).index
            fh.write(name + "\t" + "\t".join(top) + "\n")

    cmap = ct.loc[keep].assign(variable=var.to_numpy())
    cmap.reset_index()[["variable", "label", "name", "supercluster_name", "n_cells"]] \
        .to_csv(out / f"celltype_map_{level}.csv", index=False)

    print(f"{level}: kept {len(keep):,} cell types (dropped {dropped:,} with < {args.min_cells} cells); "
          f"{len(spec):,} genes; top set size {n_top:,}")

print(f"\nWrote MAGMA inputs to {out}/")