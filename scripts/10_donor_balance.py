# -*- coding: utf-8 -*-

"""Phase 5b: donor balance of every cell type, and enrichment results restricted to
donor-balanced cell types.

NOTE: the threshold was chosen after inspecting the significant clusters, so this is a
post-hoc robustness check, not a pre-planned test. Report it as such.

Run from the repo root:  python scripts/10_donor_balance.py [--threshold 0.6]

Writes  results/donor_balance_{supercluster,cluster}.csv
        results/enrichment_results_donor_balanced.csv
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

ap = argparse.ArgumentParser()
ap.add_argument("--threshold", type=float, default=0.6,
                help="keep cell types whose top donor contributes at most this share of cells")
ap.add_argument("--results", default="results/enrichment_results.csv")
ap.add_argument("--taxonomy", default="data/atlas/taxonomy_map.csv")
ap.add_argument("--meta", default="data/atlas/metadata/cell_metadata.csv")
ap.add_argument("--outdir", default="results")
args = ap.parse_args()
out = Path(args.outdir)
out.mkdir(parents=True, exist_ok=True)

tax = pd.read_csv(args.taxonomy, index_col="cluster_alias")
meta = pd.read_csv(args.meta, usecols=["cluster_alias", "donor_label"],
                   dtype={"cluster_alias": "int32", "donor_label": "category"})
print("Donor shares across the whole atlas: " + ", ".join(
    f"{k} {v:.0%}" for k, v in meta["donor_label"].value_counts(normalize=True).items()))

bal = {}
for level in ["supercluster", "cluster"]:
    lab = meta["cluster_alias"].map(tax[f"{level}_label"])
    ct = pd.crosstab(lab, meta["donor_label"])
    n = ct.sum(axis=1)
    df = pd.DataFrame({"label": ct.index, "n_cells": n.to_numpy(),
                       "top_donor": ct.idxmax(axis=1).to_numpy(),
                       "top_donor_share": (ct.max(axis=1) / n).to_numpy()})
    df.to_csv(out / f"donor_balance_{level}.csv", index=False)
    bal[level] = df.set_index("label")
    print(f"\n{level}: {len(df):,} cell types; top-donor share > 0.5: "
          f"{int((df.top_donor_share > 0.5).sum())}, > 0.6: {int((df.top_donor_share > 0.6).sum())}, "
          f"> 0.8: {int((df.top_donor_share > 0.8).sum())}")

res = pd.read_csv(args.results)
res["top_donor"] = [bal[lv].top_donor.get(lb, np.nan) for lv, lb in zip(res["level"], res["label"])]
res["top_donor_share"] = [bal[lv].top_donor_share.get(lb, np.nan)
                          for lv, lb in zip(res["level"], res["label"])]

lost = res[res["significant_bonferroni"] & (res["top_donor_share"] > args.threshold)]
print(f"\nBonferroni-significant results removed as donor-unbalanced (> {args.threshold:.0%} one donor):")
for r in lost.itertuples():
    nm = r.name if r.level == "supercluster" else r.variable
    print(f"  {r.level:12s} {r.method:5s} {nm:36s} top donor share {r.top_donor_share:.0%}  P = {r.p:.2g}")

parts = []
for (level, method), g in res[res["top_donor_share"] <= args.threshold].groupby(["level", "method"], sort=False):
    g = g.copy()
    n = len(g)
    g["n_tests"] = n
    g["p_bonferroni"] = np.minimum(g["p"] * n, 1.0)
    g["fdr_bh"] = multipletests(g["p"], method="fdr_bh")[1]
    g["significant_bonferroni"] = g["p"] < 0.05 / n
    parts.append(g.sort_values("p"))
    print(f"\n{level} | {method} | {n} donor-balanced tests, Bonferroni P < {0.05 / n:.2g}: "
          f"{int(g['significant_bonferroni'].sum())} significant")
    for r in g.head(6).itertuples():
        nm = r.name if level == "supercluster" else f"{r.variable} [{r.supercluster_name}]"
        print(f"   {nm:52s} P = {r.p:.2g}   top donor {r.top_donor_share:.0%}")

pd.concat(parts).to_csv(out / "enrichment_results_donor_balanced.csv", index=False, float_format="%.6g")
print(f"\nWrote {out}/donor_balance_*.csv and {out}/enrichment_results_donor_balanced.csv")