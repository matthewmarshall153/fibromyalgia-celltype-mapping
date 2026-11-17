# -*- coding: utf-8 -*-

"""Phase 6b: compare cell-type enrichment before and after conditioning on pan-neuronal expression.

Run from the repo root:  python scripts/12_summarise_conditional.py
Reads  results/magma/enrich_covar_{supercluster,cluster}_cond.gsa.out,
       results/enrichment_results.csv, data/atlas/magma/celltype_map_{level}.csv
Writes results/enrichment_results_conditional.csv
"""
import argparse
from pathlib import Path

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--tag", default="_cond")
ap.add_argument("--magma-dir", default="results/magma")
ap.add_argument("--map-dir", default="data/atlas/magma")
args = ap.parse_args()

orig = pd.read_csv("results/enrichment_results.csv")
frames = []
for level in ["supercluster", "cluster"]:
    df = pd.read_csv(Path(args.magma_dir) / f"enrich_covar_{level}{args.tag}.gsa.out",
                     sep=r"\s+", comment="#")
    if "FULL_NAME" in df.columns:
        df["VARIABLE"] = df["FULL_NAME"]
    df = df.rename(columns={"VARIABLE": "variable", "BETA": "beta", "SE": "se", "P": "p_conditional"})
    df = df[["variable", "beta", "se", "p_conditional"]]
    n = len(df)
    df["level"], df["n_tests"] = level, n
    df["significant_bonferroni_conditional"] = df["p_conditional"] < 0.05 / n
    cmap = pd.read_csv(Path(args.map_dir) / f"celltype_map_{level}.csv")
    df = df.merge(cmap, on="variable", how="left")
    before = (orig[(orig["level"] == level) & (orig["method"] == "covar")]
              [["variable", "p", "beta"]].rename(columns={"p": "p_unconditional", "beta": "beta_unconditional"}))
    df = df.merge(before, on="variable", how="left")
    frames.append(df.sort_values("p_conditional"))

res = pd.concat(frames, ignore_index=True)
cols = ["level", "variable", "label", "name", "supercluster_name", "n_cells", "n_tests", "beta",
        "se", "p_conditional", "significant_bonferroni_conditional", "beta_unconditional", "p_unconditional"]
res[cols].to_csv("results/enrichment_results_conditional.csv", index=False, float_format="%.6g")

for level, g in res.groupby("level", sort=False):
    n = int(g["n_tests"].iloc[0])
    print(f"\n{level} | continuous specificity, conditioned on pan_neuronal | {n} tests, "
          f"Bonferroni P < {0.05 / n:.2g}: {int(g['significant_bonferroni_conditional'].sum())} significant "
          f"(unconditional: {int((g['p_unconditional'] < 0.05 / n).sum())})")
    for r in g.head(8).itertuples():
        nm = r.name if level == "supercluster" else f"{r.variable} [{r.supercluster_name}]"
        print(f"   {str(nm):50s} P cond = {r.p_conditional:.2g}   (before: {r.p_unconditional:.2g})")
print("\nWrote results/enrichment_results_conditional.csv")