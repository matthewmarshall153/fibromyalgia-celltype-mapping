# -*- coding: utf-8 -*-

"""Phase 4: annotate MAGMA enrichment results, add multiple-testing correction, draw figures.

Run from the repo root:  python scripts/08_summarise_enrichment.py [--tag _noChr3]

Reads   results/magma/enrich_{sets,covar}_{supercluster,cluster}{tag}.gsa.out
        data/atlas/magma/celltype_map_{level}.csv
Writes  results/enrichment_results{tag}.csv
        results/figures/enrichment_supercluster{tag}.png, enrichment_cluster{tag}.png
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.stats.multitest import multipletests

ap = argparse.ArgumentParser()
ap.add_argument("--tag", default="", help="suffix used in the .gsa.out file names")
ap.add_argument("--magma-dir", default="results/magma")
ap.add_argument("--map-dir", default="data/atlas/magma")
args = ap.parse_args()

METHODS = {"sets": "Top 10% specific genes (gene set)",
           "covar": "Continuous specificity (gene property)"}
LEVELS = ["supercluster", "cluster"]
Path("results/figures").mkdir(parents=True, exist_ok=True)

frames = []
for level in LEVELS:
    cmap = pd.read_csv(Path(args.map_dir) / f"celltype_map_{level}.csv")
    for method in METHODS:
        f = Path(args.magma_dir) / f"enrich_{method}_{level}{args.tag}.gsa.out"
        df = pd.read_csv(f, sep=r"\s+", comment="#")
        if "FULL_NAME" in df.columns:          
            df["VARIABLE"] = df["FULL_NAME"]
        df = df.rename(columns={"VARIABLE": "variable", "NGENES": "n_genes", "BETA": "beta",
                                "BETA_STD": "beta_std", "SE": "se", "P": "p"})
        df = df[["variable", "n_genes", "beta", "beta_std", "se", "p"]]
        n = len(df)
        df["level"], df["method"], df["n_tests"] = level, method, n
        df["p_bonferroni"] = np.minimum(df["p"] * n, 1.0)
        df["fdr_bh"] = multipletests(df["p"], method="fdr_bh")[1]
        df["significant_bonferroni"] = df["p"] < 0.05 / n
        df = df.merge(cmap, on="variable", how="left")
        frames.append(df)

res = pd.concat(frames, ignore_index=True)
cols = ["level", "method", "variable", "label", "name", "supercluster_name", "n_cells", "n_genes", "n_tests",
        "beta", "beta_std", "se", "p", "p_bonferroni", "fdr_bh", "significant_bonferroni"]
res = res[cols].sort_values(["level", "method", "p"])
out_csv = Path(f"results/enrichment_results{args.tag}.csv")
res.to_csv(out_csv, index=False, float_format="%.6g")

# console summary
for (level, method), g in res.groupby(["level", "method"], sort=False):
    n = int(g["n_tests"].iloc[0])
    print(f"\n{level} | {METHODS[method]} | {n} tests, Bonferroni P < {0.05 / n:.2g}: "
          f"{int(g['significant_bonferroni'].sum())} significant")
    show = g.head(6)
    for r in show.itertuples():
        label = r.name if isinstance(r.name, str) and level == "supercluster" else r.variable
        parent = "" if level == "supercluster" else f"  [{r.supercluster_name}, {r.n_cells:,} cells]"
        print(f"   {label:34s} P = {r.p:.2g}{parent}")

# figure 1: superclusters 
sc = res[res["level"] == "supercluster"].copy()
sc["nlp"] = -np.log10(sc["p"])
order = (sc[sc["method"] == "covar"].sort_values("p")["variable"].tolist())
fig, axes = plt.subplots(1, 2, figsize=(11, 8), sharey=True)
for ax, method in zip(axes, METHODS):
    d = sc[sc["method"] == method].set_index("variable").loc[order]
    ax.barh(range(len(d)), d["nlp"], color="#4c72b0")
    ax.axvline(-np.log10(0.05 / len(d)), color="red", ls="--", lw=1, label="Bonferroni")
    ax.set_title(METHODS[method], fontsize=10)
    ax.set_xlabel("-log10(P), one-sided")
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels(d.index.str.replace("_", " "), fontsize=8)
axes[0].invert_yaxis()  # shared y-axis: invert once so the most significant is on top
axes[0].legend(loc="lower right", fontsize=8)
fig.suptitle(f"Fibromyalgia heritability enrichment across human brain superclusters{args.tag}")
fig.tight_layout()
fig.savefig(f"results/figures/enrichment_supercluster{args.tag}.png", dpi=200)
plt.close(fig)

# figure 2: clusters grouped by supercluster
cl = res[res["level"] == "cluster"].copy()
cl["nlp"] = -np.log10(cl["p"])
order = (cl[cl["method"] == "covar"].groupby("supercluster_name")["nlp"].max()
         .sort_values(ascending=False).index.tolist())
fig, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
for ax, method in zip(axes, METHODS):
    d = cl[cl["method"] == method]
    sns.stripplot(data=d, x="supercluster_name", y="nlp", order=order, size=3, jitter=0.25,
                  color="#4c72b0", ax=ax)
    ax.axhline(-np.log10(0.05 / d["n_tests"].iloc[0]), color="red", ls="--", lw=1)
    ax.set_ylabel("-log10(P)")
    ax.set_xlabel("")
    ax.set_title(METHODS[method], fontsize=10)
plt.setp(axes[1].get_xticklabels(), rotation=90, fontsize=7)
fig.suptitle(f"Cluster-level enrichment grouped by parent supercluster{args.tag} "
             "(dashed line = Bonferroni)")
fig.tight_layout()
fig.savefig(f"results/figures/enrichment_cluster{args.tag}.png", dpi=200)
plt.close(fig)
print(f"\nWrote {out_csv} and figures in results/figures/")