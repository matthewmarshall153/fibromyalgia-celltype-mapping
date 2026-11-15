# -*- coding: utf-8 -*-


#!/usr/bin/env python
"""Phase 1: QC of the European-only fibromyalgia GWAS (GCST90838603).

Run from the repo root:
    python scripts/01_gwas_qc.py [path/to/GCST90838603.tsv]

Outputs
    data/gwas/GCST90838603.qc.parquet   cleaned table 
    results/gwas_qc_summary.json        small summary 
    results/figures/gwas_qq.png, gwas_manhattan.png
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

GWAS = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/gwas/GCST90838603.tsv")
OUT_TABLE = Path("data/gwas/GCST90838603.qc.parquet")
OUT_JSON = Path("results/gwas_qc_summary.json")
FIG_DIR = Path("results/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_JSON.parent.mkdir(parents=True, exist_ok=True)

EXPECTED = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency", "p_value",
            "rsid", "n", "n_eff"]

# Lead variants from Table 1 of Kerrebijn et al. 2026 (multi-ancestry analysis,
# so expect similar but not identical numbers in the European-only file).
LEADS = ["rs149109767", "rs9862795", "rs809955", "rs359243", "rs62100765",
         "rs6748639", "rs1788821"]

summary = {}

# 1. Load 
print(f"Reading {GWAS} ...")
df = pd.read_csv(
    GWAS, sep="\t",
    dtype={"chromosome": "string", "effect_allele": "string", "other_allele": "string",
           "rsid": "string", "base_pair_location": "int64",
           "beta": "float32", "standard_error": "float32",
           "effect_allele_frequency": "float32", "p_value": "float64",
           "n": "float64", "n_eff": "float64"},
)
missing_cols = [c for c in EXPECTED if c not in df.columns]
if missing_cols:
    sys.exit(f"Missing expected columns: {missing_cols}")
summary["n_variants_raw"] = int(len(df))
print(f"  {len(df):,} variants loaded")

# 2. Checks and filtering 
summary["missing_per_column"] = {c: int(df[c].isna().sum()) for c in EXPECTED}
bad_p = ~((df["p_value"] > 0) & (df["p_value"] <= 1))
bad_se = ~(df["standard_error"] > 0)
bad_eaf = ~((df["effect_allele_frequency"] > 0) & (df["effect_allele_frequency"] < 1))
summary["invalid_p"] = int(bad_p.sum())
summary["invalid_se"] = int(bad_se.sum())
summary["invalid_eaf"] = int(bad_eaf.sum())

keep = ~(df[EXPECTED].isna().any(axis=1) | bad_p | bad_se | bad_eaf)
df = df[keep].copy()
summary["n_variants_after_filter"] = int(len(df))

summary["duplicate_chr_pos_alleles"] = int(
    df.duplicated(["chromosome", "base_pair_location", "effect_allele", "other_allele"]).sum())
summary["duplicate_rsid"] = int(df["rsid"].duplicated().sum())
summary["chromosome_counts"] = {str(k): int(v) for k, v in df["chromosome"].value_counts().items()}

# 3. Effective sample size distribution 
qs = [0.01, 0.05, 0.25, 0.5, 0.75, 0.8, 0.95]
summary["n_eff_quantiles"] = {str(q): float(df["n_eff"].quantile(q)) for q in qs}
summary["n_quantiles"] = {str(q): float(df["n"].quantile(q)) for q in qs}

# 4. Inflation 
# lambda_GC is only a rough check for a polygenic trait at this sample size;
# the LD score regression intercept is the proper measure (later phase).
chisq = stats.chi2.isf(df["p_value"].to_numpy(), 1)
summary["lambda_gc"] = float(np.median(chisq) / stats.chi2.ppf(0.5, 1))

# 5. Genome-wide significant hits and known lead variants
sig = df[df["p_value"] < 5e-8]
summary["n_genome_wide_significant_variants"] = int(len(sig))
top = df.nsmallest(10, "p_value")[["rsid", "chromosome", "base_pair_location",
                                   "effect_allele", "other_allele", "beta", "p_value"]]
summary["top10"] = top.assign(
    beta=top["beta"].astype(float)).to_dict(orient="records")

lead_rows = df[df["rsid"].isin(LEADS)]
summary["lead_variants_found"] = [
    {"rsid": r.rsid, "effect_allele": r.effect_allele, "other_allele": r.other_allele,
     "beta": float(r.beta), "odds_ratio": float(np.exp(r.beta)),
     "p_value": float(r.p_value)}
    for r in lead_rows.itertuples()
]
summary["lead_variants_not_found"] = sorted(set(LEADS) - set(lead_rows["rsid"]))

# 6. Save cleaned table and summary 
df.to_parquet(OUT_TABLE, index=False)
with open(OUT_JSON, "w") as fh:
    json.dump(summary, fh, indent=2)

# 7. Figures
p_sorted = np.sort(df["p_value"].to_numpy())
n = len(p_sorted)
idx = np.unique(np.geomspace(1, n, 4000).astype(int)) - 1
exp_p = (idx + 0.5) / n
fig, ax = plt.subplots(figsize=(4.5, 4.5))
ax.scatter(-np.log10(exp_p), -np.log10(p_sorted[idx]), s=6)
lim = max(-np.log10(exp_p).max(), -np.log10(p_sorted[0]))
ax.plot([0, lim], [0, lim], color="grey", lw=1)
ax.set_xlabel("Expected -log10(P)")
ax.set_ylabel("Observed -log10(P)")
ax.set_title(f"QQ plot (lambda GC = {summary['lambda_gc']:.2f})")
fig.tight_layout()
fig.savefig(FIG_DIR / "gwas_qq.png", dpi=200)
plt.close(fig)

# Manhattan 
chrom_num = df["chromosome"].replace({"X": "23", "x": "23"})
chrom_num = pd.to_numeric(chrom_num, errors="coerce")
plot_df = pd.DataFrame({"chr": chrom_num, "pos": df["base_pair_location"],
                        "logp": -np.log10(df["p_value"])}).dropna(subset=["chr"])
rest = plot_df[plot_df["logp"] < 3].sample(frac=0.02, random_state=0)
plot_df = pd.concat([plot_df[plot_df["logp"] >= 3], rest])
chr_max = plot_df.groupby("chr")["pos"].max().sort_index()
offset = chr_max.cumsum().shift(fill_value=0)
plot_df["x"] = plot_df["pos"] + plot_df["chr"].map(offset)
fig, ax = plt.subplots(figsize=(11, 4))
for c, grp in plot_df.groupby("chr"):
    ax.scatter(grp["x"], grp["logp"], s=3, color="C0" if int(c) % 2 else "C1")
ax.axhline(-np.log10(5e-8), color="red", lw=1, ls="--")
ticks = offset + chr_max / 2
ax.set_xticks(ticks.to_numpy())
ax.set_xticklabels(["X" if int(c) == 23 else str(int(c)) for c in ticks.index], fontsize=7)
ax.set_xlabel("Chromosome")
ax.set_ylabel("-log10(P)")
ax.set_title("Fibromyalgia GWAS, European-only (GCST90838603)")
fig.tight_layout()
fig.savefig(FIG_DIR / "gwas_manhattan.png", dpi=200)
plt.close(fig)

# 8. Console report
print("\n=== QC summary ===")
print(f"Variants raw / kept: {summary['n_variants_raw']:,} / {summary['n_variants_after_filter']:,}")
print(f"Duplicate chr:pos:alleles: {summary['duplicate_chr_pos_alleles']:,}; duplicate rsid: {summary['duplicate_rsid']:,}")
print(f"lambda GC: {summary['lambda_gc']:.3f}")
print(f"Genome-wide significant variants: {summary['n_genome_wide_significant_variants']:,}")
print(f"n_eff quantiles: {summary['n_eff_quantiles']}")
print("Lead variants found:")
for r in summary["lead_variants_found"]:
    print(f"  {r['rsid']}  OR={r['odds_ratio']:.3f}  P={r['p_value']:.2e}")
print("Not found:", summary["lead_variants_not_found"])
print(f"\nWrote {OUT_TABLE}, {OUT_JSON}, and figures in {FIG_DIR}/")
