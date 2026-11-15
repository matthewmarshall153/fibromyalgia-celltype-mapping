# -*- coding: utf-8 -*-

"""Phase 1b: list independent genome-wide significant loci."""
import numpy as np
import pandas as pd

df = pd.read_parquet("data/gwas/GCST90838603.qc.parquet")
sig = df[df["p_value"] < 5e-8].sort_values("p_value").copy()
sig["chromosome"] = sig["chromosome"].astype(str)

WINDOW = 500_000  # +/- 500 kb around each lead variant
leads = []
remaining = sig
while len(remaining):
    lead = remaining.iloc[0]
    leads.append(lead)
    near = (remaining["chromosome"] == lead["chromosome"]) & (
        (remaining["base_pair_location"] - lead["base_pair_location"]).abs() <= WINDOW
    )
    remaining = remaining[~near]

loci = pd.DataFrame(leads)[["chromosome", "base_pair_location", "rsid", "effect_allele",
                            "other_allele", "beta", "p_value"]]
loci["odds_ratio"] = np.exp(loci["beta"].astype(float))
chrom_order = loci["chromosome"].map(lambda v: int(v) if v.isdigit() else 99)
loci = loci.assign(_c=chrom_order).sort_values(["_c", "base_pair_location"]).drop(columns="_c")

print(f"{len(loci)} loci (500 kb distance pruning)")
print(loci.to_string(index=False))
loci.to_csv("results/gwas_loci.csv", index=False)