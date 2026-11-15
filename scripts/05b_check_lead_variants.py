# -*- coding: utf-8 -*-

"""Check that the lead variants from the loci table survived the n_eff filter and
matching.

Run from the repo root:  python scripts/05b_check_lead_variants.py
"""
import json

import pandas as pd

loci = pd.read_csv("results/gwas_loci.csv", dtype={"chromosome": str})
thresh = json.load(open("results/magma_input_summary.json"))["n_eff_threshold"]

# n_eff of each lead variant, read from the cleaned GWAS (only the needed rows)
filters = [[("chromosome", "==", str(c)), ("base_pair_location", "==", int(p))]
           for c, p in zip(loci["chromosome"], loci["base_pair_location"])]
gw = pd.read_parquet("data/gwas/GCST90838603.qc.parquet",
                     columns=["chromosome", "base_pair_location", "effect_allele", "n_eff"],
                     filters=filters)
gw["chromosome"] = gw["chromosome"].astype(str)

out = loci.merge(gw, on=["chromosome", "base_pair_location", "effect_allele"], how="left")
out["found_in_cleaned_gwas"] = out["n_eff"].notna()
out["passes_n_eff_filter"] = out["n_eff"] >= thresh

# Positions present in the MAGMA input (position-level check, alleles not compared)
snploc = pd.read_csv("data/gwas/magma_input.snploc", sep="\t", usecols=["CHR", "BP"])
present = set(zip(snploc["CHR"].astype(int), snploc["BP"].astype(int)))
chr_num = out["chromosome"].replace({"X": "23"}).astype(int)
out["position_in_magma_input"] = [(c, int(p)) in present
                                  for c, p in zip(chr_num, out["base_pair_location"])]

cols = ["chromosome", "base_pair_location", "rsid", "p_value", "n_eff",
        "found_in_cleaned_gwas", "passes_n_eff_filter", "position_in_magma_input"]
print(f"n_eff threshold: {thresh:,.0f}")
print(out[cols].to_string(index=False))
print(f"\nLead variants found in the cleaned GWAS: {int(out['found_in_cleaned_gwas'].sum())}/{len(out)}")
print(f"Lead variants kept by the filter: {int(out['passes_n_eff_filter'].sum())}/{len(out)}")
print(f"Lead variant positions present in MAGMA input: {int(out['position_in_magma_input'].sum())}/{len(out)}")