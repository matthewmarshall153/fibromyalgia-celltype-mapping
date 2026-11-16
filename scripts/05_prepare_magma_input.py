# -*- coding: utf-8 -*-

"""Phase 3a: prepare the fibromyalgia GWAS for MAGMA.

Matches GWAS variants to the reference panel (.bim) by chromosome + position + allele pair,
so it works whatever ID style the panel uses, then writes MAGMA input files that use the
panel's own variant IDs.

Usage (from the repo root):
    python scripts/05_prepare_magma_input.py

Outputs:
    data/gwas/magma_input.pval     columns SNP P N  (N = per-variant effective sample size)
    data/gwas/magma_input.snploc   columns SNP CHR BP
    results/magma_input_summary.json
"""
import argparse
import json
from pathlib import Path

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--gwas", default="data/gwas/GCST90838603.qc.parquet")
ap.add_argument("--bim", default="data/reference/magma/EUR/EUR.bim")
ap.add_argument("--out-prefix", default="data/gwas/magma_input")
ap.add_argument("--neff-frac", type=float, default=0.67,
                help="drop variants with n_eff below this fraction of the 90th percentile")
ap.add_argument("--chunk", type=int, default=2_000_000, help="bim rows per chunk")
args = ap.parse_args()

summary = {}


def norm_chr(s):
    s = s.astype("string").str.replace("^chr", "", regex=True, case=False)
    s = s.replace({"X": "23", "x": "23"})
    return pd.to_numeric(s, errors="coerce")


def allele_key(a1, a2):
    a1 = a1.astype("string").str.upper()
    a2 = a2.astype("string").str.upper()
    return (a1.where(a1 < a2, a2) + "_" + a2.where(a1 < a2, a1)).astype("category")


# GWAS side 
g = pd.read_parquet(args.gwas, columns=["chromosome", "base_pair_location", "effect_allele",
                                        "other_allele", "p_value", "n_eff"])
summary["gwas_variants_in"] = int(len(g))
g["chr"] = norm_chr(g["chromosome"])
g = g[g["chr"].between(1, 23)].copy()
g["chr"] = g["chr"].astype("int8")
g["pos"] = g["base_pair_location"].astype("int64")

thresh = args.neff_frac * g["n_eff"].quantile(0.9)
summary["n_eff_threshold"] = float(thresh)
keep = g["n_eff"] >= thresh
summary["removed_low_n_eff"] = int((~keep).sum())
g = g[keep]
g["akey"] = allele_key(g["effect_allele"], g["other_allele"])
g = g[["chr", "pos", "akey", "p_value", "n_eff"]]

before = len(g)
g = g.sort_values("n_eff", ascending=False).drop_duplicates(["chr", "pos", "akey"])
summary["removed_duplicate_keys"] = int(before - len(g))
summary["gwas_variants_after_filters"] = int(len(g))
print(f"GWAS variants after filters: {len(g):,} "
      f"(removed {summary['removed_low_n_eff']:,} low n_eff, "
      f"{summary['removed_duplicate_keys']:,} duplicate keys)")

# reference side, chunked 
matched = []
n_bim = 0
reader = pd.read_csv(args.bim, sep=r"\s+", header=None, usecols=[0, 1, 3, 4, 5],
                     names=["chr", "snp", "pos", "a1", "a2"],
                     dtype={"chr": "string", "snp": "string", "pos": "int64",
                            "a1": "string", "a2": "string"},
                     chunksize=args.chunk)
for chunk in reader:
    n_bim += len(chunk)
    chunk["chr"] = norm_chr(chunk["chr"])
    chunk = chunk[chunk["chr"].between(1, 23)]
    chunk["chr"] = chunk["chr"].astype("int8")
    chunk["akey"] = allele_key(chunk["a1"], chunk["a2"])
    m = chunk[["chr", "snp", "pos", "akey"]].merge(g, on=["chr", "pos", "akey"], how="inner")
    matched.append(m)
    print(f"  reference variants read: {n_bim:,}; matched so far: "
          f"{sum(len(x) for x in matched):,}", flush=True)

m = pd.concat(matched, ignore_index=True)
summary["reference_variants"] = int(n_bim)
summary["matched_before_id_dedup"] = int(len(m))

# MAGMA needs unique SNP IDs
dup_ids = int(m["snp"].duplicated(keep=False).sum())
summary["rows_with_duplicated_reference_id"] = dup_ids
m = m.drop_duplicates("snp", keep="first")
m = m[m["snp"].notna() & (m["snp"] != ".")]
summary["final_variants"] = int(len(m))
summary["fraction_of_gwas_matched"] = float(len(m) / max(len(g), 1))
print(f"\nFinal variants for MAGMA: {len(m):,} "
      f"({summary['fraction_of_gwas_matched']:.1%} of filtered GWAS variants)")
print(f"Rows with a duplicated reference ID (deduplicated): {dup_ids:,}")

# write 
prefix = Path(args.out_prefix)
prefix.parent.mkdir(parents=True, exist_ok=True)
Path("results").mkdir(exist_ok=True)
pd.DataFrame({"SNP": m["snp"], "P": m["p_value"], "N": m["n_eff"].round().astype(int)}) \
    .to_csv(f"{prefix}.pval", sep="\t", index=False, float_format="%.6g")
pd.DataFrame({"SNP": m["snp"], "CHR": m["chr"], "BP": m["pos"]}) \
    .to_csv(f"{prefix}.snploc", sep="\t", index=False, header=False)
with open("results/magma_input_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2)
print(f"Wrote {prefix}.pval, {prefix}.snploc and results/magma_input_summary.json")