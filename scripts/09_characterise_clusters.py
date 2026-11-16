# -*- coding: utf-8 -*-

"""Phase 5: profile the significant clusters: marker genes, donor and brain-region composition.

Run from the repo root:  python scripts/09_characterise_clusters.py
Optional: --labels CS202210140_137 CS202210140_164   (default: Bonferroni-significant clusters,
          continuous-specificity method)

Writes results/top_cluster_profiles.csv
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="results/enrichment_results.csv")
ap.add_argument("--labels", nargs="*", default=None)
ap.add_argument("--npz", nargs="+",
                default=["data/atlas/WHB-10Xv3-Neurons-raw.pseudobulk.npz",
                         "data/atlas/WHB-10Xv3-Nonneurons-raw.pseudobulk.npz"])
ap.add_argument("--taxonomy", default="data/atlas/taxonomy_map.csv")
ap.add_argument("--spec", default="data/atlas/specificity_cluster.parquet")
ap.add_argument("--genes", default="data/atlas/metadata/gene.csv")
ap.add_argument("--meta", default="data/atlas/metadata/cell_metadata.csv")
ap.add_argument("--min-cpm", type=float, default=20.0)
ap.add_argument("--out", default="results/top_cluster_profiles.csv")
args = ap.parse_args()

res = pd.read_csv(args.results)
cl = res[(res["level"] == "cluster") & (res["method"] == "covar")].set_index("label")
labels = args.labels or cl.index[cl["significant_bonferroni"]].tolist()
if not labels:
    raise SystemExit("No clusters to profile (none significant and no --labels given)")
print(f"Profiling {len(labels)} clusters")

tax = pd.read_csv(args.taxonomy, index_col="cluster_alias")
alias_of = {lab: tax.index[tax["cluster_label"] == lab].to_numpy() for lab in labels}
target = np.unique(np.concatenate(list(alias_of.values())))

# pseudobulk counts for the target subclusters 
counts, total, gene_ids = {}, {}, None
for f in args.npz:
    d = np.load(f, allow_pickle=False)
    if gene_ids is None:
        gene_ids = d["gene_ids"]
    a = d["alias"]
    sel = np.flatnonzero(np.isin(a, target))
    if len(sel) == 0:
        continue
    c, t = d["counts"][sel], d["total_counts"][sel]
    for al, row, tt in zip(a[sel], c, t):
        counts[al], total[al] = row, tt

spec = pd.read_parquet(args.spec, columns=labels)
symbols = (pd.read_csv(args.genes).drop_duplicates("gene_identifier")
           .set_index("gene_identifier")["gene_symbol"])

# donor / region composition 
meta = pd.read_csv(args.meta,
                   usecols=["cluster_alias", "donor_label", "anatomical_division_label",
                            "region_of_interest_label"],
                   dtype={"cluster_alias": "int32", "donor_label": "category",
                          "anatomical_division_label": "category",
                          "region_of_interest_label": "category"})
meta["cluster_label"] = meta["cluster_alias"].map(tax["cluster_label"])
base = meta["donor_label"].value_counts(normalize=True)
print("Donor shares across the whole atlas: " +
      ", ".join(f"{k} {v:.0%}" for k, v in base.items()))


def shares(series, k):
    vc = series.value_counts(normalize=True).head(k)
    return "; ".join(f"{i} {v:.0%}" for i, v in vc.items())


rows = []
for lab in labels:
    aliases = alias_of[lab]
    c = np.sum([counts[a] for a in aliases], axis=0)
    t = float(np.sum([total[a] for a in aliases]))
    cpm = c / t * 1e6
    s = spec[lab].reindex(gene_ids).to_numpy()
    s = np.where(cpm >= args.min_cpm, s, -1)
    top = np.argsort(s)[::-1][:10]
    markers = ", ".join(str(symbols.get(gene_ids[j], gene_ids[j])) for j in top)

    sub = meta[meta["cluster_label"] == lab]
    rows.append({
        "cluster_label": lab,
        "supercluster": tax.loc[aliases[0], "supercluster_name"],
        "n_cells": int(tax.loc[aliases, "n_cells"].sum()),
        "p_covar": cl.loc[lab, "p"] if lab in cl.index else np.nan,
        "donor_shares": shares(sub["donor_label"], 3),
        "top_divisions": shares(sub["anatomical_division_label"], 3),
        "top_regions": shares(sub["region_of_interest_label"], 3),
        "top_markers": markers,
    })
    print(f"\n{lab}  [{rows[-1]['supercluster']}, {rows[-1]['n_cells']:,} cells]")
    print(f"  markers : {markers}")
    print(f"  donors  : {rows[-1]['donor_shares']}")
    print(f"  division: {rows[-1]['top_divisions']}")
    print(f"  region  : {rows[-1]['top_regions']}")

Path(args.out).parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(args.out, index=False)
print(f"\nWrote {args.out}")