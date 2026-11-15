# -*- coding: utf-8 -*-

"""Phase 2a: build a wide lookup table mapping every subcluster to its cluster,
supercluster and neurotransmitter annotation, with integrity checks.

Run from the repo root:
    python scripts/03_taxonomy_map.py
Output: data/atlas/taxonomy_map.csv (one row per subcluster, indexed by cluster_alias)
"""
from pathlib import Path

import pandas as pd

TAX = Path("data/atlas/taxonomy")
OUT = Path("data/atlas/taxonomy_map.csv")

members = pd.read_csv(TAX / "cluster_to_cluster_annotation_membership.csv")
terms = pd.read_csv(TAX / "cluster_annotation_term.csv")

SET = "cluster_annotation_term_set_name"


def collapse(values):
    # Some levels may carry more than one term per subcluster
    return ";".join(sorted(set(map(str, values))))


names = (members.groupby(["cluster_alias", SET])["cluster_annotation_term_name"]
         .agg(collapse).unstack())
labels = (members.groupby(["cluster_alias", SET])["cluster_annotation_term_label"]
          .agg(collapse).unstack())
names.columns = [f"{c}_name" for c in names.columns]
labels.columns = [f"{c}_label" for c in labels.columns]

cells = (members[members[SET] == "subcluster"]
         .set_index("cluster_alias")["number_of_cells"].rename("n_cells"))
tax = pd.concat([names, labels, cells], axis=1).sort_index()

# integrity checks 
print(f"Subclusters in lookup table: {len(tax):,}")
print(f"Total cells: {int(tax['n_cells'].sum()):,}")
per_alias = members.groupby(["cluster_alias", SET]).size()
for level in ["subcluster", "cluster", "supercluster"]:
    bad = (per_alias.xs(level, level=SET) != 1).sum()
    print(f"Subclusters without exactly one {level} term: {bad}")
print("Distinct terms per level:",
      {lvl: int(tax[f"{lvl}_label"].nunique()) for lvl in
       ["subcluster", "cluster", "supercluster"] if f"{lvl}_label" in tax})

# Cross-check supercluster cell counts against the term table
sup_from_members = tax.groupby("supercluster_label")["n_cells"].sum()
sup_terms = (terms[terms[SET] == "supercluster"].set_index("label")["number_of_cells"])
diff = (sup_from_members - sup_terms).dropna()
print(f"Superclusters with cell-count mismatch vs term table: {(diff != 0).sum()}")

# overview of the 31 superclusters 
overview = (tax.groupby(["supercluster_label", "supercluster_name"])
            .agg(n_subclusters=("n_cells", "size"), n_cells=("n_cells", "sum"))
            .reset_index().sort_values("n_cells", ascending=False))
print("\nSuperclusters:")
print(overview[["supercluster_name", "n_subclusters", "n_cells"]].to_string(index=False))

tax.to_csv(OUT)
print(f"\nWrote {OUT}")