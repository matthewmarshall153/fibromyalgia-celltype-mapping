# -*- coding: utf-8 -*-

"""Sanity check of the specificity matrix: where do well-known marker genes rank, and how
many subclusters are too small to trust?

Run from the repo root:  python scripts/06b_marker_check.py
"""
import pandas as pd

# Marker genes (from memory, so treat disagreements as questions, not proof) and the
# supercluster each is expected to be most specific to.
MARKERS = {
    "PROX1": "Hippocampal dentate gyrus",
    "AQP4": "Astrocyte",
    "GFAP": "Astrocyte",
    "MOBP": "Oligodendrocyte",
    "MBP": "Oligodendrocyte",
    "PLP1": "Oligodendrocyte",
    "PPP1R1B": "Medium spiny neuron",
    "DRD2": "Medium spiny neuron",
    "CSF1R": "Microglia",
}

spec = pd.read_parquet("data/atlas/specificity_supercluster.parquet")
ct = pd.read_csv("data/atlas/celltypes_supercluster.csv").set_index("label")
spec.columns = ct.loc[spec.columns, "name"].to_numpy()
genes = pd.read_csv("data/atlas/metadata/gene.csv")
sym2id = genes.drop_duplicates("gene_symbol").set_index("gene_symbol")["gene_identifier"]

print(f"{'gene':9s} {'expected':27s} {'rank':>4s}  top 3 superclusters (specificity)")
for g, expected in MARKERS.items():
    gid = sym2id.get(g)
    if gid is None or gid not in spec.index:
        print(f"{g:9s} not in the protein-coding gene table")
        continue
    row = spec.loc[gid].sort_values(ascending=False)
    rank = list(row.index).index(expected) + 1 if expected in row.index else None
    top3 = ", ".join(f"{n} ({v:.2f})" for n, v in row.head(3).items())
    print(f"{g:9s} {expected:27s} {str(rank):>4s}  {top3}")

sub = pd.read_csv("data/atlas/celltypes_subcluster.csv")
print(f"\nSubclusters: {len(sub):,}")
for t in (20, 50, 100, 200, 500):
    print(f"  with fewer than {t:>3d} cells: {int((sub['n_cells'] < t).sum()):,}")