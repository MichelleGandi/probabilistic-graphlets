"""Step 2 - select the cells and genes for the network."""
import numpy as np
import pandas as pd
import scipy.sparse as sp

import _setup  # noqa: F401
from p36 import config as C
from p36.data import load_batch2, lognorm

counts, genes, cells = load_batch2()
mc = C.META_COLS
keep = ((cells[mc["condition"]] == C.CONDITION) & (cells[mc["multiplet"]] == "singlet")
        & (cells[mc["cell_type"]] == C.CELL_TYPE)).to_numpy()
X = counts[keep]
cells = cells[keep].reset_index(drop=True)
n_genes_cell = np.asarray((X > 0).sum(axis=1)).ravel()
ok = n_genes_cell >= C.MIN_GENES_PER_CELL
X, cells = X[ok], cells[ok].reset_index(drop=True)
print(f"Cells after QC (>= {C.MIN_GENES_PER_CELL} genes): {X.shape[0]:,}")

# collapse duplicated symbols, drop mitochondrial/ribosomal genes, keep detected genes
sym = genes["symbol"].to_numpy()
detect = np.asarray((X > 0).mean(axis=0)).ravel()
L = lognorm(X)
mean = np.asarray(L.mean(axis=0)).ravel()
var = np.asarray(L.multiply(L).mean(axis=0)).ravel() - mean ** 2
g = pd.DataFrame({"symbol": sym, "ensembl": genes["ensembl"], "detection": detect, "mean_log": mean, "var_log": var})
g["col"] = np.arange(len(g))
g = g[~g["symbol"].str.upper().str.startswith(C.EXCLUDE_GENE_PREFIXES)]
g = g[g["detection"] >= C.GENE_MIN_DETECTION]
g = g.sort_values("var_log", ascending=False).drop_duplicates("symbol")
g = g.head(C.N_GENES).sort_values("symbol").reset_index(drop=True)
print(f"Genes detected in >= {100 * C.GENE_MIN_DETECTION:.0f}% of cells, top {C.N_GENES} by variance: {len(g)}")

Xg = X[:, g["col"].to_numpy()].tocsr()
sp.save_npz(C.INTERIM / "counts.npz", Xg)
g.drop(columns="col").to_csv(C.INTERIM / "genes.csv", index=False)
cells[["barcode", mc["donor"]]].rename(columns={mc["donor"]: "donor"}).to_csv(C.INTERIM / "cells.csv", index=False)
pd.DataFrame([{"cells": Xg.shape[0], "genes": Xg.shape[1], "donors": cells[mc["donor"]].nunique()}]).to_csv(
    C.TABLES / "02_dataset_size.csv", index=False)
