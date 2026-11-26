"""Step 1 - verify the assay, donors, conditions and cell types BEFORE any analysis.

The project brief says: do not assume the data are single-cell. This step checks that the matrix
holds integer UMI counts for individual barcoded cells, and tabulates cells per donor, condition and
cell type."""
import numpy as np
import pandas as pd

import _setup  # noqa: F401
from p36 import config as C
from p36.data import load_batch2

counts, genes, cells = load_batch2()
mc = C.META_COLS
print(f"Matrix: {counts.shape[0]:,} barcodes x {counts.shape[1]:,} genes; "
      f"{100 * counts.nnz / np.prod(counts.shape):.2f}% non-zero")
vals = counts.data[:100000]
print("Values are integer counts:", bool(np.all(vals == np.round(vals))), "| max count:", int(counts.data.max()))
umis = np.asarray(counts.sum(axis=1)).ravel()
print(f"UMIs per barcode: median {np.median(umis):.0f} (droplet single-cell range is ~500-5,000)")
matched = cells[mc["donor"]].notna().mean()
print(f"Barcodes matched to the GEO metadata table: {100 * matched:.1f}%")

tab = (cells.dropna(subset=[mc["donor"]])
       .groupby([mc["condition"], mc["multiplet"], mc["cell_type"]]).size().unstack(0).fillna(0).astype(int))
tab.to_csv(C.TABLES / "01_cells_by_type_condition.csv")
print("\nCells by type and condition (singlets and doublets):\n", tab)

sel = cells[(cells[mc["condition"]] == C.CONDITION) & (cells[mc["multiplet"]] == "singlet")
            & (cells[mc["cell_type"]] == C.CELL_TYPE)]
per_donor = sel[mc["donor"]].value_counts().sort_index()
per_donor.rename("cells").to_csv(C.TABLES / "01_selected_cells_per_donor.csv")
print(f"\nSelected: condition={C.CONDITION}, singlets, cell type={C.CELL_TYPE}: {len(sel):,} cells")
print(per_donor.to_string())
assert per_donor.size >= 4, "Need several donors for donor-level calibration"
