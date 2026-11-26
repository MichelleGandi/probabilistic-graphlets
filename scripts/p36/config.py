"""All settings for the pipeline live here. Change a value here, not inside the scripts."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]   # scripts/p36/config.py -> project folder
RAW = ROOT / "data" / "raw"            # downloaded files (GEO + Gene Ontology)
INTERIM = ROOT / "data" / "processed"  # intermediate files written by the scripts (rebuilt by code)
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"

SEED = 7

# ---- Dataset: Kang et al. 2018, GEO GSE96583, batch 2 (8 lupus donors, 10x Chromium) ----
GEO_ACCESSION = "GSE96583"
CONDITION = "ctrl"                 # unstimulated cells only (the other condition is IFN-beta "stim")
CELL_TYPE = "CD14+ Monocytes"      # cell-type label used in the GEO metadata table
META_COLS = {"donor": "ind", "condition": "stim", "cell_type": "cell", "multiplet": "multiplets"}

# ---- Cell and gene filtering ----
MIN_GENES_PER_CELL = 200           # drop cells with fewer detected genes
GENE_MIN_DETECTION = 0.10          # keep genes detected in >= 10% of the selected cells
N_GENES = 2000                     # then keep this many most variable genes
EXCLUDE_GENE_PREFIXES = ("MT-", "RPL", "RPS")   # mitochondrial and ribosomal genes

# ---- Edge confidence by resampling ----
METACELL_SIZE = 10                 # cells summed into one metacell (reduces single-cell sparsity)
N_BOOTSTRAP = 100                  # resamples of cells within each donor
TAU_QUANTILE = 0.99                # an edge is "co-expressed" in a resample if r >= the 99th percentile
                                   # of all gene-pair correlations in the full data (fixed threshold tau)

# ---- Networks ----
P_MIN = 0.01                       # probabilistic network keeps edges with calibrated p >= this
BINARY_THRESHOLDS = (0.5,)         # binary network 1: calibrated p >= 0.5 ("more likely than not to replicate")
# binary network 2, "matched": the top-k edges by p, with k = expected edges of the probabilistic network
#   (same density, so it differs from the probabilistic network ONLY in discarding uncertainty)
# binary network 3, "naive": raw correlation r >= tau, no resampling (what a typical analysis does)
MIN_EDGES = 50                     # a binary network with fewer edges is skipped (reported in step 5)

NETWORK_LABELS = {"probabilistic": "Probabilistic GDV", "binary_matched": "Binary GDV (matched density)",
                  "binary_naive_r": "Binary GDV (naive r >= tau)"}
for _t in BINARY_THRESHOLDS:
    NETWORK_LABELS[f"binary_p{_t:g}"] = f"Binary GDV (p >= {_t:g})"

# ---- Probabilistic graphlets ----
N_REALIZATIONS = 100               # Monte Carlo samples of the probabilistic network

# ---- Gene Ontology labels ----
GO_EXCLUDE_EVIDENCE = {"IEA", "ND", "IEP"}   # electronic, no data, and expression-based evidence
GO_TERM_MIN_GENES = 20
GO_TERM_MAX_GENES = 200

# ---- Prediction ----
N_FOLDS = 5
N_REPEATS = 3
N_EXPR_PCS = 20                    # expression-only baseline: 20 principal components + 3 summary stats

for _d in (RAW, INTERIM, TABLES, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)
