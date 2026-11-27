"""Settings for the pipeline.

The values come from config/config.yaml: change them there, not here. This module only
reads that file and works out the folder paths, so the scripts can keep using
`from p36 import config as C` and names like `C.N_BOOTSTRAP`.
"""
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]   # scripts/p36/config.py -> project folder
CONFIG_FILE = ROOT / "config" / "config.yaml"
with open(CONFIG_FILE) as fh:
    _c = yaml.safe_load(fh)

# ---- Folders (fixed by the project layout) ----
RAW = ROOT / "data" / "raw"            # downloaded files (GEO + Gene Ontology)
INTERIM = ROOT / "data" / "processed"  # intermediate files written by the scripts (rebuilt by code)
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"

SEED = int(_c["seed"])

# ---- Dataset ----
GEO_ACCESSION = _c["dataset"]["geo_accession"]
CONDITION = _c["dataset"]["condition"]
CELL_TYPE = _c["dataset"]["cell_type"]
META_COLS = dict(_c["dataset"]["metadata_columns"])

# ---- Cell and gene filtering ----
MIN_GENES_PER_CELL = int(_c["filtering"]["min_genes_per_cell"])
GENE_MIN_DETECTION = float(_c["filtering"]["gene_min_detection"])
N_GENES = int(_c["filtering"]["n_genes"])
EXCLUDE_GENE_PREFIXES = tuple(_c["filtering"]["exclude_gene_prefixes"])

# ---- Edge confidence by resampling ----
METACELL_SIZE = int(_c["edge_confidence"]["metacell_size"])
N_BOOTSTRAP = int(_c["edge_confidence"]["n_bootstrap"])
TAU_QUANTILE = float(_c["edge_confidence"]["tau_quantile"])

# ---- Networks ----
P_MIN = float(_c["networks"]["p_min"])
ZERO_IF_NEVER_STABLE = bool(_c["networks"]["zero_if_never_stable"])
BINARY_THRESHOLDS = tuple(float(t) for t in _c["networks"]["binary_thresholds"])
MIN_EDGES = int(_c["networks"]["min_edges"])
# Other networks, built in step 5: "matched" = the top-k pairs by p, with k = expected edges of the
# probabilistic network (same density, differs only in discarding uncertainty); "naive" = raw
# correlation r >= tau, no resampling (what a typical analysis does).

NETWORK_LABELS = {"probabilistic": "Probabilistic GDV", "binary_matched": "Binary GDV (matched density)",
                  "binary_naive_r": "Binary GDV (naive r >= tau)"}
for _t in BINARY_THRESHOLDS:
    NETWORK_LABELS[f"binary_p{_t:g}"] = f"Binary GDV (p >= {_t:g})"

# ---- Probabilistic graphlets ----
N_REALIZATIONS = int(_c["graphlets"]["n_realizations"])

# ---- Gene Ontology labels ----
GO_EXCLUDE_EVIDENCE = set(_c["gene_ontology"]["exclude_evidence"])
GO_TERM_MIN_GENES = int(_c["gene_ontology"]["term_min_genes"])
GO_TERM_MAX_GENES = int(_c["gene_ontology"]["term_max_genes"])

# ---- Prediction ----
N_FOLDS = int(_c["prediction"]["n_folds"])
N_REPEATS = int(_c["prediction"]["n_repeats"])
N_EXPR_PCS = int(_c["prediction"]["n_expression_pcs"])

for _d in (RAW, INTERIM, TABLES, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)
