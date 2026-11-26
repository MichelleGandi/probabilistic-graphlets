"""Loading the GSE96583 files and building metacells."""
import glob
import gzip
import re

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

from . import config as C


def _find(pattern):
    hits = sorted(glob.glob(str(C.RAW / "**" / pattern), recursive=True))
    if not hits:
        raise FileNotFoundError(f"No file matching {pattern} in {C.RAW}. Run scripts/00_download.py first.")
    return hits[0]


def _read_lines(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        return [l.rstrip("\n") for l in fh if l.strip()]


def load_batch2():
    """Return (counts cells x genes CSR, genes DataFrame, cell metadata DataFrame) for batch 2."""
    genes_path = _find("*batch2.genes.tsv*")
    meta_path = _find("*batch2.total.tsne.df.tsv*")
    genes = pd.DataFrame([l.split("\t")[:2] for l in _read_lines(genes_path)], columns=["ensembl", "symbol"])
    meta = pd.read_csv(meta_path, sep="\t", index_col=0)

    blocks, cells = [], []
    for sample in ("2.1", "2.2"):
        mtx = _find(f"*_{sample}.mtx*")
        gsm = re.match(r"(GSM\d+)", mtx.split("/")[-1]).group(1)
        barcodes = _read_lines(_find(f"{gsm}*barcodes*"))
        m = scipy.io.mmread(mtx).tocsr()
        if m.shape[0] == len(genes) and m.shape[1] == len(barcodes):
            m = m.T.tocsr()                                  # 10x matrices are genes x cells
        assert m.shape == (len(barcodes), len(genes)), f"unexpected matrix shape {m.shape} for {mtx}"
        blocks.append(m)
        cells.append(pd.DataFrame({"barcode": barcodes, "sample": sample, "file": mtx.split("/")[-1]}))
    counts = sp.vstack(blocks).tocsr()
    cells = pd.concat(cells, ignore_index=True)
    cells = _attach_metadata(cells, meta)
    return counts, genes, cells


def _attach_metadata(cells, meta):
    """Match barcodes to the metadata table, tolerating '-1'-style suffixes."""
    cond_col = C.META_COLS["condition"]
    strip = lambda s: re.sub(r"-\d+$", "", str(s))
    out = []
    # sample 2.1 = control, 2.2 = IFN-beta stimulated (GEO sample descriptions)
    for sample, cond in (("2.1", "ctrl"), ("2.2", "stim")):
        sub = cells[cells["sample"] == sample].copy()
        m = meta[meta[cond_col] == cond] if cond_col in meta else meta
        exact = sub["barcode"].isin(m.index)
        if exact.mean() > 0.5:
            sub = sub.join(m, on="barcode")
        else:
            key = m.copy(); key.index = [strip(i) for i in key.index]
            key = key[~key.index.duplicated()]
            sub["_k"] = sub["barcode"].map(strip)
            sub = sub.join(key, on="_k").drop(columns="_k")
        out.append(sub)
    return pd.concat(out, ignore_index=True)


def lognorm(counts):
    """Counts per 10,000 then log1p. counts: cells x genes (dense or sparse)."""
    counts = sp.csr_matrix(counts, dtype=np.float64)
    size = np.asarray(counts.sum(axis=1)).ravel()
    size[size == 0] = 1
    x = sp.diags(1e4 / size) @ counts
    x.data = np.log1p(x.data)
    return x


def metacells(counts, donors, rng, size=C.METACELL_SIZE, resample=False):
    """Sum random groups of `size` cells within each donor; return log-normalised metacells x genes.

    resample=True first draws cells with replacement within each donor (a bootstrap)."""
    rows = []
    for d in np.unique(donors):
        idx = np.flatnonzero(donors == d)
        if resample:
            idx = rng.choice(idx, size=len(idx), replace=True)
        idx = rng.permutation(idx)
        n_groups = max(1, len(idx) // size)
        for g in np.array_split(idx, n_groups):
            rows.append(np.asarray(counts[g].sum(axis=0)).ravel())
    return lognorm(np.vstack(rows)).toarray()


def correlation(x):
    """Pearson correlation between genes (columns) of x; returns genes x genes."""
    z = x - x.mean(axis=0)
    sd = z.std(axis=0)
    sd[sd == 0] = np.inf
    z = z / sd
    r = (z.T @ z) / x.shape[0]
    np.fill_diagonal(r, 0.0)
    return r
