"""Write a small MOCK dataset in the same file formats as GSE96583 + GO, to test the pipeline end to end.

Run:  python tests/make_mock_data.py      (writes to data/raw/ ; delete it before using real data)"""
import gzip
import sys
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from p36 import config as C

rng = np.random.default_rng(0)
n_mod, mod_size, n_noise = 15, 20, 300
G = n_mod * mod_size + n_noise
symbols = [f"GENE{k:04d}" for k in range(G)] + ["MT-CO1", "RPL3"]
G = len(symbols)
donors = [str(1000 + d) for d in range(8)]
types = ["CD14+ Monocytes", "CD4 T cells", "B cells"]

def sample(sample_id, cond, n_cells):
    cells, rows = [], []
    for c in range(n_cells):
        d, t = rng.choice(donors), rng.choice(types, p=[0.4, 0.4, 0.2])
        base = np.full(G, 0.3)
        f = rng.gamma(2, 1, n_mod)                       # module activity in this cell
        for m in range(n_mod):
            base[m * mod_size:(m + 1) * mod_size] *= 1 + 2 * f[m]
        lib = rng.lognormal(7.3, 0.3)
        lam = base / base.sum() * lib
        rows.append(rng.poisson(lam))
        bc = f"{''.join(rng.choice(list('ACGT'), 14))}-1"
        cells.append((bc, d, cond, t, "singlet" if rng.random() > 0.05 else "doublet"))
    return np.array(rows), cells

meta_lines = ["\ttsne1\ttsne2\tind\tstim\tcluster\tcell\tmultiplets"]
for gsm, sid, cond, n in (("GSM2560248", "2.1", "ctrl", 2400), ("GSM2560249", "2.2", "stim", 600)):
    X, cells = sample(sid, cond, n)
    scipy.io.mmwrite(str(C.RAW / f"{gsm}_{sid}.mtx"), sp.csr_matrix(X.T))   # genes x cells, like 10x
    with gzip.open(C.RAW / f"{gsm}_barcodes.tsv.gz", "wt") as fh:
        fh.write("\n".join(c[0] for c in cells) + "\n")
    for bc, d, cnd, t, mult in cells:
        meta_lines.append(f"{bc}\t0\t0\t{d}\t{cnd}\t1\t{t}\t{mult}")
with gzip.open(C.RAW / "GSE96583_batch2.genes.tsv.gz", "wt") as fh:
    fh.write("\n".join(f"ENSG{k:011d}\t{s}" for k, s in enumerate(symbols)) + "\n")
with gzip.open(C.RAW / "GSE96583_batch2.total.tsne.df.tsv.gz", "wt") as fh:
    fh.write("\n".join(meta_lines) + "\n")

# GO: one term per module, grouped under 3 parent terms; annotations are noisy
obo = ["format-version: 1.2", "", "[Term]", "id: GO:0008150", "name: biological_process", "namespace: biological_process"]
for p in range(3):
    obo += ["", "[Term]", f"id: GO:20000{p:02d}", f"name: parent process {p}", "namespace: biological_process", "is_a: GO:0008150"]
for m in range(n_mod):
    obo += ["", "[Term]", f"id: GO:10000{m:02d}", f"name: module process {m}", "namespace: biological_process",
            f"is_a: GO:20000{m % 3:02d}"]
(C.RAW / "go-basic.obo").write_text("\n".join(obo) + "\n")
gaf = ["!gaf-version: 2.2"]
for m in range(n_mod):
    members = [symbols[m * mod_size + k] for k in range(mod_size) if rng.random() < 0.8]
    members += list(rng.choice(symbols[n_mod * mod_size:-2], 5, replace=False))
    for s in members:
        gaf.append("\t".join(["UniProtKB", "X", s, "involved_in", f"GO:10000{m:02d}", "PMID:1", "IDA", "", "P",
                              "", "", "protein", "taxon:9606", "20200101", "UniProt"]))
with gzip.open(C.RAW / "goa_human.gaf.gz", "wt") as fh:
    fh.write("\n".join(gaf) + "\n")
print("Mock data written to", C.RAW)
