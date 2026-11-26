"""Step 7 - Gene Ontology biological-process labels, and their coverage of the network genes."""
import glob

import numpy as np
import pandas as pd

import _setup  # noqa: F401
from p36 import config as C
from p36.go import read_gaf, read_obo

genes = pd.read_csv(C.INTERIM / "genes.csv")
obo = glob.glob(str(C.RAW / "go" / "go-basic.obo*"))[0]
gaf = glob.glob(str(C.RAW / "go" / "*.gaf*"))[0]   # HUMAN-uniprot.gaf.gz (GO release in metadata/sources.tsv)
parents, namespace = read_obo(obo)
ann, n_lines, n_kept = read_gaf(gaf, genes.symbol, namespace, parents)

annotated = [g for g in genes.symbol if ann.get(g)]
counts = pd.Series([t for g in annotated for t in ann[g]]).value_counts()
terms = counts[(counts >= C.GO_TERM_MIN_GENES) & (counts <= C.GO_TERM_MAX_GENES)].index.tolist()
Y = pd.DataFrame(0, index=annotated, columns=terms, dtype=np.int8)
for g in annotated:
    hit = [t for t in ann[g] if t in Y.columns]
    Y.loc[g, hit] = 1
Y.to_csv(C.INTERIM / "go_labels.csv")

names = {}
cur = None
with open(obo) as fh:
    for line in fh:
        if line.startswith("id: GO:"):
            cur = line[4:].strip()
        elif line.startswith("name:") and cur:
            names[cur] = line[6:].strip()
pd.DataFrame({"term": terms, "name": [names.get(t, "") for t in terms], "genes": counts[terms].values}).to_csv(
    C.TABLES / "07_go_terms.csv", index=False)
cov = {"network_genes": len(genes), "genes_with_BP_annotation": len(annotated),
       "coverage": len(annotated) / len(genes), "terms_used": len(terms),
       "term_size_range": f"{C.GO_TERM_MIN_GENES}-{C.GO_TERM_MAX_GENES}",
       "excluded_evidence": " ".join(sorted(C.GO_EXCLUDE_EVIDENCE))}
pd.DataFrame([cov]).to_csv(C.TABLES / "07_go_coverage.csv", index=False)
print(pd.Series(cov).to_string())
