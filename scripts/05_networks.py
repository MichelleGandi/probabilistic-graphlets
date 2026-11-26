"""Step 5 - build the probabilistic network and the binary comparison networks."""
import numpy as np
import pandas as pd

import _setup  # noqa: F401
from p36 import config as C

E = np.load(C.INTERIM / "edge_scores.npz")
p = np.load(C.INTERIM / "edge_p.npy")
i, j, r, tau = E["i"], E["j"], E["r_full"], float(E["tau"])
n = len(pd.read_csv(C.INTERIM / "genes.csv"))

nets = {"probabilistic": p >= C.P_MIN}
for t in C.BINARY_THRESHOLDS:
    nets[f"binary_p{t:g}"] = p >= t
k = int(round(p[p >= C.P_MIN].sum()))      # expected number of edges in the probabilistic network
top = np.lexsort((-r, -p))[:k]             # highest p first, ties broken by correlation
matched = np.zeros(len(p), bool); matched[top] = True
nets["binary_matched"] = matched
nets["binary_naive_r"] = r >= tau          # what a typical analysis does: threshold raw correlation once
for name in list(nets):
    if name != "probabilistic" and nets[name].sum() < C.MIN_EDGES:
        print(f"WARNING: {name} has only {nets[name].sum()} edges (< {C.MIN_EDGES}); skipped")
        del nets[name]

rows = []
for name, keep in nets.items():
    df = pd.DataFrame({"i": i[keep], "j": j[keep]})
    if name == "probabilistic":
        df["p"] = p[keep]
    df.to_csv(C.INTERIM / f"net_{name}.csv", index=False)
    deg = np.bincount(np.concatenate([df.i, df.j]), minlength=n)
    m = len(df)
    rows.append({"network": name, "genes": n, "genes_with_edges": int((deg > 0).sum()), "edges": m,
                 "expected_edges": float(df["p"].sum()) if "p" in df else m,
                 "density": 2 * m / (n * (n - 1)), "mean_degree": 2 * m / n})
tab = pd.DataFrame(rows)
tab.to_csv(C.TABLES / "05_networks.csv", index=False)
print(tab.to_string(index=False))
