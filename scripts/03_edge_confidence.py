"""Step 3 - edge stability by resampling cells within donors.

For each gene pair we ask: in what fraction of bootstrap resamples is the pair co-expressed, i.e.
correlation >= tau (a fixed threshold, the 99th percentile of all correlations in the full data)?
This is computed (a) on all donors, and (b) separately on two disjoint halves of the donors, which
step 4 uses to calibrate stability into a probability of replication in independent donors."""
import numpy as np
import pandas as pd
import scipy.sparse as sp

import _setup  # noqa: F401
from p36 import config as C
from p36.data import correlation, metacells

rng = np.random.default_rng(C.SEED)
X = sp.load_npz(C.INTERIM / "counts.npz").tocsr()
donors = pd.read_csv(C.INTERIM / "cells.csv")["donor"].astype(str).to_numpy()
n = X.shape[1]
iu = np.triu_indices(n, k=1)

r_full = correlation(metacells(X, donors, rng))[iu]
tau = float(np.quantile(r_full, C.TAU_QUANTILE))
print(f"Full-data correlation threshold tau (q={C.TAU_QUANTILE}): {tau:.3f}")


def stability(mask, n_boot):
    hits = np.zeros(len(iu[0]), dtype=np.int32)
    for b in range(n_boot):
        r = correlation(metacells(X[mask], donors[mask], rng, resample=True))[iu]
        hits += r >= tau
        if (b + 1) % 25 == 0:
            print(f"   bootstrap {b + 1}/{n_boot}")
    return hits / n_boot


ud = np.array(sorted(set(donors)), dtype=object)
perm = rng.permutation(ud)
half_a, half_b = {str(x) for x in perm[: len(ud) // 2]}, {str(x) for x in perm[len(ud) // 2:]}
in_a = np.isin(donors, list(half_a))
print("Donor halves:", sorted(half_a), "|", sorted(half_b))

print("Stability, all donors"); s_full = stability(np.ones(len(donors), bool), C.N_BOOTSTRAP)
print("Stability, half A"); s_a = stability(in_a, C.N_BOOTSTRAP)
print("Stability, half B"); s_b = stability(~in_a, C.N_BOOTSTRAP)
# replication calls: is the pair co-expressed (r >= tau) in the OTHER half's own data?
called_a = correlation(metacells(X[in_a], donors[in_a], rng))[iu] >= tau
called_b = correlation(metacells(X[~in_a], donors[~in_a], rng))[iu] >= tau

np.savez_compressed(C.INTERIM / "edge_scores.npz", i=iu[0], j=iu[1], r_full=r_full.astype(np.float32),
                    s_full=s_full.astype(np.float32), s_a=s_a.astype(np.float32), s_b=s_b.astype(np.float32),
                    called_a=called_a, called_b=called_b, tau=tau)
pd.DataFrame([{"tau": tau, "donors_half_A": " ".join(sorted(half_a)), "donors_half_B": " ".join(sorted(half_b)),
               "pairs": len(r_full), "pairs_stability_gt0": int((s_full > 0).sum())}]).to_csv(
    C.TABLES / "03_edge_confidence_summary.csv", index=False)
