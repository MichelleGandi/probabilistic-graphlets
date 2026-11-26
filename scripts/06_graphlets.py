"""Step 6 - graphlet degree vectors: probabilistic (expected counts, Monte Carlo) and binary."""
import numpy as np
import pandas as pd

import _setup  # noqa: F401
from p36 import config as C
from p36.graphlets import orbit_counts, probabilistic_orbit_counts

rng = np.random.default_rng(C.SEED)
genes = pd.read_csv(C.INTERIM / "genes.csv")
n = len(genes)
cols = [f"o{k}" for k in range(15)]

net = pd.read_csv(C.INTERIM / "net_probabilistic.csv")
print(f"Probabilistic network: {len(net):,} candidate edges, {net.p.sum():.0f} expected; {C.N_REALIZATIONS} realisations")
mean, se = probabilistic_orbit_counts(n, net[["i", "j"]].to_numpy(), net["p"].to_numpy(), C.N_REALIZATIONS, rng)
pd.DataFrame(mean, columns=cols, index=genes.symbol).to_csv(C.INTERIM / "gdv_probabilistic.csv")
pd.DataFrame(se, columns=cols, index=genes.symbol).to_csv(C.INTERIM / "gdv_probabilistic_se.csv")
rel_se = np.where(mean > 0, se / np.maximum(mean, 1e-12), 0)
print(f"Monte Carlo precision: median relative SE over non-zero entries = {np.median(rel_se[mean > 0]):.3f}")

binary = [n for n in pd.read_csv(C.TABLES / "05_networks.csv").network if n != "probabilistic"]
for name in binary:
    e = pd.read_csv(C.INTERIM / f"net_{name}.csv")
    o = orbit_counts(n, e[["i", "j"]].to_numpy())
    pd.DataFrame(o, columns=cols, index=genes.symbol).to_csv(C.INTERIM / f"gdv_{name}.csv")
    print(f"{name}: {len(e):,} edges, {(o[:, 0] > 0).sum():,} genes with at least one edge")

pd.DataFrame([{"realizations": C.N_REALIZATIONS, "median_relative_se": float(np.median(rel_se[mean > 0]))}]).to_csv(
    C.TABLES / "06_monte_carlo_precision.csv", index=False)
