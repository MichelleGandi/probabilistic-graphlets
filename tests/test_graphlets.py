"""Check that Monte Carlo probabilistic orbit counts agree with exact enumeration.

Run:  python tests/test_graphlets.py"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from p36.graphlets import exact_expected_orbits, orbit_counts, probabilistic_orbit_counts

rng = np.random.default_rng(1)
n = 7
edges = [(0, 1), (1, 2), (0, 2), (2, 3), (3, 4), (4, 5), (3, 5), (5, 6), (1, 3), (4, 6)]
probs = rng.uniform(0.2, 0.95, len(edges))

# 1) probabilities of 0/1 reduce to ordinary graphlet counting
assert np.allclose(exact_expected_orbits(n, edges, np.ones(len(edges))), orbit_counts(n, edges))

# 2) Monte Carlo converges to the exact expectation
exact = exact_expected_orbits(n, edges, probs)
mc, se = probabilistic_orbit_counts(n, edges, probs, 20000, rng)
z = np.abs(mc - exact) / np.maximum(se, 1e-9)
print("max |MC - exact| =", np.abs(mc - exact).max().round(4), "| max z-score =", z.max().round(2))
assert z.max() < 5, "Monte Carlo estimate disagrees with exact enumeration"

# 3) hand check: triangle with probabilities p1,p2,p3 -> expected 2-path count at the centre node
p1, p2, p3 = 0.9, 0.8, 0.5
e = exact_expected_orbits(3, [(0, 1), (1, 2), (0, 2)], [p1, p2, p3])
# node 1 is the middle of path 0-1-2 when edges 0-1 and 1-2 exist and 0-2 does not
assert np.isclose(e[1, 2], p1 * p2 * (1 - p3)), e[1, 2]
print("All graphlet tests passed.")
