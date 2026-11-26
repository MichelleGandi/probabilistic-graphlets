"""Graphlet degree vectors (orbits 0-14) for binary and probabilistic networks.

Probabilistic graphlets (Doria-Belenguer et al. 2020) are EXPECTED graphlet counts: every edge (u,v)
is present independently with probability p_uv. We estimate the expectation by Monte Carlo:
sample many realisations of the network, count orbits in each with ORCA (orbit-count), average.
The estimate is unbiased; its standard error shrinks as 1/sqrt(number of realisations).
`exact_expected_orbits` computes the expectation exactly by enumeration, for checking on small graphs.
"""
import itertools

import networkx as nx
import numpy as np
from orbit_count import node_orbit_counts

N_ORBITS = 15


def orbit_counts(n_nodes, edges):
    """Orbit counts (n_nodes x 15) of an unweighted graph on nodes 0..n_nodes-1."""
    G = nx.Graph()
    G.add_nodes_from(range(n_nodes))
    G.add_edges_from(map(tuple, edges))
    if G.number_of_edges() == 0:
        return np.zeros((n_nodes, N_ORBITS))
    return np.asarray(node_orbit_counts(G, graphlet_size=4, node_list=list(range(n_nodes))), dtype=float)


def probabilistic_orbit_counts(n_nodes, edges, probs, n_realizations, rng):
    """Expected orbit counts by Monte Carlo. Returns (mean, standard error), each n_nodes x 15."""
    edges = np.asarray(edges)
    probs = np.asarray(probs)
    total = np.zeros((n_nodes, N_ORBITS))
    total_sq = np.zeros((n_nodes, N_ORBITS))
    for _ in range(n_realizations):
        keep = rng.random(len(probs)) < probs
        o = orbit_counts(n_nodes, edges[keep])
        total += o
        total_sq += o ** 2
    mean = total / n_realizations
    var = np.maximum(total_sq / n_realizations - mean ** 2, 0)
    se = np.sqrt(var / n_realizations)
    return mean, se


def exact_expected_orbits(n_nodes, edges, probs):
    """Exact expected orbit counts by enumerating every edge configuration (tiny graphs only)."""
    edges = [tuple(e) for e in edges]
    if len(edges) > 18:
        raise ValueError("exact enumeration is only feasible for small graphs")
    expected = np.zeros((n_nodes, N_ORBITS))
    for keep in itertools.product([0, 1], repeat=len(edges)):
        w = np.prod([p if k else 1 - p for p, k in zip(probs, keep)])
        if w == 0:
            continue
        expected += w * orbit_counts(n_nodes, [e for e, k in zip(edges, keep) if k])
    return expected


def gdv_features(orbits):
    """Log-transform orbit counts, as in GDV similarity measures (counts span orders of magnitude)."""
    return np.log1p(orbits)
