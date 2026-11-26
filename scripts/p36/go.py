"""Gene Ontology biological-process labels for the network genes."""
import gzip
from collections import defaultdict

from . import config as C


def read_obo(path):
    """Return {term: set(parents)} for is_a and part_of, plus {term: namespace}; obsolete terms dropped."""
    parents, namespace, obsolete = defaultdict(set), {}, set()
    cur = None
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                cur = None
            elif line.startswith("[") and line.endswith("]"):
                cur = "SKIP"
            elif cur == "SKIP":
                continue
            elif line.startswith("id: GO:"):
                cur = line[4:]
            elif cur and line.startswith("namespace:"):
                namespace[cur] = line.split(": ", 1)[1]
            elif cur and line.startswith("is_a:"):
                parents[cur].add(line.split()[1])
            elif cur and line.startswith("relationship: part_of"):
                parents[cur].add(line.split()[2])
            elif cur and line == "is_obsolete: true":
                obsolete.add(cur)
    for t in obsolete:
        namespace.pop(t, None)
    return parents, namespace


def ancestors(term, parents, cache):
    if term not in cache:
        out = set()
        for p in parents.get(term, ()):
            out.add(p)
            out |= ancestors(p, parents, cache)
        cache[term] = out
    return cache[term]


def read_gaf(path, genes, namespace, parents):
    """Propagated GO BP annotations {symbol: set(terms)} for the given gene symbols."""
    genes = set(genes)
    direct = defaultdict(set)
    n_lines = n_kept = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            n_lines += 1
            symbol, qualifier, term, evidence, aspect = f[2], f[3], f[4], f[6], f[8]
            if aspect != "P" or "NOT" in qualifier or evidence in C.GO_EXCLUDE_EVIDENCE:
                continue
            if symbol in genes and namespace.get(term) == "biological_process":
                direct[symbol].add(term)
                n_kept += 1
    cache = {}
    prop = {}
    for g, terms in direct.items():
        allt = set(terms)
        for t in terms:
            allt |= ancestors(t, parents, cache)
        prop[g] = {t for t in allt if namespace.get(t) == "biological_process" and t != "GO:0008150"}
    return prop, n_lines, n_kept
