"""Step 8 - predict held-out GO annotations from each feature set.

For every GO term: a gene is positive if annotated with the term. Logistic regression is trained on 4/5
of the annotated genes and predicts the held-out 1/5 (5-fold cross-validation, repeated 3 times, the
SAME folds for every feature set so comparisons are paired). Score = AUPR (average precision)."""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from joblib import Parallel, delayed
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import _setup  # noqa: F401
from p36 import config as C
from p36.data import metacells
from p36.graphlets import gdv_features

genes = pd.read_csv(C.INTERIM / "genes.csv").set_index("symbol")
Y = pd.read_csv(C.INTERIM / "go_labels.csv", index_col=0)
G = Y.index

feats = {}
pg = pd.read_csv(C.INTERIM / "gdv_probabilistic.csv", index_col=0).loc[G]
feats["Probabilistic GDV"] = gdv_features(pg.to_numpy())
for name in pd.read_csv(C.TABLES / "05_networks.csv").network:
    if name != "probabilistic":
        gdv = pd.read_csv(C.INTERIM / f"gdv_{name}.csv", index_col=0).loc[G].to_numpy()
        feats[C.NETWORK_LABELS[name]] = gdv_features(gdv)
feats["Degree only (expected)"] = gdv_features(pg[["o0"]].to_numpy())

# expression-only baseline: summary statistics + principal components of each gene's metacell profile
rng = np.random.default_rng(C.SEED)
X = sp.load_npz(C.INTERIM / "counts.npz").tocsr()
donors = pd.read_csv(C.INTERIM / "cells.csv")["donor"].astype(str).to_numpy()
M = metacells(X, donors, rng)                      # metacells x genes
Z = (M - M.mean(0)) / np.where(M.std(0) > 0, M.std(0), 1)
pcs = PCA(n_components=C.N_EXPR_PCS, random_state=C.SEED).fit_transform(Z.T)   # genes x PCs
pcs = pd.DataFrame(pcs, index=genes.index).loc[G].to_numpy()
summ = genes.loc[G, ["mean_log", "detection", "var_log"]].to_numpy()
feats["Expression only"] = np.hstack([summ, pcs])

pd.DataFrame({k: [v.shape[1]] for k, v in feats.items()}).T.rename(columns={0: "n_features"}).to_csv(
    C.TABLES / "08_feature_sets.csv")


def evaluate(term):
    y = Y[term].to_numpy()
    cv = RepeatedStratifiedKFold(n_splits=C.N_FOLDS, n_repeats=C.N_REPEATS, random_state=C.SEED)
    splits = list(cv.split(np.zeros(len(y)), y))
    out = []
    for name, F in feats.items():
        ap, auc = [], []
        for tr, te in splits:
            model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, class_weight="balanced"))
            model.fit(F[tr], y[tr])
            s = model.predict_proba(F[te])[:, 1]
            ap.append(average_precision_score(y[te], s))
            auc.append(roc_auc_score(y[te], s))
        out.append({"term": term, "features": name, "aupr": np.mean(ap), "auroc": np.mean(auc),
                    "prevalence": y.mean(), "positives": int(y.sum())})
    return out


terms = [t for t in Y.columns if Y[t].sum() >= C.N_FOLDS]
print(f"Evaluating {len(terms)} GO terms x {len(feats)} feature sets x {C.N_FOLDS * C.N_REPEATS} folds")
res = Parallel(n_jobs=-1, verbose=5)(delayed(evaluate)(t) for t in terms)
res = pd.DataFrame([r for rr in res for r in rr])
res["aupr_fold_enrichment"] = res["aupr"] / res["prevalence"]
res.to_csv(C.TABLES / "08_predictions_per_term.csv", index=False)
print(res.groupby("features")[["aupr", "auroc", "aupr_fold_enrichment"]].median().round(3))
