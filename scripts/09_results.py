"""Step 9 - statistics and figures answering the research question."""
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import _setup  # noqa: F401
from p36 import config as C
from p36.plotting import plt, BLUE, ORANGE, GREY

res = pd.read_csv(C.TABLES / "08_predictions_per_term.csv")
binary = [C.NETWORK_LABELS[n] for n in pd.read_csv(C.TABLES / "05_networks.csv").network if n != "probabilistic"]
order = ["Probabilistic GDV"] + binary + ["Degree only (expected)", "Expression only"]
wide = res.pivot(index="term", columns="features", values="aupr")[order]
fe = res.pivot(index="term", columns="features", values="aupr_fold_enrichment")[order]
auc = res.pivot(index="term", columns="features", values="auroc")[order]

summary = pd.DataFrame({
    "median_AUPR": wide.median(), "mean_AUPR": wide.mean(),
    "median_fold_enrichment_over_random": fe.median(),
    "terms_better_than_random_%": 100 * (fe > 1).mean(),
    "median_AUROC": auc.median()}).round(3)
summary.to_csv(C.TABLES / "09_summary_by_feature_set.csv")
print(summary.to_string())

rows = []
ref = "Probabilistic GDV"
for other in order[1:]:
    d = wide[ref] - wide[other]
    stat, pval = wilcoxon(wide[ref], wide[other])
    rows.append({"comparison": f"{ref} vs {other}", "terms": len(d), "median_AUPR_difference": d.median(),
                 "probabilistic_wins_%": 100 * (d > 0).mean(), "wilcoxon_p": pval})
tests = pd.DataFrame(rows)
# Benjamini-Hochberg across the comparisons
p = tests["wilcoxon_p"].to_numpy(); o = np.argsort(p); m = len(p)
q = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
tests.loc[o, "BH_q"] = np.minimum(q, 1)
tests.to_csv(C.TABLES / "09_paired_tests.csv", index=False)
print(tests.round(4).to_string(index=False))

colors = [BLUE] + [ORANGE] * len(binary) + [GREY] * 2
fig, ax = plt.subplots(figsize=(8, 4.2))
bp = ax.boxplot([np.log2(fe[c]) for c in order], orientation="horizontal", patch_artist=True, showfliers=False, widths=0.6)
for patch, col in zip(bp["boxes"], colors):
    patch.set_facecolor(col); patch.set_alpha(0.55)
ax.set_yticks(range(1, len(order) + 1), order); ax.invert_yaxis()
ax.axvline(0, color="black", lw=0.8, ls="--")
ax.set_xlabel("log2(AUPR / random), per GO term   (0 = no better than random)")
ax.set_title("Held-out GO function prediction, CD14+ monocyte co-expression network")
plt.tight_layout(); plt.savefig(C.FIGURES / "09_aupr_by_feature_set.png"); plt.close()

best_bin = max(binary, key=lambda c: wide[c].median())
fig, ax = plt.subplots(figsize=(4.6, 4.4))
ax.scatter(wide[best_bin], wide[ref], s=10, alpha=0.5, color=BLUE)
lim = [0, max(wide[ref].max(), wide[best_bin].max()) * 1.05]
ax.plot(lim, lim, color=GREY, ls="--", lw=1)
ax.set_xlabel(f"AUPR, {best_bin}"); ax.set_ylabel("AUPR, probabilistic GDV")
ax.set_title(f"Per GO term ({100 * (wide[ref] > wide[best_bin]).mean():.0f}% above the line)")
plt.tight_layout(); plt.savefig(C.FIGURES / "09_probabilistic_vs_best_binary.png"); plt.close()

nets = pd.read_csv(C.TABLES / "05_networks.csv")
fig, ax = plt.subplots(figsize=(6, 3))
ax.barh(nets.network, nets.expected_edges, color=[BLUE] + [ORANGE] * (len(nets) - 1), alpha=0.7)
ax.invert_yaxis(); ax.set_xlabel("Edges (expected edges for the probabilistic network)")
plt.tight_layout(); plt.savefig(C.FIGURES / "05_network_sizes.png"); plt.close()
print("Figures written to", C.FIGURES)
