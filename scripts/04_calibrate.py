"""Step 4 - turn resampling stability into a calibrated probability.

Raw stability (fraction of resamples) is NOT automatically a probability that the association is real.
We calibrate it against an independent check: stability in one half of the donors is mapped, by
isotonic regression, to how often the pair is co-expressed in the OTHER half of the donors.
The final edge confidence p = probability that the association replicates in independent donors."""
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression

import _setup  # noqa: F401
from p36 import config as C
from p36.plotting import plt, BLUE, ORANGE, GREY

E = np.load(C.INTERIM / "edge_scores.npz")
s_a, s_b, called_a, called_b = E["s_a"], E["s_b"], E["called_a"], E["called_b"]


def fit(s, y):
    return IsotonicRegression(y_min=0, y_max=1, out_of_bounds="clip", increasing=True).fit(s, y.astype(float))


brier = lambda p, y: float(np.mean((p - y) ** 2))
# cross-fitted evaluation: fit on one direction, test on the other
iso_ab, iso_ba = fit(s_a, called_b), fit(s_b, called_a)
rows = []
for name, s, y, iso in (("A->B", s_b, called_a, iso_ab), ("B->A", s_a, called_b, iso_ba)):
    rows.append({"tested_on": name, "brier_raw_stability": brier(s, y), "brier_calibrated": brier(iso.predict(s), y),
                 "replication_rate_overall": float(y.mean())})
ev = pd.DataFrame(rows)
print(ev.to_string(index=False))

iso = fit(np.concatenate([s_a, s_b]), np.concatenate([called_b, called_a]))
p = 0.5 * (iso.predict(s_a) + iso.predict(s_b))
np.save(C.INTERIM / "edge_p.npy", p.astype(np.float32))
grid = np.linspace(0, 1, 101)
pd.DataFrame({"stability": grid, "calibrated_p": iso.predict(grid)}).to_csv(C.TABLES / "04_calibration_map.csv", index=False)
ev.to_csv(C.TABLES / "04_calibration_brier.csv", index=False)

# reliability diagram (pooled, out-of-sample: s_half predicts the other half)
s_all, y_all = np.concatenate([s_a, s_b]), np.concatenate([called_b, called_a]).astype(float)
p_cross = np.concatenate([iso_ba.predict(s_a), iso_ab.predict(s_b)])
bins = np.linspace(0, 1, 11)
rel = []
for label, score in (("raw stability", s_all), ("calibrated p", p_cross)):
    idx = np.clip(np.digitize(score, bins) - 1, 0, 9)
    for k in range(10):
        m = idx == k
        if m.sum() >= 20:
            rel.append({"score": label, "bin": k, "mean_score": score[m].mean(), "observed": y_all[m].mean(), "n": int(m.sum())})
rel = pd.DataFrame(rel)
rel.to_csv(C.TABLES / "04_reliability.csv", index=False)

fig, axs = plt.subplots(1, 2, figsize=(9, 4))
ax = axs[0]
ax.plot([0, 1], [0, 1], ls="--", color=GREY, lw=1, label="perfect calibration")
for label, col in (("raw stability", ORANGE), ("calibrated p", BLUE)):
    d = rel[rel.score == label]
    ax.plot(d.mean_score, d.observed, "o-", color=col, label=label)
ax.set_xlabel("Predicted probability of replication"); ax.set_ylabel("Observed replication in other donors")
ax.set_title("Reliability diagram (out-of-sample)"); ax.legend(frameon=False)
ax = axs[1]
ax.plot(grid, iso.predict(grid), color=BLUE)
ax.set_xlabel("Bootstrap stability (half of donors)"); ax.set_ylabel("Calibrated p")
ax.set_title("Calibration map (isotonic regression)")
plt.tight_layout(); plt.savefig(C.FIGURES / "04_calibration.png"); plt.close()

m = s_all > 0
print(f"Pairs with stability > 0: {int((E['s_full'] > 0).sum()):,}; with calibrated p >= {C.P_MIN}: {int((p >= C.P_MIN).sum()):,}")
