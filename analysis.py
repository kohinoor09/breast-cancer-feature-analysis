"""Capstone: which tumour-cell nuclear features differ between malignant and benign breast masses?
Dataset: Wisconsin Diagnostic Breast Cancer (UCI, public; bundled with scikit-learn).
Run:  python analysis.py
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.datasets import load_breast_cancer

# ---------- Parameters ----------
ALPHA = 0.05          # FDR / significance threshold
OUT = "results"

# ---------- Load data ----------
data = load_breast_cancer(as_frame=True)
df = data.frame
df["group"] = df["target"].map({0: "malignant", 1: "benign"})
features = list(data.feature_names)
mal, ben = df[df.group == "malignant"], df[df.group == "benign"]
print(f"Samples: {len(df)} (malignant={len(mal)}, benign={len(ben)}); features tested: {len(features)}")

# ---------- Statistics: Mann-Whitney U per feature ----------
rows = []
for f in features:
    u, p = stats.mannwhitneyu(mal[f], ben[f], alternative="two-sided")
    # effect size: rank-biserial correlation
    rbc = 2 * u / (len(mal) * len(ben)) - 1
    rows.append((f, mal[f].median(), ben[f].median(), rbc, p))
res = pd.DataFrame(rows, columns=["feature", "median_malignant", "median_benign", "rank_biserial", "p_value"])

# ---------- Multiple-testing control ----------
m = len(res)
res["p_bonferroni"] = np.minimum(res.p_value * m, 1.0)
order = res.p_value.argsort().values
ranked = res.p_value.values[order] * m / (np.arange(m) + 1)
bh = np.minimum.accumulate(ranked[::-1])[::-1]       # Benjamini-Hochberg
res.loc[order, "q_BH"] = np.minimum(bh, 1.0)
res["significant_FDR"] = res.q_BH < ALPHA
res = res.sort_values("p_value")
res.to_csv(f"{OUT}/feature_tests.csv", index=False)
print(f"Raw p<{ALPHA}: {(res.p_value<ALPHA).sum()} | Bonferroni: {(res.p_bonferroni<ALPHA).sum()} | BH-FDR: {res.significant_FDR.sum()}")

# ---------- Confounding check: features are highly correlated ----------
corr = df[features].corr(method="spearman")
corr.to_csv(f"{OUT}/feature_correlations.csv")
print("Median |Spearman r| between features:", round(corr.abs().values[np.triu_indices(m, 1)].mean(), 2))

# ---------- Figures ----------
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.scatter(res.rank_biserial, -np.log10(res.q_BH.clip(lower=1e-300)), c=np.where(res.significant_FDR, "crimson", "grey"))
ax.axhline(-np.log10(ALPHA), ls="--", c="k", lw=.8)
for _, r in res.head(5).iterrows():
    ax.annotate(r.feature, (r.rank_biserial, -np.log10(max(r.q_BH, 1e-300))), fontsize=6)
ax.set_xlabel("Effect size (rank-biserial; positive = higher in malignant)"); ax.set_ylabel("-log10(BH q-value)")
ax.set_title("Malignant vs benign: 30 features"); plt.tight_layout(); plt.savefig(f"{OUT}/volcano.png", dpi=150)

top = res.feature.iloc[0]
fig, ax = plt.subplots(figsize=(4, 4))
ax.boxplot([mal[top], ben[top]]); ax.set_xticklabels(["malignant", "benign"]); ax.set_title(top); plt.tight_layout()
plt.savefig(f"{OUT}/top_feature_boxplot.png", dpi=150)
print("Top feature:", top, "| outputs saved in", OUT)
