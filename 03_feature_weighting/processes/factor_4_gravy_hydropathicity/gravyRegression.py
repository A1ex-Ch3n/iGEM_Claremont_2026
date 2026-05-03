# ============================================================
# ★ 配置区 — 已根据你的文件名设置好，不需要再改 ★
# ============================================================

PHAGE_DIR    = "data/phage"
BACTERIA_DIR = "data/bacteria"
LABELS_CSV   = "data/labels.csv"
OUTPUT_DIR   = "results"

# ★ 告诉程序哪些文件是phage，哪些是bacteria（用文件名前缀判断）
PHAGE_PREFIX    = ["AB", "AP", "EU", "JN"]   # phage文件名开头
# bacteria用剩余文件，不需要额外配置

# ============================================================

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, roc_auc_score, roc_curve
import warnings
warnings.filterwarnings("ignore")

os.makedirs(OUTPUT_DIR, exist_ok=True)

def compute_gravy_list(faa_path):
    gravy_values = []
    skipped = 0
    for record in SeqIO.parse(faa_path, "fasta"):
        seq = str(record.seq).upper()
        clean_seq = ''.join(aa for aa in seq if aa in "ACDEFGHIKLMNPQRSTVWY")
        if len(clean_seq) < 10:
            skipped += 1
            continue
        try:
            gravy = ProteinAnalysis(clean_seq).gravy()
            gravy_values.append(gravy)
        except Exception:
            skipped += 1
    if skipped > 0:
        print(f"  [跳过] {os.path.basename(faa_path)}: {skipped} 条无效序列")
    return gravy_values

def get_mean_gravy(faa_path):
    vals = compute_gravy_list(faa_path)
    if not vals:
        raise ValueError(f"文件 {faa_path} 没有有效序列！")
    return np.mean(vals), np.median(vals), vals

print("=" * 55)
print("Step 1: 计算所有 .faa 文件的 GRAVY 值")
print("=" * 55)

gravy_cache = {}

for folder in [PHAGE_DIR, BACTERIA_DIR]:
    for fname in sorted(os.listdir(folder)):
        if not fname.endswith(".faa"):
            continue
        fid = fname.replace(".faa", "")
        fpath = os.path.join(folder, fname)
        print(f"  处理: {fid}")
        mean_g, med_g, vals = get_mean_gravy(fpath)
        gravy_cache[fid] = {"mean": mean_g, "median": med_g, "vals": vals}
        print(f"    mean GRAVY = {mean_g:.4f} | median = {med_g:.4f} | n = {len(vals)} 蛋白")

print("\n" + "=" * 55)
print("Step 2: 读取标签，构建特征表")
print("=" * 55)

labels_df = pd.read_csv(LABELS_CSV)
print(f"  共 {len(labels_df)} 对组合 | y=1: {labels_df['y'].sum()} | y=0: {(labels_df['y']==0).sum()}")

rows = []
for _, row in labels_df.iterrows():
    pid = str(row["phage_id"])
    bid = str(row["host_id"])
    if pid not in gravy_cache:
        print(f"  [警告] 找不到 phage: {pid}.faa，跳过")
        continue
    if bid not in gravy_cache:
        print(f"  [警告] 找不到 bacteria: {bid}.faa，跳过")
        continue
    gp = gravy_cache[pid]["mean"]
    gb = gravy_cache[bid]["mean"]
    x  = abs(gp - gb)
    rows.append({"phage_id": pid, "host_id": bid,
                 "gravy_phage": round(gp,4), "gravy_host": round(gb,4),
                 "x_gravy": round(x,4), "y": int(row["y"])})

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUTPUT_DIR, "gravy_features.csv"), index=False)
print("\n  特征表预览：")
print(df.to_string(index=False))

print("\n" + "=" * 55)
print("Step 3: Density Plot（蛋白质疏水性分布图）")
print("=" * 55)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
for fid, data in gravy_cache.items():
    if any(fid.startswith(p) for p in PHAGE_PREFIX):
        pd.Series(data["vals"]).plot.kde(ax=ax, label=fid, alpha=0.8)
ax.set_title("Phage GRAVY Distributions")
ax.set_xlabel("GRAVY value")
ax.set_ylabel("Density")
ax.legend(fontsize=8)
ax.axvline(0, color='gray', linestyle='--', alpha=0.4, linewidth=0.8)

ax = axes[1]
for fid, data in gravy_cache.items():
    if not any(fid.startswith(p) for p in PHAGE_PREFIX):
        pd.Series(data["vals"]).plot.kde(ax=ax, label=fid, alpha=0.8)
ax.set_title("Bacteria GRAVY Distributions")
ax.set_xlabel("GRAVY value")
ax.set_ylabel("Density")
ax.legend(fontsize=7)
ax.axvline(0, color='gray', linestyle='--', alpha=0.4, linewidth=0.8)

plt.suptitle("GRAVY Distribution: Phage vs Bacteria", fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "density_plot.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  已保存 → results/density_plot.png")

print("\n" + "=" * 55)
print("Step 4: 散点图 + Boxplot")
print("=" * 55)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
colors = ["tomato" if yi == 1 else "steelblue" for yi in df["y"]]
ax.scatter(df["x_gravy"], df["y"], c=colors, alpha=0.7, s=60, edgecolors="white")
ax.set_xlabel("|GRAVY_phage − GRAVY_host|")
ax.set_ylabel("y  (1 = infected,  0 = not)")
ax.set_title("GRAVY Feature vs Infection Label")
ax.set_yticks([0, 1])
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="tomato", label="Infected (y=1)"),
                   Patch(color="steelblue", label="Not infected (y=0)")])

ax = axes[1]
infected     = df[df["y"] == 1]["x_gravy"]
not_infected = df[df["y"] == 0]["x_gravy"]
bp = ax.boxplot([infected, not_infected],
                labels=["Infected\n(y=1)", "Not infected\n(y=0)"],
                patch_artist=True)
bp['boxes'][0].set_facecolor("tomato")
bp['boxes'][1].set_facecolor("steelblue")
for box in bp['boxes']:
    box.set_alpha(0.6)
ax.set_ylabel("|GRAVY_phage − GRAVY_host|")
ax.set_title("GRAVY Difference by Infection Status")

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "scatter_boxplot.png"), dpi=150, bbox_inches='tight')
plt.close()
print("  已保存 → results/scatter_boxplot.png")

print("\n" + "=" * 55)
print("Step 5: 统计检验")
print("=" * 55)

u_stat, p_mw = stats.mannwhitneyu(infected, not_infected, alternative="two-sided")
corr, p_corr = stats.pointbiserialr(df["y"], df["x_gravy"])
slope, intercept, r, p_lr, se = stats.linregress(df["x_gravy"], df["y"])

print(f"\n  感染组   x_gravy 均值 = {infected.mean():.4f}")
print(f"  非感染组 x_gravy 均值 = {not_infected.mean():.4f}")
print(f"\n  Mann-Whitney U = {u_stat:.1f}, p = {p_mw:.4e}")
print(f"  Point-biserial r = {corr:.4f}, p = {p_corr:.4e}")
print(f"  线性回归 R² = {r**2:.4f}, p = {p_lr:.4e}, slope = {slope:.4f}")

if p_mw < 0.05:
    print("\n  ★ 两组差异显著（p < 0.05）")
else:
    print("\n  两组差异不显著（p ≥ 0.05），单因子预测力有限")

print("\n" + "=" * 55)
print("Step 6: Logistic Regression")
print("=" * 55)

X = df[["x_gravy"]].values
y_arr = df["y"].values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
clf = LogisticRegression()
clf.fit(X_scaled, y_arr)
y_pred = clf.predict(X_scaled)
y_prob = clf.predict_proba(X_scaled)[:, 1]

print(f"\n  Logistic 系数 = {clf.coef_[0][0]:.4f}")
print(f"  截距        = {clf.intercept_[0]:.4f}")
print(f"\n  Classification Report:")
print(classification_report(y_arr, y_pred))

auc = roc_auc_score(y_arr, y_prob)
print(f"  AUC = {auc:.4f}")

fpr, tpr, _ = roc_curve(y_arr, y_prob)
plt.figure(figsize=(5, 5))
plt.plot(fpr, tpr, color="tomato", lw=2, label=f"AUC = {auc:.3f}")
plt.plot([0,1],[0,1],"k--", alpha=0.4)
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve — GRAVY feature")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "roc_curve.png"), dpi=150)
plt.close()
print("  已保存 → results/roc_curve.png")

print("\n" + "=" * 55)
print("Step 7: 保存文字报告")
print("=" * 55)

report = f"""GRAVY Feature Analysis Report
================================
数据集: {len(df)} 对 phage-host
感染(y=1): {int(df['y'].sum())} 对
不感染(y=0): {int((df['y']==0).sum())} 对

--- 各 phage/bacteria 的平均 GRAVY ---
"""
for fid, data in sorted(gravy_cache.items()):
    report += f"  {fid:<45} mean={data['mean']:.4f}  median={data['median']:.4f}  n={len(data['vals'])}\n"

report += f"""
--- 特征 x_gravy 统计 ---
感染组   mean={infected.mean():.4f}  std={infected.std():.4f}
非感染组 mean={not_infected.mean():.4f}  std={not_infected.std():.4f}

--- 统计检验 ---
Mann-Whitney U = {u_stat:.1f},  p = {p_mw:.4e}
Point-biserial r = {corr:.4f},  p = {p_corr:.4e}
线性回归 R² = {r**2:.4f},  p = {p_lr:.4e},  slope = {slope:.4f}
Logistic 系数 = {clf.coef_[0][0]:.4f},  AUC = {auc:.4f}

--- 解读 ---
slope 为{'负' if slope < 0 else '正'}：x_gravy{'越小感染概率越高' if slope < 0 else '越大感染概率越高'}（{'符合' if slope < 0 else '不符合'}生物学假设）
AUC {'> 0.7，单因子预测力尚可' if auc > 0.7 else ('> 0.5，有微弱预测力，建议结合其他因子' if auc > 0.5 else '≈ 0.5，接近随机，此因子单独预测力弱')}
"""

with open(os.path.join(OUTPUT_DIR, "report.txt"), "w", encoding="utf-8") as f:
    f.write(report)
print("  已保存 → results/report.txt")

print("\n" + "=" * 55)
print("✅ 全部完成！所有结果在 results/ 文件夹")
print("=" * 55)
