"""Rebuild the report figures from results/*.csv.  Usage: python scripts/make_figures.py"""
import pandas as pd, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def parse(label):
    m = re.match(r"(short|long)_(\d+)x(-rep\d+)?$", label)
    if m:
        return m.group(1) + ("-rep" if m.group(3) else ""), int(m.group(2))
    m = re.match(r"(hybridspades|polish)_s(\d+)_l(\d+)$", label)
    return m.group(1), int(m.group(3))

q = pd.read_csv("results/assembly_summary_v2.csv")
q = q[q["label"] != "long_30x"].copy()          # failed run (duplicated contigs), excluded from the curves
q["strategy"] = q["label"].apply(lambda s: parse(s)[0])
q["cov"] = q["label"].apply(lambda s: parse(s)[1])
metrics = [("N50", "N50 (bp, log scale)", True), ("# misassemblies", "# misassemblies", False),
           ("Genome fraction (%)", "Genome fraction (%)", False), ("# indels per 100 kbp", "Indels per 100 kbp (log scale)", True)]
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
for ax, (m, title, log) in zip(axes.ravel(), metrics):
    for s, g in q.groupby("strategy"):
        g = g.sort_values("cov")
        ax.plot(g["cov"], pd.to_numeric(g[m], errors="coerce"), "o-" if len(g) > 1 else "s", label=s)
    ax.set_xlabel("Coverage (x); hybrids plotted at their long-read coverage")
    ax.set_ylabel(title); ax.set_title(title)
    if log: ax.set_yscale("log")
    ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig("figures/fig1_assembly_metrics.png", dpi=200)

h = pd.read_csv("results/indel_context.csv")
cols = ["hp1", "hp2", "hp3", "hp4", "hp5", "hp6plus"]
share = h[cols].div(h[cols].sum(axis=1), axis=0) * 100
fig, ax = plt.subplots(figsize=(8, 4.5))
share.index = h["assembly"]
share.plot(kind="barh", stacked=True, ax=ax, colormap="viridis")
ax.set_xlabel("% of indels by homopolymer run length at the site (6 = 6 or more)")
ax.set_title("Indel errors vs reference: homopolymer context (about 6% of genome positions lie in runs >= 4)")
ax.title.set_fontsize(9)
fig.tight_layout(); fig.savefig("figures/fig2_indel_homopolymers.png", dpi=200)

m = pd.read_csv("results/missing_regions.csv")
m = m[m["assembly"] != "long_30x"]
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.barh(m["assembly"], m["missing_pct"])
for y, (p, g) in enumerate(zip(m["missing_pct"], m["GC_missing"])):
    ax.text(p, y, f"  GC of missing = {g}", va="center", fontsize=7)
ax.set_xlabel("% of reference not covered by the assembly (genome GC = 0.508)")
fig.tight_layout(); fig.savefig("figures/fig3_missing_regions.png", dpi=200)
print("figures written")
