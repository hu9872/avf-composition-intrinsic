"""Graphical abstract: 1200 x 1200 px, 300 dpi, 正方形（4 x 4 in）。

内容：上「成熟 vs 失败」组织示意（总细胞数相同，失败组髓系细胞更多、cDC2 数量不变）
      中 分解式（观测 log2FC = 组成项 + 内在项；组成项中位占 51%）
      下 结论句
输出：figures/Graphical_Abstract.png / .pdf
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 12

RED = "#C0392B"
GREY = "#5D6D7E"
DARK = "#1B2631"
LIGHT = "#F4F6F7"

fig = plt.figure(figsize=(4, 4), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")
ax.add_patch(plt.Rectangle((0, 0), 100, 100, facecolor="white", edgecolor="none"))


def tissue(cx, cy, r, n_mye, n_total, seed):
    ax.add_patch(Circle((cx, cy), r, facecolor=LIGHT, edgecolor="0.78", lw=0.8))
    rng = np.random.default_rng(seed)
    pts = []
    while len(pts) < n_total:
        t = rng.uniform(0, 2 * np.pi)
        rr = r * np.sqrt(rng.uniform(0, 0.84))
        x, y = cx + rr * np.cos(t), cy + rr * np.sin(t)
        if all((x - a) ** 2 + (y - b) ** 2 > 2.2 for a, b in pts):
            pts.append((x, y))
    for i, (x, y) in enumerate(pts):
        c = RED if i < n_mye else "#C5CBD1"
        ax.add_patch(Circle((x, y), 1.02, facecolor=c, edgecolor="white", lw=0.35))


tissue(26, 71, 13.0, 8, 34, 11)
tissue(74, 71, 13.0, 16, 34, 12)
ax.text(26, 88.5, "Matured", ha="center", fontsize=12, color=DARK)
ax.text(74, 88.5, "Failed", ha="center", fontsize=12, color=DARK)
ax.text(50, 55.4, "same total cells", ha="center", fontsize=11, color=GREY)
ax.text(50, 50.6, "myeloid 1.8x", ha="center", fontsize=12, color=RED, weight="bold")

ax.text(50, 43.5, "observed log2FC  =  composition  +  intrinsic",
        ha="center", fontsize=11.5, color=DARK)

bx0, bx1, by = 16, 84, 37.2
split = bx0 + (bx1 - bx0) * 0.51
ax.plot([bx0, split], [by, by], lw=6.5, color=RED, solid_capstyle="butt")
ax.plot([split, bx1], [by, by], lw=6.5, color="#AEB6BF", solid_capstyle="butt")
ax.text((bx0 + split) / 2, 31.4, "51% composition", ha="center", fontsize=11, color=RED)
ax.text((split + bx1) / 2, 31.4, "49% intrinsic", ha="center", fontsize=11, color=GREY)
ax.text(50, 26.2, "median across genes", ha="center", fontsize=9.5, color=GREY)

ax.text(50, 19.4, "myofibroblast abundance, fate bias and\nmatrisome programmes: unchanged",
        ha="center", va="center", fontsize=10.5, color="#4D5656")

ax.add_patch(FancyBboxPatch((8, 2.6), 84, 10.4, boxstyle="round,pad=0.4,rounding_size=1.5",
                            facecolor="#FCF3CF", edgecolor="#D4AC0D", lw=0.9))
ax.text(50, 7.8, "Most of the apparent signature\nis tissue composition",
        ha="center", va="center", fontsize=12, color=DARK, weight="bold")

fig.savefig("figures/Graphical_Abstract.png", dpi=300, facecolor="white")
fig.savefig("figures/Graphical_Abstract.pdf", facecolor="white")
print("saved graphical abstract")
