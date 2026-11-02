"""S30: 把握度/可检出效应曲线 + TOST 等效性检验。

目的（对应整改计划 I-05 / 3.2）：
1. 量化 n=6/组在 alpha=0.05、power=0.80 下可检出的最小标准化效应 d；
2. 对关键阴性结果做双单侧检验（TOST），区分「排除大效应」与「证明等效」；
3. 输出 FigS2。

设计要点：TOST 直接读取既有检验表中的 diff / cohens_d / welch p，
不重算上游统计量，保证与正文数字同源。

输出：
  results/s30_power_curve.csv
  results/s30_detectable_effect.csv
  results/s30_tost.csv
  figures/FigS2_power.pdf / .png
"""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"D:/AVF_GSE_2611"
RES = f"{ROOT}/results"
FIG = f"{ROOT}/figures"
N = 6            # 每组供体数
SESOI_D = 0.80   # 最小感兴趣效应（标准化）


def power_two_sample(d, n, alpha=0.05):
    """两独立样本 t 检验（等 n、等方差）的把握度：非中心 t 分布。"""
    df = 2 * n - 2
    nc = d * np.sqrt(n / 2.0)
    tcrit = stats.t.ppf(1 - alpha / 2, df)
    return stats.nct.sf(tcrit, df, nc) + stats.nct.cdf(-tcrit, df, nc)


def detectable_d(n, power=0.80, alpha=0.05):
    lo, hi = 0.01, 6.0
    for _ in range(90):
        mid = (lo + hi) / 2
        if power_two_sample(mid, n, alpha) < power:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def tost_from_summary(diff, d, n1=N, n2=N, sesoi=SESOI_D):
    """由 diff 与 Cohen's d 反推合并 SD，做 TOST（H0: |diff| >= sesoi x SD）。"""
    if not np.isfinite(diff) or not np.isfinite(d) or d == 0:
        return dict(sp=np.nan, se=np.nan, bound=np.nan, tost_p=np.nan,
                    equivalent=False, ci90_low=np.nan, ci90_high=np.nan)
    sp = abs(diff / d)
    se = sp * np.sqrt(1.0 / n1 + 1.0 / n2)
    df = n1 + n2 - 2
    bound = sesoi * sp
    t1 = (diff + bound) / se          # H0: diff <= -bound
    t2 = (diff - bound) / se          # H0: diff >= +bound
    tost_p = max(stats.t.sf(t1, df), stats.t.cdf(t2, df))
    t90 = stats.t.ppf(0.95, df)
    return dict(sp=sp, se=se, bound=bound, tost_p=tost_p,
                equivalent=bool(tost_p < 0.05),
                ci90_low=diff - t90 * se, ci90_high=diff + t90 * se)


def main():
    # --------------------------------------------------------- 1. 把握度曲线
    rows = []
    for n in range(3, 41):
        for d in np.arange(0.2, 4.01, 0.1):
            rows.append(dict(n=n, d=round(float(d), 2),
                             power=round(float(power_two_sample(d, n)), 4)))
    pc = pd.DataFrame(rows)
    pc.to_csv(f"{RES}/s30_power_curve.csv", index=False)

    det = pd.DataFrame([dict(n=n, alpha=0.05, power=0.80,
                             detectable_d=round(detectable_d(n), 3))
                        for n in [3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 26, 30, 40]])
    det.to_csv(f"{RES}/s30_detectable_effect.csv", index=False)
    print(det.to_string(index=False))
    d6 = float(det.loc[det.n == 6, "detectable_d"].iloc[0])

    # ------------------------------------------------------------- 2. TOST
    specs = []
    comp = pd.read_csv(f"{RES}/s5_composition_tests.csv").set_index("state")
    for st, lab in [("Myofibroblast", "Myofibroblast abundance (CLR)"),
                    ("AF_resident", "Resident AF abundance (CLR)"),
                    ("Pericyte", "Pericyte abundance (CLR)")]:
        specs.append((lab, float(comp.loc[st, "diff"]), float(comp.loc[st, "cohens_d"]),
                      float(comp.loc[st, "welch_p"])))

    fate = pd.read_csv(f"{RES}/s7_fate_tests.csv").set_index("fate")
    for st, lab in [("Myofibroblast", "Myofibroblast fate absorption"),
                    ("SMC_contractile", "Contractile SMC fate absorption"),
                    ("AF_resident", "Resident AF fate absorption")]:
        specs.append((lab, float(fate.loc[st, "diff"]), float(fate.loc[st, "cohens_d"]),
                      float(fate.loc[st, "welch_p"])))

    ecm = pd.read_csv(f"{RES}/s9_ecm_tests.csv").set_index("metric")
    for st in ["ECM1", "ECM2", "ECM3", "ECM4", "ECM5", "zonation"]:
        lab = "Zonation index" if st == "zonation" else f"Matrisome programme {st}"
        specs.append((lab, float(ecm.loc[st, "diff"]), float(ecm.loc[st, "cohens_d"]),
                      float(ecm.loc[st, "p"])))

    abn = pd.read_csv(f"{RES}/s5e_abundance_tests.csv").set_index("metric")
    for st, lab in [("myeloid_per1000", "Myeloid cells per 1,000 (positive control)"),
                    ("cdc2_D1_per1000", "cDC2 cells per 1,000")]:
        specs.append((lab, float(abn.loc[st, "diff"]), float(abn.loc[st, "cohens_d"]),
                      float(abn.loc[st, "welch_p"])))

    out = []
    for lab, diff, d, p in specs:
        r = tost_from_summary(diff, d)
        sp = r["sp"]
        excl = max(abs(r["ci90_low"]), abs(r["ci90_high"])) / sp if sp > 0 else np.nan
        out.append(dict(metric=lab, diff=round(diff, 4), cohens_d=round(d, 3),
                        welch_p=round(p, 4),
                        equivalence_bound=round(r["bound"], 4),
                        ci90_low=round(r["ci90_low"], 4),
                        ci90_high=round(r["ci90_high"], 4),
                        tost_p=round(r["tost_p"], 4),
                        equivalent_at_d080=r["equivalent"],
                        excluded_effect_from_d=round(excl, 3),
                        detectable_d_at_n6=round(d6, 3),
                        n_failed=N, n_mature=N))
    t = pd.DataFrame(out)
    t.to_csv(f"{RES}/s30_tost.csv", index=False)
    print("\n=== TOST (SESOI = d 0.80) ===")
    print(t.to_string(index=False))
    print(f"\nn=6 可检出 d = {d6:.2f} (80% power, alpha=0.05)；"
          f"等效性成立项数 = {int(t.equivalent_at_d080.sum())}/{len(t)}")

    # ------------------------------------------------------------- 3. 出图
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.5), dpi=300)
    ax = axes[0]
    for n, c in [(6, "#C0392B"), (12, "#E67E22"), (26, "#2E86C1")]:
        s = pc[pc.n == n]
        ax.plot(s.d, s.power, color=c, lw=1.5, label=f"n = {n} per group")
    ax.axhline(0.8, ls=":", color="0.4", lw=1)
    ax.axvline(d6, ls="--", color="#C0392B", lw=1)
    ax.set_xlabel("Standardised effect size (Cohen's d)", fontsize=8.5)
    ax.set_ylabel("Power", fontsize=8.5)
    ax.tick_params(labelsize=8)
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    ax.annotate(f"d = {d6:.1f}", xy=(d6, 0.8), xytext=(d6 + 0.3, 0.5),
                fontsize=7.5, color="#C0392B")
    ax.set_title("Detectable effect at 80% power", fontsize=9, loc="left")

    ax = axes[1]
    y = np.arange(len(t))
    lo = t["diff"] - t["ci90_low"]
    hi = t["ci90_high"] - t["diff"]
    ax.errorbar(t["diff"], y, xerr=[lo, hi], fmt="o", ms=3.2, lw=1,
                color="#1A5276", ecolor="#5D6D7E", capsize=2, zorder=3)
    for yy, b in zip(y, t["equivalence_bound"]):
        ax.plot([-b, b], [yy, yy], color="#E67E22", lw=4, alpha=0.25,
                solid_capstyle="butt", zorder=1)
    ax.axvline(0, color="0.3", lw=0.8, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels([m.replace(" (CLR)", "").replace(" (positive control)", "")
                        for m in t["metric"]], fontsize=6.2)
    ax.invert_yaxis()
    ax.tick_params(axis="x", labelsize=8)
    ax.set_xlabel("Difference, failed − matured (90% CI)", fontsize=8.5)
    ax.set_title("Effect bounds (SESOI d = 0.8)", fontsize=9, loc="left")
    for a, lab in zip(axes, ["A", "B"]):
        a.text(-0.16, 1.08, lab, transform=a.transAxes, fontsize=11,
               fontweight="bold", va="top", ha="left")
    fig.tight_layout()
    fig.savefig(f"{FIG}/FigS2_power.pdf", bbox_inches="tight")
    fig.savefig(f"{FIG}/FigS2_power.png", bbox_inches="tight", dpi=300)
    print("saved FigS2")


if __name__ == "__main__":
    main()
