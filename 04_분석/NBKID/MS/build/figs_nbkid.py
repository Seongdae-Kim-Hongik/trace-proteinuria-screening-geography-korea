# -*- coding: utf-8 -*-
# figs_nbkid.py - Draws Figures 1-3 (district maps, national trend, period maps with emerging hot spots, model-by-model rate ratios).
#  Input: numbers_nbkid.json, gi_nbkid.csv, results_std3, boundary file ST_SHP. Output: figs/Fig1-Fig3 (PNG, TIFF)
import os, json, numpy as np, pandas as pd, geopandas as gpd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch, Rectangle, Polygon
from matplotlib.lines import Line2D
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
SHP = os.environ.get("ST_SHP", os.path.join("data", "boundary", "final.shp"))
HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.join(HERE, "figs"); OUTF = os.path.join(os.path.dirname(HERE), "out", "Figures"); os.makedirs(FIG, exist_ok=True); os.makedirs(OUTF, exist_ok=True)
RES = os.path.abspath(os.path.join(HERE, "..", "..", "results_std3"))
N = json.load(open(os.path.join(HERE, "numbers_nbkid.json"))); GI = pd.read_csv(os.path.join(HERE, "gi_nbkid.csv"), dtype={"sigungu_cd": str})
MM = 1 / 25.4; W = 174 * MM; LW = 0.25; PANEL_LABEL_PT = 56
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "axes.labelsize": 12, "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 9.5, "legend.title_fontsize": 10, "axes.linewidth": 1.0, "lines.linewidth": 1.4})
G = gpd.read_file(SHP).set_crs("EPSG:5179");
G["cd"] = G.sigungu_cd.astype(str); G["sido"] = G.cd.str[:2]; SIDO = G.dissolve("sido"); XL = (G.total_bounds[0], G.total_bounds[2]); YL = (G.total_bounds[1], G.total_bounds[3])
def letter(fig, ax, ch, dx=0.0, dy=0.0):
    b = ax.get_position(); fig.text(b.x0 + dx, b.y1 + dy, ch.upper(), fontsize=PANEL_LABEL_PT, fontweight="bold", va="bottom", ha="left")
def mapax(ax, arrow=True):
    SIDO.boundary.plot(ax=ax, color="#707070", linewidth=0.5); ax.set_xlim(*XL); ax.set_ylim(*YL); ax.set_axis_off(); ax.set_aspect("equal")
    x0, x1 = XL; y0, y1 = YL; W_ = x1 - x0; H_ = y1 - y0; bx, by = x0 + 0.02 * W_, y0 + 0.11 * H_; hh = 0.012 * H_
    ax.add_patch(Rectangle((bx, by), 50000, hh, facecolor="black", edgecolor="black", linewidth=0.8)); ax.add_patch(Rectangle((bx + 50000, by), 50000, hh, facecolor="white", edgecolor="black", linewidth=0.8))
    for xx, lab in ((0, "0"), (50000, "50"), (100000, "100")): ax.text(bx + xx, by + hh + 0.008 * H_, lab, ha="center", va="bottom", fontsize=7)
    ax.text(bx + 108000, by + 0.5 * hh, "km", ha="left", va="center", fontsize=7)
    if arrow:
        nx, ny = x0 + 0.09 * W_, y0 + 0.52 * H_
        ax.add_patch(Polygon([[nx - 0.045 * W_, ny], [nx + 0.045 * W_, ny], [nx, ny + 0.24 * H_]], closed=True, facecolor="black", edgecolor="none")); ax.text(nx, ny + 0.255 * H_, "N", ha="center", va="bottom", fontsize=24, fontweight="bold")
def save(fig, k, dpi):
    p = os.path.join(FIG, f"Fig{k}_preview.png"); fig.savefig(p, dpi=dpi, facecolor="white"); plt.close(fig)
    im = Image.open(p).convert("RGB"); im.save(p, dpi=(dpi, dpi)); im.save(os.path.join(OUTF, f"Fig{k}.tif"), dpi=(dpi, dpi), compression="tiff_lzw"); return im.size
BR = [0, 0.80, 0.90, 1.00, 1.11, 1.25, 9]; LB = ["<0.80", "0.80\u20130.89", "0.90\u20130.99", "1.00\u20131.10", "1.11\u20131.24", "≥1.25"]
PAL = ["#2166AC", "#67A9CF", "#D1E5F0", "#FDDBC7", "#EF8A62", "#B2182B"]

# Fig 1
fig = plt.figure(figsize=(W, 5.6))
ax = fig.add_axes([0.0, 0.16, 0.54, 0.74]); gg = G.merge(GI[["sigungu_cd", "SIR", "cls"]], left_on="cd", right_on="sigungu_cd", how="left"); assert gg.SIR.notna().all() and len(gg) == 229
gg.plot(ax=ax, column="SIR", cmap=ListedColormap(PAL), norm=BoundaryNorm(BR, len(PAL)), edgecolor="#9a9a9a", linewidth=LW); mapax(ax)
gg[gg.cls == "cold"].boundary.plot(ax=ax, color="black", linewidth=0.9, linestyle=(0, (1.5, 1.0))); gg[gg.cls == "hot"].boundary.plot(ax=ax, color="black", linewidth=1.6); letter(fig, ax, "A")
h = [Patch(facecolor=c, edgecolor="#777777", linewidth=0.5, label=l) for c, l in zip(PAL, LB)] + [Patch(facecolor="none", edgecolor="black", linewidth=1.6, label="Hot spot"), Patch(facecolor="none", edgecolor="black", linewidth=0.9, linestyle=(0, (1.5, 1.0)), label="Cold spot")]
ax.legend(handles=h, title="Standardised ratio, 2011\u20132024", loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=4, frameon=False, columnspacing=0.7, handlelength=1.2)
b = fig.add_axes([0.66, 0.24, 0.32, 0.56]); yb = sorted(int(k) for k in N["rate"]); GT = N["gate"]
b.plot(yb, [N["rate"][str(y)] for y in yb], color="#B2182B", marker="o", ms=3.2, label="Trace judgement")
b.plot(yb, [GT["tr"][str(y)] for y in yb], color="#2166AC", marker="s", ms=3.0, ls="--", label="Trace in dipstick distribution")
b.set_ylim(0, 7); b.set_xticks([2011, 2014, 2017, 2020, 2024]); b.set_ylabel("% of examinees"); b.set_xlabel("Year")
b.legend(frameon=False, fontsize=7.5, loc="upper left", handlelength=1.8); b.spines[["top", "right"]].set_visible(False); letter(fig, b, "B", dx=-0.09)
s1 = save(fig, 1, 600)

WS = pd.read_csv(os.path.join(RES, "S4b_window_spatialRR.csv"), dtype={"sigungu_cd": str}); EH = pd.read_csv(os.path.join(RES, "S11_ehsa.csv"), dtype={"sigungu_cd": str})
fig = plt.figure(figsize=(W, 9.4)); pos = [[0.02, 0.585, 0.46, 0.355], [0.52, 0.585, 0.46, 0.355], [0.02, 0.175, 0.46, 0.355], [0.52, 0.175, 0.46, 0.355]]
RBR = BR
for i, (w, lab) in enumerate((("W1_2011-15", "A"), ("W2_2016-19", "B"), ("W3_2020-24", "C"))):
    a = fig.add_axes(pos[i]); d = G.merge(WS[WS.window == w][["sigungu_cd", "RR"]], left_on="cd", right_on="sigungu_cd", how="left"); assert d.RR.notna().all()
    d.plot(ax=a, column="RR", cmap=ListedColormap(PAL), norm=BoundaryNorm(RBR, len(PAL)), edgecolor="#9a9a9a", linewidth=LW); mapax(a, arrow=(i == 0)); letter(fig, a, lab)
fig.legend(handles=[Patch(facecolor=c, edgecolor="#777777", linewidth=0.5, label=l) for c, l in zip(PAL, LB)], title="Spatial relative risk", loc="lower center", bbox_to_anchor=(0.25, 0.005), ncol=4, frameon=False, columnspacing=0.6, handlelength=1.1)
a = fig.add_axes(pos[3]); CAT = {"Intensifying hot": "#67000D", "Persistent hot": "#CB181D", "Consecutive hot": "#FB6A4A", "Sporadic hot": "#FCBBA1", "New cold": "#C6DBEF", "Consecutive cold": "#6BAED6", "Persistent cold": "#08519C", "Sporadic cold": "#DEEBF7", "No pattern": "#F2F2F2"}
d = G.merge(EH[["sigungu_cd", "EHSA"]], left_on="cd", right_on="sigungu_cd", how="left"); assert d.EHSA.notna().all() and set(d.EHSA) <= set(CAT), set(d.EHSA) - set(CAT)
d.plot(ax=a, color=d.EHSA.map(CAT), edgecolor="#9a9a9a", linewidth=LW); mapax(a, arrow=False); letter(fig, a, "D")
fig.legend(handles=[Patch(facecolor=c, edgecolor="#777777", linewidth=0.5, label=f"{k} ({N['ehsa'][k]})") for k, c in CAT.items()], loc="lower center", bbox_to_anchor=(0.75, 0.005), ncol=2, frameon=False, columnspacing=0.6, handlelength=1.1, fontsize=7.2)
s2 = save(fig, 2, 600)

LAB = {"walking_rate": "Walking", "basic_livelihood_rate": "Public assistance\nrecipients"}
ORD = list(LAB); BW = lambda v, p: N["bw"].get(f"{v}|{p}")
MOD = [("Knorr-Held type I", lambda v: N["KH_typeI"][v], "#2166AC", "^"), ("Knorr-Held type II", lambda v: N["KH_typeII"][v], "#2166AC", "P"),
       ("Knorr-Held type III", lambda v: N["KH_typeIII"][v], "#2166AC", "*"), ("Knorr-Held type IV", lambda v: N["KH_typeIV"][v], "#2166AC", "D"),
       ("Restricted spatial regression", lambda v: N["RSR"][v], "#B2182B", "X"), ("Spatial+", lambda v: N["D1_spatial_plus"][v], "#B2182B", "8"),
       ("Between districts", lambda v: BW(v, "between"), "#4D9221", "<"), ("Within districts", lambda v: BW(v, "within"), "#4D9221", ">"),
       ("Suspected kidney disease", lambda v: N["nco"]["KIDS|" + v], "#8C6BB1", "s"), ("Borderline blood pressure", lambda v: N["nco"]["NBBP|" + v], "#8C6BB1", "v")]
from decimal import Decimal as _Dec, ROUND_HALF_UP as _HU
H3 = lambda x: str(_Dec(repr(float(x))).quantize(_Dec("0.001"), rounding=_HU))
vals = [t for v in ORD for j in range(len(MOD)) for t in [MOD[j][1](v)] if t] + [N["main"][v] for v in ORD]; lo_all = min(t[1] for t in vals) - 0.03; hi_all = max(t[2] for t in vals) + 0.03
fig, ax = plt.subplots(figsize=(W, 5.6)); fig.subplots_adjust(left=0.22, right=0.98, top=0.98, bottom=0.40)
off = np.linspace(-0.38, 0.38, len(MOD))
T1 = pd.read_csv(os.path.join(RES, "S9_typeI_IRR.csv")); T1C = dict(zip(T1["cov"], T1["cred"]))
VARS_N = {"walking_rate": "walking_practice_rate", "basic_livelihood_rate": "basic_livelihood_rate"}
for i, v in enumerate(ORD):
    for (lab, f, col, mk), o in zip(MOD, off):
        r = f(v)
        if not r: continue
        e, lo, hi = r; cr = lo > 1 or hi < 1 or (lab == "Knorr-Held type I" and bool(T1C.get(VARS_N[v]))); ax.plot([lo, hi], [i + o] * 2, color=col, lw=1.1); ax.plot(e, i + o, marker=mk, mfc=col if cr else "white", mec=col, ms=5, mew=1.0, ls="none")
        if e >= 1: ax.text(hi + 0.002, i + o, H3(e), ha="left", va="center", fontsize=7, color=col, bbox=dict(fc="white", ec="none", pad=0.4), zorder=5)
        else: ax.text(lo - 0.002, i + o, H3(e), ha="right", va="center", fontsize=7, color=col, bbox=dict(fc="white", ec="none", pad=0.4), zorder=5)
    pm = N["main"][v][0]; ax.plot([pm, pm], [i - 0.46, i + 0.46], color="#000000", lw=1.1, ls="--")
ax.axvline(1, color="#404040", lw=1.0, ls=":"); ax.set_yticks(range(len(ORD))); ax.set_yticklabels([LAB[v] for v in ORD]); ax.invert_yaxis()
a0 = np.floor((lo_all - 0.004) * 50) / 50; a1 = np.ceil((hi_all + 0.004) * 50) / 50; xt = list(np.round(np.arange(a0, a1 + 0.0001, 0.02), 2)); xt = xt if len(xt) <= 12 else list(np.round(np.arange(np.floor(a0 * 25) / 25, a1 + 0.04 - 0.0001, 0.04), 2))
ax.set_xlim(xt[0], xt[-1]); ax.set_xticks(xt); ax.set_xticklabels([f"{x:.2f}" for x in xt]); ax.minorticks_off()
ax.set_xlabel("Rate ratio per 1 standard deviation (95% credible interval)"); ax.spines[["top", "right"]].set_visible(False)
ax.legend(handles=[Line2D([0], [0], color=c_, marker=m, lw=1.1, ms=5, mfc=c_, mec=c_, label=l) for l, _, c_, m in MOD] + [Line2D([0], [0], color="#000000", lw=1.1, ls="--", label="Principal model (Table 2)"), Line2D([0], [0], color="#404040", marker="o", lw=0, ms=5, mfc="#404040", label="Filled: 95% CrI excludes 1"), Line2D([0], [0], color="#404040", marker="o", lw=0, ms=5, mfc="white", mec="#404040", label="Open: 95% CrI includes 1")], loc="upper center", bbox_to_anchor=(0.40, -0.16), ncol=3, frameon=False)
s3 = save(fig, 3, 600)
for k, s in ((1, s1), (2, s2), (3, s3)): print(f"Fig{k}: {s[0]}x{s[1]} px = {s[0]/600*25.4:.0f} x {s[1]/600*25.4:.0f} mm @ 600 dpi")

FT = pd.read_csv(os.path.join(RES, "S20_fitted_main.csv")); FT["mu"] = FT.E * FT.fitted_RR
A = FT.groupby("fbd_region").agg(cases=("cases", "sum"), E=("E", "sum"), mu=("mu", "sum")).reset_index(); A["obs"] = A.cases / A.E; A["fit"] = A.mu / A.E; A["res"] = A.cases / A.mu
GG = G.merge(GI[["sigungu_cd"]].assign(fbd_region=pd.read_csv(os.path.join(HERE, "gi_nbkid.csv")).fbd_region), left_on="cd", right_on="sigungu_cd", how="left").merge(A, on="fbd_region", how="left"); assert GG.obs.notna().all()
RBR2 = [0, 0.95, 0.98, 1.00, 1.02, 1.05, 9]; RLB = ["<0.95", "0.95\u20130.97", "0.98\u20130.99", "1.00\u20131.01", "1.02\u20131.04", "≥1.05"]
fig = plt.figure(figsize=(W, 4.6)); P3 = [[0.00, 0.20, 0.33, 0.66], [0.335, 0.20, 0.33, 0.66], [0.67, 0.20, 0.33, 0.66]]
for k, (col, br, lb, ch) in enumerate((("obs", BR, LB, "A"), ("fit", BR, LB, "B"), ("res", RBR2, RLB, "C"))):
    a = fig.add_axes(P3[k]); GG.plot(ax=a, column=col, cmap=ListedColormap(PAL), norm=BoundaryNorm(br, len(PAL)), edgecolor="#9a9a9a", linewidth=0.15); mapax(a, arrow=(k == 0)); letter(fig, a, ch)
    a.legend(handles=[Patch(facecolor=c, edgecolor="#777777", linewidth=0.5, label=l) for c, l in zip(PAL, lb)], loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=2, frameon=False, fontsize=8, handlelength=1.0, columnspacing=0.5)
sS1 = save(fig, "S1", 600)
