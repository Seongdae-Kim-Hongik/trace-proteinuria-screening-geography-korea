# 99G_nbkid_judgement_gate.py - Judgement check. Compares examinees judged 'normal B, kidney-function management' (KOSIS N102 national by sex and age; N101 by province)
#  with examinees whose urine dipstick protein was trace (KOSIS N065, N064) by year, sex, age and province, with the ratio of denominators
#  (N099/N098 versus the dipstick totals); also the year-to-year rank correlation of district rates.
#  Output: 04_분석/NBKID/results_gate/S0_judgement_national.csv, S0_judgement_age_sex.csv, S0_judgement_province.csv, S0_district_rankcorr.csv
import os, json, collections, statistics as st
import pandas as pd, numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); A = os.path.join(ROOT, "04_분석", "NBKID"); RS = os.path.join(A, "results_gate"); os.makedirs(RS, exist_ok=True)
SRC = os.path.join(ROOT, "00_원천", "C_NHIS_검진보험")
num = lambda x: float(str(x).replace(",", "")) if x not in ("", "-", None) else None
SXM = {"남자": "M", "여자": "F", "계": "T", "합계": "T"}
def nb_band(a):
    a = str(a).replace(" ", "")
    if a in ("19세이하", "85세이상", "계"): return {"19세이하": "0-19", "85세이상": "85+", "계": "계"}[a]
    if a.endswith("세") and "~" in a: lo, hi = a[:-1].split("~"); return f"{lo}-{hi}"
    return None
U = collections.defaultdict(float)
for p in [os.path.join(SRC, "KOSIS_NBKID_261121", "DT_35007_N065_2011_2011.json"), os.path.join(SRC, "KOSIS_JN_RR_261119", "DT_35007_N065_age_sex.json")]:
    for r in json.load(open(p)):
        v = num(r["DT"]); b = nb_band(r["C1_NM"])
        if v is None or b is None: continue
        U[(int(r["PRD_DE"]), SXM[r["C2_NM"]], b, r["C3_NM"])] += v
N102 = pd.read_csv(os.path.join(SRC, "DT_35007_N102_national_long.csv"), encoding="utf-8-sig", dtype={"grp_code": str}); N102["sx"] = N102.sex.map(SXM); N102["band"] = N102.age.map(nb_band)
N099 = pd.read_csv(os.path.join(SRC, "DT_35007_N099_national_long.csv"), encoding="utf-8-sig"); N099["sx"] = N099.sex.map(SXM); N099["band"] = N099.age.map(nb_band)
nk = N102[N102.grp_code == "00307"].groupby(["year", "sx", "band"]).value.sum(); jd = N099[N099.grp == "계"].groupby(["year", "sx", "band"]).value.sum()
rows = []
for (y, sx, b), f in nk.items():
    k = (int(y), sx, b)
    if k + ("계",) not in U or (y, sx, b) not in jd.index: continue
    tot = U[k + ("계",)]; tr = U[k + ("+-",)]; pos = sum(U[k + (c,)] for c in ["+1", "+2", "+3", "+4"]); ex = jd[(y, sx, b)]
    rows.append(dict(year=int(y), sex=sx, band=b, judged=int(ex), dip_measured=int(tot), denom_ratio=round(ex / tot, 3), nbkid=int(f), pct_nbkid=round(100 * f / ex, 3),
                     pct_trace=round(100 * tr / tot, 3), pct_1plus=round(100 * pos / tot, 3), ratio_nbkid_to_trace=round((f / ex) / (tr / tot), 3) if tr else None,
                     ratio_counts=round(f / tr, 3) if tr else None))
AS = pd.DataFrame(rows).sort_values(["year", "sex", "band"]); AS.to_csv(os.path.join(RS, "S0_judgement_age_sex.csv"), index=False, encoding="utf-8-sig")
NAT = AS[AS.band == "계"].drop(columns="band"); NAT.to_csv(os.path.join(RS, "S0_judgement_national.csv"), index=False, encoding="utf-8-sig")
N101 = pd.read_csv(os.path.join(SRC, "DT_35007_N101.csv"), encoding="utf-8-sig", dtype=str); N101 = N101[N101.C1.str.len() == 12].copy(); N101["v"] = N101.DT.map(num); N101["sx"] = N101.C2_NM.map(SXM)
N098 = pd.read_csv(os.path.join(SRC, "DT_35007_N098.csv"), encoding="utf-8-sig", dtype=str); N098 = N098[N098.C1.str.len() == 12].copy(); N098["v"] = N098.DT.map(num); N098["sx"] = N098.C2_NM.map(SXM)
F = N101[N101.C3_NM == "신장기능관리"].groupby(["PRD_DE", "C1", "sx"]).v.sum(); J = N098[N098.C3_NM == "계"].groupby(["PRD_DE", "C1", "sx"]).v.sum()
P64 = collections.defaultdict(float); NM = {}
for r in json.load(open(os.path.join(SRC, "KOSIS_NBKID_261121", "DT_35007_N064_2011_2024.json"))):
    v = num(r["DT"])
    if v is None or r["C1_NM"] == "계": continue
    P64[(r["PRD_DE"], r["C1"], SXM[r["C2_NM"]], r["C3_NM"])] += v; NM[r["C1"]] = r["C1_NM"]
prov, pl = [], []
for (y, sx), g in F.groupby(level=[0, 2]):
    fp, tp, rr = [], [], []
    for (yy, c1, _), f in g.items():
        k = (yy, c1, sx)
        if k + ("계",) not in P64 or k not in J.index: continue
        tot = P64[k + ("계",)]; tr = P64[k + ("+-",)]; a = f / J[k]; b = tr / tot
        fp.append(100 * a); tp.append(100 * b); rr.append(a / b); pl.append(dict(year=int(y), sex=sx, province=NM.get(c1, c1), pct_nbkid=round(100 * a, 3), pct_trace=round(100 * b, 3), ratio=round(a / b, 3), denom_ratio=round(J[k] / tot, 3)))
    if len(rr) < 10: continue
    prov.append(dict(year=int(y), sex=sx, n_provinces=len(rr), r_nbkid_vs_trace=round(float(np.corrcoef(fp, tp)[0, 1]), 3), spearman=round(float(pd.Series(fp).corr(pd.Series(tp), method="spearman")), 3),
                     ratio_median=round(st.median(rr), 3), ratio_min=round(min(rr), 3), ratio_max=round(max(rr), 3), ratio_cv=round(st.pstdev(rr) / st.mean(rr), 3),
                     nbkid_pct_range=f"{min(fp):.2f}-{max(fp):.2f}", trace_pct_range=f"{min(tp):.2f}-{max(tp):.2f}"))
PRV = pd.DataFrame(prov).sort_values(["sex", "year"]); PRV.to_csv(os.path.join(RS, "S0_judgement_province.csv"), index=False, encoding="utf-8-sig")
pd.DataFrame(pl).to_csv(os.path.join(RS, "S0_judgement_province_long.csv"), index=False, encoding="utf-8-sig")
E = pd.read_csv(os.path.join(ROOT, "02_패널", "NBKID_Estd3.csv"), encoding="utf-8-sig"); E["r"] = E.cases / E.examinees
W = E.pivot(index="fbd_region", columns="year", values="r")
rc = pd.DataFrame([dict(pair=f"{a}-{b}", spearman=round(W[a].corr(W[b], method="spearman"), 3)) for a, b in zip(W.columns[:-1], W.columns[1:])] + [dict(pair="2011-2024", spearman=round(W[2011].corr(W[2024], method="spearman"), 3))])
rc.to_csv(os.path.join(RS, "S0_district_rankcorr.csv"), index=False, encoding="utf-8-sig")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(NAT[NAT.sex == "T"].to_string()); print(PRV[PRV.sex == "T"].to_string()); print(rc.T.to_string())
a20 = AS[(AS.sex != "T") & (AS.band != "계")]; print("[연령·성 층] ratio_nbkid_to_trace 연도별 범위:"); print(a20.groupby("year").ratio_nbkid_to_trace.agg(["min", "median", "max"]).round(3).to_string())
