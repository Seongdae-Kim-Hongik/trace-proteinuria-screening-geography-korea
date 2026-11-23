#!/usr/bin/env python3
# person_axis_nbkid.py: national sex and age description of the trace judgement. Source tables as for the expected counts (KOSIS DT_35007_N099 examinees judged = denominator, DT_35007_N102 group 00307 = numerator), 2 sexes x 15 age bands.
#  (1) sex x age-group (5 groups) rates and rate ratios (reference women 45-64 years, Poisson ratio within period, Wald 95% CI); (2) rates per 100,000 examinees by 15 age bands and period; (3) crude and directly age-standardised rates.
#  Standard = 2015 national examinee structure (2 sexes x 15 age bands). No model fitting. Output: checks/out/person_axis_*.csv, MS/build/person_axis_nbkid.json
import os, json, numpy as np, pandas as pd
def keep(x, nd=None):   # values are stored unrounded; rounding happens once (Decimal half-up) when the manuscript is built
    return float(x)
HERE = os.path.dirname(os.path.abspath(__file__)); NB = os.path.dirname(HERE); SK = os.path.dirname(os.path.dirname(NB))
SRC = os.path.join(SK, "00_원천", "C_NHIS_검진보험"); OUT = os.path.join(HERE, "out")
YEARS = list(range(2011, 2025)); BANDS = ["0-19"] + [f"{lo}-{lo+4}" for lo in range(20, 85, 5)] + ["85+"]; SEXES = {"남자": "Men", "여자": "Women"}
PER = {"2011-2015": (2011, 2015), "2016-2019": (2016, 2019), "2020-2024": (2020, 2024), "2011-2024": (2011, 2024)}
GRP = {"0-19": "<30", "20-24": "<30", "25-29": "<30", "30-34": "30-44", "35-39": "30-44", "40-44": "30-44", "45-49": "45-64", "50-54": "45-64",
       "55-59": "45-64", "60-64": "45-64", "65-69": "65-74", "70-74": "65-74", "75-79": "75+", "80-84": "75+", "85+": "75+"}
GORD = ["<30", "30-44", "45-64", "65-74", "75+"]; REF = ("Women", "45-64")
def nb_band(a):
    a = str(a).replace(" ", "")
    if a == "19세이하": return "0-19"
    if a == "85세이상": return "85+"
    if a.endswith("세") and "~" in a: lo, hi = a[:-1].split("~"); return f"{lo}-{hi}"
A = pd.read_csv(os.path.join(SRC, "DT_35007_N099_national_long.csv"), encoding="utf-8-sig"); A = A[(A.grp == "계") & A.year.isin(YEARS)]
B = pd.read_csv(os.path.join(SRC, "DT_35007_N102_national_long.csv"), encoding="utf-8-sig", dtype={"grp_code": str}); B = B[(B.grp_code == "00307") & B.year.isin(YEARS)]
for X in (A, B): X["band"] = X.age.map(nb_band)
EX = A[A.band.notna() & A.sex.isin(SEXES)].groupby(["year", "sex", "band"]).value.sum()
NU = B[B.band.notna() & B.sex.isin(SEXES)].groupby(["year", "sex", "band"]).value.sum().reindex(EX.index).fillna(0)
# 원표 합계 대조
for y in YEARS:
    assert abs(EX.loc[y].sum() - A[(A.year == y) & (A.age == "계") & (A.sex == "합계")].value.sum()) <= 2, y
    assert abs(NU.loc[y].sum() - B[(B.year == y) & (B.age == "계") & (B.sex == "합계")].value.sum()) <= 2, y
D = pd.DataFrame({"exam": EX, "case": NU}).reset_index(); D["sex_en"] = D.sex.map(SEXES); D["agegrp"] = D.band.map(GRP)
def per(d, p): lo, hi = PER[p]; return d[(d.year >= lo) & (d.year <= hi)]
# ① 성×연령군 판정률·IRR
r1 = []
for p in PER:
    g = per(D, p).groupby(["sex_en", "agegrp"])[["case", "exam"]].sum()
    ref = g.loc[REF]
    for (s, a), v in g.iterrows():
        irr = (v.case / v.exam) / (ref.case / ref.exam); se = np.sqrt(1 / v.case + 1 / ref.case) if (s, a) != REF else 0.0
        r1.append(dict(period=p, sex=s, agegrp=a, cases=int(v.case), examinees=int(v.exam), rate_pct=keep(100 * v.case / v.exam, 3),
                       rate_per100k=keep(1e5 * v.case / v.exam, 1), IRR=keep(irr, 3), lo=keep(irr * np.exp(-1.96 * se), 3), hi=keep(irr * np.exp(1.96 * se), 3),
                       reference=f"{REF[0]} {REF[1]}"))
T1 = pd.DataFrame(r1); T1["agegrp"] = pd.Categorical(T1.agegrp, GORD); T1 = T1.sort_values(["period", "sex", "agegrp"])
# 성 비(남/여), 기간별
r1b = []
for p in PER:
    g = per(D, p).groupby("sex_en")[["case", "exam"]].sum(); m, w = g.loc["Men"], g.loc["Women"]
    irr = (m.case / m.exam) / (w.case / w.exam); se = np.sqrt(1 / m.case + 1 / w.case)
    r1b.append(dict(period=p, men_pct=keep(100 * m.case / m.exam, 3), women_pct=keep(100 * w.case / w.exam, 3), IRR_men_vs_women=keep(irr, 3),
                    lo=keep(irr * np.exp(-1.96 * se), 3), hi=keep(irr * np.exp(1.96 * se), 3)))
T1b = pd.DataFrame(r1b)
# ② 15연령(성 합) 시기별 /10만
r2 = []
for p in PER:
    g = per(D, p).groupby("band")[["case", "exam"]].sum().reindex(BANDS)
    for b, v in g.iterrows(): r2.append(dict(period=p, band=b, cases=int(v.case), examinees=int(v.exam), rate_per100k=keep(1e5 * v.case / v.exam, 1)))
T2 = pd.DataFrame(r2)
# ③ 조율·직접표준화(표준 = 2015 전국 판정 수검자 성×15연령 구성)
W = D[D.year == 2015].set_index(["sex", "band"]).exam; W = W / W.sum()
r3 = []
for p in PER:
    g = per(D, p).groupby(["sex", "band"])[["case", "exam"]].sum(); rt = (g.case / g.exam).reindex(W.index)
    asr = float((rt * W).sum()); var = float(((W ** 2) * g.case.reindex(W.index) / g.exam.reindex(W.index) ** 2).sum())
    r3.append(dict(period=p, cases=int(g.case.sum()), examinees=int(g.exam.sum()), crude_per100k=keep(1e5 * g.case.sum() / g.exam.sum(), 1),
                   asr_per100k=keep(1e5 * asr, 1), asr_lo=keep(1e5 * (asr - 1.96 * np.sqrt(var)), 1), asr_hi=keep(1e5 * (asr + 1.96 * np.sqrt(var)), 1),
                   standard="2015 national examinee population, sex x 15 age bands (KOSIS DT_35007_N099)"))
for y in YEARS:
    g = D[D.year == y].set_index(["sex", "band"]); rt = (g.case / g.exam).reindex(W.index)
    r3.append(dict(period=str(y), cases=int(g.case.sum()), examinees=int(g.exam.sum()), crude_per100k=keep(1e5 * g.case.sum() / g.exam.sum(), 1),
                   asr_per100k=keep(1e5 * float((rt * W).sum()), 1), asr_lo=np.nan, asr_hi=np.nan, standard="2015 national examinee population, sex x 15 age bands (KOSIS DT_35007_N099)"))
T3 = pd.DataFrame(r3)
T1.to_csv(os.path.join(OUT, "person_axis_sex_age_IRR.csv"), index=False, encoding="utf-8-sig")
T1b.to_csv(os.path.join(OUT, "person_axis_sex_IRR.csv"), index=False, encoding="utf-8-sig")
T2.to_csv(os.path.join(OUT, "person_axis_age_period_per100k.csv"), index=False, encoding="utf-8-sig")
T3.to_csv(os.path.join(OUT, "person_axis_crude_ASR.csv"), index=False, encoding="utf-8-sig")
# json
J = {"source": "KOSIS DT_35007_N099 (examinees, '계') and DT_35007_N102 grp 00307 (renal-function management), national, sex x 15 age bands, 2011-2024",
     "standard": "2015 national general health-screening examinee population, sex x 15 age bands (direct standardisation)",
     "reference_group": f"{REF[0]} aged {REF[1]} years", "periods": list(PER),
     "sex_age_IRR": {f"{r.period}|{r.sex}|{r.agegrp}": [r.rate_per100k, r.IRR, r.lo, r.hi] for r in T1.itertuples()},
     "sex_IRR": {r.period: dict(men_pct=r.men_pct, women_pct=r.women_pct, IRR=[r.IRR_men_vs_women, r.lo, r.hi]) for r in T1b.itertuples()},
     "age_per100k": {f"{r.period}|{r.band}": r.rate_per100k for r in T2.itertuples()},
     "crude_asr_per100k": {r.period: dict(crude=r.crude_per100k, asr=r.asr_per100k, lo=None if pd.isna(r.asr_lo) else r.asr_lo, hi=None if pd.isna(r.asr_hi) else r.asr_hi) for r in T3.itertuples()}}
json.dump(J, open(os.path.join(NB, "MS", "build", "person_axis_nbkid.json"), "w"), ensure_ascii=False, indent=1)
pd.set_option("display.width", 250)
print(T1[T1.period == "2011-2024"].to_string(index=False)); print(T1b.to_string(index=False)); print(T3.to_string(index=False))
print(T2.pivot(index="band", columns="period", values="rate_per100k").reindex(BANDS).to_string())
