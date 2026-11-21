#!/usr/bin/env python3
# 99NK2_build_Estd3_nbkid.py - District counts and expected counts for the outcome.
#  cases = examinees judged normal B, kidney-function management (KOSIS DT_35007_N101); denominator = examinees judged (DT_35007_N098).
#  Sejong is taken from its province rows. Expected counts by indirect standardisation in 2 sexes x 15 age groups using national rates
#  (DT_35007_N099, DT_35007_N102); district examinee age structure = resident population (ST_FBDZIP) x national uptake by sex and age,
#  rescaled to the observed examinees by sex. Expected counts use no added constant.
#  Output: 02_패널/NBKID_Estd3.csv (cases, examinees, E_std, E_std_age, E_pop, E_crude, pop80_share_raw)
import os, sys, json, importlib.util, unicodedata
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); AN = os.path.join(ROOT, "04_분석"); PAN = os.path.join(ROOT, "02_패널")
SRC = os.path.join(ROOT, "00_원천", "C_NHIS_검진보험")
def _load(nm, f):
    s = importlib.util.spec_from_file_location(nm, os.path.join(HERE, f)); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
FBDZIP = os.environ["ST_FBDZIP"]
cw = _load("cw", "40_region_crosswalk.py"); m50 = _load("m50", "50_build_master_panel.py")
NFC = lambda s: unicodedata.normalize("NFC", str(s))
TAG = "NBKID"; YEARS = list(range(2011, 2025))
BANDS = ["0-19"] + [f"{lo}-{lo+4}" for lo in range(20, 85, 5)] + ["85+"]
SEXES = ["남자", "여자"]; OUTS = {"": "00307"}
assert os.environ.get("ST_FBDZIP"), "ST_FBDZIP(창원 2023 교정 인구 폴더)를 지정"
def nb_band(a):
    a = str(a).replace(" ", "")
    if a == "19세이하": return "0-19"
    if a == "85세이상": return "85+"
    if a.endswith("세") and "~" in a:
        lo, hi = a[:-1].split("~"); return f"{lo}-{hi}"
    return None
A = pd.read_csv(os.path.join(SRC, "DT_35007_N099_national_long.csv"), encoding="utf-8-sig"); A = A[(A.grp == "계") & A.year.isin(YEARS)]
B = pd.read_csv(os.path.join(SRC, "DT_35007_N102_national_long.csv"), encoding="utf-8-sig", dtype={"grp_code": str}); B = B[B.grp_code.isin(OUTS.values()) & B.year.isin(YEARS)]
for X in (A, B): X["band"] = X.age.map(nb_band)
AT = A[A.age == "계"].set_index(["year", "sex"]).value
EX = A[A.band.notna() & A.sex.isin(SEXES)].set_index(["year", "sex", "band"]).value.unstack()[BANDS]
for s in SEXES: assert (abs(EX.xs(s, level="sex").sum(1) - AT.xs(s, level="sex")) <= 1).all(), s
RATE, NUM = {}, {}
for suf, code in OUTS.items():
    b = B[B.grp_code == code]; bt = b[b.age == "계"].set_index(["year", "sex"]).value
    NUM[suf] = b[b.band.notna() & b.sex.isin(SEXES)].set_index(["year", "sex", "band"]).value.unstack()[BANDS].fillna(0)
    for s in SEXES: assert (abs(NUM[suf].xs(s, level="sex").sum(1) - bt.xs(s, level="sex")) <= 1).all(), (suf, s)
    RATE[suf] = NUM[suf] / EX
print("[99NK2] 전국 의심률(%, 성 합) 연도별:", (100 * NUM[""].groupby(level="year").sum().sum(1) / EX.groupby(level="year").sum().sum(1)).round(2).to_dict())
base = pd.read_csv(os.path.join(AN, "PERIO", "panel_PERIO.csv"), encoding="utf-8-sig"); base["fbd_region"] = base.fbd_region.map(NFC); base["key"] = base.key.map(NFC)
base = base[["year", "key", "fbd_region"]].drop_duplicates()
PKEYS = set(base.key); SIDO_HDR = {k for k in cw.SIDO_STD if len(k) >= 3 and (k.endswith("시") or k.endswith("도"))}
FIVE = [f"{lo}-{lo+4}세" for lo in range(0, 85, 5)]; OLD85 = ["85-89세", "90-94세", "95-99세", "100세 이상"]
EXTRA = {("충청북도", "청원군"): "충청북도|청주시"}; rows = []
for y in YEARS:
    d = pd.read_csv(os.path.join(FBDZIP, f"{y}.csv"), encoding="cp949", dtype=str)
    d["행정구역별"] = d["행정구역별"].map(NFC).str.strip(); d["연령별"] = d["연령별"].map(NFC).str.strip()
    for s in ["전체"] + SEXES: d[s] = pd.to_numeric(d[s].str.replace(",", ""), errors="coerce")
    sido = None; city = None; blocks = (d["행정구역별"] != d["행정구역별"].shift()).cumsum()
    for _, g in d.groupby(blocks, sort=False):
        nm = g["행정구역별"].iat[0]
        if nm == "전국": continue
        if nm in SIDO_HDR:
            sido = nm; city = None
            if nm != "세종특별자치시": continue
            g = g.drop_duplicates("연령별")
        sg = nm.split("-")[-1]
        if sg.endswith("시"): city = sg
        elif sg.endswith("구") and city is not None and cw.sido_std(sido) not in ("서울시", "부산시", "대구시", "인천시", "광주시", "대전시", "울산시"): continue
        else: city = None
        key = EXTRA.get((cw.sido_std(sido), sg), cw.make_key(sido, sg))
        if key not in PKEYS and key.endswith("군") and key[:-1] + "시" in PKEYS: key = key[:-1] + "시"
        for s in SEXES:
            a = g.set_index("연령별")[s]; tot = a.get("계", np.nan)
            if pd.isna(tot): continue
            f = a.reindex(FIVE); o = a.reindex(OLD85).fillna(0); assert f.notna().all(), (y, nm, s)
            assert abs(tot - f.sum() - o.sum()) < 0.5, (y, nm, s)
            rows.append({"year": y, "sex": s, "key": NFC(key), "tot": tot, "0-19": f["15-19세"], **{b: f[b + "세"] for b in BANDS[1:-1]}, "85+": o.sum()})
G = pd.DataFrame(rows).groupby(["year", "key", "sex"], as_index=False)[["tot"] + BANDS].sum()
NATPOP = G.groupby(["year", "sex"])[BANDS].sum(); UPT = EX / NATPOP
L = pd.read_csv(os.path.join(PAN, "NHIS검진_시군구_연도_long.csv"), usecols=["year", "sido", "sigungu", "level", "sex", "category", "value", "tbl_id"], engine="pyarrow")
CAT = {"계": "exam", "신장기능관리": "c"}
L = L[L.tbl_id.isin(["DT_35007_N098", "DT_35007_N101"]) & ((L.level == "시군구") | ((L.level == "시도") & (L.sido == "세종"))) & L.year.isin(YEARS) & L.sex.isin(SEXES) & L.category.isin(CAT.keys())].copy()
L = L[~((L.tbl_id == "DT_35007_N101") & (L.category == "계"))]; L["v"] = L.category.map(CAT)
pairs = L[["sido", "sigungu"]].drop_duplicates(); mp = {(a, b): m50.map_region(a, b, PKEYS) for a, b in zip(pairs.sido, pairs.sigungu)}
L["key"] = [mp[(a, b)]["key"] for a, b in zip(L.sido, L.sigungu)]; L["role"] = [mp[(a, b)]["role"] for a, b in zip(L.sido, L.sigungu)]
sj = (L.level == "시도") & (L.sido == "세종"); L.loc[sj, "key"] = "세종시|세종시"; L.loc[sj, "role"] = "primary"; L = L[~((L.level == "시군구") & (L.sido == "세종"))]
L = L[L.key.notna()]; L["key"] = L.key.map(NFC)
hasprim = L[L.role == "primary"].groupby(["key", "year"]).size().index
L = L[~((L.role == "sub") & L.set_index(["key", "year"]).index.isin(hasprim))]
SX = L.pivot_table(index=["year", "key"], columns=["v", "sex"], values="value", aggfunc="sum")
P = base[base.year.isin(YEARS)].merge(SX.reset_index().droplevel(1, axis=1).iloc[:, :2], on=["year", "key"], how="inner") if False else None
k = base[base.year.isin(YEARS)].set_index(["year", "key"]).index
S = SX.reindex(k); ok = S.notna().all(axis=1) & (S[("exam", "남자")] > 0) & (S[("exam", "여자")] > 0)
P = base[base.year.isin(YEARS)].reset_index(drop=True)[ok.values].reset_index(drop=True); S = S[ok.values]
assert P.fbd_region.nunique() == 229 and len(P) == 229 * 14 and P.groupby("fbd_region").year.nunique().eq(len(YEARS)).all()
print(f"[99NK2] 시군구-연도 {len(k)} 중 관측 {len(P)}; 제외 {sorted(set(map(tuple, base[base.year.isin(YEARS)].reset_index(drop=True)[~ok.values][['year','fbd_region']].values.tolist())))[:12]}")
P["examinees"] = (S[("exam", "남자")] + S[("exam", "여자")]).values
for suf in OUTS: P["cases" + suf] = (S[("c" + suf, "남자")] + S[("c" + suf, "여자")]).values
MA = pd.read_parquet(os.path.join(PAN, "MASTER_시군구_연도_wide.parquet"), columns=["year", "fbd_region", "chk__정상B_신장기능관리"]); MA["fbd_region"] = MA.fbd_region.map(NFC)
c = P.merge(MA, on=["year", "fbd_region"], how="left", validate="one_to_one")
dd = (c.cases - c["chk__정상B_신장기능관리"]).abs(); print(f"[99NK2] MASTER 검산 chk__정상B_신장기능관리: 최대차 {dd.max():.0f}, MASTER 빈칸 {int(c['chk__정상B_신장기능관리'].isna().sum())} {c.loc[c['chk__정상B_신장기능관리'].isna(), ['fbd_region','year']].values.tolist()}"); assert dd.max() <= 2
nat_ex = EX.groupby(level="year").sum().sum(1); nat_c = NUM[""].groupby(level="year").sum().sum(1)
print("[99NK2] 시군구 합/전국 수검자:", (P.groupby("year").examinees.sum() / nat_ex).round(3).to_dict())
print("[99NK2] 시군구 합/전국 의심:", (P.groupby("year").cases.sum() / nat_c).round(3).to_dict())
cs = {suf: P.groupby("year")["cases" + suf].transform("sum").values for suf in OUTS}
rsc = lambda e, suf: e * cs[suf] / pd.Series(e).groupby(P.year.values).transform("sum").values
E = {suf: np.zeros(len(P)) for suf in OUTS}; EP = {suf: np.zeros(len(P)) for suf in OUTS}; EXA = 0; scs = []
for s in SEXES:
    X = P[["year", "key"]].merge(G[G.sex == s], on=["year", "key"], how="left", validate="one_to_one"); assert X.tot.notna().all(), X[X.tot.isna()][["year", "key"]]
    raw = X[BANDS].values * UPT.xs(s, level="sex").loc[X.year.values].values
    sc = S[("exam", s)].values / raw.sum(1); exam = raw * sc[:, None]; scs.append(sc); EXA = EXA + exam
    for suf in OUTS:
        E[suf] += (exam * RATE[suf].xs(s, level="sex").loc[X.year.values].values).sum(1)
        EP[suf] += (X[BANDS].values * (NUM[suf].xs(s, level="sex") / NATPOP.xs(s, level="sex")).loc[X.year.values].values).sum(1)
RATE_A = NUM[""].groupby(level="year").sum() / EX.groupby(level="year").sum()
e_age = (EXA * RATE_A.loc[P.year.values].values).sum(1)
for suf in OUTS:
    P["E_std" + suf] = rsc(E[suf], suf); P["E_pop" + suf] = rsc(EP[suf], suf)
    assert np.allclose(P.groupby("year")["E_std" + suf].sum(), P.groupby("year")["cases" + suf].sum())
P["E_std_age"] = rsc(e_age, "")
H = pd.read_csv(os.path.join(PAN, "HTND_Estd2.csv"), encoding="utf-8-sig"); H["fbd_region"] = H.fbd_region.map(NFC)
P = P.merge(H[["year", "fbd_region", "pop80_share_raw"]], on=["year", "fbd_region"], how="left", validate="one_to_one"); assert P.pop80_share_raw.notna().all()
nat = P.groupby("year").apply(lambda d: d.cases.sum() / d.examinees.sum()); P["E_crude"] = P.examinees * P.year.map(nat)
for suf in OUTS: P["SIR_std" + suf] = P["cases" + suf] / P["E_std" + suf]
P["SIR_pop"] = P.cases / P.E_pop; P["SIR_crude"] = P.cases / P.E_crude
cols = ["year", "fbd_region", "examinees", "cases", "E_std", "E_std_age", "E_pop", "E_crude", "pop80_share_raw", "SIR_std", "SIR_pop", "SIR_crude"]
fp = os.environ.get("ST_ESTD3_OUT") or os.path.join(PAN, f"{TAG}_Estd3.csv"); P[cols].to_csv(fp, index=False, encoding="utf-8-sig")
sc_all = np.concatenate(scs); r = P.E_std / P.E_crude
print(f"[99NK2] {fp}: 행 {len(P)}, 연도 {P.year.min()}-{P.year.max()}, 시군구 {P.fbd_region.nunique()}, 의심 합 {int(P.cases.sum()):,}, 수검자 합 {int(P.examinees.sum()):,}")
print(f"[99NK2] E_std/E_crude {r.min():.3f}~{r.max():.3f} · E_pop/E_std {(P.E_pop/P.E_std).min():.3f}~{(P.E_pop/P.E_std).max():.3f} · Spearman(SIR_std,SIR_pop) {P.SIR_std.corr(P.SIR_pop, method='spearman'):.3f} · Spearman(SIR_std,SIR_crude) {P.SIR_std.corr(P.SIR_crude, method='spearman'):.3f}")
print(f"[99NK2] 재척도 배수 {sc_all.min():.3f}~{sc_all.max():.3f} (P5 {np.quantile(sc_all,.05):.3f}, P95 {np.quantile(sc_all,.95):.3f})")
