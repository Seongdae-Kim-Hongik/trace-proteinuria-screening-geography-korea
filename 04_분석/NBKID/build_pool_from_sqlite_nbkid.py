#!/usr/bin/env python3
# build_pool_from_sqlite_nbkid.py - Builds the district-year covariate pool from the reference district-year dataset (SQLite, read-only; path in ST_GPT_DATASET).
#  Reads the raw (*_raw) columns for 2011-2024, checks that values agree across the tags that carry them, and adds the outcome and
#  denominator from 02_패널/NBKID_Estd3.csv. Corrections: Gunwi 2023 internal-medicine specialists (0 in the year of transfer to Daegu)
#  replaced by the mean of 2022 and 2024; screening obesity-management share z-standardised within 2011-2018 and 2019-2024 (judgement change).
#  Output: 02_패널/EXT_canonical_pool_NBKID_261121.csv, results_gate/S0_outcome_vs_sqlite.csv, S0_base_hash.txt
import os, sys, sqlite3, hashlib, unicodedata
import numpy as np, pandas as pd
NFC = lambda s: unicodedata.normalize("NFC", str(s))
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE)); PAN = os.path.join(ROOT, "02_패널"); RG = os.path.join(HERE, "results_gate")
DB = os.environ["ST_GPT_DATASET"]
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
h0 = sha(DB); m0 = os.stat(DB).st_mtime
con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
M = {"walking_rate": "walking_practice_rate", "moderate_pa_rate": "moderate_pa_rate", "public_sport_fac_per100k": "public_sports_facilities_per_100k", "village_sport_fac_per100k": "village_sports_facilities_per_100k",
     "obesity_rate": "obesity_rate", "high_risk_drinking": "high_risk_drinking_rate", "smoking_rate": "current_smoking_rate", "screening_uptake": "screening_uptake", "im_specialists_per100k": "im_specialists_per_100k",
     "basic_livelihood_rate": "basic_livelihood_rate", "urban_population_share": "urban_population_share", "physicians_per_1000": "physicians_per_1000", "sex_ratio": "sex_ratio", "fiscal_independence": "fiscal_independence",
     "pop80_share": "pop80_share", "hit_and_run_rate": "hit_and_run_rate", "traffic_culture_index": "traffic_culture_index", "fires_per_10000": "fires_per_10000", "tertiary_educated_share": "tertiary_educated_share",
     "living_alone_elderly": "older_single_household_share",
     "susp_htn_share": "screening_suspected_hypertension_share", "susp_dm_share": "screening_suspected_diabetes_share",
     "htn_treat_rate": "hypertension_treatment_rate", "dm_treat_rate": "diabetes_treatment_rate",
     "annual_mean_temp": "annual_mean_temperature", "farm_pop_share": "farm_population_share_ward0",
     "scr_obesity_share": "screening_normalB_obesity_share"}
cols = ",".join(f'"{v}_raw"' for v in M.values())
D = pd.read_sql(f'select tag, year, fbd_region, sigungu_cd, population, {cols} from disease_panel where year between 2011 and 2024', con); D["fbd_region"] = D.fbd_region.map(NFC)
key = D[D.tag == "ALZ"][["year", "fbd_region", "sigungu_cd", "population"]].sort_values(["year", "fbd_region"]).reset_index(drop=True)
assert len(key) == 229 * 14 and key.fbd_region.nunique() == 229
P = key.copy(); log = []
for k, v in M.items():
    c = v + "_raw"; ok = D.groupby("tag")[c].apply(lambda s: s.notna().sum() == 3206); tags = [t for t in ok[ok].index if not (k == "public_sport_fac_per100k" and t == "BRCA")]
    W = D[D.tag.isin(tags)].pivot_table(index=["year", "fbd_region"], columns="tag", values=c)
    spread = float((W.max(1) - W.min(1)).abs().max()); assert spread < 1e-9, (k, spread)
    s = W.iloc[:, 0].rename(k); P = P.merge(s.reset_index(), on=["year", "fbd_region"], how="left", validate="one_to_one"); log.append((k, c, len(tags), tags[0], spread))
RY = pd.read_sql('select year, fbd_region, "eloc__DT_1YL202104E__미충족의료율(표준화율)" as unmet_need_rate, "chk__정상B_신장기능관리" as sq_nbkid, "chk__수검률" as sq_uptake from region_year where year between 2011 and 2024', con); RY["fbd_region"] = RY.fbd_region.map(NFC)
P = P.merge(RY[["year", "fbd_region", "unmet_need_rate"]], on=["year", "fbd_region"], how="left", validate="one_to_one")
print("미충족의료 결측 칸:", int(P.unmet_need_rate.isna().sum()))
E = pd.read_csv(os.path.join(PAN, "NBKID_Estd3.csv"), encoding="utf-8-sig"); E["fbd_region"] = E.fbd_region.map(NFC)
P = P.merge(E[["year", "fbd_region", "cases", "examinees"]].rename(columns={"cases": "nbkid_cases", "examinees": "nbkid_examinees"}), on=["year", "fbd_region"], how="left", validate="one_to_one")
assert P.nbkid_cases.notna().all() and (P.nbkid_examinees > 0).all()
d80 = (P.set_index(["year", "fbd_region"]).pop80_share - E.set_index(["year", "fbd_region"]).pop80_share_raw.reindex(P.set_index(["year", "fbd_region"]).index)).abs().max(); print("pop80 SQLite vs E 파일 최대차:", d80)
O = E.merge(RY, on=["year", "fbd_region"], how="left", validate="one_to_one"); O["diff"] = O.cases - O.sq_nbkid
O[["year", "fbd_region", "cases", "sq_nbkid", "diff"]].to_csv(os.path.join(RG, "S0_outcome_vs_sqlite.csv"), index=False, encoding="utf-8-sig")
print("결과변수 SQLite 대조: 겹치는 칸", int(O.sq_nbkid.notna().sum()), "최대차", float(O["diff"].abs().max()), "SQLite 빈칸", O.loc[O.sq_nbkid.isna(), ["fbd_region", "year"]].values.tolist())
miss = P.drop(columns=["unmet_need_rate"]).isna().sum(); assert miss.sum() == 0, miss[miss > 0]
_g = P.fbd_region.str.contains("군위"); _v = P.loc[_g].set_index("year").im_specialists_per100k
assert float(_v.loc[2023]) == 0.0, _v.loc[2023]
P.loc[_g & (P.year == 2023), "im_specialists_per100k"] = (float(_v.loc[2022]) + float(_v.loc[2024])) / 2
P["scr_obesity_pstd"] = np.nan
for a_, b_ in ((2011, 2018), (2019, 2024)):
    m_ = P.year.between(a_, b_); x_ = P.loc[m_, "scr_obesity_share"]; P.loc[m_, "scr_obesity_pstd"] = (x_ - x_.mean()) / x_.std()
P.to_csv(os.path.join(PAN, "EXT_canonical_pool_NBKID_before_corrections.csv"), index=False, encoding="utf-8-sig"); print("풀 v3", P.shape, "군위 2023 im =", float(P.loc[_g & (P.year == 2023), "im_specialists_per100k"].iloc[0]))
import glob as _glob
S3P = os.environ.get("ST_IM_PARQUET") or (_glob.glob(os.path.join(os.path.expanduser("~"), "Library", "CloudStorage", "GoogleDrive-*", "*", "S.K", "G.Downloads", "독립변수 2", "국민건강결과_최종.parquet")) or [""])[0]
S3 = pd.read_parquet(S3P, columns=["year", "region", "내과전문의"]); S3["fbd_region"] = S3.region.map(NFC)
Q = P.merge(S3[["year", "fbd_region", "내과전문의"]], on=["year", "fbd_region"], how="left", validate="one_to_one"); assert Q["내과전문의"].notna().all()
cnt = Q.set_index(["fbd_region", "year"])["내과전문의"].astype(float).copy()
for r_ in cnt.index.get_level_values(0).unique():
    cnt.loc[(r_, 2018)] = (cnt.loc[(r_, 2017)] + cnt.loc[(r_, 2019)]) / 2
_gw = [r_ for r_ in cnt.index.get_level_values(0).unique() if "군위" in r_][0]; assert cnt.loc[(_gw, 2023)] == 0.0
cnt.loc[(_gw, 2023)] = (cnt.loc[(_gw, 2022)] + cnt.loc[(_gw, 2024)]) / 2
P["im_specialists_per100k_v3_double"] = P["im_specialists_per100k"]
P["im_specialists_per100k"] = [cnt.loc[(r_, y_)] / pop_ * 1e5 for r_, y_, pop_ in zip(P.fbd_region, P.year, P.population)]
_sj = P.fbd_region.str.contains("세종"); _f13 = float(P.loc[_sj & (P.year == 2013), "farm_pop_share"].iloc[0])
P["farm_pop_share_v3_interp"] = P["farm_pop_share"]; P.loc[_sj & (P.year >= 2014), "farm_pop_share"] = _f13
_gn = P.fbd_region.str.contains("강남구") & (P.year == 2011)
print("v4 내과(강남 2011) =", round(float(P.loc[_gn, "im_specialists_per100k"].iloc[0]), 2), "/ 세종 농가 2014–24 =", round(_f13, 2), "/ 전국 중앙값", P.groupby("year").im_specialists_per100k.median().round(1).to_dict())
P.drop(columns=["im_specialists_per100k_v3_double", "farm_pop_share_v3_interp"]).to_csv(os.path.join(PAN, "EXT_canonical_pool_NBKID.csv"), index=False, encoding="utf-8-sig")
P[["year", "fbd_region", "population", "im_specialists_per100k_v3_double", "im_specialists_per100k", "farm_pop_share_v3_interp", "farm_pop_share"]].to_csv(os.path.join(HERE, "checks", "out", "v4_corrections_im_farm.csv"), index=False, encoding="utf-8-sig")
con.close(); h1 = sha(DB); m1 = os.stat(DB).st_mtime
open(os.path.join(RG, "S0_base_hash.txt"), "w").write(f"기준 SQLite {DB}\n전 sha256={h0} mtime={m0}\n후 sha256={h1} mtime={m1}\n같음={h0 == h1 and m0 == m1}\n" + "\n".join(f"{a}<-{b} 태그 {c}개(대표 {d}) 태그 간 최대차 {e}" for a, b, c, d, e in log) + "\n")
print("기준 파일 해시 전후 같음:", h0 == h1 and m0 == m1)
