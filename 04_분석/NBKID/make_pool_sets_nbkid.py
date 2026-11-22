#!/usr/bin/env python3
# make_pool_sets_nbkid.py - Writes the model-settings files sets_NBKID_<version>.json (std3, s577, s2, s8, sc, sp), fixed before any covariate model was fitted,
#  and their SHA-256 hashes (sets_NBKID_hash_before_fit.txt). Principal set S1_narrative_full = 9 hypothesised district characteristics;
#  5 fixed adjustments; 3 negative controls; sensitivity-only variables (not in the principal model) listed in sensitivity_only_extra.
#  Input: 02_패널/EXT_canonical_pool_NBKID_261121.csv; the two comparison-outcome panels (panel_<TAG>.csv, 02_패널/<TAG>_Estd.csv) are only read for a check.
import os, json, sqlite3, unicodedata, hashlib
import numpy as np, pandas as pd
NFC = lambda s: unicodedata.normalize("NFC", str(s))
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE)); PAN = os.path.join(ROOT, "02_패널"); AN = os.path.join(ROOT, "04_분석")
POOL = "EXT_canonical_pool_NBKID.csv"
P = pd.read_csv(os.path.join(PAN, POOL), encoding="utf-8-sig"); P["fbd_region"] = P.fbd_region.map(NFC)
assert len(P) == 229 * 14 and P.drop(columns=["unmet_need_rate"]).notna().all().all()
print("풀", P.shape, "결측 0")
for t2 in ["KIDS", "NBBP"]:
    g = pd.read_csv(os.path.join(AN, t2, f"panel_{t2}.csv"), usecols=["year", "fbd_region", "cases"]); e = pd.read_csv(os.path.join(PAN, f"{t2}_Estd.csv"))
    print(t2, len(g), g.year.min(), g.year.max(), "O/E", round(g.cases.sum() / e.E_std.sum(), 4))
FORCED = ["urban_population_share", "physicians_per_1000", "sex_ratio", "fiscal_independence", "pop80_share"]
NEG = ["hit_and_run_rate", "traffic_culture_index", "fires_per_10000"]
SENS = ["tertiary_educated_share", "living_alone_elderly", "moderate_pa_rate", "public_sport_fac_per100k", "annual_mean_temp", "susp_htn_share", "obesity_rate", "scr_obesity_pstd", "im_specialists_per100k"]
HYP = ["susp_dm_share", "smoking_rate", "high_risk_drinking",
       "htn_treat_rate", "dm_treat_rate",
       "farm_pop_share", "basic_livelihood_rate", "screening_uptake", "walking_rate"]
OBES = "ext__obesity_rate"
X = lambda L: ["ext__" + v for v in L]
XK = lambda v: "ext__" + v
EXPECT = {
 "susp_htn_share": "+",
 "susp_dm_share": "+",
 "obesity_rate": "+",
 "smoking_rate": "+",
 "high_risk_drinking": "+",
 "htn_treat_rate": "-",
 "dm_treat_rate": "-",
 "im_specialists_per100k": "-",
 "annual_mean_temp": "+",
 "farm_pop_share": "+",
 "basic_livelihood_rate": "+",
 "screening_uptake": "-",
 "walking_rate": "-",
}
SETS = {"S1_narrative_full": X(HYP),
        "S2_core": X(["susp_dm_share", "htn_treat_rate", "basic_livelihood_rate", "screening_uptake", "walking_rate"]),
        "S3_metabolic_risk": X(["susp_dm_share", "smoking_rate", "high_risk_drinking"]),
        "S4_care_reach": X(["htn_treat_rate", "dm_treat_rate", "screening_uptake"]),
        "S5_work_deprivation": X(["farm_pop_share", "basic_livelihood_rate"]),
        "S6_risk_and_care": X(["susp_dm_share", "htn_treat_rate", "dm_treat_rate", "im_specialists_per100k"])}
pov = [{"year": 2023, "fbd_region": "경상남도창원시", "population": 1009038, "source": "KOSIS DT_1YL20651E, 2023 total for the integrated Changwon city"}]
allv = ["nbkid_cases", "nbkid_examinees"] + FORCED + NEG + SENS + HYP
LAB_M = "sex-and-age-standardised expected trace-proteinuria (normal B kidney-function management) judgements, 2 sexes x 15 age bands (<=19, 20-24 to 80-84, 85+) (national rates DT_35007_N102 '신장기능관리'/N099; district examinee age structure approximated from resident population x national uptake, sex totals observed from DT_35007_N098), 2011-2024"
WIN = {"W1_2011-15": [2011, 2015], "W2_2016-19": [2016, 2019], "W3_2020-24": [2020, 2024]}
def make():
    return {
     "tag": "NBKID", "outcome_column": "ext__nbkid_cases",
     "outcome_label_en": "Examinees judged 'normal B (borderline) - kidney-function management' at the National General Health Screening (urine dipstick protein trace (±); creatinine and eGFR have no normal-B band; NHIS health screening statistics DT_35007_N101, annual)",
     "disease_en": "Trace proteinuria (dipstick ±) flagged at health screening",
     "disease_ko": "검진 정상B(경계) 신장기능관리 판정(요단백 약양성 ±)",
     "pathway": "trace-proteinuria epidemiology set (measured metabolic risk, treatment reach, agricultural work, deprivation, screening selection, walking)",
     "selection_basis": "Outcome selected by prespecified spatial screening criteria (Moran's I, persistence, urban-rural R2, BYM2 mixing parameter)",
     "offset": "ext__nbkid_examinees", "exclude_years": [],
     "forced": X(FORCED), "negative_controls": X(NEG),
     "offset_E_column": "E_std", "offset_E_file": "NBKID_Estd3.csv", "offset_E_label": LAB_M,
     "graph": "593-edge contact graph, ST_GRAPH=04_분석/NBKID/graph_593.graph",
     "population_override": pov, "sensitivity_only_extra": X(SENS),
     "version": "std3 trace-proteinuria epidemiology set (measured glucose risk, smoking and drinking, treatment reach, agricultural work, deprivation, screening selection, walking), 2011-2024; temperature, suspected hypertension and obesity sensitivity only",
     "fixed_on": "fixed before any covariate model was fitted",
     "sets": SETS, "transform": {XK(v): "none" for v in FORCED + NEG + SENS + HYP},
     "expected_direction": {XK(k): v for k, v in EXPECT.items()},
     "rules": {"likelihood": "Poisson or negative binomial, whichever forced-adjustment-only model has the lower WAIC", "set_choice": "Lowest WAIC; within 4 WAIC points, the set with fewer variables", "credible": "95% CrI excludes 1", "vif_drop": 10, "multiverse_specs": 200,
               "main_model_note": "Principal model = S1_narrative_full (9 hypothesised characteristics + 5 fixed adjustments); no variable with VIF > 10 is kept.",
               "years": [2011, 2024], "windows": WIN, "covid_post_from": 2020, "windows_note": "Three periods fixed before results were seen."},
     "note": "Outcome: examinees judged normal B (borderline), kidney-function management (urine dipstick protein trace) at the general health screening (NHIS, KOSIS DT_35007_N101), 229 districts, 2011-2024; denominator = examinees judged; expected counts standardised for sex and 15 age groups.",
     "extra_columns": {"ext__" + v: [POOL, v] for v in allv}}
OUT = {}
for ver in ["std3", "s577", "s2", "sc", "sp", "s8"]:
    j = make()
    if ver == "s577": j["sens_note"] = "sensitivity: graph: queen contiguity (577 edges)"
    if ver == "s2": j["sens_note"] = "sensitivity: no fixed adjustment for the share aged 80 or older"; j["forced"] = X(FORCED[:4])
    if ver == "sp": j["sens_note"] = "sensitivity: population-based sex-age standardised expected counts"; j["offset_E_column"] = "E_pop"; j["offset_E_label"] = "population-based sex-age standardised E (resident population x national trace-proteinuria judgement per resident rate)"
    if ver == "sc": j["sens_note"] = "sensitivity: unstandardised expected counts"; j["offset_E_column"] = "E_crude"; j["offset_E_label"] = "crude E (examinees x national annual trace-proteinuria judgement rate)"
    if ver == "s8": j["sens_note"] = "sensitivity: age-only standardised expected counts"; j["offset_E_column"] = "E_std_age"; j["offset_E_label"] = "age-only standardised E"
    OUT[ver] = j
H = []
for ver, j in OUT.items():
    fp = os.path.join(HERE, f"sets_NBKID_{ver}.json"); json.dump(j, open(fp, "w"), ensure_ascii=False, indent=1)
    H.append(f"sets_NBKID_{ver}.json sha256={hashlib.sha256(open(fp, 'rb').read()).hexdigest()}")
open(os.path.join(HERE, "sets_NBKID_hash_before_fit.txt"), "w").write("# NBKID sets json SHA-256, written before any covariate model was fitted\n" + "\n".join(H) + "\n")
print("\n".join(H))
