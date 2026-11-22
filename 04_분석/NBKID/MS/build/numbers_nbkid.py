# numbers_nbkid.py - Computes every number reported in the manuscript from the result tables.
#  Input: 04_분석/NBKID/results_std3 (principal; includes S60_extra.csv from 88V_nbkid_extra.R), results_s577, results_s2, results_s8, results_sc,
#  results_sp, results_gate, panel_NBKID_std3.csv, ST_GRAPH, ST_SHP. Output: numbers_nbkid.json, gi_nbkid.csv
import os, json, unicodedata, numpy as np, pandas as pd
from decimal import Decimal, ROUND_HALF_UP
def rnd(x, nd=0):
    q = Decimal(1).scaleb(-nd); v = float(Decimal(repr(float(x))).quantize(q, rounding=ROUND_HALF_UP)); return int(v) if nd == 0 else v
from scipy.stats import spearmanr, norm
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.abspath(os.path.join(HERE, "..", "..")); AN = os.path.dirname(B); ROOT = os.path.dirname(AN); PAN = os.path.join(ROOT, "02_패널")
RS = lambda v="std3": os.path.join(B, f"results_{v}"); NFC = lambda s: unicodedata.normalize("NFC", str(s))
SHP = os.environ.get("ST_SHP") or os.path.join(os.path.expanduser("~"), "Library/CloudStorage", next((d for d in os.listdir(os.path.expanduser("~/Library/CloudStorage")) if d.startswith("GoogleDrive")), ""), "내 드라이브/S.K/기타/R_Studio/FBD_DATA_ZIP/final.shp")
VARS = {"susp_dm": "screening_suspected_diabetes_share_v", "smoking_rate": "current_smoking_rate", "high_risk_drinking": "high_risk_drinking_rate",
        "htn_treat": "hypertension_treatment_rate_v", "dm_treat": "diabetes_treatment_rate_v",         "farm_share": "farm_population_share_v", "basic_livelihood_rate": "basic_livelihood_rate", "screening_uptake": "screening_uptake", "walking_rate": "walking_practice_rate",
        "urban_population_share": "urban_population_share", "physicians_per_1000": "physicians_per_1000", "sex_ratio": "sex_ratio", "fiscal_independence": "fiscal_independence", "pop80_share": "pop80_share"}
SVARS = {"temperature": "annual_mean_temperature_v", "susp_htn": "screening_suspected_hypertension_share_v", "obesity_rate": "obesity_rate_sr", "scr_obesity": "ext__scr_obesity_pstd"}
HYP = list(VARS)[:9]; ROB = ["walking_rate", "basic_livelihood_rate"]
def tab(f, v="std3"):
    d = pd.read_csv(os.path.join(RS(v), f)); d["cov"] = d["cov"].astype(str); return d
def get(d, c):
    r = d[d["cov"] == c].iloc[0]; return [float(float(r.IRR)), float(float(r.lo)), float(float(r.hi))]
N = {}
P = pd.read_csv(os.path.join(B, "panel_NBKID_std3.csv")); P["fbd_region"] = P.fbd_region.map(NFC)
N["n_rows"] = len(P); N["n_districts"] = P.fbd_region.nunique(); N["cases_total_m"] = rnd(P.cases.sum() / 1e6, 2); N["exam_total_m"] = rnd(P.offset.sum() / 1e6, 1)
g = P.groupby("year"); N["rate"] = {int(k): rnd(100 * v, 2) for k, v in (g.cases.sum() / g.offset.sum()).items()}
N["exam"] = {int(k): rnd(v / 1e6, 2) for k, v in g.offset.sum().items()}; N["cases"] = {int(k): rnd(v / 1e6, 2) for k, v in g.cases.sum().items()}
L = P.groupby("fbd_region").apply(lambda d: d.cases.sum() / d.E.sum(), include_groups=False)
N["sir_min"], N["sir_max"] = rnd(L.min(), 2), rnd(L.max(), 2); N["sir_p10"], N["sir_p90"] = rnd(L.quantile(.1), 2), rnd(L.quantile(.9), 2); N["p90p10"] = rnd(L.quantile(.9) / L.quantile(.1), 2)
N["sir_top"] = [[r, rnd(v, 2)] for r, v in L.sort_values(ascending=False).head(8).items()]; N["sir_bottom"] = [[r, rnd(v, 2)] for r, v in L.sort_values().head(8).items()]
P["dens"] = P.population / (P.area_m2 / 1e6)
D = P.groupby("fbd_region").agg(wk=("walking_practice_rate_raw", "mean"), up=("screening_uptake_raw", "mean"), dens=("dens", "mean"), p80=("pop80_share_raw", "mean"), bl=("basic_livelihood_rate_raw", "mean"), ht=("screening_suspected_hypertension_share_v_raw", "mean"), tp=("annual_mean_temperature_v_raw", "mean"), urb=("urban_population_share_raw", "mean"))
N["rho_bl_dens"] = rnd(spearmanr(D.bl, D.dens).correlation, 2); N["rho_bl_up"] = rnd(spearmanr(D.bl, D.up).correlation, 2); N["rho_bl_p80"] = rnd(spearmanr(D.bl, D.p80).correlation, 2)
N["rho_sir_bl"] = rnd(spearmanr(L.reindex(D.index), D.bl).correlation, 2); N["rho_sir_up"] = rnd(spearmanr(L.reindex(D.index), D.up).correlation, 2); N["rho_sir_urb"] = rnd(spearmanr(L.reindex(D.index), D.urb).correlation, 2); _rng = np.random.default_rng(20261019); _x = L.reindex(D.index).values; _u = D.urb.values; _bs = [spearmanr(_x[_i], _u[_i]).correlation for _i in (_rng.integers(0, len(_x), len(_x)) for _ in range(2000))]; N["rho_sir_urb_ci"] = [rnd(float(np.percentile(_bs, 2.5)), 2), rnd(float(np.percentile(_bs, 97.5)), 2)]; N["rho_sir_ht"] = rnd(spearmanr(L.reindex(D.index), D.ht).correlation, 2)
WIN = {"2011-2015": (2011, 2015), "2016-2019": (2016, 2019), "2020-2024": (2020, 2024), "2011-2024": (2011, 2024)}
DESC = {**VARS, **{k: v for k, v in SVARS.items() if k != "scr_obesity"}}
raw = {k: (DESC[k] + "_raw") for k in DESC}; assert all(c in P for c in raw.values()), [c for c in raw.values() if c not in P]
T1 = {}
for w, (a, b) in WIN.items():
    s = P[P.year.between(a, b)]; sw = s.groupby("fbd_region").apply(lambda d: d.cases.sum() / d.E.sum(), include_groups=False)
    T1[w] = dict(n=len(s), cases=float(s.cases.sum()), exam=float(s.offset.sum()), rate=float(s.cases.sum() / s.offset.sum()), sir_med=float(sw.median()), sir_p10=float(sw.quantile(.1)), sir_p90=float(sw.quantile(.9)),
                 cov={k: [float(s[raw[k]].mean()), float(s[raw[k]].std())] for k in DESC})
N["T1"] = T1
mo = pd.read_csv(os.path.join(RS(), "S11b_moran_yearly_EB.csv")); N["moran_min"], N["moran_max"] = rnd(mo.moran_I.min(), 2), rnd(mo.moran_I.max(), 2); N["moran_mk_tau"], N["moran_mk_p"] = float(mo.trend_tau.iloc[0]), float(mo.trend_p.iloc[0])
_hy = pd.read_csv(os.path.join(RS(), "S3_main_hyper.csv")); _ph = _hy[_hy.hyperpar.str.startswith("Phi")].iloc[0]; N["phi_ci"] = [rnd(_ph["mean"], 2), rnd(_ph["X0.025quant"], 2), rnd(_ph["X0.975quant"], 2)]
_cv = pd.read_csv(os.path.join(RS(), "S5_covid.csv")); _cp = _cv[_cv["cov"] == "post"].iloc[0]; N["covid_post"] = [float(_cp.IRR), float(_cp.lo), float(_cp.hi)]
_rf = pd.read_csv(os.path.join(RS(), "S1b_refit_twice.csv")).iloc[0]; N["refit_waic"] = [float(_rf.WAIC1), float(_rf.WAIC2)]
KP = pd.read_csv(os.path.join(RS(), "S12_persistence.csv")).iloc[0]; N["kendall_w"] = rnd(float(KP.kendall_W), 2); N["kendall_perm_p"] = float(KP.perm_p); N["kendall_B"] = int(KP.B); N["kendall_null975"] = float(KP["perm_null_97.5"])
_P = pd.read_csv(os.path.join(B, "panel_NBKID_std3.csv")); _X = _P[list(VARS.values())].astype(float); _X = (_X - _X.mean()) / _X.std(); _vif = {}
for _c in _X.columns:
    _y = _X[_c].values; _Z = np.column_stack([np.ones(len(_X)), _X.drop(columns=_c).values]); _b = np.linalg.lstsq(_Z, _y, rcond=None)[0]; _vif[_c] = 1 / (((_y - _Z @ _b) ** 2).sum() / ((_y - _y.mean()) ** 2).sum())
N["vif_max"] = float(max(_vif.values())); N["vif_max_var"] = max(_vif, key=_vif.get); N["n_hyp"] = len(HYP); N["n_main"] = len(VARS)
Y = P.pivot_table(index="fbd_region", columns="year", values="SIR"); N["rho_2011_2024"] = rnd(spearmanr(Y[2011], Y[2024]).correlation, 2); _rg = np.random.default_rng(20261019); _yy = Y[[2011, 2024]].dropna().values; _bb = [spearmanr(_yy[_i, 0], _yy[_i, 1]).correlation for _i in (_rg.integers(0, len(_yy), len(_yy)) for _ in range(2000))]; N["rho_2011_2024_ci"] = [rnd(float(np.percentile(_bb, 2.5)), 2), rnd(float(np.percentile(_bb, 97.5)), 2)]
import geopandas as gpd
G = open(os.environ["ST_GRAPH"]).read().split("\n"); nreg = int(G[0]); nb = {}
for line in G[1:]:
    t = line.split()
    if len(t) >= 2: nb[int(t[0])] = [int(x) for x in t[2:]]
area = pd.DataFrame({"fbd_region": gpd.read_file(SHP)["region"].values}); area["fbd_region"] = area.fbd_region.map(NFC); assert len(area) == nreg
id2r = dict(zip(range(1, nreg + 1), area.fbd_region))
x = np.log(L.reindex([id2r[i] for i in range(1, nreg + 1)]).values); assert not np.isnan(x).any(); xb, sd = x.mean(), x.std(); gi = []
for i in range(1, nreg + 1):
    w = np.zeros(nreg); w[i - 1] = 1; w[[j - 1 for j in nb.get(i, [])]] = 1; W = w.sum()
    gi.append(((w * x).sum() - xb * W) / (sd * np.sqrt((nreg * (w ** 2).sum() - W ** 2) / (nreg - 1))))
GI = pd.DataFrame({"fbd_region": [id2r[i] for i in range(1, nreg + 1)], "SIR": np.exp(x), "z": gi}); GI["p"] = 2 * norm.sf(GI.z.abs())
GI["cls"] = np.where((GI.p < 0.05) & (GI.z > 0), "hot", np.where((GI.p < 0.05) & (GI.z < 0), "cold", "other"))
S6 = pd.read_csv(os.path.join(RS(), "S6_spatial_RR.csv")); S6["fbd_region"] = S6.fbd_region.map(NFC); GI = GI.merge(S6[["fbd_region", "sigungu_cd"]], on="fbd_region", validate="one_to_one")
GI.to_csv(os.path.join(HERE, "gi_nbkid.csv"), index=False, encoding="utf-8-sig")
SIDO = {"서울": "Seoul", "부산": "Busan", "대구": "Daegu", "인천": "Incheon", "광주": "Gwangju", "대전": "Daejeon", "울산": "Ulsan", "세종": "Sejong", "경기": "Gyeonggi", "강원": "Gangwon", "충청북": "Chungbuk",
        "충청남": "Chungnam", "전라북": "Jeonbuk", "전라남": "Jeonnam", "경상북": "Gyeongbuk", "경상남": "Gyeongnam", "제주": "Jeju"}
GI["sido"] = GI.fbd_region.map(lambda r: next(v for k, v in SIDO.items() if r.startswith(k)))
for c in ("hot", "cold"):
    s = GI[GI.cls == c]; N[f"gi_{c}"] = len(s); N[f"gi_{c}_sido"] = s.sido.value_counts().to_dict(); N[f"gi_{c}_sir_med"] = rnd(s.SIR.median(), 2)
    N[f"gi_{c}_gun"] = int(s.fbd_region.str.endswith("군").sum()); N[f"gi_{c}_gu"] = int(s.fbd_region.str.endswith("구").sum()); N[f"gi_{c}_si"] = len(s) - N[f"gi_{c}_gun"] - N[f"gi_{c}_gu"]
from statsmodels.stats.multitest import multipletests
_rj = multipletests(GI.p, 0.05, "fdr_bh")[0]; N["gi_hot_fdr"] = int(((GI.z > 0) & _rj).sum()); N["gi_cold_fdr"] = int(((GI.z < 0) & _rj).sum())
DD = D.join(GI.set_index("fbd_region").cls); N["hc"] = {c: {k: rnd(float(DD[DD.cls == c][k].median()), 2) for k in ("wk", "up", "dens", "p80", "bl", "ht", "tp", "urb")} for c in ("hot", "other", "cold")}
M = tab("S3_main_IRR.csv"); N["main"] = {k: get(M, VARS[k]) for k in VARS}
F = pd.read_csv(os.path.join(RS(), "S3_main_fit.csv")).iloc[0]; N["waic_main"] = float(F.WAIC); N["phi"] = float(F.phi); N["n_hot_bym"] = int(F.n_hotspot); N["n_cold_bym"] = int(F.n_coldspot)
N["resid_moran"] = float(F.resid_moran_I); N["resid_moran_p"] = float(F.resid_moran_p); N["obs_moran"] = float(F.obs_moran_I)
LK = pd.read_csv(os.path.join(RS(), "S1_likelihood.csv")); N["waic_forced"] = float(LK[LK.family == "nbinomial"].WAIC.iloc[0]); N["waic_forced_pois"] = float(LK[LK.family == "poisson"].WAIC.iloc[0])
hy = pd.read_csv(os.path.join(RS(), "S3_main_hyper.csv")); N["hyper"] = hy.round(3).values.tolist()
NC = tab("S3_negcontrol.csv"); N["nc"] = {c.replace("ext__", ""): get(NC, c) for c in NC["cov"]}
MV = pd.read_csv(os.path.join(RS(), "S7_multiverse_summary.csv")); MV["cov"] = MV["cov"].astype(str)
N["mv"] = {k: [rnd(float(MV[MV["cov"] == VARS[k]].cred_frac.iloc[0]), 3), rnd(float(MV[MV["cov"] == VARS[k]].sign_agree.iloc[0]), 3), int(MV[MV["cov"] == VARS[k]].n_specs.iloc[0])] for k in HYP}
for t in ("KH_typeI", "KH_typeII", "KH_typeIII", "KH_typeIV"):
    f = "S9_typeI_IRR.csv" if t == "KH_typeI" else f"S43_part_{t}.csv"
    d = tab(f); N[t] = {k: get(d, VARS[k]) for k in VARS}; N[t + "_waic"] = float(d.WAIC.iloc[0])
d = pd.read_csv(os.path.join(RS(), "S50_resid_moran_by_spec.csv")).set_index("spec"); N["rm_waic"] = float(d.loc["principal", "WAIC"]); N["rm_spec"] = {s: rnd(float(d.loc[s, "resid_moran_I"]), 2) for s in d.index}
for a in ("E1_lag1", "E2_lag2", "D1_spatial_plus", "C4_sido_year_iid", "C1_poisson_lognormal", "G1_kchs_ma3"):
    d = tab(f"S41_part_{a}.csv"); N[a] = {k: get(d, VARS[k]) for k in VARS if (d["cov"] == VARS[k]).any()}
for a in ("no_spatial", "iid_only", "icar_besag", "RSR"):
    d = tab(f"S43_part_{a}.csv"); N[a] = {k: get(d, VARS[k]) for k in VARS if (d["cov"] == VARS[k]).any()}
d = tab("S41_part_bw_twoway.csv"); N["bw"] = {k + "|" + p: get(d, VARS[k] + "__" + p) for k in ROB for p in ("between", "within") if (d["cov"] == VARS[k] + "__" + p).any()}
def _wsd(v):
    x = P[v].astype(float); dev = x - P.groupby("fbd_region")[v].transform("mean") - P.groupby("year")[v].transform("mean") + x.mean(); return float(dev.std() / x.std())
N["bw_wsd"] = {k: [_wsd(VARS[k])] + [float(N["bw"][k + "|within"][j]) ** _wsd(VARS[k]) for j in range(3)] for k in ROB if k + "|within" in N["bw"]}
N["nco"] = {}
for a in ("KIDS", "NBBP"):
    d = tab(f"S41_part_H1_nco_{a}.csv"); N["nco"].update({f"{a}|{k}": get(d, VARS[k]) for k in VARS})
W = tab("S4_window_IRR.csv"); N["win"] = {f"{k}|{w}": get(W[W.window == w], VARS[k]) for k in VARS for w in W.window.unique()}
RET = pd.read_csv(os.path.join(RS(), "S41_retention_by_cov.csv")); N["ret"] = {k: [int(RET[RET["cov"] == VARS[k]].kept.iloc[0]), int(RET[RET["cov"] == VARS[k]].n_specs.iloc[0])] for k in ROB if (RET["cov"] == VARS[k]).any()}
EV = pd.read_csv(os.path.join(RS(), "S41_evalues.csv")); N["evalue"] = {k: [float(EV[EV["cov"] == VARS[k]].E_point.iloc[0]), float(EV[EV["cov"] == VARS[k]].E_limit.iloc[0])] for k in ROB if (EV["cov"] == VARS[k]).any()}
VP = pd.read_csv(os.path.join(RS(), "S42_variance_partition.csv")).set_index("component"); N["vp"] = {k: rnd(100 * float(v), 1) for k, v in VP.share_of_linear_predictor_var.items()}
PPC = pd.read_csv(os.path.join(RS(), "S42_ppc.csv")).iloc[0]; N["ppc95"] = rnd(100 * PPC.coverage95, 1)
OOS = pd.read_csv(os.path.join(RS(), "S41_oos2024.csv")).iloc[0]; N["oos_cov"] = rnd(100 * OOS.coverage95, 1); N["oos_mape"] = rnd(OOS.median_abs_pct_error, 1)
RG = pd.read_csv(os.path.join(RS(), "S50_region_slopes.csv")); N["reg"] = {f"{r.covariate}|{r.region}": [r.IRR, r.lo, r.hi] for r in RG.itertuples()}
_ph = os.path.join(RS(), "S46_nbkid_posthoc.csv")
if os.path.exists(_ph): PH = pd.read_csv(_ph); N["ph"] = {f"{r.analysis}|{r.cov}": [float(r.IRR), float(r.lo), float(r.hi)] for r in PH.itertuples()}
LK5 = pd.read_csv(os.path.join(RS(), "S1c_likelihood5.csv")); N["lik5"] = {r.family: float(r.WAIC) for r in LK5.itertuples()}
LC = pd.read_csv(os.path.join(RS(), "S1c_likelihood5_coef.csv")); N["lik5c"] = {f"{r.model}|{r.cov}": [r.est, r.lo, r.hi] for r in LC.itertuples()}
RF = pd.read_csv(os.path.join(RS(), "S1b_refit_twice.csv")).iloc[0]; N["refit"] = dict(dwaic=float(RF.dWAIC), dcoef=float(RF.max_abs_dlogcoef))
RM = pd.read_csv(os.path.join(RS(), "S8_residual_moran.csv")); r2 = RM[RM.quantity == "pearson_residual_area_mean"].iloc[0]; N["resid_moran2"] = [float(r2.moran_I), float(r2.expected_I), float(r2.p_perm)]
GN = pd.read_csv(os.path.join(B, "results_gate", "S0_judgement_national.csv")); GP = pd.read_csv(os.path.join(B, "results_gate", "S0_judgement_province.csv")); GR = pd.read_csv(os.path.join(B, "results_gate", "S0_district_rankcorr.csv"))
gt = GN[GN.sex == "T"].set_index("year"); gp = GP[GP.sex == "T"].set_index("year")
AS = pd.read_csv(os.path.join(B, "results_gate", "S0_judgement_age_sex.csv")); AS = AS[(AS.sex != "T") & (AS.band != "계")]; asm = AS.groupby("year").ratio_nbkid_to_trace.median()
_rc = GR.set_index("pair").spearman
N["gate"] = dict(pct_2011=float(gt.loc[2011, "pct_nbkid"]), pct_2024=float(gt.loc[2024, "pct_nbkid"]), tr_2011=float(gt.loc[2011, "pct_trace"]), tr_2024=float(gt.loc[2024, "pct_trace"]),
                 p1_2011=float(gt.loc[2011, "pct_1plus"]), p1_2024=float(gt.loc[2024, "pct_1plus"]), ratio_min=float(gt.ratio_nbkid_to_trace.min()), ratio_max=float(gt.ratio_nbkid_to_trace.max()),
                 as_med_min=float(asm.min()), as_med_max=float(asm.max()), r_min=float(gp.r_nbkid_vs_trace.min()), r_max=float(gp.r_nbkid_vs_trace.max()),
                 rank_min=float(_rc.drop("2011-2024").min()), rank_max=float(_rc.drop("2011-2024").max()), rank_1124=float(_rc["2011-2024"]),
                 pct=dict((int(y), float(v)) for y, v in gt.pct_nbkid.items()), tr=dict((int(y), float(v)) for y, v in gt.pct_trace.items()))
for v in ("s577", "s2", "s8", "sc", "sp"):
    if not os.path.exists(os.path.join(RS(v), "S3_main_IRR.csv")): continue
    d = tab("S3_main_IRR.csv", v); N[v] = {k: get(d, VARS[k]) for k in VARS if (d["cov"] == VARS[k]).any()}
    f = pd.read_csv(os.path.join(RS(v), "S3_main_fit.csv")).iloc[0]; N[v + "_fit"] = dict(phi=float(f.phi), obs_moran=float(f.obs_moran_I), n=int(f.n_obs), waic=float(f.WAIC))
for t2, key in (("KIDS", "rho_kids"), ("NBBP", "rho_nbbp")):
    c2 = pd.read_csv(os.path.join(AN, t2, f"panel_{t2}.csv")).merge(pd.read_csv(os.path.join(PAN, f"{t2}_Estd.csv")), on=["year", "fbd_region"]); c2["fbd_region"] = c2.fbd_region.map(NFC)
    S2 = c2.groupby("fbd_region").apply(lambda d: d.cases.sum() / d.E_std.sum(), include_groups=False); j = pd.concat([L, S2], axis=1).dropna(); N[key] = [rnd(spearmanr(j[0], j[1]).correlation, 2), len(j)]
S6r = pd.read_csv(os.path.join(RS(), "S6_spatial_RR.csv")); N["rr_s"] = [rnd(S6r.RR.min(), 2), rnd(S6r.RR.max(), 2), rnd(S6r.RR.quantile(.1), 2), rnd(S6r.RR.quantile(.9), 2)]
EH = pd.read_csv(os.path.join(RS(), "S11_ehsa_summary.csv")); N["ehsa"] = {r.category: int(r.n) for r in EH.itertuples()}
WF = pd.read_csv(os.path.join(RS(), "S4b_window_fit.csv")); N["win_phi"] = {r.window: [rnd(r.phi, 2), rnd(r.phi_lo, 2), rnd(r.phi_hi, 2)] for r in WF.itertuples()}; N["win_hot"] = {r.window: int(r.n_hotspot) for r in WF.itertuples()}; N["win_cold"] = {r.window: int(r.n_coldspot) for r in WF.itertuples()}
N["chance_n"] = rnd(N["n_hyp"] * 0.05, 2); N["mv_n"] = 200
IT = pd.read_csv(os.path.join(RS(), "S4_interaction.csv")); IT["cov"] = IT["cov"].astype(str)
N["int"] = {f"{r.cov}|{r.window}": [float(float(r.IRR_ratio)), float(float(r.lo)), float(float(r.hi))] for r in IT.itertuples()}
S60 = pd.read_csv(os.path.join(RS(), "S60_extra.csv")); INV = {**{v: k for k, v in VARS.items()}, **{v: k for k, v in SVARS.items()}}
N["s60"] = {a: {INV.get(r.cov, r.cov): [float(float(r.IRR)), float(float(r.lo)), float(float(r.hi))] for r in g.itertuples()} for a, g in S60.groupby("analysis")}
N["s60_waic"] = {a: rnd(float(g.WAIC.iloc[0]), 1) for a, g in S60.groupby("analysis")}
def _pairs(g):
    L = open(g).read().split("\n")[1:]; return sum(int(l.split()[1]) for l in L if l.strip()) // 2
N["graph_pairs"] = _pairs(os.path.join(B, "graph_593.graph")); N["graph_pairs_q"] = _pairs(os.path.join(B, "results_s577", "nb_NBKID.graph"))
import re as _re
N["splus_k"] = int(_re.search(r"K (\d+)", open(os.path.join(B, "run_std3_88gprep.log")).read()).group(1))
import glob as _g
_lo = [pd.read_csv(f) for f in _g.glob(os.path.join(RS(), "S41_part_F[12]_*.csv"))]; N["loo_n"] = len(_lo); N["loo_walk_cred"] = int(sum(bool(d[d["cov"] == "walking_practice_rate"].cred.iloc[0]) for d in _lo))
json.dump(N, open(os.path.join(HERE, "numbers_nbkid.json"), "w"), ensure_ascii=False, indent=1, default=float)
if __name__ == "__main__":
    print({k: N[k] for k in ("cases_total_m", "exam_total_m", "rate", "sir_p10", "sir_p90", "p90p10", "gi_hot", "gi_cold", "gi_hot_sido", "gi_cold_sido", "hc", "kendall_w", "rho_2011_2024", "moran_min", "moran_max", "rm_spec", "vp", "ppc95", "oos_cov", "evalue", "rho_bl_dens", "rho_sir_urb", "rho_sir_bl")})
