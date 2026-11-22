# Do kidney risk factors explain why trace dipstick proteinuria clusters in the central urban districts of Korea?

Analysis code for a study of where the "normal B (borderline), kidney-function management" judgement (urine dipstick protein trace) is concentrated across the 229 districts of South Korea at the national general health screening, 2011-2024, and whether district kidney risk factors explain the pattern.

This code is released under the MIT licence. Version 1.0.8. An archived copy of each version is deposited on Zenodo (version 1.0.8: https://doi.org/10.5281/zenodo.24450374; all versions: concept record https://doi.org/10.5281/zenodo.24419401).

## Model

The principal model (`S1_narrative_full`, version `std3`) contains 9 hypothesised district characteristics (screening-measured suspected diabetes share, current smoking, high-risk drinking, treatment rates among adults with diagnosed hypertension and diabetes, farm-household population share, public assistance recipient rate, screening uptake, walking) and 5 fixed adjustments (urban population share, physicians, sex ratio, fiscal independence, share aged 80 or older), with a BYM2 spatial term, a first-order random walk in time and an unstructured year term (negative binomial likelihood, offset = log expected count).

Set choice. `81_fit_inla.R` (step S2) fits every candidate set in `sets_NBKID_std3.json` with the forced adjustments and applies the rule coded there and recorded in the settings file (`rules.set_choice`): lowest WAIC; within 4 WAIC points, the set with fewer hypothesised characteristics (ties by LCPO). Under this rule the lowest WAIC was that of `S2_core` (5 characteristics: suspected diabetes share, hypertension treatment rate, public assistance recipient rate, screening uptake, walking; WAIC 45,004.2); `S5_work_deprivation` (2 characteristics: farm-household population share, public assistance recipient rate; WAIC 45,006.3) lies within 4 WAIC points with fewer characteristics, so the rule would have chosen `S5_work_deprivation`. The principal model is the full 9-characteristic set `S1_narrative_full` (WAIC 45,009.9), hashed (SHA-256) and recorded in the analysis log (earlier sets were revised after collinearity, series-break and data-quality checks; `rules.main_model_note`, read as `main_set` by `_st_config.py`; `81_fit_inla.R` then overrides the rule-selected set and writes both, `best_set` and `rule_best_set`, to `HM_best_NBKID_std3.txt`). `S5_work_deprivation`, `S2_core` and the other candidate sets are reported as sensitivity sets.

Sensitivity analyses only (never in the principal model): annual mean temperature, screening-measured suspected hypertension share, survey obesity rate and the period-standardised screening obesity-management share, each added one at a time (`88V_nbkid_extra.R`); refits without Sejong, without the districts with boundary changes (Sejong, Cheongju, Gunwi), without 2020-2021 and from 2015; BYM2 mixing-parameter priors Beta(1,1) and Jeffreys Beta(0.5,0.5) in place of the PC prior. Further sensitivity versions and the full sensitivity suite are listed below.

Data corrections in `build_pool_from_sqlite_nbkid.py`: the Sejong farm population share, observed for 2011-2013 only, is held at its 2013 value. Internal-medicine specialists are not used in any reported model: the reference value had been divided by population twice, and the source head counts are averages rather than totals for cities with general districts, with a 2018 break and 13 cities recorded as 0 in 2024 (the recomputed column remains in the pool code for transparency).

## Contents

This repository contains code only. No data, adjacency graph, boundary file, result tables or CSV files are included. The processed district-year data are available from the author on reasonable request; all source statistics are public (see `DATA_DICTIONARY.md`).

| Path | What it is |
|---|---|
| `01_스크립트/` | Analysis scripts and `_var_labels.json` (variable labels) |
| `04_분석/NBKID/` | Covariate pool builder, settings writer (`make_pool_sets_nbkid.py`), model settings hashed (SHA-256) and recorded in the analysis log (`sets_NBKID_std3.json` principal; `s577` queen contiguity graph, `s2` without the share aged 80 or older, `s8` age-only standardised expected counts, `sc` unstandardised expected counts, `sp` population-based expected counts) and the two run scripts |
| `04_분석/NBKID/MS/build/` | `numbers_nbkid.py` (every number in the manuscript) and `figs_nbkid.py` (Figures 1-3) |
| `DATA_DICTIONARY.md` | Variables, units and sources |
| `OUTPUT_MAP.md` | Which script writes which output table |

Folder names and some column names are in Korean because the source tables are Korean; every script starts with an English header comment (purpose, inputs, outputs).

## Environment variables

| Variable | Meaning |
|---|---|
| `ST_ROOT` | Working root that holds `00_원천/` (downloaded source tables) and `02_패널/` (district-year tables); defaults to the repository root |
| `ST_DATA_ROOT` | Folder holding the district-year population panel `식중독최종.csv` (229 districts) |
| `ST_GPT_DATASET` | Path of the reference district-year dataset (SQLite file, opened read-only) |
| `ST_FBDZIP` | Folder of resident population by sex and 5-year age group, one file per year (`<year>.csv`) |
| `ST_SHP` | District boundary shapefile (229 districts, 2024 boundaries); default `data/boundary/final.shp` |
| `ST_GRAPH` | Principal adjacency graph in INLA graph format (229 districts, boundary-file order) |
| `KOSIS_API_KEY` | The user's own KOSIS Open API key (only for the download script; no key is included) |
| `ST_GATE_DIR` | Optional folder with external data-check scripts (`data_gate.py`, `panel_gate.py`); the checks are skipped when it is not set |
| `ST_VER`, `TAG` | Model version (`std3`, `s577`, `s2`, `s8`, `sc`, `sp`) and outcome tag (`NBKID`), set by the run scripts |

## Pipeline

Run from `ST_ROOT` after setting the environment variables above.

| Step | Script | Output |
|---|---|---|
| 1 | `01_스크립트/99G_nbkid_fetch_N064.py`, `99G_nbkid_judgement_gate.py` | Judgement check: published judgement versus the published urine dipstick distribution by year, sex, age and province; year-to-year district rank correlation |
| 2 | `01_스크립트/99NK2_build_Estd3_nbkid.py` (district key list read from the master panel `02_패널/MASTER_시군구_연도_wide.parquet` written by `50_build_master_panel.py`, with codes from `40_region_crosswalk.py`; resident population from `ST_FBDZIP`; share aged 80 or older from `02_패널/HTND_Estd2.csv`, see Inputs) | District counts (code-based linkage; Sejong from province rows) and expected counts by indirect standardisation in 2 sexes x 15 age groups |
| 3 | `04_분석/NBKID/build_pool_from_sqlite_nbkid.py`, `make_pool_sets_nbkid.py` | Covariate pool read from the reference district-year dataset; model settings `sets_NBKID_*.json`, hashed (SHA-256) and recorded in the analysis log; earlier sets were revised after collinearity, series-break and data-quality checks |
| 4 | `04_분석/NBKID/run_A_gate.sh` (`80_build_disease_panel_v.py`) | Analytic panel and optional data checks |
| 5 | `04_분석/NBKID/run_std3_chain.sh` (`81_fit_inla.R`, `81R_nbkid_refit_binom.R`, `85_window_maps.R`, `86_yearly_maps_ehsa.py`, `87_window_grid.py`, `88g_sens_prep.py`, `88g_sens_suite.R`, `88h_ppc_varpart.R`, `88k_struct_kh.R`, `88k_struct_kh_anem.R`, `89g_sens_combine_nbkid.py`) | Principal model, 200-subset multiverse, period models, Knorr-Held models, restricted spatial regression, Spatial+, alternative specifications, sensitivity versions (s577, s2, s8, sc, sp) |
| 6 | `01_스크립트/88S_nbkid_reverse_posthoc.R` | Post hoc analyses of the association opposite to the hypothesised direction |
| 7 | `01_스크립트/88V_nbkid_extra.R` | Additional sensitivity analyses: sensitivity-only variables added one at a time, district and year exclusions, mixing-parameter priors (`S60_extra.csv`) |
| 8 | `04_분석/NBKID/MS/build/numbers_nbkid.py`, `figs_nbkid.py` | Every number reported in the manuscript (`numbers_nbkid.json`) and Figures 1-3 |

## Inputs not written by this repository

| File (under `ST_ROOT`) | Read by | Content |
|---|---|---|
| `02_패널/HTND_Estd2.csv` | `99NK2_build_Estd3_nbkid.py` | `pop80_share_raw`, district share of residents aged 80 or older by year (resident registration statistics), used as a fixed adjustment |
| `04_분석/KIDS/panel_KIDS.csv`, `02_패널/KIDS_Estd.csv` | `make_pool_sets_nbkid.py` (check only), `88g_sens_prep.py`, `numbers_nbkid.py` | Comparison outcome: examinees judged suspected kidney disease, district counts and expected counts (KOSIS NHIS screening tables) |
| `04_분석/NBBP/panel_NBBP.csv`, `02_패널/NBBP_Estd.csv` | as above | Comparison outcome: examinees judged normal B, blood-pressure management, district counts and expected counts |

The comparison outcomes enter only the negative-control-outcome sensitivity analysis and the district SIR correlations reported in the manuscript; they are built in the same way as the NBKID counts and expected counts (step 2) from their own KOSIS tables.

Expected counts and offsets use no added constant (offset = log E; the panel has no zero counts).

R jobs are wrapped in `rslot.sh` and run with one INLA thread (`ST_THREADS=1`).

## Software

R 4.6.0 with R-INLA 25.10.19, sf, spdep, jsonlite, ggplot2 and dplyr; Python 3.9 with pandas, numpy, scipy, statsmodels, geopandas, matplotlib, pyarrow and Pillow. Newer INLA versions may shift DIC/WAIC by a few points without changing the rate ratios.

## Licence

MIT (see `LICENSE`).
