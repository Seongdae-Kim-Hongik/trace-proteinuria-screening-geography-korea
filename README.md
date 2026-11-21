# Do kidney risk factors explain why trace dipstick proteinuria clusters in the central urban districts of Korea?

Analysis code for a study of where the "normal B (borderline), kidney-function management" judgement (urine dipstick protein trace) is concentrated across the 229 districts of South Korea at the national general health screening, 2011-2024, and whether district kidney risk factors explain the pattern.

This code is released under the MIT licence. Version 1.0.3. An archived copy of each version is deposited on Zenodo (version 1.0.3: https://doi.org/10.5281/zenodo.24430432; all versions: concept record https://doi.org/10.5281/zenodo.24419401).

## Model

The principal model (`S1_narrative_full`, version `std3`) contains 10 hypothesised district characteristics (screening-measured suspected diabetes share, current smoking, high-risk drinking, treatment rates among adults with diagnosed hypertension and diabetes, internal-medicine specialists, farm-household population share, public assistance recipient rate, screening uptake, walking) and 5 fixed adjustments (urban population share, physicians, sex ratio, fiscal independence, share aged 80 or older), with a BYM2 spatial term, a first-order random walk in time and an unstructured year term (negative binomial likelihood, offset = log expected count).

Set choice. `81_fit_inla.R` (step S2) fits every candidate set in `sets_NBKID_std3.json` with the forced adjustments and applies the rule coded there and recorded in the settings file (`rules.set_choice`): lowest WAIC; within 4 WAIC points, the set with fewer hypothesised characteristics (ties by LCPO). Under this rule the 5-characteristic core set `S2_core` (suspected diabetes share, hypertension treatment rate, public assistance recipient rate, screening uptake, walking; WAIC 45,000.9) would have been chosen over the full set (WAIC 45,009.7). The principal model was fixed before fitting to the full 10-characteristic set `S1_narrative_full` (`rules.main_model_note`, read as `main_set` by `_st_config.py`; `81_fit_inla.R` then overrides the rule-selected set and writes both, `best_set` and `rule_best_set`, to `HM_best_NBKID_std3.txt`). `S2_core` and the other candidate sets are reported as sensitivity sets.

Sensitivity analyses only (never in the principal model): annual mean temperature, screening-measured suspected hypertension share, survey obesity rate and the period-standardised screening obesity-management share, each added one at a time (`88V_nbkid_extra.R`); refits without Sejong, without the districts with boundary changes (Sejong, Cheongju, Gunwi), without 2020-2021 and from 2015; BYM2 mixing-parameter priors Beta(1,1) and Jeffreys Beta(0.5,0.5) in place of the PC prior. Further sensitivity versions and the full sensitivity suite are listed below.

## Contents

This repository contains code only. No data, adjacency graph, boundary file, result tables or CSV files are included. The processed district-year data are available from the author on reasonable request; all source statistics are public (see `DATA_DICTIONARY.md`).

| Path | What it is |
|---|---|
| `01_스크립트/` | Analysis scripts and `_var_labels.json` (variable labels) |
| `04_분석/NBKID/` | Covariate pool builder, settings writer (`make_pool_sets_nbkid.py`), model settings fixed before fitting (`sets_NBKID_std3.json` principal; `s577` queen contiguity graph, `s2` without the share aged 80 or older, `s8` age-only standardised expected counts, `sc` unstandardised expected counts, `sp` population-based expected counts) and the two run scripts |
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
| 2 | `01_스크립트/99NK2_build_Estd3_nbkid.py` (district key and master panel from `40_region_crosswalk.py`, `50_build_master_panel.py`; resident population from `ST_FBDZIP`) | District counts (code-based linkage; Sejong from province rows) and expected counts by indirect standardisation in 2 sexes x 15 age groups |
| 3 | `04_분석/NBKID/build_pool_from_sqlite_nbkid.py`, `make_pool_sets_nbkid.py` | Covariate pool read from the reference district-year dataset; model settings `sets_NBKID_*.json`, fixed (and hashed) before any covariate model was fitted |
| 4 | `04_분석/NBKID/run_A_gate.sh` (`80_build_disease_panel_v.py`) | Analytic panel and optional data checks |
| 5 | `04_분석/NBKID/run_std3_chain.sh` (`81_fit_inla.R`, `81R_nbkid_refit_binom.R`, `85_window_maps.R`, `86_yearly_maps_ehsa.py`, `87_window_grid.py`, `88g_sens_prep.py`, `88g_sens_suite.R`, `88h_ppc_varpart.R`, `88k_struct_kh.R`, `88k_struct_kh_anem.R`, `89g_sens_combine_nbkid.py`) | Principal model, 200-subset multiverse, period models, Knorr-Held models, restricted spatial regression, Spatial+, alternative specifications, sensitivity versions (s577, s2, s8, sc, sp) |
| 6 | `01_스크립트/88S_nbkid_reverse_posthoc.R` | Post hoc analyses of the association opposite to the prespecified direction |
| 7 | `01_스크립트/88V_nbkid_extra.R` | Additional sensitivity analyses: sensitivity-only variables added one at a time, district and year exclusions, mixing-parameter priors (`S60_extra.csv`) |
| 8 | `04_분석/NBKID/MS/build/numbers_nbkid.py`, `figs_nbkid.py` | Every number reported in the manuscript (`numbers_nbkid.json`) and Figures 1-3 |

Expected counts and offsets use no added constant (offset = log E; the panel has no zero counts).

R jobs are wrapped in `rslot.sh` and run with one INLA thread (`ST_THREADS=1`).

## Software

R 4.6.0 with R-INLA 25.10.19, sf, spdep, jsonlite, ggplot2 and dplyr; Python 3.9 with pandas, numpy, scipy, statsmodels, geopandas, matplotlib, pyarrow and Pillow. Newer INLA versions may shift DIC/WAIC by a few points without changing the rate ratios.

## Licence

MIT (see `LICENSE`).
