# Data dictionary

Unit of analysis: district (si-gun-gu, 229 districts, 2024 boundaries) x year, 2011-2024 (3206 district-years). In the analytic panel every covariate is z-standardised over all district-years (rate ratios are per 1 SD); the model columns carry the prefix `ext__`.

| Variable | Role | Description | Unit | Source |
|---|---|---|---|---|
| `nbkid_cases` | outcome | Examinees judged normal B, kidney-function management (urine dipstick protein trace) | count | NHIS general health screening statistics, KOSIS DT_35007_N101 |
| `nbkid_examinees` | denominator | Examinees judged at the general health screening | count | KOSIS DT_35007_N098 |
| `E_std` | offset | Expected count, indirect standardisation by 2 sexes x 15 age groups (offset = log E) | count | Derived: KOSIS DT_35007_N099, DT_35007_N102; resident registration population |
| `E_std_age / E_crude / E_pop` | offset (sensitivity) | Age-only standardised / unstandardised / population-based expected counts | count | Derived as above |
| `susp_htn_share` | sensitivity only | Screening-measured suspected hypertension share | % of examinees | NHIS screening statistics |
| `susp_dm_share` | covariate | Screening-measured suspected diabetes share | % of examinees | NHIS screening statistics |
| `obesity_rate` | sensitivity only | Obesity rate (self-reported) | % | Community Health Survey |
| `smoking_rate` | covariate | Current smoking rate | % | Community Health Survey |
| `high_risk_drinking` | covariate | High-risk drinking rate | % | Community Health Survey |
| `htn_treat_rate` | covariate | Treatment rate among adults with diagnosed hypertension | % | Community Health Survey |
| `dm_treat_rate` | covariate | Treatment rate among adults with diagnosed diabetes (age 30+) | % | Community Health Survey |
| `im_specialists_per100k` | not used | Internal-medicine specialists: excluded from all reported models because of source-table defects (double division by population in the reference dataset, city averages for cities with general districts, 2018 break, 2024 zeros) | per 100 000 residents | Health-resource statistics (KOSIS) |
| `annual_mean_temp` | sensitivity only | Annual mean temperature (inverse-distance interpolation of station data) | deg C | Korea Meteorological Administration |
| `farm_pop_share` | covariate | Farm-household population share (urban wards set to 0) | % | Agricultural census (KOSIS) |
| `basic_livelihood_rate` | covariate | Basic livelihood security (public assistance) recipient rate | % | Public assistance statistics (KOSIS) |
| `screening_uptake` | covariate | General health screening uptake | % | NHIS screening statistics |
| `walking_rate` | covariate | Walking practice rate | % | Community Health Survey |
| `urban_population_share` | fixed adjustment | Urban-zone population share | % | KOSIS local indicators |
| `physicians_per_1000` | fixed adjustment | Physicians | per 1000 residents | Health-resource statistics (KOSIS) |
| `sex_ratio` | fixed adjustment | Sex ratio | men per 100 women | Resident registration statistics |
| `fiscal_independence` | fixed adjustment | Fiscal autonomy (own revenue and general grants as a share of the budget) | % | Local finance statistics (KOSIS) |
| `pop80_share` | fixed adjustment | Population aged 80 or older | % | Resident registration statistics |
| `hit_and_run_rate` | negative control | Hit-and-run traffic accident rate | % | KOSIS local indicators |
| `traffic_culture_index` | negative control | Traffic culture index | index | KOSIS local indicators |
| `fires_per_10000` | negative control | Fires | per 10 000 residents | Fire statistics (KOSIS) |
| `tertiary_educated_share` | sensitivity only | Tertiary-educated share | % | Population and housing census (KOSIS) |
| `living_alone_elderly` | sensitivity only | Older single-person household share | % | Population and housing census (KOSIS) |
| `moderate_pa_rate` | sensitivity only | Moderate-to-vigorous physical activity rate | % | Community Health Survey |
| `public_sport_fac_per100k` | sensitivity only | Public sports facilities | per 100 000 residents | Public sports facility register |
| `scr_obesity_pstd` | sensitivity only | Screening obesity-management (normal B) share, z-standardised within 2011-2018 and 2019-2024 (judgement change between 2018 and 2019) | z | NHIS screening statistics |
| `unmet_need_rate` | sensitivity only | Unmet medical need (standardised rate) | % | Community Health Survey (KOSIS local indicators) |

Principal covariate set `S1_narrative_full` = the 9 covariates marked `covariate` above; fixed adjustments are in every model. Comparison outcomes in the sensitivity suite: examinees judged suspected kidney disease (`KIDS`) and normal B, blood-pressure management (`NBBP`), read from their own district-year panels.

Source tables for the outcome (KOSIS): DT_35007_N101 (district, by sex), DT_35007_N098 (district examinees judged), DT_35007_N099 and DT_35007_N102 (national, by sex and age), DT_35007_N064 and DT_35007_N065 (urine dipstick protein distribution, judgement check). Judgement criteria: Ministry of Health and Welfare notice on the general health screening, annexed tables (2014, 2017 and 2024 editions); urine dipstick trace is the only component of this judgement.

Labels (English and Korean) for the 27 source columns used by this analysis: `01_스크립트/_var_labels.json`.
