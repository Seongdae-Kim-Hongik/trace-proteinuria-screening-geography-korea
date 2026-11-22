# Output map

Result tables are written to `04_분석/NBKID/results_<version>/` (principal version `std3`; sensitivity versions `s577`, `s2`, `s8`, `sc`, `sp`) unless another folder is shown.

| Script | Output | Content |
|---|---|---|
| `99G_nbkid_judgement_gate.py` | results_gate/S0_judgement_national.csv, S0_judgement_age_sex.csv, S0_judgement_province.csv, S0_district_rankcorr.csv | Judgement check (Supplementary) |
| `99NK2_build_Estd3_nbkid.py` | 02_패널/NBKID_Estd3.csv | Counts, denominators and expected counts |
| `build_pool_from_sqlite_nbkid.py` | 02_패널/EXT_canonical_pool_NBKID_261121.csv | Covariate pool |
| `make_pool_sets_nbkid.py` | 04_분석/NBKID/sets_NBKID_<version>.json, sets_NBKID_hash_before_fit.txt | Model settings, fixed and hashed before the reported models were fitted |
| `80_build_disease_panel_v.py` | panel_NBKID_<version>.csv, panel_META_*, config_* | Analytic panel, imputation counts, VIF |
| `81_fit_inla.R` | results_<version>/S1_likelihood.csv, S3_main_IRR.csv, S3_main_fit.csv, S3_negcontrol.csv, S6_spatial_RR.csv, S8_residual_moran.csv, S4_window_IRR.csv, S5_covid.csv, S6_moran_year.csv, S7_multiverse_*.csv, S9_typeI_IRR.csv, S10_offset_sensitivity_IRR.csv | Principal model, period models, multiverse |
| `81R_nbkid_refit_binom.R` | S1b_refit_twice.csv, S1c_likelihood5.csv, S1c_likelihood5_coef.csv, S0_repro_meta.csv | Reproducibility and five likelihoods |
| `85_window_maps.R` | S4b_window_spatialRR.csv, S4b_window_fit.csv, S4b_window_moran.csv | Period hot and cold districts |
| `86_yearly_maps_ehsa.py` | S11_ehsa.csv, S11_ehsa_summary.csv, S11b_moran_yearly_EB.csv, S12_persistence.csv, S12b_persistence_district.csv, S12c_yearly_EB_SIR.csv | Emerging hot spots, yearly Moran's I, persistence |
| `87_window_grid.py` | S4c_window_grid_table.csv | Period rate-ratio grid |
| `88g_sens_prep.py / 88g_sens_suite.R` | S41_panel_*.csv, S41_part_<analysis>.csv, S41_hot_<analysis>.csv, S41_oos2024.csv | Sensitivity suite |
| `88h_ppc_varpart.R` | S42_ppc.csv, S42_pit_hist.csv, S42_variance_partition.csv | Posterior predictive checks, variance partition |
| `88k_struct_kh.R / 88k_struct_kh_anem.R` | S43_part_<analysis>.csv, S43_sensitivity_IRR.csv, S50_*.csv | Structural sensitivity, Knorr-Held, restricted spatial regression |
| `88S_nbkid_reverse_posthoc.R` | S46_nbkid_posthoc.csv | Post hoc analyses (opposite direction) |
| `88V_nbkid_extra.R` | S60_extra.csv | Sensitivity-only variables, district and year exclusions, mixing-parameter priors |
| `89g_sens_combine_nbkid.py` | S41_retention_by_cov.csv, S41_noise_rubin.csv, S41_nco_detail.csv, S41_evalues.csv, S41_leaveout_ranges.csv | Sensitivity summary |
| `numbers_nbkid.py` | MS/build/numbers_nbkid.json, gi_nbkid.csv | Every number in the manuscript and supplementary tables |
| `figs_nbkid.py` | MS/build/figs/Fig1-Fig3 | Figures 1-3 |
