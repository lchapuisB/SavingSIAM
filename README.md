# Should I Save *and* Should I Go? A Spatial Model of Climate Change Adaptation with Saving and Migration

**Author:** Léon Chapuis (Paris School of Economics & École normale supérieure – PSL)  
**Degree:** Master's Thesis — M2 *Analysis and Policy in Economics* (APE), Paris School of Economics (2026)  
**Supervisor:** Katheline Schubert · **Referee:** Marc Fleurbaey

---

## Abstract

I build a Dynamic Quantitative Spatial Model featuring heterogeneous agents capable of both saving and migrating across all 51 U.S. states to examine adaptation to a spatially heterogeneous climate change. This paper contributes to an emerging literature in climate macroeconomics by providing a tractable structural framework to jointly evaluate saving and migration at the individual worker level within a dynamic spatial model. I build on an existing framework which incorporates a CARA preference structure that separates workers' migration and saving choices, ensuring tractability despite the curse of dimensionality. The main result demonstrates that climate damages induce significant spatial reallocation and precautionary saving, which smooth the drop in consumption and restrict welfare losses to 0.5%, while altering the spatial distribution of wealth. Anticipation helps to speed the transition and limits welfare drop on the transition path. Counterfactual scenarios reveal both substitutability and complementarity between saving and migration; suppressing either of these adaptive margins amplifies welfare losses, showing the interaction between wealth accumulation and spatial mobility in adapting to climate damages.

**JEL classification:** Q54, R12, R23, E21, J61  
**Keywords:** Climate Change, Migration, Saving, Spatial Economics, Dynamic Model

---

## Repository Structure

- **`Final Versions/For the Master/`**  
  Final compiled Master's thesis PDF: [`Chapuis L. - Should I Save and Should I go?.pdf`](Final%20Versions/For%20the%20Master/Chapuis%20L.%20-%20Should%20I%20Save%20and%20Should%20I%20go%3F.pdf).
- **`Draft/`**  
  LaTeX manuscript source (`main.tex`, `appendix_data.tex`, `references.bib`, `main.bbl`), compiled draft (`main.pdf`), mathematical derivations (`mathematical_documentation.md`), figures (`figures/`), and tables (`tables/`).
- **`Presentation/`**  
  Defense presentation slides (`presentation.pdf`, `presentation.tex`, `make_presentation.py`, `hsiang_map.png`).
- **`Script/`**  
  Replication code in R and MATLAB:
  - `data_cleaning.R`: Cleans raw state- and county-level datasets and computes climate damage shocks across wages, amenities, and housing costs.
  - `main.m`: Main MATLAB orchestrator calling calibration, steady-state computation, transition-path solvers, counterfactuals, and figure generation.
  - `static/`: Baseline 2010 steady-state model and parameter calibration (`model_static.m`, `calibrate_param.m`, `calibrate_param_nu.m`, `compute_Xi.m`, `compute_SS_moments.m`).
  - `dynamics/`: Transition path solvers under perfect foresight and myopia (`solve_transition_backward.m`, `solve_transition_myopia.m`, `simulate_Lambda_forward.m`, `compute_atkinson_welfare.m`, `compute_consumption_equivalents.m`).
  - `scenarios/`: Baseline and counterfactual scenario runners (`run_scenario.m`, `run_did.m`, `compare_scenarios.m`, `compare_all_audit.m`).
  - `no_saving/`: Hand-to-mouth structural counterfactual model without saving.
  - `plots/` & `helpers/`: Plotting scripts and numerical utilities.
- **`Data/`**  
  Raw state- and county-level datasets (BEA state GDP & personal income, ACS, IRS interstate migration flows `stateinflow1516.csv`, county-level climate damages `county_damages_by_sector.csv`, natural and urban amenities in `amenities/`, state wealth inequality, and border matrix).  
  > **Note on omitted large files:** `Data/cps.csv.gz` (164 MB IPUMS CPS microdata extract) is excluded to comply with GitHub's 100 MB file-size limit. Its processed outputs (including `Output/data/national_variances_2080.csv` and the state-level shock/calibration CSVs) are included directly in `Output/data/`.
- **`Output/`**  
  Cleaned calibration inputs (`Output/data/*.csv`) and generated simulation figures across baseline and counterfactual scenarios. *(Large pre-computed `.mat` simulation caches are excluded from version control and can be regenerated via `Script/main.m`.)*

---

## Replication Instructions

1. **Data preparation (optional — processed inputs are already included in `Output/data/`):**  
   To regenerate the cleaned CSV inputs from raw files, obtain the IPUMS CPS extract (`Data/cps.csv.gz`), update the local file paths in `Script/data_cleaning.R`, and run it in R.
2. **Model calibration, simulations, and figures (MATLAB):**  
   Open `Script/main.m` in MATLAB, set `RUN_SCENARIO = true;` on line 13 (so that the model calibrates and computes all steady states, transition paths, and counterfactuals rather than loading excluded `.mat` caches), and run `Script/main.m`.
