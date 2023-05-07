<a href="https://zenodo.org/badge/latestdoi/330802321"><img src="https://zenodo.org/badge/330802321.svg" alt="DOI"></a>

# spectral_slopes_aging
This repository contains data and scripts for:  
Pathania, A., Euler, M.J., Clark. M., Cowan, R., Duff, K., & Lohse, K.R. (in press). Resting EEG spectral slopes are associated with age-related differences in information processing speed. Biological Psychology.

1. Aging_study_updated_20211022: R code consistent with resubmission to Biological Psychology. Primary changes from the pre-print are: (1) using spectral parameterization rather than regression based methods to extract the spectral slope and (2) analyzing 5 RBANS domains rather than all 12 subtets.
2. MASTER_EO_and_EC_EEG_KRL: Raw power spectra that are used for several plots and reduced into MASTER_FOOOF.
3. MASTER_FOOF: power spectra for eyes-closed resting data aggregated across region. These data are fed into Python for spectral parameterization.
4. Aging_temp.ipynb: Jupyter Notebook file with Python code to conduct spectral parameterization on MASTER_FOOF
5. Aging_study_FOOOF_var: output of spectral parameterization showing exponents, offsets, and model fit statistics for each subject in each region. 
6. RBANS_aging_study_08262021: spreadsheet containing demographic and behavioral data for all participants
7. descriptive_group_stats: descriptive statistics for each group that are summarized in Table 1 of the paper.
