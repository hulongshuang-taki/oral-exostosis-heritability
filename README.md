# Oral exostosis heritability analysis

This repository contains the R code used for the parent-offspring analysis of
oral exostosis in the Tohoku Medical Megabank Organization cohort.

## Scripts

- `00_prepare_wave_data.R`: prepares the Wave 2-priority sensitivity data set
- `01_descriptive_statistics.R`: descriptive characteristics of the study sample
- `02_prevalence_analysis.R`: age- and sex-specific prevalence analyses
- `03_heritability_analysis.R`: tetrachoric and Falconer heritability estimates
- `04_multivariable_models.R`: logistic regression and stratified analyses

The scripts were run with R 4.4.3. The heritability analysis uses `polycor`
version 0.8-1.

## Data availability

Individual-level data are not included because they contain sensitive
participant information and are governed by the data access policies of the
Tohoku Medical Megabank Organization (ToMMo). Qualified researchers may apply
for access through ToMMo's established data distribution procedures:
https://www.megabank.tohoku.ac.jp/researchers/utilization/dist

## Contact

Longshuang Hu  
Tohoku University  
GitHub: [hulongshuang-taki](https://github.com/hulongshuang-taki)
