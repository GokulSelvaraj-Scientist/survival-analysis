# Survival Analysis: Lung Cancer Clinical Data

## Overview
This project performs survival analysis on a publicly available lung cancer clinical dataset using R. It includes Kaplan-Meier survival curves, log-rank tests, and Cox proportional hazards regression to identify clinical and demographic factors associated with patient survival outcomes.

## Why This Matters
Survival analysis is one of the most critical analytical tools in clinical drug development. Every Phase III oncology trial relies on survival endpoints — overall survival (OS) and progression-free survival (PFS) — as primary measures of drug efficacy. Understanding which patient factors predict survival is essential for:

- **Clinical trial design** — identifying prognostic factors helps define eligibility criteria and stratification variables that reduce trial variability and increase statistical power
- **Patient stratification** — Cox regression identifies independent predictors that can be used to enrich trials with patients most likely to benefit from treatment
- **Endpoint selection** — understanding baseline survival patterns informs whether OS or PFS is the more appropriate primary endpoint for a given patient population
- **Regulatory submissions** — survival analyses are required by FDA and EMA to demonstrate clinical benefit in oncology drug applications
- **Biomarker development** — performance status scores like ECOG — shown here to be a strong independent predictor — are routinely used as inclusion/exclusion criteria and stratification factors in cancer trials

In a real drug development setting, this type of analysis would be applied to Phase II/III trial data to characterize treatment benefit across patient subgroups and support labeling claims.

## Dataset
- **Source:** `lung` dataset from the `survival` R package (NCCTG Lung Cancer Study)
- **Type:** Clinical survival data
- **Samples:** 228 advanced lung cancer patients
- **Variables:** Survival time, vital status, age, sex, ECOG performance score, Karnofsky scores, weight loss

## Methods
- Kaplan-Meier survival estimation
- Log-rank test for group comparisons
- Cox proportional hazards regression for multivariate analysis
- Forest plot visualization of hazard ratios

## Outputs
| File | Description |
|---|---|
| `km_overall.png` | Overall Kaplan-Meier survival curve with confidence intervals |
| `km_by_sex.png` | Survival curves stratified by sex with p-value |
| `km_by_ecog.png` | Survival curves stratified by ECOG performance score |
| `cox_forest_plot.png` | Forest plot of Cox model hazard ratios |
| `survival_summary_by_sex.csv` | Median survival times by sex with confidence intervals |

## Key Findings
- Significant survival difference between male and female patients (p = 0.0016) — females survive longer, consistent with published literature on sex differences in lung cancer biology
- ECOG performance score is the strongest independent predictor of survival (HR > 1.67 per unit increase, p < 0.0001) — patients with limited activity have dramatically worse outcomes
- Age and weight loss show weaker independent effects after adjusting for ECOG and sex
- Median survival is approximately 259 days across the full cohort, consistent with historical data for advanced NSCLC prior to targeted therapy era

## Biological and Clinical Interpretation
The strong effect of ECOG performance status is clinically significant — it reflects the patient's overall physiological reserve and ability to tolerate treatment. In modern oncology trials, ECOG 0-1 is typically required for enrollment in aggressive treatment arms. The sex difference in survival likely reflects both biological differences (higher rates of EGFR mutations in female never-smokers) and behavioral differences in smoking history. These findings are directly applicable to designing stratification schemes in clinical trials to ensure balanced treatment arms.

## How to Run
1. Install R and RStudio
2. Install required packages:
```r
install.packages(c("survival", "survminer", "ggplot2", "dplyr", "readr"))
```
3. Run `survival_analysis.R` in RStudio

## Requirements
- R >= 4.0
- survival, survminer, ggplot2, dplyr, readr

## Author
**Gokul Selvaraj, PhD**
GitHub: [GokulSelvaraj-Scientist](https://github.com/GokulSelvaraj-Scientist)
