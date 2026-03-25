# ============================================================
# Survival Analysis: TCGA Lung Adenocarcinoma (LUAD)
# Author: Gokul Selvaraj
# GitHub: GokulSelvaraj-Scientist
# Description: Kaplan-Meier survival analysis and Cox proportional
#              hazards regression using TCGA clinical data
# ============================================================

# --- Load Libraries ---
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(readr)

# --- Install if needed ---
# install.packages(c("survival", "survminer", "ggplot2", "dplyr", "readr"))

# ============================================================
# We use the built-in lung cancer dataset from the survival
# package — publicly available, no download needed
# lung: NCCTG Lung Cancer dataset
# Variables:
#   time    - survival time in days
#   status  - 1=censored, 2=dead
#   age     - age in years
#   sex     - 1=male, 2=female
#   ph.ecog - ECOG performance score (0=good, 4=bedridden)
#   ph.karno, pat.karno - Karnofsky performance scores
#   meal.cal - calories consumed at meals
#   wt.loss  - weight loss in last 6 months
# ============================================================

data(lung)

# --- Data Cleaning ---
lung_clean <- lung %>%
  filter(!is.na(time), !is.na(status)) %>%
  mutate(
    sex_label  = ifelse(sex == 1, "Male", "Female"),
    age_group  = ifelse(age >= median(age, na.rm = TRUE), "Older", "Younger"),
    ecog_group = case_when(
      ph.ecog == 0 ~ "Fully active",
      ph.ecog == 1 ~ "Restricted",
      ph.ecog >= 2 ~ "Limited/Bedridden",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(ecog_group))

cat("Sample size after cleaning:", nrow(lung_clean), "\n")
cat("Median survival time:", median(lung_clean$time), "days\n")

# --- Plot 1: Overall Kaplan-Meier Survival Curve ---
surv_obj <- Surv(time = lung_clean$time, event = lung_clean$status == 2)
fit_overall <- survfit(surv_obj ~ 1, data = lung_clean)

overall_km <- ggsurvplot(
  fit_overall,
  data          = lung_clean,
  conf.int      = TRUE,
  risk.table    = TRUE,
  xlab          = "Time (days)",
  ylab          = "Survival Probability",
  title         = "Overall Kaplan-Meier Survival Curve",
  subtitle      = "NCCTG Lung Cancer Dataset",
  palette       = "#2C7BB6",
  ggtheme       = theme_classic(base_size = 13),
  risk.table.height = 0.25
)

ggsave("km_overall.png",
       print(overall_km, newpage = FALSE),
       width = 9, height = 7, dpi = 300)
cat("Saved: km_overall.png\n")

# --- Plot 2: Survival by Sex ---
fit_sex <- survfit(surv_obj ~ sex_label, data = lung_clean)

sex_km <- ggsurvplot(
  fit_sex,
  data          = lung_clean,
  conf.int      = TRUE,
  pval          = TRUE,
  risk.table    = TRUE,
  xlab          = "Time (days)",
  ylab          = "Survival Probability",
  title         = "Survival by Sex",
  subtitle      = "NCCTG Lung Cancer Dataset",
  legend.title  = "Sex",
  palette       = c("#E63946", "#457B9D"),
  ggtheme       = theme_classic(base_size = 13),
  risk.table.height = 0.25
)

ggsave("km_by_sex.png",
       print(sex_km, newpage = FALSE),
       width = 9, height = 7, dpi = 300)
cat("Saved: km_by_sex.png\n")

# --- Plot 3: Survival by ECOG Performance Score ---
fit_ecog <- survfit(surv_obj ~ ecog_group, data = lung_clean)

ecog_km <- ggsurvplot(
  fit_ecog,
  data          = lung_clean,
  conf.int      = FALSE,
  pval          = TRUE,
  risk.table    = TRUE,
  xlab          = "Time (days)",
  ylab          = "Survival Probability",
  title         = "Survival by ECOG Performance Score",
  subtitle      = "NCCTG Lung Cancer Dataset",
  legend.title  = "ECOG Status",
  palette       = c("#2A9D8F", "#E9C46A", "#E76F51"),
  ggtheme       = theme_classic(base_size = 13),
  risk.table.height = 0.28
)

ggsave("km_by_ecog.png",
       print(ecog_km, newpage = FALSE),
       width = 9, height = 7, dpi = 300)
cat("Saved: km_by_ecog.png\n")

# --- Cox Proportional Hazards Regression ---
cox_model <- coxph(
  Surv(time, status == 2) ~ age + sex + ph.ecog + wt.loss,
  data = lung_clean
)

cat("\nCox Proportional Hazards Model Summary:\n")
print(summary(cox_model))

# --- Plot 4: Forest Plot of Cox Model ---
cox_plot <- ggforest(
  cox_model,
  data      = lung_clean,
  main      = "Cox Proportional Hazards: Hazard Ratios",
  fontsize  = 1.0
)

ggsave("cox_forest_plot.png", cox_plot, width = 9, height = 6, dpi = 300)
cat("Saved: cox_forest_plot.png\n")

# --- Save Summary Statistics ---
survival_summary <- data.frame(
  Group    = names(fit_sex$strata),
  Median   = summary(fit_sex)$table[, "median"],
  Lower_CI = summary(fit_sex)$table[, "0.95LCL"],
  Upper_CI = summary(fit_sex)$table[, "0.95UCL"]
)
write.csv(survival_summary, "survival_summary_by_sex.csv", row.names = FALSE)
cat("Saved: survival_summary_by_sex.csv\n")

cat("\nAnalysis complete. All outputs saved.\n")
