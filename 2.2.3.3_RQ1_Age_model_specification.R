# =========================
# run_model_test3.R
# Long-running brms model for tmux
# =========================

options(mc.cores = parallel::detectCores())

cat("=====================================\n")
cat("Start time:", as.character(Sys.time()), "\n")
cat("=====================================\n")

# ---- libraries ----
suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(ggplot2)
  library(cmdstanr)
  library(bayesplot)
  library(tidybayes)
  library(modelr)
  library(posterior)
})

cat("Packages loaded successfully.\n")
cat("CmdStan path:", cmdstanr::cmdstan_path(), "\n")

# ---- paths ----
data_path <- "/work/Thesis/Data/RL_all_annotations_updated.csv"
model_dir <- "/work/Thesis/R_scripts/model_test3"

model_file  <- file.path(model_dir, "model_test3_backup")
final_rds   <- file.path(model_dir, "model_test3_final.rds")
summary_txt <- file.path(model_dir, "model_test3_summary.txt")
newdata_rds <- file.path(model_dir, "model_test3_newdata_grid.rds")
epred_rds   <- file.path(model_dir, "model_test3_epred_grid.rds")
pred_rds    <- file.path(model_dir, "model_test3_pred_grid.rds")

# ---- create output directory if needed ----
if (!dir.exists(model_dir)) {
  dir.create(model_dir, recursive = TRUE)
}

# ---- load data ----
cat("Reading data...\n")
df <- read.csv(data_path)

# ---- filter data ----
df_analysis <- df %>%
  filter(RL_s >= -4, RL_s <= 7)

cat("Filtered data dimension:\n")
print(dim(df_analysis))

# ---- preprocessing ----
df_analysis$Year_z <- as.numeric(scale(df_analysis$Year))

df_analysis$Annotator <- factor(
  df_analysis$Annotator,
  levels = c("manual", "vtc2")
)

df_analysis$Speaker <- factor(
  df_analysis$Speaker,
  levels = c("FEM", "KCHI", "MAL")
)

df_analysis$File <- factor(df_analysis$File)
df_analysis$Child <- factor(df_analysis$Child)

cat("Post-processing data dimension:\n")
print(dim(df_analysis))

# ---- NA check ----
required_vars <- c("RL_s", "Annotator", "Speaker", "Year_z", "File", "Child")
na_check <- colSums(is.na(df_analysis[, required_vars]))
cat("NA check:\n")
print(na_check)

if (any(na_check > 0)) {
  stop("There are missing values in required columns. Please fix them before fitting.")
}

# ---- model formula ----
latency_f1 <- bf(
  RL_s ~ 0 + Annotator:Speaker + Annotator:Speaker:Year_z +
    (0 + Annotator:Speaker || File) +
    (0 + Annotator:Speaker + Annotator:Speaker:Year_z | Child)
)

cat("Model formula defined.\n")

# ---- priors ----
prior_1 <- c(
  prior(normal(0, 0.5), class = "b", coef = "Annotatormanual:SpeakerFEM"),
  prior(normal(0.3, 0.5), class = "b", coef = "Annotatormanual:SpeakerKCHI"),
  prior(normal(0, 0.5), class = "b", coef = "Annotatormanual:SpeakerMAL"),
  prior(normal(0, 0.5), class = "b", coef = "Annotatorvtc2:SpeakerFEM"),
  prior(normal(0.3, 0.5), class = "b", coef = "Annotatorvtc2:SpeakerKCHI"),
  prior(normal(0, 0.5), class = "b", coef = "Annotatorvtc2:SpeakerMAL"),
  
  prior(normal(0, 0.25), class = "b", coef = "Annotatormanual:SpeakerFEM:Year_z"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatormanual:SpeakerKCHI:Year_z"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatormanual:SpeakerMAL:Year_z"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatorvtc2:SpeakerFEM:Year_z"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatorvtc2:SpeakerKCHI:Year_z"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatorvtc2:SpeakerMAL:Year_z"),
  
  prior(student_t(3, 0, 0.2), class = "sd", lb = 0),
  prior(student_t(3, 0, 0.3), class = "sigma", lb = 0),
  prior(gamma(5, 4), class = "nu"),
  prior(lkj(2), class = "cor")
)

cat("Priors defined.\n")

# ---- fit model ----
cat("Starting model fitting...\n")
cat("Model cache file:", model_file, "\n")

model_test3 <- brm(
  formula = latency_f1,
  data = df_analysis,
  family = student(),
  prior = prior_1,
  save_pars = save_pars(all = TRUE),
  file = model_file,
  file_refit = "on_change",
  sample_prior = TRUE,
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  iter = 3000,
  warmup = 1500,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 14
  ),
  threads = threading(4)
)

cat("Model fitting completed.\n")

# ---- save full model object ----
saveRDS(model_test3, final_rds)
cat("Full model object saved to:", final_rds, "\n")

# ---- save compact summary ----
sink(summary_txt)
cat("=====================================\n")
cat("Model summary\n")
cat("Generated at:", as.character(Sys.time()), "\n")
cat("=====================================\n\n")
print(summary(model_test3))
cat("\n\n")
cat("=====================================\n")
cat("Session info\n")
cat("=====================================\n\n")
print(sessionInfo())
sink()

cat("Summary saved to:", summary_txt, "\n")

# ---- create compact prediction grid ----
# IMPORTANT:
# Do NOT run posterior_predict/posterior_epred on all observed rows.
# The dataset is too large for that.
cat("Creating compact prediction grid...\n")

year_seq <- seq(
  from = min(df_analysis$Year_z, na.rm = TRUE),
  to   = max(df_analysis$Year_z, na.rm = TRUE),
  length.out = 100
)

newdata_grid <- expand.grid(
  Annotator = factor(c("manual", "vtc2"), levels = c("manual", "vtc2")),
  Speaker   = factor(c("FEM", "KCHI", "MAL"), levels = c("FEM", "KCHI", "MAL")),
  Year_z    = year_seq
)

saveRDS(newdata_grid, newdata_rds)
cat("Prediction grid saved to:", newdata_rds, "\n")

# ---- posterior expected predictions on compact grid ----
cat("Running posterior_epred on compact grid...\n")
epred_grid <- posterior_epred(
  model_test3,
  newdata = newdata_grid,
  re_formula = NA
)
saveRDS(epred_grid, epred_rds)
rm(epred_grid)
gc()

cat("epred grid saved to:", epred_rds, "\n")

# ---- posterior predictive draws on compact grid ----
cat("Running posterior_predict on compact grid...\n")
pred_grid <- posterior_predict(
  model_test3,
  newdata = newdata_grid,
  re_formula = NA
)
saveRDS(pred_grid, pred_rds)
rm(pred_grid)
gc()

cat("pred grid saved to:", pred_rds, "\n")

cat("=====================================\n")
cat("End time:", as.character(Sys.time()), "\n")
cat("=====================================\n")

