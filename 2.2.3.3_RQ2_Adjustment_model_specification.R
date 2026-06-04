library(dplyr)
library(brms)
library(posterior)
library(cmdstanr)
library(tidyr)

# cmdstanr 

install.packages(
  "cmdstanr",
  repos = c("https://mc-stan.org/r-packages/", "https://cloud.r-project.org")
)

library(cmdstanr)
check_cmdstan_toolchain() 
install_cmdstan(overwrite = TRUE) 


data_path <- "RL_all_annotations_updated.csv"
model_dir <- "Pathfinder"
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)


# ---- read data ----
df <- read.csv(data_path)

# ---- filter data ----
df_analysis <- df %>%
  filter(RL_s >= -4, RL_s <= 7)

cat("Filtered data dimension:\n")
print(dim(df_analysis))

# ---- preprocessing ----
df_analysis <- df_analysis %>%
  mutate(
    Year_z = as.numeric(scale(Year)),
    Annotator = factor(Annotator, levels = c("manual", "vtc2")),
    Speaker = factor(Speaker, levels = c("FEM", "KCHI", "MAL")),
    File = factor(File),
    Child = factor(Child)
  )

cat("Post-processing data dimension:\n")
print(dim(df_analysis))

# ---- NA check ----
required_vars <- c(
  "RL_s", "Annotator", "Speaker", "Year_z",
  "File", "Child", "PrevInterLatency", "PrevSelfLatency"
)

na_check <- colSums(is.na(df_analysis[, required_vars]))
cat("NA check:\n")
print(na_check)


df_analysis <- df_analysis %>%
  filter(
    !is.na(RL_s),
    !is.na(Annotator),
    !is.na(Speaker),
    !is.na(Year_z),
    !is.na(File),
    !is.na(Child),
    !is.na(PrevInterLatency),
    !is.na(PrevSelfLatency)
  )

cat("NA excluded. Final data dimension:\n")
print(dim(df_analysis))


# ---- save cleaned data ----
saveRDS(
  df_analysis,
  file.path(model_dir, "pathfinder4_df_analysis.rds")
)

write.csv(
  df_analysis,
  file.path(model_dir, "pathfinder4_df_analysis.csv"),
  row.names = FALSE
)


# ---- model formula ----
latency_f2 <- bf(
  RL_s ~ 0 +
    Annotator:Speaker +
    Annotator:Speaker:Year_z +
    Annotator:Speaker:PrevInterLatency +
    Annotator:Speaker:PrevSelfLatency +
    Annotator:Speaker:PrevInterLatency:Year_z +
    Annotator:Speaker:PrevSelfLatency:Year_z +
    (0 + Annotator:Speaker || File) +
    (0 + Annotator:Speaker +
       Annotator:Speaker:Year_z +
       Annotator:Speaker:PrevInterLatency +
       Annotator:Speaker:PrevSelfLatency +
       Annotator:Speaker:PrevInterLatency:Year_z +
       Annotator:Speaker:PrevSelfLatency:Year_z || Child)
)


# ---- priors ----
prior_1 <- c(
  # global fallback prior for fixed effects
  prior(normal(0, 0.3), class = "b"),
  
  # baseline cell means
  prior(normal(0, 0.5), class = "b", coef = "Annotatormanual:SpeakerFEM"),
  prior(normal(0.3, 0.5), class = "b", coef = "Annotatormanual:SpeakerKCHI"),
  prior(normal(0, 0.5), class = "b", coef = "Annotatormanual:SpeakerMAL"),
  prior(normal(0, 0.5), class = "b", coef = "Annotatorvtc2:SpeakerFEM"),
  prior(normal(0.3, 0.5), class = "b", coef = "Annotatorvtc2:SpeakerKCHI"),
  prior(normal(0, 0.5), class = "b", coef = "Annotatorvtc2:SpeakerMAL"),
  
  # Year effects
  prior(normal(0, 0.25), class = "b", coef = "Annotatormanual:SpeakerFEM:Year_z"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatormanual:SpeakerKCHI:Year_z"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatormanual:SpeakerMAL:Year_z"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatorvtc2:SpeakerFEM:Year_z"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatorvtc2:SpeakerKCHI:Year_z"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatorvtc2:SpeakerMAL:Year_z"),
  
  # PrevInterLatency effects
  prior(normal(0, 0.25), class = "b", coef = "Annotatormanual:SpeakerFEM:PrevInterLatency"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatormanual:SpeakerKCHI:PrevInterLatency"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatormanual:SpeakerMAL:PrevInterLatency"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatorvtc2:SpeakerFEM:PrevInterLatency"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatorvtc2:SpeakerKCHI:PrevInterLatency"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatorvtc2:SpeakerMAL:PrevInterLatency"),
  
  # PrevSelfLatency effects
  prior(normal(0, 0.25), class = "b", coef = "Annotatormanual:SpeakerFEM:PrevSelfLatency"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatormanual:SpeakerKCHI:PrevSelfLatency"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatormanual:SpeakerMAL:PrevSelfLatency"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatorvtc2:SpeakerFEM:PrevSelfLatency"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatorvtc2:SpeakerKCHI:PrevSelfLatency"),
  prior(normal(0, 0.25), class = "b", coef = "Annotatorvtc2:SpeakerMAL:PrevSelfLatency"),
  
  # Year x PrevInterLatency
  prior(normal(0, 0.1), class = "b", coef = "Annotatormanual:SpeakerFEM:Year_z:PrevInterLatency"),
  prior(normal(0, 0.1), class = "b", coef = "Annotatormanual:SpeakerKCHI:Year_z:PrevInterLatency"),
  prior(normal(0, 0.1), class = "b", coef = "Annotatormanual:SpeakerMAL:Year_z:PrevInterLatency"),
  prior(normal(0, 0.1), class = "b", coef = "Annotatorvtc2:SpeakerFEM:Year_z:PrevInterLatency"),
  prior(normal(0, 0.1), class = "b", coef = "Annotatorvtc2:SpeakerKCHI:Year_z:PrevInterLatency"),
  prior(normal(0, 0.1), class = "b", coef = "Annotatorvtc2:SpeakerMAL:Year_z:PrevInterLatency"),
  
  # Year x PrevSelfLatency
  prior(normal(0, 0.1), class = "b", coef = "Annotatormanual:SpeakerFEM:Year_z:PrevSelfLatency"),
  prior(normal(0, 0.1), class = "b", coef = "Annotatormanual:SpeakerKCHI:Year_z:PrevSelfLatency"),
  prior(normal(0, 0.1), class = "b", coef = "Annotatormanual:SpeakerMAL:Year_z:PrevSelfLatency"),
  prior(normal(0, 0.1), class = "b", coef = "Annotatorvtc2:SpeakerFEM:Year_z:PrevSelfLatency"),
  prior(normal(0, 0.1), class = "b", coef = "Annotatorvtc2:SpeakerKCHI:Year_z:PrevSelfLatency"),
  prior(normal(0, 0.1), class = "b", coef = "Annotatorvtc2:SpeakerMAL:Year_z:PrevSelfLatency"),
  
  # hierarchical priors
  prior(student_t(3, 0, 0.2), class = "sd", lb = 0),
  prior(student_t(3, 0, 0.3), class = "sigma", lb = 0),
  prior(gamma(5, 4), class = "nu")
  # prior(lkj(2), class = "cor")
)


# ---- save formula and priors ----
saveRDS(
  latency_f2,
  file.path(model_dir, "pathfinder4_formula.rds")
)

saveRDS(
  prior_1,
  file.path(model_dir, "pathfinder4_priors.rds")
)

# ---- check prior names ----
prior_check <- get_prior(
  formula = latency_f2,
  data = df_analysis,
  family = student()
)

write.csv(
  prior_check,
  file.path(model_dir, "pathfinder4_prior_check.csv"),
  row.names = FALSE
)


# ---- run pathfinder with cmdstanr ----

stan_file <- file.path(model_dir, "pathfinder4.stan")
stan_data_file <- file.path(model_dir, "pathfinder4_standata.json")

# 1. generate Stan code
stan_code <- make_stancode(
  formula = latency_f2,
  data = df_analysis,
  family = student(),
  prior = prior_1
)

writeLines(stan_code, stan_file)

# 2. generate standata
standata4 <- make_standata(
  formula = latency_f2,
  data = df_analysis,
  family = student(),
  prior = prior_1
)

saveRDS(
  standata4,
  file.path(model_dir, "pathfinder4_standata.rds")
)

cmdstanr::write_stan_json(
  standata4,
  stan_data_file
)

# 3. compile model
mod <- cmdstan_model(stan_file)

# 4. run pathfinder
pf <- mod$pathfinder(
  data = stan_data_file,
  seed = 2026,
  output_dir = model_dir,
  output_basename = "pathfinder4",
  
  num_paths = 12,
  single_path_draws = 3000,
  
  psis_resample = FALSE,
  calculate_lp = TRUE,
  
  max_lbfgs_iters = 1000,
  num_elbo_draws = 100,
  history_size = 50,
  
  save_single_paths = TRUE,
  save_cmdstan_config = TRUE,
  refresh = 100
)

saveRDS(
  pf,
  file.path(model_dir, "pathfinder4_cmdstanr_pathfinder.rds")
)

pf$save_output_files(dir = model_dir)


class(pf)
pf

# ============================================================
# Extract and check raw draws
# ============================================================

draws_raw <- pf$draws(format = "df")

saveRDS(
  draws_raw,
  file.path(model_dir, "pathfinder4_draws_raw.rds")
)

write.csv(
  draws_raw,
  file.path(model_dir, "pathfinder4_draws_raw.csv"),
  row.names = FALSE
)

b_raw_cols <- grep("^b\\[", names(draws_raw), value = TRUE)

b_sd_check_raw <- draws_raw |>
  dplyr::select(dplyr::all_of(b_raw_cols)) |>
  dplyr::summarise(dplyr::across(everything(), sd, na.rm = TRUE))

cat("Raw b-parameter SD check:\n")
print(t(b_sd_check_raw))

if (all(as.numeric(b_sd_check_raw[1, ]) == 0)) {
  stop("All raw fixed-effect Pathfinder draws have sd = 0. Do not continue.")
}

pf$return_codes()


draws_raw <- pf$draws(format = "df")
dim(draws_raw)

pf$output_files()
cat(tail(strsplit(pf$output(), "\n")[[1]], 80), sep = "\n")



# quality check
# pf$summary()
all(as.matrix(pf$return_codes()) == 0)

# ============================================================
# Rename b[1], b[2], ... into readable brms names
# ============================================================

b_names <- colnames(standata4$X)

b_map <- data.frame(
  raw_name = paste0("b[", seq_along(b_names), "]"),
  full_name = paste0("b_", b_names)
)

write.csv(
  b_map,
  file.path(model_dir, "pathfinder4_b_name_map.csv"),
  row.names = FALSE
)

draws_named <- draws_raw

for (i in seq_len(nrow(b_map))) {
  old <- b_map$raw_name[i]
  new <- b_map$full_name[i]
  
  if (old %in% names(draws_named)) {
    names(draws_named)[names(draws_named) == old] <- new
  }
}

saveRDS(
  draws_named,
  file.path(model_dir, "pathfinder4_draws_named.rds")
)

write.csv(
  draws_named,
  file.path(model_dir, "pathfinder4_draws_named.csv"),
  row.names = FALSE
)

# ---- check renamed draws ----

b_cols <- grep("^b_", names(draws_named), value = TRUE)

b_sd_check_named <- draws_named |>
  dplyr::select(dplyr::all_of(b_cols)) |>
  dplyr::summarise(dplyr::across(everything(), sd, na.rm = TRUE))

cat("Renamed b-parameter SD check:\n")
print(t(b_sd_check_named))

if (all(as.numeric(b_sd_check_named[1, ]) == 0)) {
  stop("All renamed fixed-effect Pathfinder draws have sd = 0. Do not continue.")
}

# ============================================================
# Fixed-effect summary
# ============================================================

b_summary_named <- draws_named |>
  dplyr::select(dplyr::all_of(b_cols)) |>
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "term",
    values_to = "value"
  ) |>
  dplyr::group_by(term) |>
  dplyr::summarise(
    median = median(value),
    lower_ci_95 = quantile(value, 0.025),
    upper_ci_95 = quantile(value, 0.975),
    mean = mean(value),
    sd = sd(value),
    .groups = "drop"
  )

saveRDS(
  b_summary_named,
  file.path(model_dir, "pathfinder4_b_summary_named.rds")
)

write.csv(
  b_summary_named,
  file.path(model_dir, "pathfinder4_b_summary_named.csv"),
  row.names = FALSE
)

cat("Named Pathfinder draws and fixed-effect summaries saved.\n")


