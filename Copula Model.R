# ==============================================================================
#  CONSTRUCTION OF A NOVELTY-MATCHED EPIC COHORT USING A COPULA MODEL
# ==============================================================================
# PURPOSE: 
#   1. Runs the EPIC simulation multiple times to generate a large pool of simulated patients.
#   2. Uses a Copula to score how well EPIC agents match the correlation relationship between variables in NOVELTY.
#   3. Uses Iterative Proportional Fitting (IPF) to align EPIC to the NOVELTY marginal distributions (e.g., exact % of smokers).
#   4. Selects the best-fitting agents to create a comparison cohort.
# ==============================================================================

suppressPackageStartupMessages({
  library(epicUS)
  library(dplyr)
  library(copula)
  library(qs) 
})

# ==============================================================================
#  1. USER-CONFIGURABLE PARAMETERS
# ==============================================================================
SEED_GLOBAL             <- 136L     # Seed
K_runs                  <- 1500L    # Number of runs of the EPIC model
N_match                 <- 10000L   # The final desired sample size of the matched cohort
HORIZON_YEARS           <- 6L       # Time horizon
N_AGENTS                <- 3.5e6    # Agents per run 

## OUTPUT CONTROLS
SAVE_DIR                <- "~/External Validation Novelty/epic_histories"
SAVE_FORMAT             <- "rds"   
SAVE_INDIVIDUAL         <- TRUE     # Save one file per patient
SAVE_COMBINED           <- TRUE     # Save one large file 
SAVE_MANIFEST           <- TRUE     # Save a CSV linking IDs to their calculated weights
CLEAN_TMP_EVENTDUMPS    <- TRUE     # Delete temporary raw simulation outputs to save disk space

set.seed(SEED_GLOBAL)

# ==============================================================================
#  2. LOAD TARGET BUNDLE (NOVELTY DATA)
# ==============================================================================
# PURPOSE:
#  This section imports the relationship between patient characteristics fron NOVELTY. 
#  Matching characteristics include: Gender, Age, mMRC, prior moderate exacerbations, 
#  prior severe exacerbations, GOLD stage, and smoking status.
# ==============================================================================

target <- readRDS("novelty_COPDLLN.rds")
var_order <- target$var_order
Sigma <- target$Sigma
dimnames(Sigma) <- list(var_order, var_order)

# Continuous marginal parameters
ml_age  <- target$meanlog_ageCOPD_baseline;  sl_age  <- target$sdlog_ageCOPD_baseline

# Binary / categorical / counts NOVELTY baseline characteristics 
p_sex           <- target$p_sex_baseline
p_mmrc          <- target$p_mMRC_baseline
lev_smoking     <- as.integer(target$smoking_levels)
prob_smoking    <- as.numeric(target$smoking_probs)
lev_gold        <- as.integer(target$GOLD_levels)
prob_gold       <- as.numeric(target$GOLD_probs)
lambda_exac_mod <- target$lambda_exacerbations_moderate_baseline
lambda_exac_sev <- target$lambda_exacerbations_severeplus_baseline


# ==============================================================================
#  3. HELPER FUNCTIONS
# ==============================================================================

# --- A. Copula Math Helpers (PIT Transformations) ---
# PURPOSE:
#   Copula models cannot accept raw clinical units (e.g., "Age 65" or "GOLD 3").
#   They require all data to be transformed into Uniform Distributions [0, 1].
#   These functions perform that transformation (The "Probability Integral Transform" or PIT).
#
# EXAMPLE (GOLD STAGE):
#   - Context: GOLD 1 represents the "healthiest" 20% of the population.
#   - Input: An EPIC agent has GOLD 1.
#   - Logic: The `pit_categorical_from_probs` function assigns a random number between 0.00 and 0.20.
#   - Output: u = 0.15 (This value is now ready for the Copula).

clip01 <- function(U, eps = 1e-9) {pmin(pmax(U, eps), 1 - eps)}

pit_categorical_from_probs <- function(values, level_order, level_probs) {
  idx <- match(as.integer(values), as.integer(level_order))
  if (anyNA(idx)) stop("Found NAs or values not in level_order.")
  p <- as.numeric(level_probs); p <- p / sum(p)
  cum_all <- cumsum(p); lower_all <- c(0, head(cum_all, -1))
  lower <- lower_all[idx]; upper <- cum_all[idx]
  clip01(stats::runif(length(values), min = lower, max = upper))
}

pit_bernoulli <- function(y01, p) {
  out <- rep(NA_real_, length(y01))
  i0 <- y01 == 0; i1 <- y01 == 1
  out[i0] <- runif(sum(i0, na.rm = TRUE), 0, 1 - p)
  out[i1] <- runif(sum(i1, na.rm = TRUE), 1 - p, 1)
  clip01(out)
}

pit_poisson <- function(k, lambda) {
  if (any(lambda <= 0 | !is.finite(lambda))) stop("Lambda must be positive.")
  upper <- ppois(k, lambda = lambda)
  lower <- ppois(k - 1, lambda = lambda) 
  u <- runif(length(k), min = lower, max = upper)
  clip01(u)
}

log_pmf_cat <- function(x, levels, probs) {
  p <- probs[match(as.integer(x), as.integer(levels))]
  p <- p + 1e-9; p <- p / sum(p)
  log(p)
}

TAU_CAT <- 1.00

# --- B. Alignment & Weighting Helpers ---

# ALIGNMENT: Ensures the Simulation Matrix (U) has columns in the EXACT order 
# required by the Target Correlation Matrix (Sigma).
strict_align_R <- function(U_epic, Sigma, var_order) {
  U <- U_epic[, var_order, drop = FALSE]
  R <- cov2cor(as.matrix(Sigma[var_order, var_order, drop = FALSE]))
  R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
  list(U = U, R = R)
}

# WEIGHTING: Iterative Proportional Fitting (IPF)
# PURPOSE: 
#   Calibrates the sample weights so that the marginal distributions of the simulation 
#   exactly match the target population parameters (e.g., 48% Female, 40% Current Smokers), 
#   while preserving the correlation structure between variables established by the Copula.
#
# DESCRIPTION:
#   The following iteratively adjusts weights to satisfy one marginal constraint at a time.
#   1. Adjustment A (e.g., Sex): Weights are scaled so the sex ratio matches the target (e.g., 48% Female).
#      Side Effect: This scaling may inadvertently shift the distribution of correlated variables (e.g., Smoking)
#   2. Adjustment B (e.g., Smoking): Weights are re-scaled to match the smoking target.
#      Side Effect: This may slightly disrupt the previously established Sex ratio.
#   3. Iteration: The process cycles through all variables repeatedly. With each pass, 
#      the error (discrepancy) decreases until the weights converge to a stable solution 
#      that satisfies all marginal targets simultaneously.

ipf_weights <- function(w, factors, levels_list, targets_list, max_iter = 10000, tol = 1e-8) {
  stopifnot(length(factors) == length(levels_list))
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  w <- w / sum(w)
  
  for (iter in seq_len(max_iter)) {
    max_delta <- 0
    for (j in seq_along(factors)) {
      f <- factors[[j]]
      lev <- as.integer(levels_list[[j]])
      tgt <- as.numeric(targets_list[[j]])
      
      emp <- as.numeric(tapply(w, factor(f, levels = lev), sum))
      emp[is.na(emp)] <- 0
      emp <- emp + 1e-9; emp <- emp / sum(emp)
      
      ratio <- tgt / emp
      adj <- ratio[match(f, lev)]
      
      w <- w * adj
      w[!is.finite(w)] <- 0
      w <- w / sum(w) # Renormalize immediately
      
      max_delta <- max(max_delta, max(abs(ratio - 1)))
    }
    if (max_delta < tol) break
  }
  w
}

# Helper to split exacerbation counts into bins for IPF 
poisson_bin <- function(lambda, tail_p = 0.1, hard_cap = 15) {
  if (!is.finite(lambda) || lambda <= 0) stop("ERROR")
  K <- as.integer(ceiling(qpois(1 - tail_p, lambda)))
  K <- max(1L, min(K, hard_cap))
  levels <- 0L:K
  probs  <- dpois(0:(K - 1L), lambda = lambda)
  tail   <- 1 - sum(probs)
  probs  <- c(probs, tail)
  probs  <- pmax(probs, 1e-9); probs <- probs / sum(probs)
  list(levels = levels, probs = probs, K = K)
}

# --- C. I/O Helpers (Saving/Loading) ---
# Helper functions for saving temporary raw data (Used to retrieve full histories later)

.pkg_avail <- function(pkg) requireNamespace(pkg, quietly = TRUE)

.save_events_run <- function(df, run_idx) {
  path_qs  <- file.path(TMP_DIR, sprintf("events_run_%03d.qs", run_idx))
  path_rds <- file.path(TMP_DIR, sprintf("events_run_%03d.rds", run_idx))
  if (.pkg_avail("qs")) {
    qs::qsave(df, path_qs, preset = "high")
    return(path_qs)
  } else {
    saveRDS(df, path_rds)
    return(path_rds)
  }
}

.load_events_run <- function(run_idx) {
  path_qs  <- file.path(TMP_DIR, sprintf("events_run_%03d.qs", run_idx))
  path_rds <- file.path(TMP_DIR, sprintf("events_run_%03d.rds", run_idx))
  if (file.exists(path_qs) && .pkg_avail("qs")) return(qs::qread(path_qs))
  if (file.exists(path_rds)) return(readRDS(path_rds))
  stop(sprintf("Missing cached events for run %d.", run_idx))
}

.write_one <- function(df, filepath, format = SAVE_FORMAT) {
  format <- tolower(format)
  if (format == "rds") {
    saveRDS(df, filepath)
  } else if (format == "csv") {
    utils::write.csv(df, filepath, row.names = FALSE)
  } else if (format == "parquet") {
    if (!.pkg_avail("arrow")) stop("arrow not installed.")
    arrow::write_parquet(df, filepath)
  } else {
    stop("Unsupported SAVE_FORMAT")
  }
}

# ==============================================================================
# 4. COPULA MATCHING: JOINT LIKELIHOOD CALCULATION
# ==============================================================================
# PURPOSE: 
#   Takes raw event data from one EPIC run and calculates the "Match Score" (Likelihood) for every agent.
#
# PROCESS STEPS:
#   1. Filter: Keep only diagnosed patients within the 2016-2018 index window.
#   2. Define: Create variables (Age at Index, Smoking Status) to match NOVELTY definitions.
#   3. Transform (PIT): Convert all clinical variables into U-values [0,1] using the Helper functions.
#   4. Construct Matrix: Combine U-values into the input matrix for the Copula.
#   5. Score: Calculate 'll_joint' (Log Likelihood) indicating how well the agent matches NOVELTY.
#
# Total Score Calculation:
#   - Combines the Correlation Score (Copula) with the Individual Trait Scores (Marginals).
#   - Result: 'll_joint' (Total Log Likelihood).
#
# OUTPUT:
#  - 'll_joint': A single number for each agent. High numbers = Good Match. Low numbers = Bad Match.
#  - These scores are converted to weights (w = exp(ll_joint)) 

compute_pool_candidates <- function(all_events2) {
  
  # --- Filter & Select Candidates ---
  cand <- all_events2 %>%
    mutate(calendar_year = floor(2015 + time_at_creation + local_time)) %>%
    filter(diagnosis > 0, gold > 0, FEV1_baseline > 0) %>%
    group_by(id) %>%
    slice_min(local_time, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(index_time = calendar_year) %>%
    filter(
      index_time >= 2016 & index_time < 2018,
      local_time >= 1 
    ) %>%
    mutate(age_at_index = age_at_creation + local_time) %>%
    mutate(
      smoking_raw = as.integer(as.character(smoking_status)),
      pack_years  = as.numeric(pack_years),
      gold        = as.integer(gold),
      exac_mod    = as.integer(exac_history_n_moderate),
      exac_sev    = as.integer(exac_history_n_severe_plus),
      smoking_status = case_when(
        smoking_raw == 1L ~ 1L,
        pack_years <= 0.1 ~ 0L,
        TRUE ~ 2L
      )
    )
  
  # --- Bin for IPF ---
  bin_mod <- poisson_bin(lambda_exac_mod)
  bin_sev <- poisson_bin(lambda_exac_sev)
  cand$exac_mod_bin <- pmin(pmax(cand$exac_mod, 0L), bin_mod$K)
  cand$exac_sev_bin <- pmin(pmax(cand$exac_sev, 0L), bin_sev$K)
  
  # --- PIT Transformation (Clinical -> Uniform) ---
  u_sex  <- pit_bernoulli(cand$female, p = p_sex)
  u_age  <- clip01(plnorm(pmax(cand$age_at_index, 1e-6), meanlog = ml_age, sdlog = sl_age))
  
  u_smoking <- pit_categorical_from_probs(cand$smoking_status, level_order = lev_smoking, level_probs = prob_smoking)
  u_gold    <- pit_categorical_from_probs(cand$gold,           level_order = lev_gold,    level_probs = prob_gold)
  
  u_mmrc <- pit_bernoulli(cand$dyspnea, p = p_mmrc)
  u_exac_mod <- pit_poisson(cand$exac_mod, lambda = lambda_exac_mod)
  u_exac_sev <- pit_poisson(cand$exac_sev, lambda = lambda_exac_sev)
  
  # --- Matrix Construction ---
  U_epic <- cbind(
    sex_baseline                      = u_sex,
    ageCOPD_baseline                  = u_age,
    smoking_baseline                  = u_smoking,
    GOLDgrade_baseline                = u_gold,
    mMRC_baseline                     = u_mmrc,
    exacerbations_moderate_baseline   = u_exac_mod,
    exacerbations_severeplus_baseline = u_exac_sev
  )
  
  # --- Copula Calculation ---
  aligned <- strict_align_R(U_epic, Sigma, var_order)
  U <- aligned$U
  R <- aligned$R
  
  gcop   <- copula::normalCopula(param = copula::P2p(R), dim = ncol(R), dispstr = "un")
  ll_cop <- copula::dCopula(U, gcop, log = TRUE)
  
  # --- Total Likelihood Scoring ---
  ll_cat <-
    log_pmf_cat(cand$smoking_status, levels = lev_smoking, probs = prob_smoking) +
    log_pmf_cat(cand$gold,           levels = lev_gold,    probs = prob_gold)
  
  ll_marg <-
    dlnorm(pmax(cand$age_at_COPD, 1e-6), meanlog = ml_age,  sdlog = sl_age,  log = TRUE) +
    dbinom(cand$female, size = 1, prob = p_sex,  log = TRUE) +
    dbinom(cand$dyspnea, size = 1, prob = p_mmrc, log = TRUE) +
    dpois(cand$exac_mod, lambda = lambda_exac_mod, log = TRUE) +
    dpois(cand$exac_sev, lambda = lambda_exac_sev, log = TRUE) +
    TAU_CAT * ll_cat
  
  cand$ll_joint <- ll_cop + ll_marg
  cand
}


# ==============================================================================
#  5. EXECUTION: SIMULATION LOOP
# ==============================================================================
# PURPOSE:
#   A single run of EPIC (even with 3.5 million agents) may not produce enough simulations that matches
#   NOVELTY characteristics
#   To solve this, we run the simulation multiple times in a loop (represented by K_runs).
#   - K_runs: How many separate batches we run.
#   - compute_pool_candidates: Processes just that one batch and saves the potential matches.
#   
#   The result is a pooled list that will be merged later.

# Create Directories
dir.create(SAVE_DIR, recursive = TRUE, showWarnings = FALSE)
TMP_DIR <- file.path(SAVE_DIR, "tmp_events")
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)

all_pools <- vector("list", K_runs)

for (k in seq_len(K_runs)) {
  message(sprintf("=== EPIC run %d/%d ===", k, K_runs))
  
  settings <- get_default_settings()
  settings$record_mode   <- 2
  settings$n_base_agents <- N_AGENTS
  init_session(settings = settings)
  on.exit(try(terminate_session(), silent = TRUE), add = TRUE)
  
  input <- get_input()
  input$values$global_parameters$time_horizon <- HORIZON_YEARS
  
  run(input = input$values)
  all_events2 <- as.data.frame(Cget_all_events_matrix())
  
  invisible(.save_events_run(all_events2, k))
  
  pool_k <- compute_pool_candidates(all_events2)
  pool_k$.run_id <- k
  all_pools[[k]] <- pool_k
  
  terminate_session()
}

# ==============================================================================
#  6. EXECUTION: WEIGHTING & SELECTION
# ==============================================================================
# PURPOSE:
#   Now that the loop is finished, we have K separate pools of candidates.
#   
#   Next steps are to:
#   1. Combine: Merge them into one large dataframe (cand_all).
#   2. Weight (Copula): Calculate base weights from the Joint Likelihood scores.
#   3. Weight (IPF): Refine those weights to ensure exact marginal matches (e.g. 48% Female).
#   4. Sample: Randomly draw the final N_match patients based on these weights.

cand_all <- dplyr::bind_rows(all_pools)
if (!nrow(cand_all)) stop("No candidates available after pooling.")

cand_all <- cand_all %>% distinct(id, .run_id, .keep_all = TRUE)
message(sprintf("Aggregated pool size: %d", nrow(cand_all)))

# Base Weights
w <- exp(cand_all$ll_joint - max(cand_all$ll_joint, na.rm = TRUE))
w[!is.finite(w)] <- 0
w <- w / sum(w)

# IPF Targets
bin_mod <- poisson_bin(lambda_exac_mod)
bin_sev <- poisson_bin(lambda_exac_sev)
prob_sex_vec  <- c(1 - p_sex, p_sex)
prob_mmrc_vec <- c(1 - p_mmrc, p_mmrc)
lev_binary    <- c(0L, 1L)

w_ipf <- ipf_weights(
  w = w,
  factors = list(
    gold    = cand_all$gold,
    smoke   = cand_all$smoking_status,
    ex_mod  = cand_all$exac_mod_bin,
    ex_sev  = cand_all$exac_sev_bin,
    sex     = cand_all$female,        
    dyspnea = cand_all$dyspnea        
  ),
  levels_list = list(
    gold    = lev_gold,
    smoke   = lev_smoking,
    ex_mod  = bin_mod$levels,
    ex_sev  = bin_sev$levels,
    sex     = lev_binary,             
    dyspnea = lev_binary             
  ),
  targets_list = list(
    gold    = prob_gold,
    smoke   = prob_smoking,
    ex_mod  = bin_mod$probs,
    ex_sev  = bin_sev$probs,
    sex     = prob_sex_vec,           
    dyspnea = prob_mmrc_vec           
  ),
  max_iter = 10000,
  tol = 1e-8
)

# Selection
cand_all$w_final <- w_ipf
cand_all <- cand_all[cand_all$w_final > 0, ]
w_norm <- cand_all$w_final / sum(cand_all$w_final)

set.seed(SEED_GLOBAL + 57L)
size_final <- min(N_match, nrow(cand_all))
idx_final <- sample.int(n = nrow(cand_all), size = size_final, replace = FALSE, prob = w_norm)
epic_selected <- cand_all[idx_final, , drop = FALSE]

cat("Selected:", nrow(epic_selected), "records.\n")

# Save baseline cohort details for SMD calculation 
saveRDS(epic_selected, "epic_selected.rds")

# ==============================================================================
#  7. EXECUTION: EXPORT HISTORIES
# ==============================================================================
# PURPOSE:
#   We currently only have the "Index Date" info for the selected patients.
#   We need their FULL event history (drugs, exacerbations over time).
#   
#   Instead of keeping storing all data, we:
#   1. Look at our list of selected patients.
#   2. Load the temporary file for Run 1.
#   3. Extract ONLY the patients we need from Run 1.
#   4. Repeat for Run 2, Run 3, etc.

hist_dir <- file.path(SAVE_DIR, "per_id")
if (SAVE_INDIVIDUAL) dir.create(hist_dir, showWarnings = FALSE)

# Helper to create unique IDs (e.g., "Agent10_Run2")
.make_global_id <- function(id_vec, run_i) paste0(as.character(id_vec), "_r", sprintf("%03d", run_i))

combined_list <- list()
manifest <- data.frame(
  run_id = integer(), id = character(), global_id = character(),
  weight = numeric(), file_path = character(), stringsAsFactors = FALSE
)

# Split the selected list by Run ID so we can process one batch at a time
by_run <- split(epic_selected, epic_selected$.run_id, drop = TRUE)

# Load the large raw file for this run
for (run_key in names(by_run)) {
  run_i <- as.integer(run_key)
  sel_i <- by_run[[run_key]]
  
  # Create global IDs to match the selection
  ev_i  <- .load_events_run(run_i)
  ev_i$id  <- as.character(ev_i$id)
  sel_i$id <- as.character(sel_i$id)
  sel_i$global_id <- .make_global_id(sel_i$id, run_i)
  
  # Keep only the events for our selected patients
  sub_i <- ev_i[ev_i$id %in% sel_i$id, , drop = FALSE]
  if (!nrow(sub_i)) next
  
  sub_i$.run_id   <- run_i
  sub_i$global_id <- .make_global_id(sub_i$id, run_i)
  
  # Save Individual Files (Per Patient)
  if (SAVE_INDIVIDUAL) {
    seen <- character(0)
    for (gid in unique(sub_i$global_id)) {
      rows_id <- sub_i[sub_i$global_id == gid, , drop = FALSE]
      
      # File naming logic
      base <- paste0("history_", gid)
      ext  <- switch(tolower(SAVE_FORMAT), rds="rds", csv="csv", parquet="parquet")
      fpath <- file.path(hist_dir, paste0(base, ".", ext))
      
      # Handle potential duplicates
      kdup <- 1L
      while (fpath %in% seen || file.exists(fpath)) {
        fpath <- file.path(hist_dir, paste0(base, "_", kdup, ".", ext))
        kdup <- kdup + 1L
      }
      seen <- c(seen, fpath)
      
      # Write file
      .write_one(rows_id, fpath, SAVE_FORMAT)
      
      # Update Manifest (Tracking Sheet)
      wt <- tryCatch(sel_i$w_final[match(gid, sel_i$global_id)], error = function(e) NA_real_)
      manifest <- rbind(manifest, data.frame(
        run_id = run_i, id = unique(rows_id$id)[1], global_id = gid,
        weight = wt, file_path = fpath, stringsAsFactors = FALSE
      ))
    }
  }
  # Accumulate for the combined file
  if (SAVE_COMBINED) combined_list[[run_key]] <- sub_i
}

# Write the combined file (all patients in one large data table)
if (SAVE_COMBINED && length(combined_list)) {
  combined_df <- dplyr::bind_rows(combined_list)
  comb_path <- file.path(SAVE_DIR, paste0("combined_histories.", switch(tolower(SAVE_FORMAT), rds="rds", csv="csv", parquet="parquet")))
  .write_one(combined_df, comb_path, SAVE_FORMAT)
}

# Write the file (Key to weights and file paths)
if (SAVE_MANIFEST) {
  man_path <- file.path(SAVE_DIR, "manifest_per_id.csv")
  utils::write.csv(manifest, man_path, row.names = FALSE)
}

if (CLEAN_TMP_EVENTDUMPS) {
  unlink(TMP_DIR, recursive = TRUE, force = TRUE)
}

# ==============================================================================
#  8. EPIC MULTIPLE-RUN TEST
# ==============================================================================
# PURPOSE: 
#   Verifies that the 'combined_histories' file matches the 'epic_selected' cohort
#   It prevents ID mismatches between runs.

if (exists("combined_df") && nrow(combined_df) > 0) {
  
  message("\n============================================================")
  message(" MULTIPLE EPIC RUN TEST ")
  message("============================================================")
  
  # --- Create Composite Keys (ID + Run) ---
  # We use paste() to create a unique ID because IDs repeat across runs
  target_ids <- paste(epic_selected$id, epic_selected$.run_id, sep = "_")
  actual_ids <- unique(paste(combined_df$id, combined_df$.run_id, sep = "_"))
  
  # --- TEST: Check if all simulated patients are captured---
  missing_patients <- setdiff(target_ids, actual_ids)
  n_missing <- length(missing_patients)
  
  if (n_missing == 0) {
    message("[PASS] All ", length(target_ids), " selected patients are present in the history file.")
  } else {
    warning(sprintf("[FAIL]: Missing %d patients! (e.g., %s)", 
                    n_missing, paste(head(missing_patients, 3), collapse=", ")))
  }
  
  # --- TEST: Did we include any extra simulated patients that weren't meant to be captured? ---
  extra_patients <- setdiff(actual_ids, target_ids)
  n_extra <- length(extra_patients)
  
  if (n_extra == 0) {
    message("[PASS]: No simulated patients that weren't meant to be included were identified.")
  } else {
    warning(sprintf("[FAIL]: Found %d extra patients who should not be there! (e.g., %s)", 
                    n_extra, paste(head(extra_patients, 3), collapse=", ")))
  }
  
  # --- TEST: Check IDs (Spot Check) ---
  # Verify that the Run IDs in the history file make sense.
  
  # Sample 100 random patients to check
  sample_keys <- sample(actual_ids, min(100, length(actual_ids)))
  
  message("\n--- Spot Check (Verifying Run IDs) ---")
  pass_count <- 0
  
  for (key in sample_keys) {
    # Parse the expected Run ID from our key string "AgentID_RunID"
    parts <- strsplit(key, "_")[[1]]
    expected_id <- parts[1]
    expected_run <- as.integer(parts[2])
    
    # Check the actual data
    match_row <- combined_df %>% 
      filter(id == expected_id, .run_id == expected_run) %>% 
      slice(1)
    
    if (nrow(match_row) > 0) {
      pass_count <- pass_count + 1
    } else {
      warning(sprintf("Mismatch found for %s (Run %d)", expected_id, expected_run))
    }
  }
  
  if (pass_count == length(sample_keys)) {
    message(sprintf("[PASS] Spot check passed for %d random patients.", pass_count))
  } else {
    warning("[FAIL] Spot check failed. Data may be misaligned.")
  }
  
  message("\n============================================================")
  
} else {
  message("Skipping Multi-Run Test: 'combined_df' not found.")
}