# =============================================================================
# Recurring Operations: Group Comparison Analysis
# =============================================================================
# Methods covered:
#   1. Kaplan-Meier + log-rank (time to FIRST re-operation)
#   2. Cox Proportional Hazards (time to first re-operation)
#   3. Andersen-Gill model (all re-operations, independent increments)
#   4. Prentice-Williams-Peterson (PWP) model (ordered recurring events)
#   5. Negative Binomial regression (total count of re-operations)
#   6. Mean Cumulative Function (MCF) plot
# =============================================================================

# ---- 0. Install / load packages ----
pkgs <- c("survival", "survminer", "MASS", "dplyr", "ggplot2", "tibble")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(42)

# =============================================================================
# 1. SIMULATE DATA
# =============================================================================
# We'll generate:
#   - n patients per group
#   - Individual follow-up times (right-censoring)
#   - Recurring re-operations during follow-up
# Group A (control): higher re-operation rate, shorter time to first event
# Group B (treatment): lower re-operation rate, longer time to first event
# =============================================================================

simulate_patients <- function(n_per_group = 100) {

  sim_group <- function(n, group_label, base_rate, max_followup_mean) {

    records <- vector("list", n)

    for (i in seq_len(n)) {
      pid        <- paste0(group_label, "_", i)
      followup   <- runif(1, min = max_followup_mean * 0.5,
                              max = max_followup_mean * 1.5)

      # Simulate event times from an exponential process
      event_times <- c()
      t_current   <- 0
      repeat {
        gap     <- rexp(1, rate = base_rate)
        t_next  <- t_current + gap
        if (t_next >= followup) break
        event_times <- c(event_times, t_next)
        t_current   <- t_next
      }

      records[[i]] <- list(
        id         = pid,
        group      = group_label,
        followup   = followup,
        n_reops    = length(event_times),
        event_times = event_times
      )
    }
    records
  }

  group_a <- sim_group(n_per_group, "A", base_rate = 0.30, max_followup_mean = 12)
  group_b <- sim_group(n_per_group, "B", base_rate = 0.12, max_followup_mean = 18)
  c(group_a, group_b)
}

patients <- simulate_patients(n_per_group = 120)

# ---- Build a tidy patient-level summary ----
patient_df <- tibble(
  id       = sapply(patients, `[[`, "id"),
  group    = factor(sapply(patients, `[[`, "group")),
  followup = sapply(patients, `[[`, "followup"),
  n_reops  = sapply(patients, `[[`, "n_reops")
)

cat("=== Patient summary ===\n")
patient_df %>%
  group_by(group) %>%
  summarise(
    n            = n(),
    median_fup   = round(median(followup), 1),
    mean_reops   = round(mean(n_reops), 2),
    pct_any_reop = round(mean(n_reops > 0) * 100, 1)
  ) %>%
  print()

# =============================================================================
# 2. DATA FORMATS
# =============================================================================

# ---- 2a. Wide format: one row per patient, time to FIRST re-op ----
first_reop_df <- patient_df %>%
  mutate(
    time_first = mapply(function(rec, fup) {
      et <- rec$event_times
      if (length(et) == 0) fup else min(et)
    }, patients, followup),
    event_first = as.integer(n_reops > 0)
  )

# ---- 2b. Long (counting process) format: one row per event interval ----
# Format:  id | group | tstart | tstop | status | event_num
# Used for: Andersen-Gill, PWP, MCF
make_long_format <- function(patients, patient_df) {

  rows <- vector("list", nrow(patient_df))

  for (i in seq_along(patients)) {
    rec  <- patients[[i]]
    fup  <- rec$followup
    et   <- sort(rec$event_times)
    grp  <- rec$group
    pid  <- rec$id

    if (length(et) == 0) {
      # No events: single interval [0, followup), status = 0
      rows[[i]] <- tibble(
        id        = pid,
        group     = grp,
        tstart    = 0,
        tstop     = fup,
        status    = 0L,
        event_num = 1L
      )
    } else {
      times  <- c(0, et, fup)
      n_int  <- length(times) - 1
      rows[[i]] <- tibble(
        id        = pid,
        group     = grp,
        tstart    = times[1:n_int],
        tstop     = times[2:(n_int + 1)],
        status    = c(rep(1L, length(et)), 0L),
        event_num = 1:n_int
      )
    }
  }

  bind_rows(rows) %>%
    mutate(group = factor(group))
}

long_df <- make_long_format(patients, patient_df)

cat("\n=== Long format (first 10 rows) ===\n")
print(head(long_df, 10))

# =============================================================================
# 3. KAPLAN-MEIER + LOG-RANK  (time to FIRST re-operation)
# =============================================================================
cat("\n\n========== 1. KAPLAN-MEIER + LOG-RANK ==========\n")

km_fit <- survfit(Surv(time_first, event_first) ~ group, data = first_reop_df)
print(km_fit)

# Log-rank test
logrank <- survdiff(Surv(time_first, event_first) ~ group, data = first_reop_df)
cat("\nLog-rank test:\n")
print(logrank)

# P-value
p_logrank <- pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
cat(sprintf("Log-rank p-value: %.4f\n", p_logrank))

# Plot
km_plot <- ggsurvplot(
  km_fit,
  data        = first_reop_df,
  pval        = TRUE,
  conf.int    = TRUE,
  risk.table  = TRUE,
  palette     = c("#E74C3C", "#2ECC71"),
  title       = "Time to First Re-operation (Kaplan-Meier)",
  xlab        = "Time (months)",
  ylab        = "Re-operation-free Probability",
  legend.labs = c("Group A (control)", "Group B (treatment)")
)
print(km_plot)

# =============================================================================
# 4. COX PROPORTIONAL HAZARDS  (time to FIRST re-operation)
# =============================================================================
cat("\n\n========== 2. COX PROPORTIONAL HAZARDS ==========\n")

cox_fit <- coxph(Surv(time_first, event_first) ~ group, data = first_reop_df)
cat("\nCox PH model summary:\n")
print(summary(cox_fit))

# Test PH assumption (Schoenfeld residuals)
cat("\nTest of PH assumption (Schoenfeld residuals):\n")
ph_test <- cox.zph(cox_fit)
print(ph_test)
# p > 0.05 → PH assumption not violated

# =============================================================================
# 5. ANDERSEN-GILL MODEL  (all re-operations, independent increments)
# =============================================================================
cat("\n\n========== 3. ANDERSEN-GILL MODEL ==========\n")
# Assumes each event interval is independent (like multiple independent subjects)
# Uses (tstart, tstop] counting process notation

ag_fit <- coxph(
  Surv(tstart, tstop, status) ~ group + cluster(id),
  data   = long_df,
  method = "efron"
)
cat("\nAndersen-Gill model summary:\n")
print(summary(ag_fit))
cat("  HR > 1 means Group B has HIGHER re-op rate than Group A\n")
cat("  HR < 1 means Group B has LOWER  re-op rate (protective)\n")

# =============================================================================
# 6. PRENTICE-WILLIAMS-PETERSON (PWP) MODEL  (ordered recurring events)
# =============================================================================
cat("\n\n========== 4. PWP MODEL (ordered recurring events) ==========\n")
# Conditions on prior event — models k-th re-operation separately.
# stratify by event number so each stratum gets its own baseline hazard.

# Limit to events 1-4 (small strata beyond that)
pwp_df <- long_df %>%
  filter(event_num <= 4) %>%
  mutate(strata = factor(event_num))

pwp_fit <- coxph(
  Surv(tstart, tstop, status) ~ group + strata(strata) + cluster(id),
  data   = pwp_df,
  method = "efron"
)
cat("\nPWP model summary:\n")
print(summary(pwp_fit))

# Per-stratum (per event number) analysis
cat("\nPer-event-number HRs:\n")
for (k in 1:4) {
  df_k  <- long_df %>% filter(event_num == k)
  if (sum(df_k$status) < 5) next   # skip if too few events
  fit_k <- coxph(Surv(tstart, tstop, status) ~ group + cluster(id), data = df_k)
  s     <- summary(fit_k)
  hr    <- round(s$conf.int[1, 1], 3)
  lo    <- round(s$conf.int[1, 3], 3)
  hi    <- round(s$conf.int[1, 4], 3)
  pv    <- round(s$coefficients[1, 5], 4)
  cat(sprintf("  Event %d: HR = %.3f (95%% CI %.3f–%.3f), p = %.4f\n",
              k, hr, lo, hi, pv))
}

# =============================================================================
# 7. NEGATIVE BINOMIAL REGRESSION  (total re-operation count)
# =============================================================================
cat("\n\n========== 5. NEGATIVE BINOMIAL REGRESSION ==========\n")
# offset(log(followup)) accounts for different observation lengths

nb_fit <- glm.nb(
  n_reops ~ group + offset(log(followup)),
  data = patient_df
)
cat("\nNegative Binomial model summary:\n")
print(summary(nb_fit))

# Rate ratio with CIs
rr     <- exp(coef(nb_fit))
rr_ci  <- exp(confint(nb_fit))
cat("\nRate Ratios (relative to Group A):\n")
print(cbind(RR = round(rr, 3), round(rr_ci, 3)))

cat(sprintf(
  "\nGroup B has %.1f%% %s re-operations per unit time compared to Group A\n",
  abs(1 - rr["groupB"]) * 100,
  ifelse(rr["groupB"] < 1, "fewer", "more")
))

# =============================================================================
# 8. MEAN CUMULATIVE FUNCTION (MCF)
# =============================================================================
cat("\n\n========== 6. MEAN CUMULATIVE FUNCTION (MCF) ==========\n")
# Uses the Nelson-Aalen estimator applied to recurring events.
# survfit() with id= argument computes the MCF.

mcf_fit <- survfit(
  Surv(tstart, tstop, status) ~ group,
  data = long_df,
  id   = id       # identifies subject; cumhaz gives the MCF for recurring events
)

# Extract MCF data for plotting
mcf_data <- bind_rows(
  tibble(
    time  = mcf_fit[1]$time,
    mcf   = mcf_fit[1]$cumhaz,
    lower = mcf_fit[1]$cumhaz - 1.96 * mcf_fit[1]$std.err,
    upper = mcf_fit[1]$cumhaz + 1.96 * mcf_fit[1]$std.err,
    group = "A"
  ),
  tibble(
    time  = mcf_fit[2]$time,
    mcf   = mcf_fit[2]$cumhaz,
    lower = mcf_fit[2]$cumhaz - 1.96 * mcf_fit[2]$std.err,
    upper = mcf_fit[2]$cumhaz + 1.96 * mcf_fit[2]$std.err,
    group = "B"
  )
)

mcf_plot <- ggplot(mcf_data, aes(x = time, y = mcf, color = group, fill = group)) +
  geom_step(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  scale_color_manual(
    values = c("A" = "#E74C3C", "B" = "#2ECC71"),
    labels = c("A" = "Group A (control)", "B" = "Group B (treatment)")
  ) +
  scale_fill_manual(
    values = c("A" = "#E74C3C", "B" = "#2ECC71"),
    labels = c("A" = "Group A (control)", "B" = "Group B (treatment)")
  ) +
  labs(
    title    = "Mean Cumulative Function (MCF) of Re-operations",
    subtitle = "Expected cumulative number of re-operations per patient over time",
    x        = "Time (months)",
    y        = "Mean Cumulative Number of Re-operations",
    color    = "Group",
    fill     = "Group"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

print(mcf_plot)

# =============================================================================
# 9. FOREST PLOT — all effect estimates in one figure
# =============================================================================
# All models report a ratio (HR or RR) of Group B relative to Group A.
#   ratio < 1  →  Group B has FEWER / LATER re-operations  (better)
#   ratio > 1  →  Group B has MORE  / EARLIER re-operations (worse)

# Collect all estimates
extract_hr <- function(fit, coef_row = 1) {
  s  <- summary(fit)
  ci <- exp(confint(fit))
  if ("conf.int" %in% names(s)) {
    hr <- s$conf.int[coef_row, 1]
    lo <- s$conf.int[coef_row, 3]
    hi <- s$conf.int[coef_row, 4]
    pv <- s$coefficients[coef_row, ncol(s$coefficients)]
  } else {
    hr <- exp(s$coefficients[coef_row, 1])
    lo <- ci[coef_row, 1]
    hi <- ci[coef_row, 2]
    pv <- s$coefficients[coef_row, 4]
  }
  list(est = hr, lo = lo, hi = hi, p = pv)
}

cox_e <- extract_hr(cox_fit)
ag_e  <- extract_hr(ag_fit)
pwp_e <- extract_hr(pwp_fit)

# Negative binomial rate ratio
nb_s  <- summary(nb_fit)
nb_ci <- exp(confint(nb_fit))
nb_e  <- list(
  est = exp(nb_s$coefficients["groupB", 1]),
  lo  = nb_ci["groupB", 1],
  hi  = nb_ci["groupB", 2],
  p   = nb_s$coefficients["groupB", 4]
)

forest_df <- tibble(
  Model  = factor(
    c("Cox PH\n(time to 1st re-op)",
      "Andersen-Gill\n(all re-ops)",
      "PWP\n(ordered re-ops)",
      "Neg. Binomial\n(total count)"),
    levels = rev(c(
      "Cox PH\n(time to 1st re-op)",
      "Andersen-Gill\n(all re-ops)",
      "PWP\n(ordered re-ops)",
      "Neg. Binomial\n(total count)"))
  ),
  Estimate = c(cox_e$est, ag_e$est, pwp_e$est, nb_e$est),
  Lower    = c(cox_e$lo,  ag_e$lo,  pwp_e$lo,  nb_e$lo),
  Upper    = c(cox_e$hi,  ag_e$hi,  pwp_e$hi,  nb_e$hi),
  P        = c(cox_e$p,   ag_e$p,   pwp_e$p,   nb_e$p),
  Type     = c("Time-to-event", "Time-to-event", "Time-to-event", "Count")
)

forest_df <- forest_df %>%
  mutate(
    label  = sprintf("%.2f (%.2f–%.2f)\np=%s",
                     Estimate, Lower, Upper,
                     ifelse(P < 0.001, "<0.001", sprintf("%.3f", P))),
    sig    = ifelse(P < 0.05, "Significant", "Non-significant"),
    colour = ifelse(Estimate < 1, "#2ECC71", "#E74C3C")
  )

forest_plot <- ggplot(forest_df,
    aes(x = Estimate, y = Model, xmin = Lower, xmax = Upper)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(colour = colour), height = 0.2, linewidth = 0.8) +
  geom_point(aes(colour = colour, shape = sig), size = 3.5) +
  geom_text(aes(x = max(Upper) * 1.05, label = label),
            hjust = 0, size = 3, lineheight = 0.9) +
  scale_colour_identity() +
  scale_shape_manual(values = c("Significant" = 16, "Non-significant" = 1)) +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.45)),
    trans  = "log10"
  ) +
  labs(
    title    = "Group B vs Group A — Effect Estimates (Forest Plot)",
    subtitle = "Ratio < 1 (green): Group B has fewer/later re-operations (BETTER)\nRatio > 1 (red): Group B has more/earlier re-operations (WORSE)",
    x        = "Hazard Ratio / Rate Ratio (log scale)",
    y        = NULL,
    shape    = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

print(forest_plot)

# =============================================================================
# 10. PLAIN-LANGUAGE COMPARISON SUMMARY
# =============================================================================
cat("\n\n========== GROUP COMPARISON SUMMARY ==========\n")

# Helper: interpret ratio with direction
interpret <- function(label, est, lo, hi, p, unit = "re-ops") {
  direction <- if (est < 1) "FEWER" else "MORE"
  colour    <- if (est < 1) "BETTER" else "WORSE"
  sig_str   <- if (p < 0.05) sprintf("(p=%.4f, statistically significant)", p)
               else           sprintf("(p=%.4f, NOT significant)", p)
  pct       <- abs(1 - est) * 100
  cat(sprintf(
    "  %-35s  HR/RR = %.2f (95%% CI %.2f-%.2f)\n    → Group B has %.1f%% %s %s than Group A  [%s]  %s\n\n",
    label, est, lo, hi, pct, direction, unit, colour, sig_str
  ))
}

# KM medians
km_med <- surv_median(km_fit)
cat("--- Time to First Re-operation (Kaplan-Meier medians) ---\n")
for (i in seq_len(nrow(km_med))) {
  cat(sprintf("  %s:  median re-op-free time = %.1f months  (95%% CI %.1f-%.1f)\n",
              km_med$strata[i],
              km_med$median[i],
              km_med$lower[i],
              km_med$upper[i]))
}
p_logrank <- pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
cat(sprintf("  Log-rank p-value: %.4f  %s\n\n",
            p_logrank,
            ifelse(p_logrank < 0.05,
                   "→ Groups differ significantly in time to first re-operation.",
                   "→ No significant difference in time to first re-operation.")))

cat("--- Effect Estimates (Group B relative to Group A) ---\n")
interpret("Cox PH (time to 1st re-op)",  cox_e$est, cox_e$lo, cox_e$hi, cox_e$p)
interpret("Andersen-Gill (all re-ops)",  ag_e$est,  ag_e$lo,  ag_e$hi,  ag_e$p)
interpret("PWP (ordered re-ops)",        pwp_e$est, pwp_e$lo, pwp_e$hi, pwp_e$p)
interpret("Neg. Binomial (total count)", nb_e$est,  nb_e$lo,  nb_e$hi,  nb_e$p,
          unit = "re-ops per unit time")

# Overall verdict
all_ests <- c(cox_e$est, ag_e$est, pwp_e$est, nb_e$est)
any_sig  <- any(c(cox_e$p, ag_e$p, pwp_e$p, nb_e$p) < 0.05)
all_same_dir <- all(all_ests < 1) || all(all_ests > 1)

cat("--- OVERALL VERDICT ---\n")
if (all_same_dir && any_sig) {
  better <- if (all(all_ests < 1)) "B" else "A"
  worse  <- if (better == "B") "A" else "B"
  cat(sprintf(
    "  All models consistently show Group %s has BETTER outcomes\n  (fewer re-operations and/or longer time to re-operation)\n  compared to Group %s, and at least one result is statistically\n  significant.\n",
    better, worse
  ))
} else if (all_same_dir && !any_sig) {
  better <- if (all(all_ests < 1)) "B" else "A"
  cat(sprintf(
    "  All models point in the same direction (Group %s appears better)\n  but no result reaches statistical significance.\n  Consider increasing sample size.\n", better
  ))
} else {
  cat("  Results are inconsistent across models.\n")
  cat("  Check model assumptions and data quality.\n")
}

cat("\nDone.\n")
