#!/usr/bin/env Rscript
# Home-tank activity during recovery from sleep restriction: statistics and figures.
#
# Input
#   results/activity_binned.csv  (written by run_pipeline.py); one row per
#   recording, tank and 10-minute bin.
#
# Output (all in results/analysis/)
#   session_values.csv        tank-session values used in every test
#   tank_means.csv            per-tank means, by time of day and overall
#   stats_by_day.csv          group comparison on each recovery day
#   stats_by_time_of_day.csv  group effect within each recording time
#   stats_overall_model.csv   mixed model over all sessions
#   ratio_by_tank_day.csv     dark/light ratio, one row per tank and day
#   stats_ratio.csv           tests on the dark/light ratio
#   power.csv                 observed effect sizes and the n they would need
#   sensitivity.csv           same tests with uncertain tank-days removed
#   figure_activity.pdf/.png  two-panel figure, panel_a.pdf, panel_b.pdf
#   session_info.txt          R and package versions of the run
#
# Usage
#   Rscript analysis.R                      # run from the folder holding results/
#   Rscript analysis.R path/to/activity_binned.csv
#   SLEEP_ROOT=/path/to/data Rscript analysis.R
#
# Design notes that the tests follow
#   The tank is the experimental unit: fish within a tank share water, feeding
#   and handling, so they are not independent. There are four tanks per
#   condition, two in each cohort, and each tank is recorded repeatedly, which
#   is why tank enters the mixed models as a random intercept.
#   Condition is confounded with camera and rack position (camera 003 recorded
#   sleep-restricted tanks, camera 004 control tanks, in both cohorts), so
#   between-condition differences in absolute activity cannot be separated from
#   position effects. The dark/light ratio is formed inside one tank on one day
#   and is free of that confound.
#   The light and dark phases use different motion thresholds, so the ratio is
#   partly method-dependent and is read qualitatively; comparisons between
#   conditions within one time of day are unaffected, both conditions being
#   processed identically.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(lme4)
  library(lmerTest)
})

# ---------------------------------------------------------------- settings ---

ROOT <- Sys.getenv("SLEEP_ROOT", unset = ".")

# The binned file can be given as the first argument; otherwise it is looked for
# under ROOT, first at the path the pipeline writes and then at the one used
# while these data were first processed.
cli_args <- commandArgs(trailingOnly = TRUE)
candidates <- c(cli_args[1],
                file.path(ROOT, "results", "activity_binned.csv"),
                file.path(ROOT, "results_v3", "activity_binned_v3.csv"))
candidates <- candidates[!is.na(candidates)]
INPUT <- candidates[file.exists(candidates)][1]
if (is.na(INPUT)) {
  stop("No input file found. Looked for:\n  ", paste(candidates, collapse = "\n  "))
}
OUTDIR <- Sys.getenv("SLEEP_OUT", unset = file.path(dirname(INPUT), "analysis"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Colours and labels used in the manuscript figures.
COND_COLOURS <- c("CTRL/7DL" = "#4DAF4A", "SR/7DL" = "#984EA3")
TIME_LABELS  <- c("06" = "06:00 (dark)", "14" = "14:00 (light)", "21" = "21:00 (lights off)")

# First recording day of each cohort; recovery day is counted from these dates.
COHORT_START <- c(b1 = "2026-04-29", b2 = "2026-08-26")

# Tank-days at or after a death date that was recorded from memory rather than
# from the daily records. Used only in the sensitivity analysis.
UNCERTAIN_DEATHS <- tribble(
  ~batch, ~camera, ~tank,   ~from_date,
  "b2",   "004",   "tank1", "2026-08-27",
  "b2",   "003",   "tank2", "2026-08-31"
)

ALPHA <- 0.05

# ----------------------------------------------------------------- helpers ---

#' Significance marks from a raw p-value.
stars <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < ALPHA, "*", ""))))
}

#' Two-sample test, chosen by a Shapiro-Wilk check on each sample.
#' Welch t-test when both samples pass, Wilcoxon rank-sum test otherwise.
two_sample_test <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) {
    return(list(p = NA_real_, test = "n too small"))
  }
  normal_enough <- function(v) {
    if (length(v) >= 3 && sd(v) > 0) shapiro.test(v)$p.value > ALPHA else TRUE
  }
  if (normal_enough(x) && normal_enough(y)) {
    list(p = t.test(x, y)$p.value, test = "Welch t")
  } else {
    list(p = suppressWarnings(wilcox.test(x, y)$p.value), test = "Wilcoxon")
  }
}

#' Cohen's d for two independent samples (pooled standard deviation).
cohens_d <- function(x, y) {
  s <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) /
              (length(x) + length(y) - 2))
  (mean(y) - mean(x)) / s
}

write_table <- function(x, name) {
  path <- file.path(OUTDIR, name)
  write.csv(x, path, row.names = FALSE)
  path
}

# -------------------------------------------------------------------- data ---

# camera, date and time are identifiers, not numbers: read them as text so that
# leading zeros survive (camera 003, date 01052026, time 060000).
raw <- read.csv(INPUT, stringsAsFactors = FALSE,
                colClasses = c(batch = "character", camera = "character",
                               date = "character", time = "character",
                               tank = "character"))

sessions <- raw %>%
  filter(tank %in% c("tank1", "tank2"), time_bin != 0) %>%   # first bin: handling
  mutate(
    date_r    = as.Date(date, format = "%d%m%Y"),
    hour      = substr(time, 1, 2),
    condition = factor(ifelse(group == "Control", "CTRL/7DL", "SR/7DL"),
                       levels = names(COND_COLOURS)),
    tank_id   = paste(batch, camera, tank, sep = "_"),
    cohort    = factor(batch)
  ) %>%
  filter(hour %in% names(TIME_LABELS)) %>%
  group_by(cohort, condition, tank_id, date_r, hour) %>%
  summarise(activity = mean(mean_motion),      # mean over the 10-minute bins
            fish_detected = mean(mean_fish),   # detector output, quality control
            n_bins = n(), .groups = "drop") %>%
  mutate(
    recovery_day = as.integer(date_r - as.Date(COHORT_START[as.character(cohort)])),
    time_of_day  = factor(TIME_LABELS[hour], levels = unname(TIME_LABELS)),
    phase        = ifelse(hour == "14", "light", "dark")
  )

stopifnot(nrow(sessions) > 0, n_distinct(sessions$tank_id) == 8)

cat(sprintf("Input: %s\nTank-sessions: %d from %d tanks (%d per condition)\n\n",
            INPUT, nrow(sessions), n_distinct(sessions$tank_id),
            n_distinct(sessions$tank_id) / 2))

write_table(sessions, "session_values.csv")

tank_means <- bind_rows(
  sessions %>%
    group_by(condition, cohort, tank_id, time_of_day) %>%
    summarise(activity = mean(activity), n_sessions = n(), .groups = "drop"),
  sessions %>%
    group_by(condition, cohort, tank_id) %>%
    summarise(time_of_day = "all", activity = mean(activity),
              n_sessions = n(), .groups = "drop")
)
write_table(tank_means, "tank_means.csv")

# ------------------------------------- group comparison on each recovery day ---
# Exploratory, and shown in the figure with raw p-values; the false discovery
# rate across days within one time of day is reported alongside.

stats_by_day <- sessions %>%
  group_by(time_of_day, recovery_day) %>%
  group_modify(function(d, key) {
    ctrl <- d$activity[d$condition == "CTRL/7DL"]
    sr   <- d$activity[d$condition == "SR/7DL"]
    res  <- two_sample_test(ctrl, sr)
    tibble(n_ctrl = length(ctrl), n_sr = length(sr),
           mean_ctrl = mean(ctrl), mean_sr = mean(sr),
           diff_pct = 100 * (mean(sr) / mean(ctrl) - 1),
           test = res$test, p_raw = res$p)
  }) %>%
  group_by(time_of_day) %>%
  mutate(p_fdr = p.adjust(p_raw, method = "fdr"), mark = stars(p_raw)) %>%
  ungroup()

write_table(stats_by_day, "stats_by_day.csv")

# --------------------------------- group effect within each time of day ------
# One mixed model per recording time, over all days at once: condition and
# recovery day as fixed effects, cohort as a fixed effect because the two
# cohorts differ in fish and in video compression, tank as a random intercept.

stats_by_time <- lapply(levels(sessions$time_of_day), function(tod) {
  d <- filter(sessions, time_of_day == tod)
  model <- lmer(activity ~ condition + recovery_day + cohort + (1 | tank_id), data = d)
  p <- anova(model)["condition", "Pr(>F)"]
  tank_level <- d %>% group_by(condition, tank_id) %>%
    summarise(a = mean(activity), .groups = "drop")
  data.frame(
    time_of_day = tod,
    mean_ctrl = mean(d$activity[d$condition == "CTRL/7DL"]),
    mean_sr   = mean(d$activity[d$condition == "SR/7DL"]),
    diff_pct  = 100 * (mean(d$activity[d$condition == "SR/7DL"]) /
                         mean(d$activity[d$condition == "CTRL/7DL"]) - 1),
    p_condition = p,
    n_tank_sessions = nrow(d),
    n_tanks = n_distinct(d$tank_id),
    d_cohen = cohens_d(tank_level$a[tank_level$condition == "CTRL/7DL"],
                       tank_level$a[tank_level$condition == "SR/7DL"])
  )
}) %>% bind_rows() %>%
  mutate(p_fdr = p.adjust(p_condition, "fdr"), mark = stars(p_condition))

write_table(stats_by_time, "stats_by_time_of_day.csv")

# ------------------------------------------------- model over all sessions ---
# Does activity change across recovery days, and differently between conditions
# (the condition by day interaction)?

overall <- lmer(activity ~ recovery_day * condition + time_of_day + cohort +
                  (1 | tank_id), data = sessions)
overall_table <- as.data.frame(anova(overall))
overall_table$term <- rownames(overall_table)
overall_table <- overall_table[, c("term", setdiff(names(overall_table), "term"))]
write_table(overall_table, "stats_overall_model.csv")

# ---------------------------------------------------- dark/light ratio -------
# Within a tank and day: dark activity (mean of the 06:00 and 21:00 recordings)
# divided by light activity (14:00). The tank's own level - camera, reflections,
# compression - cancels out. Ratios are multiplicative, so they are averaged and
# tested on the log scale and back-transformed for reporting (geometric means).

ratio_by_tank_day <- sessions %>%
  select(cohort, condition, tank_id, date_r, recovery_day, hour, activity) %>%
  pivot_wider(names_from = hour, values_from = activity,
              names_prefix = "h") %>%
  mutate(dark = rowMeans(cbind(h06, h21), na.rm = TRUE),
         light = h14,
         ratio = dark / light,
         log_ratio = log(ratio)) %>%
  filter(is.finite(log_ratio))

write_table(ratio_by_tank_day, "ratio_by_tank_day.csv")

ratio_by_tank <- ratio_by_tank_day %>%
  group_by(condition, cohort, tank_id) %>%
  summarise(ratio = exp(mean(log_ratio)), mean_log = mean(log_ratio),
            n_days = n(), .groups = "drop")

# Does each condition keep a day-night difference (ratio different from 1)?
within_condition <- ratio_by_tank_day %>%
  group_by(condition) %>%
  summarise(n_tank_days = n(),
            ratio = exp(mean(log_ratio)),
            se_low = exp(mean(log_ratio) - sd(log_ratio) / sqrt(n())),
            se_high = exp(mean(log_ratio) + sd(log_ratio) / sqrt(n())),
            p = t.test(log_ratio)$p.value, .groups = "drop") %>%
  mutate(comparison = paste(condition, "ratio vs 1 (one-sample t on log ratio)"),
         mark = stars(p))

# Between conditions, with the tank as the unit (n = 4 per condition) and, for
# comparison, a mixed model over all tank-days.
tank_ctrl <- ratio_by_tank$mean_log[ratio_by_tank$condition == "CTRL/7DL"]
tank_sr   <- ratio_by_tank$mean_log[ratio_by_tank$condition == "SR/7DL"]
between_tank <- two_sample_test(tank_ctrl, tank_sr)

ratio_model <- lmer(log_ratio ~ condition + recovery_day + cohort + (1 | tank_id),
                    data = ratio_by_tank_day)
p_model <- anova(ratio_model)["condition", "Pr(>F)"]
effect_model <- exp(fixef(ratio_model)[["conditionSR/7DL"]])

stats_ratio <- bind_rows(
  within_condition %>%
    transmute(comparison, n = n_tank_days, effect = ratio, p, mark),
  tibble(
    comparison = sprintf("SR/7DL vs CTRL/7DL, per tank (%s)", between_tank$test),
    n = nrow(ratio_by_tank),
    effect = exp(mean(tank_sr)) / exp(mean(tank_ctrl)),
    p = between_tank$p,
    mark = stars(between_tank$p)),
  tibble(
    comparison = "SR/7DL vs CTRL/7DL, mixed model over tank-days",
    n = nrow(ratio_by_tank_day),
    effect = effect_model,
    p = p_model,
    mark = stars(p_model))
)
write_table(stats_ratio, "stats_ratio.csv")

# ------------------------------------------------------------------ power ----
# What the study could detect with four tanks per condition, and the number of
# tanks per condition that 80% power would have required for the effects seen.

power_rows <- bind_rows(
  stats_by_time %>%
    transmute(measure = paste("activity at", time_of_day), d = d_cohen),
  tibble(measure = "dark/light ratio (log)",
         d = cohens_d(tank_ctrl, tank_sr))
) %>%
  rowwise() %>%
  mutate(
    n_per_group_for_80pct = ceiling(power.t.test(delta = abs(d), sd = 1,
                                                 power = 0.8)$n),
    power_with_4_tanks = power.t.test(n = 4, delta = abs(d), sd = 1)$power
  ) %>%
  ungroup()
write_table(power_rows, "power.csv")

# ------------------------------------------------------------ sensitivity ----
# One death date in cohort 2 comes from memory rather than from the records.
# Repeat the two headline tests without the tank-days at or after those dates.

excluded <- UNCERTAIN_DEATHS %>%
  mutate(tank_id = paste(batch, camera, tank, sep = "_"),
         from_date = as.Date(from_date))

sessions_sens <- sessions %>%
  anti_join(excluded, by = "tank_id") %>%
  bind_rows(
    sessions %>%
      inner_join(excluded, by = "tank_id") %>%
      filter(date_r < from_date) %>%
      select(names(sessions))
  )

ratio_sens <- ratio_by_tank_day %>%
  anti_join(excluded, by = "tank_id") %>%
  bind_rows(
    ratio_by_tank_day %>%
      inner_join(excluded, by = "tank_id") %>%
      filter(date_r < from_date) %>%
      select(names(ratio_by_tank_day))
  )

sens_overall <- lmer(activity ~ recovery_day * condition + time_of_day + cohort +
                       (1 | tank_id), data = sessions_sens)
sens_ratio <- lmer(log_ratio ~ condition + recovery_day + cohort + (1 | tank_id),
                   data = ratio_sens)

sensitivity <- tibble(
  test = c("condition x recovery day (all sessions)",
           "condition x recovery day (uncertain tank-days removed)",
           "condition effect on dark/light ratio (all tank-days)",
           "condition effect on dark/light ratio (uncertain tank-days removed)"),
  n = c(nrow(sessions), nrow(sessions_sens),
        nrow(ratio_by_tank_day), nrow(ratio_sens)),
  p = c(anova(overall)["recovery_day:condition", "Pr(>F)"],
        anova(sens_overall)["recovery_day:condition", "Pr(>F)"],
        p_model,
        anova(sens_ratio)["condition", "Pr(>F)"])
)
write_table(sensitivity, "sensitivity.csv")

# ----------------------------------------------------------------- figures ---

fig_theme <- theme_classic(base_size = 11, base_family = "sans") +
  theme(axis.line = element_line(colour = "black", linewidth = 0.4),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        axis.text = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(10, "pt"))

panel_box <- theme(plot.background = element_rect(colour = "black", fill = NA,
                                                  linewidth = 0.6),
                   plot.margin = margin(10, 10, 8, 10))

# Panel A: mean activity per recovery day, one facet per recording time.
timecourse <- sessions %>%
  group_by(condition, time_of_day, recovery_day) %>%
  summarise(mean = mean(activity), se = sd(activity) / sqrt(n()), n = n(),
            .groups = "drop")
write_table(timecourse, "timecourse_by_day.csv")

marks <- stats_by_day %>% filter(mark != "")

panel_a <- ggplot(timecourse, aes(recovery_day, mean,
                                  colour = condition, fill = condition)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.22,
              colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.9, show.legend = FALSE) +
  geom_point(size = 3, shape = 21, colour = "black", stroke = 0.5) +
  facet_wrap(~ time_of_day, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = COND_COLOURS, name = NULL) +
  scale_fill_manual(values = COND_COLOURS, name = NULL) +
  scale_x_continuous(breaks = sort(unique(timecourse$recovery_day))) +
  labs(x = "Recovery day", y = "Activity index  (a.u.)") +
  fig_theme +
  theme(legend.position = c(0.995, 0.995), legend.justification = c(1, 1)) +
  guides(colour = "none",
         fill = guide_legend(override.aes = list(shape = 21, size = 3.4,
                                                 colour = "black", stroke = 0.5)))

if (nrow(marks) > 0) {
  panel_a <- panel_a +
    geom_text(data = marks,
              aes(recovery_day, max(timecourse$mean + timecourse$se) * 1.04,
                  label = mark),
              inherit.aes = FALSE, size = 4)
}

# Panel B: dark/light ratio, bar at the geometric mean with the four tank values.
ratio_summary <- within_condition %>% select(condition, ratio, se_high, p, mark) %>%
  mutate(mark = ifelse(mark == "", "n.s.", mark))
p_between_tank <- between_tank$p

y_max  <- max(c(ratio_by_tank$ratio, ratio_summary$se_high))
y_mark <- y_max * 1.06
y_bar  <- y_max * 1.20
y_text <- y_max * 1.26

panel_b <- ggplot(ratio_summary, aes(condition, ratio, fill = condition)) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.5, alpha = 0.55) +
  geom_errorbar(aes(ymin = ratio, ymax = se_high), width = 0.14,
                linewidth = 0.5, colour = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey35",
             linewidth = 0.45) +
  geom_point(data = ratio_by_tank, aes(y = ratio),
             position = position_jitter(width = 0.13, height = 0),
             size = 2.6, shape = 21, colour = "black", stroke = 0.4) +
  geom_text(aes(y = y_mark, label = mark), size = 4.2, colour = "black") +
  annotate("segment", x = 1, xend = 2, y = y_bar, yend = y_bar, linewidth = 0.45) +
  annotate("text", x = 1.5, y = y_text, size = 3.4,
           label = sprintf("%s (p = %.2f)",
                           ifelse(p_between_tank < ALPHA, stars(p_between_tank), "n.s."),
                           p_between_tank)) +
  scale_fill_manual(values = COND_COLOURS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04)),
                     limits = c(0, y_max * 1.33)) +
  labs(x = NULL, y = "Dark / light activity ratio") +
  fig_theme

figure <- ((panel_a + panel_box) /
             ((panel_b + panel_box) + plot_spacer() + plot_layout(widths = c(1, 1.35)))) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16, family = "sans"))

save_pdf <- function(plot, file, width, height) {
  pdf(file.path(OUTDIR, file), width = width, height = height,
      useDingbats = FALSE, family = "sans")   # text stays editable in the PDF
  print(plot)
  invisible(dev.off())
}

save_pdf(figure, "figure_activity.pdf", 9.0, 7.0)
save_pdf(panel_a, "panel_a.pdf", 9.0, 3.4)
save_pdf(panel_b, "panel_b.pdf", 4.0, 3.8)
ggsave(file.path(OUTDIR, "figure_activity.png"), figure,
       width = 9.0, height = 7.0, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "session_info.txt"))

# ------------------------------------------------------------------ report ---

cat("Group comparison within each time of day\n")
print(stats_by_time %>%
        mutate(across(where(is.numeric), ~signif(.x, 3))) %>% as.data.frame(),
      row.names = FALSE)

cat("\nDark/light ratio\n")
print(stats_ratio %>%
        mutate(across(where(is.numeric), ~signif(.x, 3))) %>% as.data.frame(),
      row.names = FALSE)

cat("\nMixed model over all sessions\n")
print(signif(as.data.frame(anova(overall)), 3))

cat("\nSensitivity\n")
print(sensitivity %>% mutate(p = signif(p, 3)) %>% as.data.frame(), row.names = FALSE)

cat(sprintf("\nTables and figures written to %s\n", normalizePath(OUTDIR)))
