# ---------------------------------------------------------------------------
# Fig 2  (CTRL vs SR)  and  Fig 5  (CTRL/7DL vs SR/7DL)
# NTT 6-min session totals: kinematics (row A) + ethogram (row B)
#
#   run from a clean session:  source("01_Fig2_Fig5_NTT_totals.R")
#   outputs ->  Fig2.png/.pdf, Fig5.png/.pdf, Stats_Fig2.csv, Stats_Fig5.csv
# ---------------------------------------------------------------------------

## ---- settings -------------------------------------------------------------

BASE_DIR <- "C:/Users/Barba/Desktop/Sleep_Deprivation_Study"
OUT_DIR  <- file.path(BASE_DIR, "Figures")

ID_FILE  <- "id.xlsx"
NTT_FILE <- "NTT_6min_with_behaviors.csv"

# group colours used throughout the paper
COL <- c(A = "#1F77B4", B = "#F28E2B", C = "#2E9E4F", D = "#8E5FA8")

# panel titles for the two figures
LAB <- c(A = "CTRL", B = "SR", C = "CTRL/7DL", D = "SR/7DL")

BAR_TINT <- 0.55   # 0 = full colour bar, 1 = white
ALPHA   <- 0.05
DOT_SZ  <- 1.5
BAR_W   <- 0.6
FIG_W   <- 13.5    # inches
FIG_H   <- 7.0

setwd(BASE_DIR)
dir.create(OUT_DIR, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# row letters: drawn as bold titles instead of patchwork tags, because the
# `&` operator for tag styling is not available in every patchwork/ggplot2
# combination. blank titles keep the panels in a row the same height.
# bars are drawn in a pale tint of the group colour with a dark outline, and
# the individual points carry the full-strength colour - matching the figure
# style used through the paper.
lighten <- function(col, amount = BAR_TINT) {
  v <- grDevices::col2rgb(col)
  grDevices::rgb(t(255 - (255 - v) * (1 - amount)), maxColorValue = 255)
}

row_tag <- function(p, lab) {
  p + ggtitle(lab) +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0,
                                    margin = margin(b = 2)))
}

## ---- column lookup --------------------------------------------------------
# tracking export has changed names between versions, so resolve by pattern
# instead of hard-coding. first match wins; stops loudly if nothing matches.

grab <- function(d, patterns, what) {
  for (p in patterns) {
    hit <- grep(p, names(d), ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  stop("could not find a column for '", what, "' in the NTT file.\n",
       "  columns present: ", paste(names(d), collapse = ", "))
}

## ---- data -----------------------------------------------------------------

ids <- read_excel(ID_FILE)
names(ids) <- tolower(names(ids))
names(ids)[names(ids) == "id"]    <- "ID"
names(ids)[names(ids) == "group"] <- "Group"
names(ids)[names(ids) == "sex"]   <- "Sex"
ids$ID    <- as.character(ids$ID)
ids$Group <- toupper(trimws(ids$Group))

ntt <- read.csv(NTT_FILE, stringsAsFactors = FALSE, check.names = TRUE)
names(ntt)[grep("^(ID|Fish|Animal|Subject)", names(ntt), ignore.case = TRUE)[1]] <- "ID"
ntt$ID <- as.character(ntt$ID)

cn <- list(
  dist    = grab(ntt, c("^Distance_m$", "^Distance_travelled", "^Total_distance", "Distance.*_m$"), "distance travelled"),
  speed   = grab(ntt, c("^Mean_speed", "Average_speed"), "mean speed"),
  meander = grab(ntt, c("^Meander", "Meander_deg"), "meander"),
  bottom  = grab(ntt, c("Mean_dist.*bottom", "Distance_from_bottom", "Mean_height"), "mean distance from bottom"),
  nb      = grab(ntt, c("^NB_time_s$", "^NB_time", "Normal.*time"), "NB time"),
  im      = grab(ntt, c("^IM_time_s$", "^IM_time", "Immobil.*_time_s$"), "IM time"),
  hs      = grab(ntt, c("^HS_time_s$", "^HS_time", "Hyperactiv.*time"), "HS time"),
  tr      = grab(ntt, c("^TR_time_s$", "^TR_time", "Thrash.*time"), "TR time")
)

dat <- ntt %>%
  transmute(ID,
            Distance   = .data[[cn$dist]],
            Speed      = .data[[cn$speed]]   * 100,   # m/s -> cm/s
            Meander    = .data[[cn$meander]],
            Bottom     = .data[[cn$bottom]]  * 100,   # m   -> cm
            NB         = .data[[cn$nb]],
            IM         = .data[[cn$im]],
            HS         = .data[[cn$hs]],
            TR         = .data[[cn$tr]]) %>%
  inner_join(ids[, c("ID", "Group", "Sex")], by = "ID")

# speed/height may already be in cm in some exports - sanity check
if (median(dat$Speed, na.rm = TRUE) > 200) {
  message("mean speed looks like it was already cm/s - undoing the x100")
  dat$Speed <- dat$Speed / 100
}
if (median(dat$Bottom, na.rm = TRUE) > 60) {
  message("distance from bottom looks like it was already cm - undoing the x100")
  dat$Bottom <- dat$Bottom / 100
}

## ---- panel definitions (order = panel order in the figure) ----------------

PANELS <- list(
  list(var = "Distance", lab = "Distance travelled (m)"),
  list(var = "Speed",    lab = "Mean speed (cm/s)"),
  list(var = "Meander",  lab = "Meander (\u00b0/m)"),
  list(var = "Bottom",   lab = "Mean distance from bottom (cm)"),
  list(var = "NB",       lab = NA),                 # label depends on figure
  list(var = "IM",       lab = "Immobility (IM) time (s)"),
  list(var = "HS",       lab = "Hyperactive swimming (HS) time (s)"),
  list(var = "TR",       lab = "Thrashing (TR) time (s)")
)

# the ethogram's "normal" state is labelled differently in the two figures
NB_LAB <- c("AB" = "Baseline Exploration (BE) time (s)",
            "CD" = "Normal swimming (NB) time (s)")

## ---- statistics -----------------------------------------------------------
# Shapiro-Wilk per group -> F-test on variances -> Student's or Welch's t;
# Mann-Whitney U if either group departs from normality.

compare_two <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
  if (n1 < 3 || n2 < 3) stop("fewer than 3 usable values in a group")

  sw1 <- shapiro.test(x)$p.value
  sw2 <- shapiro.test(y)$p.value

  if (sw1 > ALPHA && sw2 > ALPHA) {
    vp <- var.test(x, y)$p.value
    eq <- vp > ALPHA
    tt <- t.test(x, y, var.equal = eq)
    sp <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
    out <- data.frame(Test      = if (eq) "Student's t" else "Welch's t",
                      Statistic = unname(tt$statistic),
                      df        = unname(tt$parameter),
                      p         = tt$p.value,
                      Effect    = "Cohen's d",
                      Effect_size = (mean(x) - mean(y)) / sp,
                      F_var_p   = vp)
  } else {
    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    out <- data.frame(Test      = "Mann-Whitney U",
                      Statistic = unname(wt$statistic),
                      df        = NA_real_,
                      p         = wt$p.value,
                      Effect    = "rank-biserial r",
                      Effect_size = 2 * unname(wt$statistic) / (n1 * n2) - 1,
                      F_var_p   = NA_real_)
  }
  cbind(data.frame(n1 = n1, n2 = n2,
                   Mean1 = mean(x), SEM1 = sd(x) / sqrt(n1),
                   Mean2 = mean(y), SEM2 = sd(y) / sqrt(n2),
                   Shapiro_p1 = sw1, Shapiro_p2 = sw2),
        out)
}

stars <- function(p) {
  if (is.na(p))      ""    else
  if (p < 0.001)     "***" else
  if (p < 0.01)      "**"  else
  if (p < 0.05)      "*"   else ""
}

## ---- plotting -------------------------------------------------------------

bar_panel <- function(d, var, ylab, g1, g2, star) {

  d <- d[d$Group %in% c(g1, g2) & is.finite(d[[var]]), ]
  d$Group <- factor(d$Group, levels = c(g1, g2))

  summ <- d %>%
    group_by(Group) %>%
    summarise(m = mean(.data[[var]]), sem = sd(.data[[var]]) / sqrt(n()),
              .groups = "drop")

  top  <- max(c(d[[var]], summ$m + summ$sem), na.rm = TRUE)
  head <- if (nzchar(star)) 1.18 else 1.06

  p <- ggplot(summ, aes(Group, m, fill = Group)) +
    geom_col(width = BAR_W, colour = "black", linewidth = 0.4) +
    geom_errorbar(aes(ymin = m, ymax = m + sem), width = 0.18, linewidth = 0.5) +
    geom_jitter(data = d, aes(y = .data[[var]]), width = 0.13, height = 0,
                size = DOT_SZ, shape = 21, fill = COL[as.character(d$Group)],
                colour = "grey20", stroke = 0.4) +
    scale_fill_manual(values = lighten(COL[c(g1, g2)]), guide = "none") +
    scale_x_discrete(labels = LAB[c(g1, g2)]) +
    scale_y_continuous(limits = c(0, top * head), expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 11) +
    theme(axis.text  = element_text(colour = "black"),
          axis.title.y = element_text(size = 10),
          axis.line  = element_line(linewidth = 0.4),
          plot.margin = margin(4, 6, 2, 4))

  if (nzchar(star))
    p <- p + annotate("text", x = 1.5, y = top * 1.08, label = star,
                      size = 5.2, vjust = 0.4)
  p
}

build_figure <- function(g1, g2, nb_lab, file_stem) {

  tab <- list()
  plots <- list()

  for (i in seq_along(PANELS)) {
    v   <- PANELS[[i]]$var
    ylb <- if (is.na(PANELS[[i]]$lab)) nb_lab else PANELS[[i]]$lab

    d <- dat[dat$Group %in% c(g1, g2), ]
    st <- compare_two(d[[v]][d$Group == g1], d[[v]][d$Group == g2])

    tab[[i]] <- cbind(
      data.frame(Panel      = paste0(if (i <= 4) "A" else "B", ((i - 1) %% 4) + 1),
                 Measure    = ylb,
                 Comparison = paste(LAB[g1], "vs", LAB[g2])),
      st)

    plots[[i]] <- bar_panel(dat, v, ylb, g1, g2, stars(st$p))
  }

  tab <- do.call(rbind, tab)
  tab$p_BH  <- p.adjust(tab$p, "BH")
  tab$Signif <- vapply(tab$p, stars, character(1))

  for (i in 1:8)
    plots[[i]] <- row_tag(plots[[i]], if (i == 1) "A" else if (i == 5) "B" else " ")

  fig <- wrap_plots(plots[1:4], nrow = 1) / wrap_plots(plots[5:8], nrow = 1)

  ggsave(file.path(OUT_DIR, paste0(file_stem, ".png")), fig,
         width = FIG_W, height = FIG_H, dpi = 600, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(file_stem, ".pdf")), fig,
         width = FIG_W, height = FIG_H, device = cairo_pdf)

  num <- sapply(tab, is.numeric)
  tab[num] <- lapply(tab[num], function(z) round(z, 4))
  write.csv(tab, file.path(OUT_DIR, paste0("Stats_", file_stem, ".csv")), row.names = FALSE)

  tab
}

## ---- run ------------------------------------------------------------------

s2 <- build_figure("A", "B", NB_LAB["AB"], "Fig2")
s5 <- build_figure("C", "D", NB_LAB["CD"], "Fig5")

cat("\n--- Fig 2: CTRL vs SR ---\n")
print(s2[, c("Panel", "Measure", "n1", "n2", "Mean1", "Mean2", "Test", "p", "p_BH", "Signif")],
      row.names = FALSE)
cat("\n--- Fig 5: CTRL/7DL vs SR/7DL ---\n")
print(s5[, c("Panel", "Measure", "n1", "n2", "Mean1", "Mean2", "Test", "p", "p_BH", "Signif")],
      row.names = FALSE)
cat("\nfigures and stats tables written to", OUT_DIR, "\n")
