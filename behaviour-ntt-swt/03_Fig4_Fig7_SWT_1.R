# ---------------------------------------------------------------------------
# Fig 4  (CTRL vs SR)  and  Fig 7  (CTRL/7DL vs SR/7DL)
# Shallow water test, 5 min
#
#   A  session totals: distance travelled while active (m), time immobile (s)
#   B  30-s time course of the same two measures, mean +/- SEM as a ribbon,
#      per-bin asterisks
#
#   outputs -> Fig4.png/.pdf, Fig7.png/.pdf,
#              Stats_Fig4.csv, Stats_Fig7.csv  (totals + per-bin)
# ---------------------------------------------------------------------------

## ---- settings -------------------------------------------------------------

BASE_DIR <- "C:/Users/Barba/Desktop/Sleep_Deprivation_Study"
OUT_DIR  <- file.path(BASE_DIR, "Figures")

ID_FILE    <- "id.xlsx"
TOTAL_FILE <- "Summary_Trunk_5min.csv"
BIN_FILE   <- "Summary_Trunk_30s.csv"

BIN_WIDTH <- 30                 # s
N_BINS    <- 10

COL <- c(A = "#1F77B4", B = "#F28E2B", C = "#2E9E4F", D = "#8E5FA8")
LAB <- c(A = "CTRL", B = "SR", C = "CTRL/7DL", D = "SR/7DL")

RIBBON_ALPHA <- 0.22
BAR_TINT     <- 0.55   # 0 = full colour bar, 1 = white
ALPHA        <- 0.05
STAR_ON_FDR  <- FALSE           # asterisks follow raw p; BH also reported

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
# inline legend placement moved from `legend.position = c(x, y)` to
# `legend.position.inside` in ggplot2 3.5; keep both routes so the script
# runs on either.
inside_legend <- function(x, y) {
  if (utils::packageVersion("ggplot2") >= "3.5.0")
    theme(legend.position = "inside", legend.position.inside = c(x, y))
  else
    theme(legend.position = c(x, y))
}

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

## ---- helpers --------------------------------------------------------------

grab <- function(d, patterns, what) {
  for (p in patterns) {
    hit <- grep(p, names(d), ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  stop("no column matching '", what, "'. present: ", paste(names(d), collapse = ", "))
}

fix_id <- function(d) {
  i <- grep("^(ID|Fish|Animal|Subject)", names(d), ignore.case = TRUE)[1]
  if (is.na(i)) stop("no ID column")
  names(d)[i] <- "ID"
  d$ID <- as.character(d$ID)
  d
}

compare_two <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
  if (n1 < 3 || n2 < 3) return(NULL)
  sw1 <- shapiro.test(x)$p.value
  sw2 <- shapiro.test(y)$p.value
  if (sw1 > ALPHA && sw2 > ALPHA) {
    vp <- var.test(x, y)$p.value
    tt <- t.test(x, y, var.equal = vp > ALPHA)
    sp <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
    res <- data.frame(Test = if (vp > ALPHA) "Student's t" else "Welch's t",
                      Statistic = unname(tt$statistic), df = unname(tt$parameter),
                      p = tt$p.value, Effect = "Cohen's d",
                      Effect_size = (mean(x) - mean(y)) / sp)
  } else {
    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
    res <- data.frame(Test = "Mann-Whitney U", Statistic = unname(wt$statistic),
                      df = NA_real_, p = wt$p.value, Effect = "rank-biserial r",
                      Effect_size = 2 * unname(wt$statistic) / (n1 * n2) - 1)
  }
  cbind(data.frame(n1 = n1, n2 = n2,
                   Mean1 = mean(x), SEM1 = sd(x) / sqrt(n1),
                   Mean2 = mean(y), SEM2 = sd(y) / sqrt(n2),
                   Shapiro_p1 = sw1, Shapiro_p2 = sw2), res)
}

stars <- function(p) {
  if (is.na(p)) "" else if (p < 0.001) "***" else
  if (p < 0.01) "**" else if (p < 0.05) "*" else ""
}

## ---- data -----------------------------------------------------------------

ids <- read_excel(ID_FILE)
names(ids) <- tolower(names(ids))
names(ids)[names(ids) == "id"]    <- "ID"
names(ids)[names(ids) == "group"] <- "Group"
names(ids)[names(ids) == "sex"]   <- "Sex"
ids$ID    <- as.character(ids$ID)
ids$Group <- toupper(trimws(ids$Group))

tot <- fix_id(read.csv(TOTAL_FILE, stringsAsFactors = FALSE, check.names = TRUE))
bin <- fix_id(read.csv(BIN_FILE,   stringsAsFactors = FALSE, check.names = TRUE))

c_act_t <- grab(tot, c("Distance_active", "Active_distance", "Distance.*active"), "active distance (totals)")
c_imm_t <- grab(tot, c("Time_immobile", "Immobile_time", "Immobility_time"),      "immobility time (totals)")
c_act_b <- grab(bin, c("Distance_active", "Active_distance", "Distance.*active"), "active distance (bins)")
c_imm_b <- grab(bin, c("Time_immobile", "Immobile_time", "Immobility_time"),      "immobility time (bins)")
c_bin   <- grab(bin, c("^Bin$", "^Bin_", "^Interval", "^Window", "^Time_bin"),    "bin index")

totals <- tot %>%
  transmute(ID, Active = .data[[c_act_t]], Immobile = .data[[c_imm_t]]) %>%
  inner_join(ids[, c("ID", "Group", "Sex")], by = "ID")

bins <- bin %>%
  transmute(ID, Bin = .data[[c_bin]],
            Active = .data[[c_act_b]], Immobile = .data[[c_imm_b]]) %>%
  inner_join(ids[, c("ID", "Group")], by = "ID")

# bin index may be stored as 1..10 or as the end time in seconds
if (max(bins$Bin, na.rm = TRUE) <= N_BINS) {
  bins$Time <- bins$Bin * BIN_WIDTH
} else {
  bins$Time <- bins$Bin
}

MEASURES <- list(
  list(var = "Active",   lab = "Distance travelled while active (m)"),
  list(var = "Immobile", lab = "Time immobile (s)")
)

## ---- panels ---------------------------------------------------------------

bar_panel <- function(var, ylab, g1, g2, star) {
  d <- totals[totals$Group %in% c(g1, g2) & is.finite(totals[[var]]), ]
  d$Group <- factor(d$Group, levels = c(g1, g2))
  s <- d %>% group_by(Group) %>%
    summarise(m = mean(.data[[var]]), sem = sd(.data[[var]]) / sqrt(n()), .groups = "drop")

  top  <- max(c(d[[var]], s$m + s$sem))
  head <- if (nzchar(star)) 1.18 else 1.06

  p <- ggplot(s, aes(Group, m, fill = Group)) +
    geom_col(width = 0.6, colour = "black", linewidth = 0.4) +
    geom_errorbar(aes(ymin = m, ymax = m + sem), width = 0.18, linewidth = 0.5) +
    geom_jitter(data = d, aes(y = .data[[var]]), width = 0.13, height = 0,
                shape = 21, size = 1.8, fill = COL[as.character(d$Group)],
                colour = "grey20", stroke = 0.4) +
    scale_fill_manual(values = lighten(COL[c(g1, g2)]), guide = "none") +
    scale_x_discrete(labels = LAB[c(g1, g2)]) +
    scale_y_continuous(limits = c(0, top * head), expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 11) +
    theme(axis.text = element_text(colour = "black"),
          axis.line = element_line(linewidth = 0.4))

  if (nzchar(star))
    p <- p + annotate("text", x = 1.5, y = top * 1.08, label = star, size = 5.2)
  p
}

tc_panel <- function(var, ylab, g1, g2, bt) {
  d <- bins[bins$Group %in% c(g1, g2) & is.finite(bins[[var]]), ]
  d$Group <- factor(d$Group, levels = c(g1, g2))

  s <- d %>% group_by(Group, Time) %>%
    summarise(m = mean(.data[[var]]), sem = sd(.data[[var]]) / sqrt(n()), .groups = "drop")

  ann <- bt %>%
    filter(Measure == ylab) %>%
    mutate(star = vapply(if (STAR_ON_FDR) p_BH else p, stars, character(1))) %>%
    filter(nzchar(star)) %>%
    left_join(s %>% group_by(Time) %>% summarise(y = max(m + sem), .groups = "drop"),
              by = "Time")

  ytop <- max(s$m + s$sem) * 1.14

  p <- ggplot(s, aes(Time, m, colour = Group, fill = Group)) +
    geom_ribbon(aes(ymin = m - sem, ymax = m + sem), alpha = RIBBON_ALPHA,
                colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.6) +
    scale_colour_manual(values = COL[c(g1, g2)], labels = LAB[c(g1, g2)], name = NULL) +
    scale_fill_manual(values = COL[c(g1, g2)], guide = "none") +
    scale_x_continuous(breaks = seq(BIN_WIDTH, N_BINS * BIN_WIDTH, by = BIN_WIDTH)) +
    scale_y_continuous(limits = c(0, ytop), expand = expansion(mult = c(0, 0))) +
    labs(x = "Time (s)", y = ylab) +
    theme_classic(base_size = 11) +
    theme(axis.text = element_text(colour = "black"),
          axis.line = element_line(linewidth = 0.4),
          legend.background = element_blank(),
          legend.key.height = unit(0.38, "cm")) +
    inside_legend(0.86, 0.92)

  if (nrow(ann))
    p <- p + geom_text(data = ann, aes(x = Time, y = y * 1.06, label = star),
                       inherit.aes = FALSE, size = 4.2, colour = "black")
  p
}

## ---- assemble -------------------------------------------------------------

build_figure <- function(g1, g2, file_stem) {

  # session totals
  tt <- lapply(MEASURES, function(m) {
    r <- compare_two(totals[[m$var]][totals$Group == g1],
                     totals[[m$var]][totals$Group == g2])
    cbind(data.frame(Panel = "A", Measure = m$lab, Time = NA_real_), r)
  })
  tt <- do.call(rbind, tt)
  tt$p_BH <- p.adjust(tt$p, "BH")

  # per bin, FDR corrected across the ten bins within each measure
  bb <- lapply(MEASURES, function(m) {
    out <- lapply(sort(unique(bins$Time)), function(tm) {
      sub <- bins[bins$Time == tm, ]
      r <- compare_two(sub[[m$var]][sub$Group == g1], sub[[m$var]][sub$Group == g2])
      if (is.null(r)) return(NULL)
      cbind(data.frame(Panel = "B", Measure = m$lab, Time = tm), r)
    })
    out <- do.call(rbind, out)
    out$p_BH <- p.adjust(out$p, "BH")
    out
  })
  bb <- do.call(rbind, bb)

  rowA <- row_tag(bar_panel(MEASURES[[1]]$var, MEASURES[[1]]$lab, g1, g2,
                            stars(tt$p[tt$Measure == MEASURES[[1]]$lab])), "A") |
          row_tag(bar_panel(MEASURES[[2]]$var, MEASURES[[2]]$lab, g1, g2,
                            stars(tt$p[tt$Measure == MEASURES[[2]]$lab])), " ")

  rowB <- row_tag(tc_panel(MEASURES[[1]]$var, MEASURES[[1]]$lab, g1, g2, bb), "B") |
          row_tag(tc_panel(MEASURES[[2]]$var, MEASURES[[2]]$lab, g1, g2, bb), " ")

  fig <- rowA / rowB

  ggsave(file.path(OUT_DIR, paste0(file_stem, ".png")), fig,
         width = 8.5, height = 7.2, dpi = 600, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(file_stem, ".pdf")), fig,
         width = 8.5, height = 7.2, device = cairo_pdf)

  st <- rbind(tt, bb) %>%
    mutate(Comparison = paste(LAB[g1], "vs", LAB[g2]),
           Signif = vapply(p, stars, character(1))) %>%
    relocate(Comparison, .after = Measure)
  num <- sapply(st, is.numeric)
  st[num] <- lapply(st[num], function(z) round(z, 4))
  write.csv(st, file.path(OUT_DIR, paste0("Stats_", file_stem, ".csv")), row.names = FALSE)
  st
}

s4 <- build_figure("A", "B", "Fig4")
s7 <- build_figure("C", "D", "Fig7")

cat("\n--- Fig 4: CTRL vs SR ---\n")
print(s4[, c("Panel", "Measure", "Time", "Mean1", "Mean2", "Test", "p", "p_BH", "Signif")],
      row.names = FALSE)
cat("\n--- Fig 7: CTRL/7DL vs SR/7DL ---\n")
print(s7[, c("Panel", "Measure", "Time", "Mean1", "Mean2", "Test", "p", "p_BH", "Signif")],
      row.names = FALSE)
cat("\nwritten to", OUT_DIR, "\n")
