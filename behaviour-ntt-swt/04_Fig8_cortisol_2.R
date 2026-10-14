# ---------------------------------------------------------------------------
# Fig 8   whole-body cortisol
#
#   A  CTRL vs SR            B  CTRL/7DL vs SR/7DL
#
# input: Cortisol.csv with columns  ID, Cortisol_ng_g   (Group is taken from
# id.xlsx; a Group column in the csv is used only if the ID is not in id.xlsx)
#
#   outputs -> Fig8.png/.pdf, Stats_Fig8.csv
# ---------------------------------------------------------------------------

## ---- settings -------------------------------------------------------------

BASE_DIR <- "C:/Users/Barba/Desktop/Sleep_Deprivation_Study"
OUT_DIR  <- file.path(BASE_DIR, "Figures")

ID_FILE   <- "id.xlsx"
CORT_FILE <- "Cortisol.csv"

COL <- c(A = "#1F77B4", B = "#F28E2B", C = "#2E9E4F", D = "#8E5FA8")
LAB <- c(A = "CTRL", B = "SR", C = "CTRL/7DL", D = "SR/7DL")

YLAB    <- "Cortisol (ng/g body weight)"
BAR_TINT <- 0.55   # 0 = full colour bar, 1 = white
ALPHA   <- 0.05
SHOW_NS <- FALSE           # bracket drawn only where the comparison is significant

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

## ---- data -----------------------------------------------------------------

if (!file.exists(CORT_FILE)) {
  tmpl <- data.frame(ID = integer(0), Group = character(0), Cortisol_ng_g = numeric(0))
  write.csv(tmpl, CORT_FILE, row.names = FALSE)
  stop(CORT_FILE, " was not found - an empty template has been written to\n  ",
       file.path(BASE_DIR, CORT_FILE), "\nfill in ID and Cortisol_ng_g and re-run.")
}

ids <- read_excel(ID_FILE)
names(ids) <- tolower(names(ids))
names(ids)[names(ids) == "id"]    <- "ID"
names(ids)[names(ids) == "group"] <- "Group"
names(ids)[names(ids) == "sex"]   <- "Sex"
ids$ID    <- as.character(ids$ID)
ids$Group <- toupper(trimws(ids$Group))

cort <- read.csv(CORT_FILE, stringsAsFactors = FALSE, check.names = TRUE)
names(cort)[grep("^(ID|Fish|Animal|Sample)", names(cort), ignore.case = TRUE)[1]] <- "ID"
vc <- grep("cortisol|ng", names(cort), ignore.case = TRUE, value = TRUE)[1]
if (is.na(vc)) stop("no cortisol column in ", CORT_FILE)

cort$ID <- as.character(cort$ID)
cort$Value <- suppressWarnings(as.numeric(cort[[vc]]))

gi <- grep("^group$", names(cort), ignore.case = TRUE)[1]
if (!is.na(gi)) names(cort)[gi] <- "Group_csv"

dat <- left_join(cort, ids[, c("ID", "Group")], by = "ID")

# fall back to a Group column in the csv for any sample not listed in id.xlsx
if (!is.na(gi)) {
  from_csv <- toupper(trimws(as.character(dat$Group_csv)))
  dat$Group[is.na(dat$Group)] <- from_csv[is.na(dat$Group)]
}

dat <- dat[!is.na(dat$Group) & is.finite(dat$Value), ]

if (any(!cort$ID %in% dat$ID))
  warning("dropped sample(s) with no group or no value: ",
          paste(setdiff(cort$ID, dat$ID), collapse = ", "))

cat("n per group:\n"); print(table(dat$Group))

## ---- statistics -----------------------------------------------------------
# n = 4 per group, so normality is tested but the non-parametric route is the
# likely one; both branches reported in the table.

compare_two <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
  sw1 <- if (n1 >= 3) shapiro.test(x)$p.value else NA_real_
  sw2 <- if (n2 >= 3) shapiro.test(y)$p.value else NA_real_
  norm <- !is.na(sw1) && !is.na(sw2) && sw1 > ALPHA && sw2 > ALPHA

  if (norm) {
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
  if (p < 0.01) "**" else if (p < 0.05) "*" else if (SHOW_NS) "ns" else ""
}

## ---- panel ----------------------------------------------------------------

# both panels share one y scale, as in the published figure
Y_TOP <- max(dat$Value, na.rm = TRUE) * 1.20

cort_panel <- function(g1, g2, star) {
  d <- dat[dat$Group %in% c(g1, g2), ]
  d$Group <- factor(d$Group, levels = c(g1, g2))
  s <- d %>% group_by(Group) %>%
    summarise(m = mean(Value), sem = sd(Value) / sqrt(n()), .groups = "drop")

  top <- max(c(d$Value, s$m + s$sem))
  br  <- Y_TOP * 0.90                   # bracket height, same in both panels
  ylim_top <- Y_TOP

  p <- ggplot(s, aes(Group, m, fill = Group)) +
    geom_col(width = 0.6, colour = "black", linewidth = 0.4) +
    geom_errorbar(aes(ymin = m, ymax = m + sem), width = 0.18, linewidth = 0.5) +
    geom_jitter(data = d, aes(y = Value), width = 0.11, height = 0,
                shape = 21, size = 2.4, fill = COL[as.character(d$Group)],
                colour = "grey20", stroke = 0.5) +
    scale_fill_manual(values = lighten(COL[c(g1, g2)]), guide = "none") +
    scale_x_discrete(labels = LAB[c(g1, g2)]) +
    scale_y_continuous(limits = c(0, ylim_top), expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = YLAB) +
    theme_classic(base_size = 11) +
    theme(axis.text = element_text(colour = "black"),
          axis.line = element_line(linewidth = 0.4))

  if (nzchar(star)) {
    tick <- Y_TOP * 0.02
    p <- p +
      annotate("segment", x = 1, xend = 2, y = br, yend = br, linewidth = 0.4) +
      annotate("segment", x = 1, xend = 1, y = br - tick, yend = br, linewidth = 0.4) +
      annotate("segment", x = 2, xend = 2, y = br - tick, yend = br, linewidth = 0.4) +
      annotate("text", x = 1.5, y = br + tick,
               label = star, size = if (star == "ns") 3.4 else 5.2,
               vjust = if (star == "ns") -0.1 else 0.2)
  }
  p
}

## ---- assemble -------------------------------------------------------------

st <- rbind(
  cbind(data.frame(Panel = "A", Measure = YLAB, Comparison = "CTRL vs SR"),
        compare_two(dat$Value[dat$Group == "A"], dat$Value[dat$Group == "B"])),
  cbind(data.frame(Panel = "B", Measure = YLAB, Comparison = "CTRL/7DL vs SR/7DL"),
        compare_two(dat$Value[dat$Group == "C"], dat$Value[dat$Group == "D"])))
st$Signif <- vapply(st$p, stars, character(1))

fig <- row_tag(cort_panel("A", "B", st$Signif[1]), "A") |
       row_tag(cort_panel("C", "D", st$Signif[2]), "B")

ggsave(file.path(OUT_DIR, "Fig8.png"), fig, width = 6.6, height = 3.4,
       dpi = 600, bg = "white")
ggsave(file.path(OUT_DIR, "Fig8.pdf"), fig, width = 6.6, height = 3.4,
       device = cairo_pdf)

num <- sapply(st, is.numeric)
st[num] <- lapply(st[num], function(z) round(z, 4))
write.csv(st, file.path(OUT_DIR, "Stats_Fig8.csv"), row.names = FALSE)

cat("\n--- Fig 8 ---\n")
print(st[, c("Panel", "Comparison", "n1", "n2", "Mean1", "SEM1", "Mean2", "SEM2",
             "Test", "p", "Signif")], row.names = FALSE)
cat("\nwritten to", OUT_DIR, "\n")
