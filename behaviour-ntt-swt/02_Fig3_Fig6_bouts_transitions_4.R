# ---------------------------------------------------------------------------
# Fig 3  (CTRL vs SR)  and  Fig 6  (CTRL/7DL vs SR/7DL)
#
#   A  left   state changes per minute
#      right  mean bout duration per state (log10 y axis)
#   B         state-transition diagrams, one per group
#
# panel A comes from the 6-min ethogram columns (time + episode count per
# state).  panel B needs frame-level states: either Predicted_Behaviors.rds
# or a pre-computed Transition_matrices.csv (Group, From, To, Prob).
#
#   outputs -> Fig3.png/.pdf, Fig6.png/.pdf, Stats_Fig3.csv, Stats_Fig6.csv,
#              Transition_matrices.csv
# ---------------------------------------------------------------------------

## ---- settings -------------------------------------------------------------

BASE_DIR <- "C:/Users/Barba/Desktop/Sleep_Deprivation_Study"
OUT_DIR  <- file.path(BASE_DIR, "Figures")

ID_FILE    <- "id.xlsx"
NTT_FILE   <- "NTT_6min_with_behaviors.csv"
FRAME_FILE <- "Predicted_Behaviors.rds"          # frame-level classifier output
TRANS_FILE <- "Transition_matrices.csv"          # cached / supplied alternative

SESSION_MIN <- 6                                 # NTT duration, minutes
FPS         <- 30                                # frame rate of the recordings

COL <- c(A = "#1F77B4", B = "#F28E2B", C = "#2E9E4F", D = "#8E5FA8")
LAB <- c(A = "CTRL", B = "SR", C = "CTRL/7DL", D = "SR/7DL")

# state order used on the bout axis and in the diagrams
STATES     <- c("HS", "IM", "BE", "TR")
STATE_FULL <- c(HS = "HS", IM = "IM", BE = "BE", TR = "TR")

# the classifier writes NB where the figures say BE
STATE_RAW   <- c("IM", "NB", "HS", "TR")       # as stored in the rds
STATE_ALIAS <- c(NB = "BE", BE = "BE", IM = "IM", HS = "HS", TR = "TR")

MIN_PROB   <- 0.08   # transitions rarer than this are not drawn
NODE_R_MIN <- 0.19   # node radius in data units: smallest and largest state
NODE_R_MAX <- 0.34
BOW        <- 0.30   # how far the arcs bow away from the straight line
LABEL_PAD  <- 0.17   # gap between an arc and its probability label
BAR_TINT  <- 0.55   # 0 = full colour bar, 1 = white
ALPHA     <- 0.05

setwd(BASE_DIR)
dir.create(OUT_DIR, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
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

## ---- data: session-level ethogram ----------------------------------------

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

# map each figure state to its time / episode columns in the csv.
# patterns are tried in order, most specific first - "^Immobil" would
# otherwise pick up the tracking-based immobility columns instead of the
# ethogram's IM state.
col_for <- function(state, kind) {
  pats <- switch(state,
    BE = c("^NB_", "^BE_", "^Normal"),
    IM = c("^IM_", "^Immobil"),
    HS = c("^HS_", "^Hyperactiv"),
    TR = c("^TR_", "^Thrash"))
  tail_pat <- if (kind == "time") "time" else "episode"
  ok <- grepl(tail_pat, names(ntt), ignore.case = TRUE)
  for (p in pats) {
    hit <- names(ntt)[ok & grepl(p, names(ntt), ignore.case = TRUE)]
    if (length(hit)) return(hit[1])
  }
  stop("no ", kind, " column found for state ", state,
       ". columns present: ", paste(names(ntt), collapse = ", "))
}

eth <- data.frame(ID = ntt$ID)
for (s in STATES) {
  tc <- col_for(s, "time"); ec <- col_for(s, "episodes")
  eth[[paste0(s, "_time")]] <- ntt[[tc]]
  eth[[paste0(s, "_eps")]]  <- ntt[[ec]]
  eth[[paste0(s, "_bout")]] <- ifelse(ntt[[ec]] > 0, ntt[[tc]] / ntt[[ec]], NA_real_)
}
eth$Changes_per_min <- rowSums(eth[paste0(STATES, "_eps")], na.rm = TRUE) / SESSION_MIN
eth$Total_time      <- rowSums(eth[paste0(STATES, "_time")], na.rm = TRUE)

eth <- inner_join(eth, ids[, c("ID", "Group", "Sex")], by = "ID")

## ---- statistics -----------------------------------------------------------

compare_two <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
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

## ---- panel A left: state changes per minute -------------------------------

panel_changes <- function(g1, g2, star) {
  d <- eth[eth$Group %in% c(g1, g2), ]
  d$Group <- factor(d$Group, levels = c(g1, g2))
  s <- d %>% group_by(Group) %>%
    summarise(m = mean(Changes_per_min), sem = sd(Changes_per_min) / sqrt(n()), .groups = "drop")

  top  <- max(c(d$Changes_per_min, s$m + s$sem))
  head <- if (nzchar(star)) 1.18 else 1.06

  p <- ggplot(s, aes(Group, m, fill = Group)) +
    geom_col(width = 0.6, colour = "black", linewidth = 0.4) +
    geom_errorbar(aes(ymin = m, ymax = m + sem), width = 0.18, linewidth = 0.5) +
    geom_jitter(data = d, aes(y = Changes_per_min), width = 0.13, height = 0,
                shape = 21, size = 1.8, fill = COL[as.character(d$Group)],
                colour = "grey20", stroke = 0.4) +
    scale_fill_manual(values = lighten(COL[c(g1, g2)]), guide = "none") +
    scale_x_discrete(labels = LAB[c(g1, g2)]) +
    scale_y_continuous(limits = c(0, top * head), expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = "State changes per minute") +
    theme_classic(base_size = 11) +
    theme(axis.text = element_text(colour = "black"),
          axis.line = element_line(linewidth = 0.4))

  if (nzchar(star))
    p <- p + annotate("text", x = 1.5, y = top * 1.08, label = star, size = 5.2)
  p
}

## ---- panel A right: mean bout duration per state --------------------------

panel_bouts <- function(g1, g2, bt) {
  d <- eth %>%
    filter(Group %in% c(g1, g2)) %>%
    select(ID, Group, all_of(paste0(STATES, "_bout"))) %>%
    pivot_longer(-c(ID, Group), names_to = "State", values_to = "Bout") %>%
    mutate(State = factor(sub("_bout$", "", State), levels = STATES),
           Group = factor(Group, levels = c(g1, g2))) %>%
    filter(is.finite(Bout), Bout > 0)

  s <- d %>% group_by(State, Group) %>%
    summarise(m = mean(Bout), sem = sd(Bout) / sqrt(n()), .groups = "drop")

  ann <- bt %>%
    transmute(State = factor(State, levels = STATES), star = vapply(p, stars, character(1))) %>%
    filter(nzchar(star)) %>%
    left_join(s %>% group_by(State) %>% summarise(y = max(m + sem), .groups = "drop"),
              by = "State")

  pd <- position_dodge(width = 0.55)

  p <- ggplot(s, aes(State, m, colour = Group)) +
    geom_errorbar(aes(ymin = m - sem, ymax = m + sem), width = 0.16,
                  linewidth = 0.6, position = pd) +
    geom_point(size = 2.8, position = pd) +
    scale_colour_manual(values = COL[c(g1, g2)], labels = LAB[c(g1, g2)], name = NULL) +
    scale_y_log10() +
    annotation_logticks(sides = "l", outside = FALSE,
                        short = unit(0.05, "cm"), mid = unit(0.08, "cm"),
                        long = unit(0.12, "cm")) +
    labs(x = NULL, y = "Mean bout duration (s)") +
    theme_classic(base_size = 11) +
    theme(axis.text = element_text(colour = "black"),
          axis.line = element_line(linewidth = 0.4),
          legend.background = element_blank(),
          legend.key.height = unit(0.38, "cm")) +
    inside_legend(0.86, 0.92)

  if (nrow(ann))
    p <- p + geom_text(data = ann, aes(x = State, y = y * 1.45, label = star),
                       inherit.aes = FALSE, size = 5, colour = "black")
  p
}

## ---- transition matrices --------------------------------------------------
# conventions taken from Script_NTT_Transitions_1.R, which reproduces
# NTT_6min_with_behaviors.csv exactly:
#   * frames 0 .. (360*30 - 1) only
#   * an episode is a run of >= MINRUN consecutive frames of the same state;
#     shorter runs are absorbed into the preceding run before anything is
#     counted, which removes classifier flicker
#   * transitions are taken BETWEEN BOUTS, not between frames, so the diagonal
#     is empty by construction and each row is normalised over the three
#     possible destinations
# computed per fish, then averaged within group so no one recording dominates.

FPS    <- 30
NFRAME <- 360 * FPS
MINRUN <- 10          # frames; shorter runs are absorbed

# runs shorter than minrun are absorbed into the preceding run, then re-rle'd
absorb_short <- function(v, minrun = MINRUN) {
  r <- rle(v)
  keep <- r$lengths >= minrun
  if (!any(keep)) return(rep(r$values[which.max(r$lengths)], length(v)))
  out <- rep(NA_character_, length(r$values))
  out[keep] <- r$values[keep]
  first_keep <- which(keep)[1]
  if (first_keep > 1) out[1:(first_keep - 1)] <- r$values[first_keep]
  for (i in seq_along(out)) if (is.na(out[i])) out[i] <- out[i - 1]
  rep(out, r$lengths)
}

# bout-level matrix for one fish: consecutive distinct bouts, rows summing to 1
bout_matrix <- function(v) {
  b <- rle(v)$values
  if (length(b) < 2) return(NULL)   # levels are the classifier's own names
  m <- table(factor(head(b, -1), STATE_RAW), factor(tail(b, -1), STATE_RAW))
  m <- as.matrix(m); storage.mode(m) <- "double"
  rs <- rowSums(m)
  sweep(m, 1, ifelse(rs == 0, NA, rs), "/")
}

transitions_from_frames <- function(path) {
  message("reading ", path, " ...")
  pred <- readRDS(path)
  if (!is.list(pred) || is.data.frame(pred))
    stop(path, " was expected to be a list of per-fish data frames")

  # ids live in the list names, e.g. "007_something" -> 7
  fish <- as.integer(gsub("^0*(\\d+).*", "\\1", trimws(names(pred))))
  message("  ", length(pred), " recordings; ",
          sum(!is.na(fish)), " with a readable id")

  per_fish <- lapply(seq_along(pred), function(i) {
    d <- pred[[i]]
    stc <- if ("behavior.smoothed" %in% names(d)) "behavior.smoothed" else names(d)[ncol(d)]
    if ("frame" %in% names(d)) {
      d <- d[d$frame >= 0 & d$frame < NFRAME, ]
      d <- d[order(d$frame), ]
    } else {
      d <- head(d, NFRAME)
    }
    v <- toupper(trimws(as.character(d[[stc]])))
    v <- v[v %in% STATE_RAW]
    if (length(v) < 2) return(NULL)

    m <- bout_matrix(absorb_short(v))
    if (is.null(m)) return(NULL)

    expand.grid(From = STATE_RAW, To = STATE_RAW, stringsAsFactors = FALSE) %>%
      mutate(ID = as.character(fish[i]),
             Prob = mapply(function(a, b) m[a, b], From, To))
  })

  per_fish <- bind_rows(per_fish)
  if (!nrow(per_fish)) stop("no usable sequences in ", path)

  matched <- length(intersect(unique(per_fish$ID), ids$ID))
  message("  ", matched, " of ", length(unique(per_fish$ID)),
          " recordings matched an ID in ", ID_FILE)
  if (matched == 0)
    stop("no recording matched an ID in ", ID_FILE,
         ".\n  rds ids look like: ", paste(head(unique(per_fish$ID), 5), collapse = ", "))

  per_fish %>%
    inner_join(ids[, c("ID", "Group")], by = "ID") %>%
    filter(From != To) %>%
    # rename the classifier's NB to the figure's BE
    mutate(From = unname(STATE_ALIAS[From]),
           To   = unname(STATE_ALIAS[To])) %>%
    group_by(Group, From, To) %>%
    summarise(Prob = mean(Prob, na.rm = TRUE), .groups = "drop")
}

if (file.exists(FRAME_FILE)) {
  # a problem reading the frame file should not cost us panel A
  trans <- tryCatch({
    tm <- transitions_from_frames(FRAME_FILE)
    write.csv(tm, file.path(OUT_DIR, TRANS_FILE), row.names = FALSE)
    tm
  }, error = function(e) {
    warning("could not build transition matrices: ", conditionMessage(e),
            "\n  panel B will be skipped")
    NULL
  })
} else if (file.exists(TRANS_FILE) || file.exists(file.path(OUT_DIR, TRANS_FILE))) {
  f <- if (file.exists(TRANS_FILE)) TRANS_FILE else file.path(OUT_DIR, TRANS_FILE)
  message("using cached transition matrices from ", f)
  trans <- read.csv(f, stringsAsFactors = FALSE)
} else {
  warning("neither ", FRAME_FILE, " nor ", TRANS_FILE,
          " found - panel B will be skipped")
  trans <- NULL
}

## ---- panel B: transition diagram ------------------------------------------

# square layout, same node positions in every panel
NODE <- data.frame(
  State = c("IM", "TR", "BE", "HS"),
  x     = c(-1,    1,   -1,    1),
  y     = c( 1,    1,   -1,   -1),
  stringsAsFactors = FALSE)

occupancy <- function(g) {
  d <- eth[eth$Group == g, ]
  tot <- rowSums(d[paste0(STATES, "_time")], na.rm = TRUE)
  sapply(STATES, function(s) mean(d[[paste0(s, "_time")]] / tot, na.rm = TRUE))
}

# quadratic bezier from (x1,y1) to (x2,y2), bowed sideways by `bow`. drawing the
# arcs ourselves rather than with geom_curve keeps the geometry in data units,
# so arrow heads meet the circles and the labels sit on the arc at any figure
# size - geom_curve works in device space and drifts when the figure is resized.
arc_points <- function(x1, y1, x2, y2, bow, n = 64) {
  dx <- x2 - x1; dy <- y2 - y1
  cx <- (x1 + x2) / 2 - dy * bow
  cy <- (y1 + y2) / 2 + dx * bow
  tt <- seq(0, 1, length.out = n)
  data.frame(x = (1 - tt)^2 * x1 + 2 * (1 - tt) * tt * cx + tt^2 * x2,
             y = (1 - tt)^2 * y1 + 2 * (1 - tt) * tt * cy + tt^2 * y2)
}

arc_apex <- function(x1, y1, x2, y2, bow) {
  dx <- x2 - x1; dy <- y2 - y1
  c(x = (x1 + x2) / 2 - dy * bow / 2,
    y = (y1 + y2) / 2 + dx * bow / 2)
}

network_panel <- function(g) {

  occ <- occupancy(g)
  rad <- NODE_R_MIN + (NODE_R_MAX - NODE_R_MIN) * sqrt(occ / max(occ))

  nd <- NODE
  nd$r   <- rad[nd$State]
  nd$pct <- paste0(round(100 * occ[nd$State]), "%")
  # the percentage sits just outside the circle, pushed away from the centre
  len <- sqrt(nd$x^2 + nd$y^2)
  nd$px <- nd$x + nd$x / len * (nd$r + 0.20)
  nd$py <- nd$y + nd$y / len * (nd$r + 0.20)

  ed <- trans %>%
    filter(Group == g, Prob >= MIN_PROB,
           From %in% NODE$State, To %in% NODE$State) %>%   # drop unknown states
    left_join(NODE, by = c("From" = "State")) %>% rename(x1 = x, y1 = y) %>%
    left_join(NODE, by = c("To"   = "State")) %>% rename(x2 = x, y2 = y)

  if (!nrow(ed)) {
    warning("no transitions above MIN_PROB for group ", g)
    return(ggplot() + theme_void())
  }

  # the two directions of a pair bow opposite ways, forming a lens
  # both directions take the SAME bow: the control point is offset perpendicular
  # to the direction of travel, and reversing the direction already flips that
  # perpendicular. opposite bow values cancel the flip and stack the two arcs on
  # top of each other instead of opening them into a lens.
  ed$bow <- BOW

  paths <- list(); labs <- list()
  for (i in seq_len(nrow(ed))) {
    x1 <- ed$x1[i]; y1 <- ed$y1[i]; x2 <- ed$x2[i]; y2 <- ed$y2[i]
    dx <- x2 - x1;  dy <- y2 - y1;  L <- sqrt(dx^2 + dy^2)
    # trim to the circle edges before building the arc
    r1 <- rad[ed$From[i]]; r2 <- rad[ed$To[i]] + 0.04
    ax1 <- x1 + dx / L * r1; ay1 <- y1 + dy / L * r1
    ax2 <- x2 - dx / L * r2; ay2 <- y2 - dy / L * r2

    p <- arc_points(ax1, ay1, ax2, ay2, ed$bow[i])
    p$grp <- i; p$Prob <- ed$Prob[i]
    paths[[i]] <- p

    # unname first: a named value would be taken as a row name by data.frame
    ap <- unname(arc_apex(ax1, ay1, ax2, ay2, ed$bow[i]))
    # step further along the same normal the arc was bowed by, so the label
    # ends up on the outside of its own arc
    aL <- sqrt((ax2 - ax1)^2 + (ay2 - ay1)^2)
    nx <- -(ay2 - ay1) / aL; ny <- (ax2 - ax1) / aL
    labs[[i]] <- data.frame(x   = ap[1] + nx * LABEL_PAD,
                            y   = ap[2] + ny * LABEL_PAD,
                            lab = sprintf("%.2f", ed$Prob[i]),
                            row.names = NULL, stringsAsFactors = FALSE)
  }
  paths <- do.call(rbind, paths)
  labs  <- do.call(rbind, labs)

  # node circles as polygons, so their radius is in the same units as the arcs
  circles <- do.call(rbind, lapply(seq_len(nrow(nd)), function(i) {
    th <- seq(0, 2 * pi, length.out = 90)
    data.frame(x = nd$x[i] + nd$r[i] * cos(th),
               y = nd$y[i] + nd$r[i] * sin(th),
               State = nd$State[i])
  }))

  # arrow head drawn only on the last segment of each arc
  heads <- paths %>% group_by(grp) %>% slice_tail(n = 2) %>% ungroup()

  ggplot() +
    geom_path(data = paths, aes(x, y, group = grp, linewidth = Prob),
              colour = COL[[g]], lineend = "round") +
    geom_path(data = heads, aes(x, y, group = grp, linewidth = Prob),
              colour = COL[[g]],
              arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +
    geom_polygon(data = circles, aes(x, y, group = State),
                 fill = COL[[g]], colour = "black", linewidth = 0.5) +
    geom_text(data = nd, aes(x, y, label = State), colour = "white",
              fontface = "bold", size = 4.0) +
    geom_text(data = nd, aes(px, py, label = pct), colour = "grey30",
              size = 3.4) +
    geom_text(data = labs, aes(x, y, label = lab), colour = "grey15",
              size = 3.3) +
    scale_linewidth(range = c(0.5, 3.2), guide = "none") +
    coord_equal(xlim = c(-1.75, 1.75), ylim = c(-1.75, 1.75), clip = "off") +
    labs(subtitle = LAB[[g]]) +
    theme_void(base_size = 11) +
    theme(plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"),
          plot.margin = margin(6, 10, 6, 10))
}

## ---- assemble -------------------------------------------------------------

build_figure <- function(g1, g2, file_stem) {

  tab <- list()

  d1 <- eth$Changes_per_min[eth$Group == g1]
  d2 <- eth$Changes_per_min[eth$Group == g2]
  tab[["changes"]] <- cbind(data.frame(Panel = "A1", Measure = "State changes per minute"),
                            compare_two(d1, d2))

  bt <- lapply(STATES, function(s) {
    v <- paste0(s, "_bout")
    cbind(data.frame(Panel = "A2", State = s,
                     Measure = paste0(s, " mean bout duration (s)")),
          compare_two(eth[[v]][eth$Group == g1], eth[[v]][eth$Group == g2]))
  })
  bt <- do.call(rbind, bt)
  bt$p_BH <- p.adjust(bt$p, "BH")          # FDR within the four states

  pA <- (row_tag(panel_changes(g1, g2, stars(tab[["changes"]]$p)), "A") |
         row_tag(panel_bouts(g1, g2, bt), " ")) +
        plot_layout(widths = c(1, 1.7))

  if (!is.null(trans)) {
    pB  <- row_tag(network_panel(g1), "B") | row_tag(network_panel(g2), " ")
    fig <- (pA / pB) + plot_layout(heights = c(1, 1.15))
    H <- 9.0
  } else {
    fig <- pA
    H <- 4.2
  }

  ggsave(file.path(OUT_DIR, paste0(file_stem, ".png")), fig,
         width = 9.5, height = H, dpi = 600, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(file_stem, ".pdf")), fig,
         width = 9.5, height = H, device = cairo_pdf)

  st <- bind_rows(
    tab[["changes"]] %>% mutate(State = NA_character_, p_BH = NA_real_),
    bt) %>%
    mutate(Comparison = paste(LAB[g1], "vs", LAB[g2]),
           Signif = vapply(p, stars, character(1))) %>%
    relocate(Comparison, .after = Measure)

  # occupancy printed alongside, since it sets the node sizes in panel B
  occ <- bind_rows(
    data.frame(Group = LAB[[g1]], State = STATES, Occupancy_pct = round(100 * occupancy(g1), 1)),
    data.frame(Group = LAB[[g2]], State = STATES, Occupancy_pct = round(100 * occupancy(g2), 1)))
  write.csv(occ, file.path(OUT_DIR, paste0("Occupancy_", file_stem, ".csv")), row.names = FALSE)

  num <- sapply(st, is.numeric)
  st[num] <- lapply(st[num], function(z) round(z, 4))
  write.csv(st, file.path(OUT_DIR, paste0("Stats_", file_stem, ".csv")), row.names = FALSE)
  st
}

s3 <- build_figure("A", "B", "Fig3")
s6 <- build_figure("C", "D", "Fig6")

cat("\n--- Fig 3: CTRL vs SR ---\n")
print(s3[, c("Panel", "Measure", "Mean1", "Mean2", "Test", "p", "p_BH", "Signif")], row.names = FALSE)
cat("\n--- Fig 6: CTRL/7DL vs SR/7DL ---\n")
print(s6[, c("Panel", "Measure", "Mean1", "Mean2", "Test", "p", "p_BH", "Signif")], row.names = FALSE)
cat("\nwritten to", OUT_DIR, "\n")
