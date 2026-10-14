# Behavioural tests and whole-body cortisol

Part of the repository for [paper reference]; the home-tank video analysis is in
`../home-tank-activity`.

Four R scripts, one per figure pair. Each builds the figure and writes the
statistics table behind it.

| Script | Figures | Content |
|---|---|---|
| `01_ntt_totals.R` | Fig 2, Fig 5 | Novel tank test, 6-minute session totals: kinematics (row A) and ethogram state times (row B) |
| `02_bouts_transitions.R` | Fig 3, Fig 6 | State changes per minute and mean bout duration per state (row A), state-transition diagrams per group (row B) |
| `03_swt.R` | Fig 4, Fig 7 | Shallow water test, 5 minutes: session totals (row A) and 30-second time course (row B) |
| `04_cortisol.R` | Fig 8 | Whole-body cortisol, both comparisons side by side |

Figures 2, 3, 4 compare CTRL with SR; figures 5, 6, 7 compare CTRL/7DL with
SR/7DL; figure 8 shows both comparisons.

## Groups

| Code in the data | Group | Colour |
|---|---|---|
| A | CTRL, tested immediately after the protocol | blue, `#1F77B4` |
| B | SR, tested immediately after the protocol | orange, `#F28E2B` |
| C | CTRL/7DL, tested seven days later | green, `#2E9E4F` |
| D | SR/7DL, tested seven days later | purple, `#8E5FA8` |

## Input files

All of them live in `data/`:

```
id.xlsx                       ID, Group (A-D), Sex
NTT_6min_with_behaviors.csv   one row per fish: distance, mean speed, meander,
                              mean distance from bottom, and time and episode
                              count for each ethogram state (NB, IM, HS, TR)
Summary_Trunk_5min.csv        shallow water test totals: distance travelled
                              while active, time immobile
Summary_Trunk_30s.csv         the same two measures per 30-second bin
Cortisol.csv                  ID and cortisol in ng per g body weight
Predicted_Behaviors.rds       frame-level states from the classifier, used for
                              the transition diagrams
```

Column names are resolved by pattern rather than by exact match, because the
tracking export has changed names between versions; a script stops and lists the
columns it did find when nothing matches. Fish are joined to groups through
`id.xlsx` by ID.

If `Predicted_Behaviors.rds` is absent, `02_bouts_transitions.R` falls back to a
cached `Transition_matrices.csv` (Group, From, To, Prob), and if that is missing
too it writes the figure without the transition diagrams rather than failing.

## Running

```bash
Rscript 01_ntt_totals.R
Rscript 02_bouts_transitions.R
Rscript 03_swt.R
Rscript 04_cortisol.R
```

Outputs go to `results/`: each figure as PNG at 600 dpi and as PDF with live
text, plus `Stats_FigN.csv`. `02_bouts_transitions.R` also writes
`Occupancy_FigN.csv` (the share of session time in each state, which sets the
node sizes in the diagrams) and `Transition_matrices.csv`.

The data folder and the output folder can be moved without editing the scripts:

```bash
SLEEP_STUDY_DIR=/path/to/data SLEEP_STUDY_OUT=/path/to/out Rscript 03_swt.R
```

R 4.x with readxl, dplyr, tidyr, ggplot2 and patchwork.

## Statistics

The animal is the unit. Each comparison runs the same route: Shapiro-Wilk on
both groups; if both pass, an F test on the variances decides between Student's
and Welch's t-test; otherwise a Mann-Whitney U test. The test that was used is
named in every row of the statistics table, together with the group means and
standard errors, the Shapiro-Wilk p-values, the test statistic, the degrees of
freedom, and an effect size, Cohen's d for the t-tests and the rank-biserial
correlation for Mann-Whitney.

Asterisks on the figures follow the raw p-value (`*` below 0.05, `**` below
0.01, `***` below 0.001). Benjamini-Hochberg corrected p-values are reported
alongside in the tables: across the eight measures of the novel tank test,
across the four states for bout duration, and across the ten bins within each
measure of the time course. In `03_swt.R`, `STAR_ON_FDR` switches the asterisks
to the corrected values.

## Figure conventions

Bars show the group mean with the standard error, drawn in a pale tint of the
group colour with a black outline, and every animal appears as a point in the
full group colour, so the sample size is visible in the figure itself. Time
courses show the mean with a standard error ribbon. Bout durations are drawn on
a log scale, since they span two orders of magnitude.

In the transition diagrams, node size is the share of session time spent in that
state, arrow width is the transition probability, and transitions below
`MIN_PROB` (0.08) are left out to keep the diagram readable. Transitions are
counted between bouts rather than between frames, so the diagonal is empty by
construction and each row is normalised over the three possible destinations.
Runs shorter than `MINRUN` (10 frames, a third of a second) are absorbed into
the preceding bout first, which removes classifier flicker.

The ethogram state that the classifier stores as NB is shown as BE, baseline
exploration, in the figures of the first comparison; the alias is set in
`02_bouts_transitions.R` and the label in `01_ntt_totals.R`.
