# Zebrafish sleep restriction: behaviour and home-tank activity

Code and processed data for [paper reference], a study of adult zebrafish
exposed to seven days of sleep restriction, tested immediately after the
protocol and again seven days later.

## The experiment

Sleep was restricted by shortening the dark phase from ten to six hours per day
for seven days. Four groups were compared:

| Code | Group |
|---|---|
| CTRL | control, tested immediately after the protocol |
| SR | sleep restricted, tested immediately after the protocol |
| CTRL/7DL | control, tested seven days later |
| SR/7DL | sleep restricted, tested seven days later |

Three measurements are covered by this repository: the novel tank test, the
shallow water test, cortisol and video monitoring of activity in the home tanks during
the recovery period.

## What is in this repository

```
behaviour-ntt-swt/        novel tank test, shallow water test and cortisol: R analysis
home-tank-activity/       video analysis of activity in the home tanks
```

Each folder has its own README with the measures, the files, and the commands
that reproduce its outputs. Both folders ship the processed data they analyse,
so every table and figure can be reproduced without the raw videos.

**behaviour-ntt-swt** holds one standalone R script per figure. Both produce
whole-session measures and 30-second time courses, for the two comparisons the
paper reports (CTRL against SR, and CTRL/7DL against SR/7DL), with individual
animals shown as points, normality checked before each test, raw p-values marked
on the figures.

**home-tank-activity** holds the Python pipeline that measures activity from the
home-tank recordings and the R script that produces the statistics and the
figure. Activity is measured as the fraction of pixels that change between
consecutive frames inside each tank region, a measure that does not depend on
tracking individual fish.

## Reproducing the published results

Both folders run on the processed data included here.

```bash
cd behaviour-ntt-swt   
cd home-tank-activity  
```

Reprocessing the raw recordings needs the videos and is described in
`home-tank-activity/README.md`.

Requirements: R 4.x with dplyr, tidyr, ggplot2, patchwork, lme4 and lmerTest;
for the video pipeline, Python 3.10 with the packages in
`home-tank-activity/requirements.txt` and FFmpeg on the path.

## Colours used in the figures

Blue for CTRL and orange for SR; green for CTRL/7DL and purple for SR/7DL. The
same palette is used in both folders.

## Data

The processed data are included: per-animal measures for the two behavioural
tests, and per-tank, per-recording, per-10-minute-bin activity for the home-tank
recordings.

## How to read the statistics

The behavioural tests use the animal as the unit. The home-tank recordings use
the tank as the unit: fish within a tank share water, feeding and handling, so
they are not independent replicates, which leaves four tanks per condition and
makes that part of the study exploratory.

## Citation

If you use this code or the data, please cite:
> Fontana et al. (in preparation). [paper title]
