# Home-tank activity in group-housed zebrafish

Video analysis and statistics for the home-tank activity measurements reported
in [paper reference]. Everything needed to reproduce the published figure and
tables from the raw recordings is here: the code, the parameters that were used,
the validation data, and the limitations of the measurement.

## What is measured

**Activity index (primary measure).** The fraction of pixels inside a tank
region that change between consecutive processed frames, at 2 frames per second.
It reflects the movement of all fish in a tank together. It is measured from the
images alone and does not depend on the fish detector.

**Fish count per frame (quality control only).** The detector reports how many
fish it can see, which is used to check that a tank is being observed properly.
Its accuracy is not sufficient for measures of individual fish, so no per-fish
endpoint is computed or reported. Individual tracking, swimming speed, immobility
and sleep bouts are not part of this repository: they were tested during
development, failed validation, and were removed.

## Validation of the detector

An observer counted fish in 24 frames (48 tank-frames) from four recordings, on
clean images, without seeing the detector output.

| Subset | Tank-frames | Mean absolute deviation | Mean signed deviation | Within 1 fish |
|---|---|---|---|---|
| All | 48 | 0.79 fish | -0.08 | 81% |
| Recordings not used for parameter selection | 24 | 0.67 fish | -0.17 | 88% |

Manual counts averaged 5.04 of 6 fish visible per tank-frame, because fish
occlude one another in the side view, so six detections is not the correct target
in every frame. This is why per-fish measures were dropped.

Raw validation data: `validation/manual_counts.csv` and `validation/frames/`.

## Files

```
sleep_core.py          shared configuration and image analysis
run_pipeline.py        main run over all recordings  -> results/
tune_parameters.py     parameter selection by grid search
validate.py            counting sheets, detector score, checks on the output
analysis.R             statistics, tables and figures -> results/analysis/
tuned_params.json      detection parameters per cohort, camera and light phase
fish_alive.csv         number of fish alive per tank from a given date
```

Recordings are named `CAMERA R DDMMYYYY HHMMSS.mp4`. The first cohort used
underscores and the second spaces; both are read.

## Reproducing the analysis

Python 3.10 with OpenCV, NumPy, Matplotlib, and FFmpeg on the path
(`pip install -r requirements.txt`). R 4.x with dplyr, tidyr, ggplot2,
patchwork, lme4 and lmerTest.

```bash
# 1. process every recording (one to two hours; resumable)
python run_pipeline.py --folder <folder with the recordings> --workers 2

# 2. statistics, tables and figures
Rscript analysis.R          # run from the folder that holds results/
```

Optional, not needed to reproduce the published numbers:

```bash
python tune_parameters.py --folder <folder>    # re-select detection parameters
python validate.py sheet    --folder <folder>  # frames for a blinded observer
python validate.py score    --folder <folder>  # detector against manual counts
python validate.py activity --folder <folder>  # checks on the pipeline output
```

`tune_parameters.py` overwrites `tuned_params.json`, including the values used
in the paper; keep a copy before re-running it.

## Output

`results/activity_binned.csv`, one row per recording, tank and 10-minute bin:

| Column | Meaning |
|---|---|
| video, batch, camera, group, date, time | recording identifiers |
| day_night | light phase of the bin |
| tank, time_bin | tank and 10-minute bin index |
| mean_motion | activity index, the measure analysed |
| mean_activity | fraction of the region differing from the background |
| mean_fish, pct_frames_6, mean_blobs | detector output, quality control |
| n_fish | fish alive in that tank on that date, from `fish_alive.csv` |
| n_frames | frames contributing to the bin |

`results/qc_summary.csv` holds one row per recording. `analysis.R` writes its
tables, the figure and `session_info.txt` to `results/analysis/`.

## Processing parameters, as used

Two frames per second processed; 3 x 3 median filter; the mean intensity of each
frame matched to the background to remove camera gain drift; 3 x 3 elliptical
opening; 10-minute bins, with the first bin of each recording discarded as
handling-related.

Motion threshold: 16 grey levels in the light phase, 12 in the dark phase.

Tank regions (x0, y0, x1, y1):

| Cohort | Camera | Tank 1 (left) | Tank 2 (right) |
|---|---|---|---|
| 1 | 003 | 45, 350, 600, 710 | 640, 350, 1238, 710 |
| 1 | 004 | 42, 300, 610, 690 | 645, 300, 1240, 690 |
| 2 | 003 | 15, 350, 610, 715 | 655, 350, 1240, 715 |
| 2 | 004 | 20, 305, 622, 715 | 672, 305, 1250, 715 |

Detection parameters are in `tuned_params.json`. They were selected per cohort,
camera and light phase by grid search over the background percentile (15, 50,
85), the difference threshold (6 to 28), the closing kernel (5, 9, 13), the
fragment merge distance (0 to 70 px) and the minimum blob area (20 to 250 px),
minimising the deviation of the per-frame count from the number of fish in the
tank. Because that criterion uses the fish count itself, the detector was also
checked against independent manual counts, as described above.

Light phase: taken from the hour each recording started (06:00 and 21:00 dark,
14:00 light). A recording is split when the mean frame brightness changes
abruptly inside it, measured as the median of a 20-second window after the step
against the 20 seconds before, by at least 8 grey levels, in the direction the
hour allows, with a minimum segment of 60 s. Six of 78 recordings were split
this way. A simple brightness cutoff was not usable: camera auto-exposure raises
the picture back to daytime levels within minutes of lights-off.

## Design and limitations

- Condition is confounded with camera and rack position: camera 003 recorded
  sleep-restricted tanks and camera 004 control tanks, in both cohorts.
  Differences between conditions in absolute activity cannot be separated from
  differences between recording positions. Comparisons within a tank, such as
  the dark/light ratio, are not affected.
- The experimental unit is the tank, not the fish: four tanks per condition, two
  per cohort, six fish per tank (five in three tanks after the fish lost during
  the second cohort). Fish within a tank are not independent replicates.
- Power is low. For the dark-phase difference observed here, 80% power would
  require about ten tanks per condition; for the dark/light ratio, about twenty.
  These analyses are exploratory.
- The light and dark phases use different motion thresholds, so dark/light
  ratios are partly method-dependent and should be read qualitatively.
  Comparisons between conditions within one time of day are unaffected, since
  both conditions are processed identically.
- The second cohort was recorded with stronger compression (about 1.1 against
  2.2 Mbit/s) and a variable frame rate. Cohort is a fixed effect in every model.
- One death date (second cohort, camera 004, tank 1) is not in the records.
  `analysis.R` repeats the main tests without the affected tank-days and writes
  the comparison to `sensitivity.csv`.
- No individual fish tracking and no sleep measures.

## Suggested changes for future recordings

1. Swap cameras between conditions across cohorts, or alternate tanks between
   cameras, so that condition is not confounded with position.
2. Fix camera exposure and gain, and record at a constant frame rate and bitrate.
3. Record baseline days before treatment, so each tank serves as its own control.
4. Use one motion threshold for both light phases, or calibrate the noise floor
   separately for each.

## Citation and licence

If you use this code, please cite [paper reference]. Released under [licence].
The raw recordings are [available at / available on request] because of their
size.
