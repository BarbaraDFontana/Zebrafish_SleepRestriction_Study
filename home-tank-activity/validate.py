#!/usr/bin/env python3
"""Validation tools for the home-tank recordings.

Three modes.

sheet
    Draw random frames per cohort, camera and light phase, each saved as one
    image with the clean frame on top and the detections below, together with a
    CSV to fill in. An observer counts the fish in each tank on the top image
    and types the number into the manual_count column; the detected count
    appears only in the CSV, so the counting is blind.
        python validate.py sheet --folder "D:/Sleep/mp4"
    Output: validation/frames/*.png and validation/manual_counts.csv

score
    Compare the detector with the manual counts: mean absolute deviation, bias,
    and the share of frames that match exactly or within one fish, per cohort,
    camera and light phase.
        python validate.py score --folder "D:/Sleep/mp4"
    Output: validation/validation_score.csv

activity
    Checks on the pipeline output: whether the detected fish count sits near the
    number of fish, whether activity is lower in the dark than in the light in
    every cohort, whether control tanks differ between cohorts, and whether the
    light-phase labels match the recording hours.
        python validate.py activity --folder "D:/Sleep/mp4"
    Output: validation/activity_validation.csv and activity_validation.png
"""

import argparse
import csv
import glob
import os
import random
from concurrent.futures import ProcessPoolExecutor

import cv2
import numpy as np

import sleep_core as sc

SHEET_SAMPLE_SEC = 10.0
RESERVOIR_EVERY_SEC = 30.0
ACTIVITY_COLUMNS = ["mean_fish", "pct_frames_6", "mean_motion", "mean_activity"]


# --- sheet -------------------------------------------------------------------


def sheet_clip(job):
    """Draw the sampled frames of one recording and return their CSV rows."""
    clip, batch, camera, light, n_frames, outdir, seed, sample_sec = job
    rng = random.Random(seed)
    tuned = sc.load_tuned()
    rois = sc.rois_for(batch, camera)

    brights, frames = [], []
    for frame in sc.iter_frames(clip, 1.0 / sample_sec):
        brights.append(float(frame.mean()))
        frames.append(frame.copy())
    if len(frames) < 10:
        return []

    def percentile_for(is_night):
        phase = "night" if is_night else "day"
        return sc.params_for(tuned, batch, camera, phase)["bg_pct"]

    # The background is built from one frame every 30 s, as in run_pipeline.py,
    # so that the detections shown here are the ones the pipeline would make.
    every = max(1, int(round(RESERVOIR_EVERY_SEC / sample_sec)))
    reservoir = [(i, f) for i, f in enumerate(frames) if i % every == 0]
    segments, backgrounds = sc.make_segments(
        brights, reservoir, sample_sec, pct_fn=percentile_for,
        cutoff=sc.cutoff_for(batch, camera), hour=sc.parse_name(clip)[2][:2])

    want_night = light == "night"
    candidates = [i for start, end, is_night in segments if is_night == want_night
                  for i in range(start, end)]
    if not candidates:
        return []

    params = sc.params_for(tuned, batch, camera, light)
    y_low = max(0, min(r[1] for r in rois) - 20)
    y_high = min(frames[0].shape[0], max(r[3] for r in rois) + 20)

    rows = []
    for frame_index in sorted(rng.sample(candidates, min(n_frames, len(candidates)))):
        segment_index = next(k for k, (start, end, _) in enumerate(segments)
                             if start <= frame_index < end)
        frame = frames[frame_index]
        clean = cv2.cvtColor(frame, cv2.COLOR_GRAY2BGR)
        annotated = clean.copy()

        counts = []
        for x0, y0, x1, y1 in rois:
            bg_prepped, bg_mean = sc.prep_bg(backgrounds[segment_index], (x0, y0, x1, y1))
            _, fish = sc.detect(frame[y0:y1, x0:x1], bg_prepped, bg_mean, params)
            counts.append(sum(f[3] for f in fish))
            cv2.rectangle(annotated, (x0, y0), (x1, y1), (0, 255, 0), 1)
            for _, cx, cy, n_fish_blob in fish:
                centre = (int(cx + x0), int(cy + y0))
                cv2.circle(annotated, centre, 14, (0, 0, 255), 2)
                if n_fish_blob > 1:
                    cv2.putText(annotated, f"x{n_fish_blob}",
                                (centre[0] + 12, centre[1] - 10),
                                cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 255, 255), 2)

        name = os.path.splitext(os.path.basename(clip))[0].replace(" ", "_")
        image_id = f"{batch}_{camera}_{light}_{name}_f{frame_index}"
        top = clean[y_low:y_high]
        for tank_index, (x0, y0, _, _) in enumerate(rois):
            cv2.putText(top, f"tank{tank_index + 1}", (x0 + 5, y0 - y_low - 4),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 255, 0), 2)
        separator = np.full((6, top.shape[1], 3), 255, np.uint8)
        cv2.imwrite(os.path.join(outdir, image_id + ".png"),
                    np.vstack([top, separator, annotated[y_low:y_high]]))

        for tank_index in (0, 1):
            rows.append(dict(image=image_id + ".png", batch=batch, camera=camera,
                             light=light, clip=os.path.basename(clip),
                             time_s=int(frame_index * sample_sec),
                             tank=f"tank{tank_index + 1}",
                             detected=counts[tank_index], manual_count=""))
    return rows


def do_sheet(args):
    validation_dir = os.path.join(args.folder, "validation")
    frames_dir = os.path.join(validation_dir, "frames")
    os.makedirs(frames_dir, exist_ok=True)

    clips = sorted(glob.glob(os.path.join(args.folder, "*.mp4")))
    rng = random.Random(args.seed)
    jobs = []
    for batch in sc.BATCHES:
        for camera in sc.GROUPS:
            mine = [c for c in clips
                    if sc.batch_of(c) == batch and sc.parse_name(c)[0] == camera]
            for light, times in (("day", ("140000",)), ("night", ("210000", "060000"))):
                pool = [c for c in mine if sc.parse_name(c)[2] in times]
                for clip in rng.sample(pool, min(args.clips_per_group, len(pool))):
                    per_clip = int(np.ceil(args.frames_per_group / args.clips_per_group))
                    jobs.append((clip, batch, camera, light, per_clip, frames_dir,
                                 rng.randint(0, 10 ** 6), SHEET_SAMPLE_SEC))

    print(f"Building the counting sheet from {len(jobs)} recordings")
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(sheet_clip, jobs):
            rows.extend(result)
    if not rows:
        print("No frames could be sampled.")
        return

    out_path = os.path.join(validation_dir, "manual_counts.csv")
    if os.path.exists(out_path):
        out_path = out_path.replace(".csv", "_new.csv")
        print("manual_counts.csv exists already, writing manual_counts_new.csv")
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows) // 2} frames written to {frames_dir}")
    print(f"Fill in the manual_count column of {out_path}")


# --- score -------------------------------------------------------------------


def do_score(args):
    validation_dir = os.path.join(args.folder, "validation")
    with open(os.path.join(validation_dir, "manual_counts.csv")) as handle:
        rows = [r for r in csv.DictReader(handle) if r["manual_count"].strip() != ""]
    if not rows:
        print("No manual counts have been filled in yet.")
        return

    groups = {}
    for row in rows:
        key = (row["batch"], row["camera"], row["light"])
        groups.setdefault(key, []).append(
            (int(row["detected"]), int(float(row["manual_count"]))))
    groups[("all", "", "")] = [pair for values in list(groups.values()) for pair in values]

    summary = []
    for key, pairs in sorted(groups.items()):
        detected = np.array([p[0] for p in pairs])
        manual = np.array([p[1] for p in pairs])
        if len(pairs) > 2 and detected.std() and manual.std():
            correlation = float(np.corrcoef(detected, manual)[0, 1])
        else:
            correlation = float("nan")
        summary.append(dict(
            batch=key[0], camera=key[1], light=key[2], n=len(pairs),
            mae=round(float(np.mean(np.abs(detected - manual))), 3),
            bias_detected_minus_manual=round(float(np.mean(detected - manual)), 3),
            pct_exact=round(100 * float(np.mean(detected == manual)), 1),
            pct_within1=round(100 * float(np.mean(np.abs(detected - manual) <= 1)), 1),
            manual_mean=round(float(manual.mean()), 2), r=round(correlation, 3)))
        last = summary[-1]
        print(f"{key[0]:4s} {key[1]:4s} {key[2]:6s} n={len(pairs):3d}  "
              f"MAE {last['mae']:.2f}  bias {last['bias_detected_minus_manual']:+.2f}  "
              f"exact {last['pct_exact']:.0f}%  within 1 {last['pct_within1']:.0f}%  "
              f"manual mean {last['manual_mean']}")

    out_path = os.path.join(validation_dir, "validation_score.csv")
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)
    print(f"Written to {out_path}")


# --- activity ----------------------------------------------------------------


def do_activity(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    validation_dir = os.path.join(args.folder, "validation")
    os.makedirs(validation_dir, exist_ok=True)
    source = os.path.join(args.folder, "results", "activity_binned.csv")
    with open(source) as handle:
        rows = list(csv.DictReader(handle))

    # Average the bins of each tank-recording so that every recording counts once.
    per_recording = {}
    for row in rows:
        key = (row["batch"], row["group"], row["day_night"], row["video"],
               row["tank"], row["time"])
        per_recording.setdefault(key, []).append(
            [float(row[column]) for column in ACTIVITY_COLUMNS])
    records = [dict(batch=key[0], group=key[1], light=key[2], time=key[5],
                    **dict(zip(ACTIVITY_COLUMNS, np.mean(values, axis=0))))
               for key, values in per_recording.items()]

    summary = []
    for key in sorted({(r["batch"], r["group"], r["light"]) for r in records}):
        subset = [r for r in records if (r["batch"], r["group"], r["light"]) == key]
        entry = dict(batch=key[0], group=key[1], light=key[2],
                     n_tank_recordings=len(subset))
        for column in ACTIVITY_COLUMNS:
            values = np.array([s[column] for s in subset])
            entry[column + "_mean"] = round(float(values.mean()), 5)
            entry[column + "_sd"] = round(
                float(values.std(ddof=1)) if len(values) > 1 else 0.0, 5)
        flags = []
        if abs(entry["mean_fish_mean"] - sc.EXPECTED_FISH) > 0.5:
            flags.append("fish_count_off_target")
        entry["flags"] = ";".join(flags)
        summary.append(entry)

    print("\ncohort group             phase  n   fish   % at 6   activity index")
    for entry in summary:
        print(f"{entry['batch']:6s} {entry['group']:17s} {entry['light']:5s} "
              f"{entry['n_tank_recordings']:3d}  {entry['mean_fish_mean']:.2f}  "
              f"{entry['pct_frames_6_mean']:6.0f}   {entry['mean_motion_mean']:.5f}  "
              f"{entry['flags']}")

    lookup = {(r["batch"], r["group"], r["light"]): r for r in summary}

    print("\nDark/light activity ratio (a working pipeline gives a ratio below 1 "
          "in control tanks):")
    for batch in sorted({r["batch"] for r in summary}):
        for group in sorted({r["group"] for r in summary}):
            day = lookup.get((batch, group, "day"))
            night = lookup.get((batch, group, "night"))
            if day and night and day["mean_motion_mean"]:
                ratio = night["mean_motion_mean"] / day["mean_motion_mean"]
                note = "   check" if ratio >= 1 else ""
                print(f"  {batch} {group:17s} {ratio:.2f}{note}")

    print("\nControl tanks, second cohort against the first (a large difference "
          "means a cohort effect to keep in the model):")
    for light in ("day", "night"):
        first = lookup.get(("b1", "Control", light))
        second = lookup.get(("b2", "Control", light))
        if first and second:
            for column in ("mean_motion", "mean_activity", "mean_fish"):
                a = first[column + "_mean"]
                b = second[column + "_mean"]
                change = 100 * (b - a) / a if a else float("nan")
                note = "   check" if abs(change) > 25 else ""
                print(f"  {light:5s} {column:15s} b1 {a:.4f}  b2 {b:.4f}  "
                      f"({change:+.0f}%){note}")

    qc_path = os.path.join(args.folder, "results", "qc_summary.csv")
    if os.path.exists(qc_path):
        print("\nLight phase assigned per recording hour (14h light, 06h and 21h "
              "dark; 'mixed' means a lights switch was found inside the recording):")
        with open(qc_path) as handle:
            qc_rows = list(csv.DictReader(handle))
        table = {}
        for row in qc_rows:
            key = (row["batch"], row["camera"], row["time"][:2] + "h")
            table.setdefault(key, {}).setdefault(row["day_night"], []).append(row["date"])
        for key in sorted(table):
            counts = ", ".join(f"{phase} {len(dates)}"
                               for phase, dates in sorted(table[key].items()))
            unexpected = [f"{phase}: {' '.join(dates)}"
                          for phase, dates in table[key].items()
                          if (key[2] == "14h" and phase != "day")
                          or (key[2] == "06h" and phase != "night")]
            note = f"   check {unexpected}" if unexpected else ""
            print(f"  {key[0]} {key[1]} {key[2]}  {counts}{note}")

    out_csv = os.path.join(validation_dir, "activity_validation.csv")
    with open(out_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)

    metrics = ["mean_fish", "mean_motion", "mean_activity"]
    colours = {"Control": "#4DAF4A", "SleepDeprivation": "#984EA3"}
    figure, axes = plt.subplots(len(metrics), 2, figsize=(10, 3.2 * len(metrics)))
    for i, column in enumerate(metrics):
        for j, light in enumerate(("day", "night")):
            axis = axes[i, j]
            labels, data, box_colours = [], [], []
            for batch in sorted({r["batch"] for r in records}):
                for group in ("Control", "SleepDeprivation"):
                    values = [r[column] for r in records
                              if r["batch"] == batch and r["group"] == group
                              and r["light"] == light]
                    if values:
                        labels.append(f"{batch}\n{'Ctrl' if group == 'Control' else 'SR'}")
                        data.append(values)
                        box_colours.append(colours[group])
            if data:
                boxes = axis.boxplot(data, patch_artist=True, widths=0.6)
                axis.set_xticks(range(1, len(labels) + 1))
                axis.set_xticklabels(labels)
                for patch, colour in zip(boxes["boxes"], box_colours):
                    patch.set_facecolor(colour)
                    patch.set_alpha(0.5)
                for position, values in enumerate(data):
                    axis.scatter(np.random.normal(position + 1, 0.05, len(values)),
                                 values, s=8, c="k", alpha=0.6)
            if column == "mean_fish":
                axis.axhline(sc.EXPECTED_FISH, ls="--", c="grey")
            axis.set_title(f"{column} - {light}", fontsize=10)
            axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    out_png = os.path.join(validation_dir, "activity_validation.png")
    figure.savefig(out_png, dpi=130)
    print(f"\nWritten to {out_csv} and {out_png}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validation tools for the home-tank recordings.")
    parser.add_argument("mode", choices=["sheet", "score", "activity"],
                        help="sheet: build frames to count; score: detector against "
                             "the manual counts; activity: checks on the output")
    parser.add_argument("--folder", default=".",
                        help="folder holding the recordings and the results folder")
    parser.add_argument("--frames-per-group", type=int, default=12,
                        help="frames per cohort, camera and light phase (default: 12)")
    parser.add_argument("--clips-per-group", type=int, default=2,
                        help="recordings sampled per cohort, camera and light phase")
    parser.add_argument("--seed", type=int, default=42,
                        help="random seed for the frame selection (default: 42)")
    parser.add_argument("--workers", type=int,
                        default=max(1, min(4, (os.cpu_count() or 2) // 2)),
                        help="worker processes to run in parallel")
    return parser.parse_args()


def main():
    args = parse_args()
    {"sheet": do_sheet, "score": do_score, "activity": do_activity}[args.mode](args)


if __name__ == "__main__":
    main()
