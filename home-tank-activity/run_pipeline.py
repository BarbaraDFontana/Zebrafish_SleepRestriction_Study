#!/usr/bin/env python3
"""Measure home-tank activity in every recording of a folder.

Each recording is read twice. The first pass collects the mean brightness of
every processed frame and a reservoir of frames spread across the recording;
these give the light-phase segments and one background image per segment. The
second pass measures, for each tank region and each frame, the fraction of
pixels that changed since the previous frame (the activity index analysed in
the paper) and the number of fish the detector finds (quality control only).
Values are averaged within 10-minute bins.

Inputs
    Recordings named CAMERA R DDMMYYYY HHMMSS.mp4 in the folder given with
    --folder; sleep_core.py for the layout and detection settings;
    tuned_params.json and fish_alive.csv next to sleep_core.py.

Outputs, written to <folder>/results/
    activity_binned.csv  one row per recording, tank and 10-minute bin
    qc_summary.csv       one row per recording

Usage
    python run_pipeline.py --folder "D:/Sleep/mp4"
    python run_pipeline.py --folder "D:/Sleep/mp4" --batch b2 --workers 4

The run is resumable: recordings already listed in qc_summary.csv are skipped.
Delete the results folder for a fresh run, for example after re-selecting
parameters.
"""

import argparse
import csv
import glob
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

import sleep_core as sc

BIN_SEC_DEFAULT = 600.0
RESERVOIR_EVERY_SEC = 30.0

BIN_HEADER = ["video", "batch", "camera", "group", "date", "time", "day_night",
              "tank", "time_bin", "mean_fish", "pct_frames_6", "mean_blobs",
              "mean_activity", "mean_motion", "n_fish", "n_frames"]
QC_HEADER = ["clip", "batch", "camera", "group", "date", "time", "day_night",
             "n_light_segments", "window_min", "proc_fps", "n_frames",
             "bg_brightness"] + \
            [f"tank{t}_{k}" for t in (1, 2)
             for k in ("fish_median", "fish_mean", "pct6", "pct_within1")]


def process_clip(clip, tuned, window_min, bin_sec):
    """Analyse one recording. Returns (qc row, binned rows) or None if skipped."""
    batch = sc.batch_of(clip)
    camera, day, clock = sc.parse_name(clip)
    rois = sc.rois_for(batch, camera)
    if batch is None or rois is None:
        return None

    max_sec = int(window_min * 60) if window_min else 0
    dt = 1.0 / sc.PROC_FPS
    reservoir_every = max(1, int(RESERVOIR_EVERY_SEC / dt))

    # Pass one: brightness trace and background reservoir.
    brights, reservoir = [], []
    for index, frame in enumerate(sc.iter_frames(clip, sc.PROC_FPS, max_sec)):
        brights.append(float(frame.mean()))
        if index % reservoir_every == 0:
            reservoir.append((index, frame.copy()))
    if len(brights) < 5:
        return None

    def percentile_for(is_night):
        light = "night" if is_night else "day"
        return sc.params_for(tuned, batch, camera, light)["bg_pct"]

    segments, backgrounds = sc.make_segments(
        brights, reservoir, dt, pct_fn=percentile_for,
        cutoff=sc.cutoff_for(batch, camera), hour=clock[:2])
    del reservoir

    segment_of_frame = np.zeros(len(brights), int)
    for index, (start, end, _) in enumerate(segments):
        segment_of_frame[start:end] = index

    segment_config = []
    for index, (_, _, is_night) in enumerate(segments):
        light = "night" if is_night else "day"
        params = sc.params_for(tuned, batch, camera, light)
        segment_config.append(dict(
            params=params, light=light,
            tanks=[sc.prep_bg(backgrounds[index], box) for box in rois]))

    # Pass two: detection and frame-to-frame pixel change.
    bins = [dict(), dict()]
    fish_counts = [[], []]
    previous_frame = [None, None]
    previous_segment = [None, None]
    n_frames = 0

    for frame_index, frame in enumerate(sc.iter_frames(clip, sc.PROC_FPS, max_sec)):
        if frame_index >= len(segment_of_frame):
            break
        n_frames += 1
        seconds = frame_index * dt
        bin_index = int(seconds // bin_sec)
        segment_index = segment_of_frame[frame_index]
        config = segment_config[segment_index]
        params = config["params"]

        for tank_index, (x0, y0, x1, y1) in enumerate(rois):
            roi = frame[y0:y1, x0:x1]
            bg_prepped, bg_mean = config["tanks"][tank_index]
            occupancy, fish = sc.detect(roi, bg_prepped, bg_mean, params)
            n_fish_frame = sum(f[3] for f in fish)
            fish_counts[tank_index].append(n_fish_frame)

            prepared = sc.prep(roi, bg_mean)
            motion = None
            # Only compare frames inside the same light segment: the background
            # and the motion threshold change at a lights switch.
            if (previous_frame[tank_index] is not None
                    and previous_segment[tank_index] == segment_index):
                motion = sc.motion_fraction(prepared, previous_frame[tank_index],
                                            params["motion_thr"])
            previous_frame[tank_index] = prepared
            previous_segment[tank_index] = segment_index

            totals = bins[tank_index].setdefault(
                bin_index, dict(occupancy=0.0, blobs=0, fish=0, exactly_six=0,
                                frames=0, night=0, motion=0.0, motion_frames=0))
            totals["occupancy"] += occupancy
            totals["blobs"] += len(fish)
            totals["fish"] += n_fish_frame
            totals["exactly_six"] += (n_fish_frame == sc.EXPECTED_FISH)
            totals["frames"] += 1
            totals["night"] += config["light"] == "night"
            if motion is not None:
                totals["motion"] += motion
                totals["motion_frames"] += 1

    group = sc.GROUPS.get(camera, "")
    name = os.path.basename(clip)
    alive = sc.load_fish_alive()

    rows = []
    for tank_index in (0, 1):
        n_alive = sc.n_fish_for(alive, batch, camera, f"tank{tank_index + 1}", day)
        for bin_index in sorted(bins[tank_index]):
            totals = bins[tank_index][bin_index]
            frames = max(totals["frames"], 1)
            rows.append(dict(
                video=name, batch=batch, camera=camera, group=group,
                date=day, time=clock,
                day_night="night" if totals["night"] > totals["frames"] / 2 else "day",
                tank=f"tank{tank_index + 1}", time_bin=bin_index,
                mean_fish=round(totals["fish"] / frames, 3),
                pct_frames_6=round(100.0 * totals["exactly_six"] / frames, 1),
                mean_blobs=round(totals["blobs"] / frames, 3),
                mean_activity=round(totals["occupancy"] / frames, 6),
                mean_motion=round(totals["motion"] / max(totals["motion_frames"], 1), 6),
                n_fish=n_alive,
                n_frames=totals["frames"]))

    lights = {config["light"] for config in segment_config}
    qc = dict(clip=name, batch=batch, camera=camera, group=group, date=day, time=clock,
              day_night=(next(iter(lights)) if len(segments) == 1 else "mixed"),
              n_light_segments=len(segments), window_min=window_min,
              proc_fps=sc.PROC_FPS, n_frames=n_frames,
              bg_brightness=round(float(np.mean(brights)), 2))
    for tank_index in (0, 1):
        counts = np.array(fish_counts[tank_index]) if fish_counts[tank_index] else np.zeros(1)
        tank = f"tank{tank_index + 1}"
        qc[f"{tank}_fish_median"] = int(np.median(counts))
        qc[f"{tank}_fish_mean"] = round(float(counts.mean()), 3)
        qc[f"{tank}_pct6"] = round(100.0 * float(np.mean(counts == sc.EXPECTED_FISH)), 1)
        qc[f"{tank}_pct_within1"] = round(
            100.0 * float(np.mean(np.abs(counts - sc.EXPECTED_FISH) <= 1)), 1)
    return qc, rows


def append_rows(path, rows, header):
    """Append rows to a CSV, writing the header when the file is new."""
    is_new = not os.path.exists(path)
    with open(path, "a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        if is_new:
            writer.writeheader()
        writer.writerows(rows)


def _worker(job):
    """Run one recording in a worker process, reporting errors instead of raising."""
    clip, tuned, window_min, bin_sec = job
    try:
        return clip, process_clip(clip, tuned, window_min, bin_sec), None
    except Exception as error:
        return clip, None, repr(error)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Measure home-tank activity in every recording of a folder.")
    parser.add_argument("--folder", default=".",
                        help="folder holding the .mp4 recordings (default: current)")
    parser.add_argument("--batch", default=None,
                        help="process only this cohort, as named in sleep_core.BATCHES")
    parser.add_argument("--clip", default=None,
                        help="process only this recording, given as a file name")
    parser.add_argument("--window-min", type=float, default=0.0,
                        help="analyse only the first N minutes of each recording "
                             "(0 = whole recording, the default)")
    parser.add_argument("--bin-min", type=float, default=BIN_SEC_DEFAULT / 60,
                        help="length of the output bins in minutes (default: 10)")
    parser.add_argument("--workers", type=int,
                        default=max(1, min(4, (os.cpu_count() or 2) // 2)),
                        help="worker processes to run in parallel")
    return parser.parse_args()


def main():
    args = parse_args()

    outdir = os.path.join(args.folder, "results")
    os.makedirs(outdir, exist_ok=True)
    qc_csv = os.path.join(outdir, "qc_summary.csv")
    bin_csv = os.path.join(outdir, "activity_binned.csv")

    tuned = sc.load_tuned()
    if not tuned:
        print("tuned_params.json not found: falling back to the default parameters")

    done = set()
    if os.path.exists(qc_csv):
        with open(qc_csv) as handle:
            done = {row["clip"] for row in csv.DictReader(handle)}

    clips = sorted(glob.glob(os.path.join(args.folder, "*.mp4")))
    clips = [c for c in clips if os.path.basename(c) not in done]
    if args.batch:
        clips = [c for c in clips if sc.batch_of(c) == args.batch]
    if args.clip:
        clips = [c for c in clips if os.path.basename(c) == args.clip]

    undated = [c for c in clips if sc.batch_of(c) is None]
    for clip in undated:
        print(f"skipped {os.path.basename(clip)}: date outside every cohort range")
    clips = [c for c in clips if c not in undated]

    print(f"{len(clips)} recordings to process with {args.workers} workers")

    started = time.time()
    processed = 0
    jobs = [(c, tuned, args.window_min, args.bin_min * 60) for c in clips]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(_worker, job) for job in jobs]
        for future in as_completed(futures):
            clip, result, error = future.result()
            name = os.path.basename(clip)
            if error:
                print(f"failed  {name}: {error}")
                continue
            if result is None:
                print(f"skipped {name}: no tank regions configured, or unreadable")
                continue
            qc, rows = result
            append_rows(qc_csv, [qc], QC_HEADER)
            append_rows(bin_csv, rows, BIN_HEADER)
            processed += 1
            print(f"done    {name} [{qc['batch']} {qc['group']}] {qc['day_night']}  "
                  f"fish tank1 {qc['tank1_fish_mean']:.2f} ({qc['tank1_pct6']:.0f}% at 6)  "
                  f"tank2 {qc['tank2_fish_mean']:.2f} ({qc['tank2_pct6']:.0f}% at 6)  "
                  f"[{processed}/{len(clips)}, {time.time() - started:.0f} s]")

    print(f"{processed} recordings processed, written to {outdir}")


if __name__ == "__main__":
    main()
