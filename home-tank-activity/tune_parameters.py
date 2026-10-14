#!/usr/bin/env python3
"""Select detection parameters per cohort, camera and light phase.

The detector is used for quality control: it reports how many fish are visible
per frame. Its parameters are chosen so that this count sits as close as
possible to the number of fish in the tank.

Procedure
    A few recordings are taken per cohort, camera and light phase, spread over
    the recording days. From each one, a frame is sampled every --sample-sec
    seconds across the whole hour, the per-segment background is built exactly
    as the main pipeline builds it, and --grid-frames frames are kept for
    scoring. Every combination of the grids below is then scored as the mean
    absolute difference between the detected count and the number of fish, over
    all scored tank-frames; the lowest score wins.

Grids
    bg_pct      background percentile (50 is the median; 85 keeps a dark fish
                resting in one place out of the background, 15 does the same for
                fish that appear bright under infrared light)
    thr         difference threshold against the background
    close_k     closing kernel, which fills gaps inside one fish
    merge_dist  distance below which blob pieces are treated as one fish
    min_area    smallest blob accepted as a fish
    area_count  whether an oversized blob counts as several fish

Outputs
    tuned_params.json           the winning parameters, read by run_pipeline.py
                                and validate.py. Running this script overwrites
                                that file, including the values published with
                                the paper; keep a copy before re-running.
    <folder>/tuning/tuning_grid.csv   every combination and its score
    <folder>/tuning/tuning_best.csv   winners, next to the default parameters
    <folder>/tuning/tuning_clips.csv  the recordings sampled, with brightness

Usage
    python tune_parameters.py --folder "D:/Sleep/mp4"
    python tune_parameters.py --folder "D:/Sleep/mp4" --batch b2

A low score can also arise from errors that cancel out, so the result is
confirmed against blinded manual counts with validate.py.
"""

import argparse
import csv
import glob
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import sleep_core as sc

BGPCT_GRID = [15, 50, 85]
THR_GRID = list(range(6, 30, 2))
CLOSE_GRID = [5, 9, 13]
MERGE_GRID = [0, 15, 25, 40, 55, 70]
MINAREA_GRID = [20, 30, 45, 60, 80, 100, 130, 160, 200, 250]
AREACOUNT_GRID = [0, 1]


def pick_clips(clips, n):
    """Up to n recordings, spread evenly across the list."""
    if len(clips) <= n:
        return clips
    index = np.linspace(0, len(clips) - 1, n).round().astype(int)
    return [clips[i] for i in index]


def sample_clip(job):
    """Sample one recording and pre-compute blob pieces for the grid.

    Returns (clip, {(bg_pct, thr, close_k): [(light, tank, pieces), ...]}, info).
    Thresholding and morphology are done once per (bg_pct, thr, close_k); the
    merge distance, minimum area and area counting are applied later, during
    scoring, because they work on the pieces alone.
    """
    clip, batch, camera, sample_sec, grid_frames = job
    rois = sc.rois_for(batch, camera)
    y_low = min(r[1] for r in rois)
    y_high = max(r[3] for r in rois)
    rois_cropped = [(x0, y0 - y_low, x1, y1 - y_low) for x0, y0, x1, y1 in rois]

    brights, frames = [], []
    for frame in sc.iter_frames(clip, 1.0 / sample_sec):
        brights.append(float(frame.mean()))
        frames.append(frame[y_low:y_high].copy())
    if len(frames) < 10:
        return clip, None, None

    dt = float(sample_sec)
    segments, backgrounds = sc.make_segments(
        brights, list(enumerate(frames)), dt, pct_list=BGPCT_GRID,
        cutoff=sc.cutoff_for(batch, camera), hour=sc.parse_name(clip)[2][:2])
    segment_of_frame = np.zeros(len(frames), int)
    for index, (start, end, _) in enumerate(segments):
        segment_of_frame[start:end] = index

    background_cache = {}
    keep = np.linspace(0, len(frames) - 1, min(grid_frames, len(frames))).astype(int)
    pieces_by_setting = {(q, t, k): []
                         for q in BGPCT_GRID for t in THR_GRID for k in CLOSE_GRID}

    for frame_index in keep:
        segment_index = segment_of_frame[frame_index]
        light = "night" if segments[segment_index][2] else "day"
        for tank_index, box in enumerate(rois_cropped):
            x0, y0, x1, y1 = box
            for percentile in BGPCT_GRID:
                cache_key = (segment_index, tank_index, percentile)
                if cache_key not in background_cache:
                    background_cache[cache_key] = sc.prep_bg(
                        backgrounds[segment_index][percentile], box)
                bg_prepped, bg_mean = background_cache[cache_key]
                diff = sc.diff_image(frames[frame_index][y0:y1, x0:x1], bg_prepped, bg_mean)
                for threshold in THR_GRID:
                    for close_k in CLOSE_GRID:
                        _, pieces = sc.raw_pieces(diff, threshold, close_k)
                        pieces_by_setting[(percentile, threshold, close_k)].append(
                            (light, tank_index, pieces))

    night_fraction = sum(end - start for start, end, is_night in segments if is_night) / len(frames)
    info = dict(clip=os.path.basename(clip), n_frames=len(frames),
                mean_brightness=round(float(np.mean(brights)), 1),
                min_brightness=round(float(np.min(brights)), 1),
                max_brightness=round(float(np.max(brights)), 1),
                cutoff=sc.cutoff_for(batch, camera),
                pct_night=round(100 * night_fraction),
                n_segments=len(segments))
    return clip, pieces_by_setting, info


def score_group(frame_lists, bg_pct, thr, close_k, merge_dist, min_area, area_count):
    """Score one parameter combination over a list of per-frame blob lists.

    frame_lists holds pieces already merged at merge_dist. When area_count is
    set, the area of a single fish is taken as the median blob area of this
    combination, and a blob above MULTI_FISH_RATIO times that area counts as
    several fish.
    """
    kept = [[b for b in frame if min_area <= b[0] <= sc.MAX_BLOB_AREA]
            for frame in frame_lists]

    single_area = None
    if area_count:
        areas = [b[0] for frame in kept for b in frame]
        single_area = float(np.median(areas)) if areas else None

    counts = []
    for frame in kept:
        total = 0
        for area, _, _ in frame:
            if single_area and area > sc.MULTI_FISH_RATIO * single_area:
                total += max(1, int(round(area / single_area)))
            else:
                total += 1
        counts.append(total)

    counted = np.array(counts, float)
    expected = sc.EXPECTED_FISH
    return dict(bg_pct=bg_pct, thr=thr, close_k=close_k, merge_dist=merge_dist,
                min_area=min_area, area_count=area_count,
                single_area=round(single_area, 1) if single_area else "",
                mae=round(float(np.mean(np.abs(counted - expected))), 4),
                bias=round(float(np.mean(counted - expected)), 3),
                pct_exact6=round(float(np.mean(counted == expected)) * 100, 1),
                pct_within1=round(float(np.mean(np.abs(counted - expected) <= 1)) * 100, 1),
                mean_count=round(float(counted.mean()), 3),
                sd_count=round(float(counted.std()), 3),
                n_tank_frames=len(counts))


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Select detection parameters per cohort, camera and light phase.")
    parser.add_argument("--folder", default=".",
                        help="folder holding the .mp4 recordings (default: current)")
    parser.add_argument("--batch", default=None,
                        help="tune only this cohort, as named in sleep_core.BATCHES")
    parser.add_argument("--clips-per-light", type=int, default=3,
                        help="recordings sampled per light phase (default: 3)")
    parser.add_argument("--sample-sec", type=float, default=6.0,
                        help="seconds between sampled frames (default: 6)")
    parser.add_argument("--grid-frames", type=int, default=200,
                        help="frames per recording kept for scoring (default: 200)")
    parser.add_argument("--workers", type=int,
                        default=max(1, min(4, (os.cpu_count() or 2) // 2)),
                        help="worker processes to run in parallel")
    return parser.parse_args()


def main():
    args = parse_args()

    outdir = os.path.join(args.folder, "tuning")
    os.makedirs(outdir, exist_ok=True)
    clips = sorted(glob.glob(os.path.join(args.folder, "*.mp4")))

    jobs = []
    for batch in sc.BATCHES:
        if args.batch and batch != args.batch:
            continue
        for camera in sc.GROUPS:
            mine = [c for c in clips
                    if sc.batch_of(c) == batch and sc.parse_name(c)[0] == camera]
            day = [c for c in mine if sc.parse_name(c)[2] == "140000"]
            night = [c for c in mine if sc.parse_name(c)[2] in ("210000", "060000")]
            for clip in (pick_clips(day, args.clips_per_light)
                         + pick_clips(night, args.clips_per_light)):
                jobs.append((clip, batch, camera, args.sample_sec, args.grid_frames))
    print(f"Sampling {len(jobs)} recordings with {args.workers} workers")

    collected = {}
    clip_rows = []
    started = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for job, (clip, sampled, info) in zip(jobs, pool.map(sample_clip, jobs)):
            batch, camera = job[1], job[2]
            if sampled is None:
                print(f"  unreadable: {os.path.basename(clip)}")
                continue
            for setting, entries in sampled.items():
                for light, _, pieces in entries:
                    collected.setdefault((batch, camera, light), {}) \
                             .setdefault(setting, []).append(pieces)
            print(f"  sampled {info['clip']}  brightness {info['mean_brightness']} "
                  f"[{info['min_brightness']}-{info['max_brightness']}] "
                  f"cutoff {info['cutoff']}  dark {info['pct_night']}%  "
                  f"segments {info['n_segments']}  ({time.time() - started:.0f} s)")
            clip_rows.append(info)

    tuned = sc.load_tuned()
    grid_rows, best_rows = [], []
    for (batch, camera, light), by_setting in sorted(collected.items()):
        results = []
        for (percentile, threshold, close_k), frame_pieces in by_setting.items():
            for merge_dist in MERGE_GRID:
                merged = [sc.merge_pieces(p, merge_dist) for p in frame_pieces]
                for min_area in MINAREA_GRID:
                    for area_count in AREACOUNT_GRID:
                        results.append(score_group(merged, percentile, threshold,
                                                   close_k, merge_dist, min_area,
                                                   area_count))
        for row in results:
            grid_rows.append(dict(batch=batch, camera=camera, light=light, **row))

        # Ties are broken towards more frames with the exact count, then the
        # median background, then the higher threshold.
        results.sort(key=lambda r: (r["mae"], -r["pct_exact6"],
                                    abs(r["bg_pct"] - 50), -r["thr"]))
        best = results[0]

        defaults = sc.DEFAULTS[light]
        default_pieces = by_setting.get((50, defaults["thr"], defaults["close_k"]))
        default_score = (score_group(default_pieces, 50, defaults["thr"],
                                     defaults["close_k"], 0, defaults["min_area"], 0)
                         if default_pieces else None)
        if default_score:
            best_rows.append(dict(batch=batch, camera=camera, light=light,
                                  kind="defaults", **default_score))
        best_rows.append(dict(batch=batch, camera=camera, light=light,
                              kind="tuned", **best))

        tuned.setdefault(batch, {}).setdefault(camera, {})[light] = dict(
            bg_pct=best["bg_pct"], thr=best["thr"], close_k=best["close_k"],
            merge_dist=best["merge_dist"], min_area=best["min_area"],
            single_area=(best["single_area"] if best["area_count"] else None),
            motion_thr=sc.DEFAULTS[light]["motion_thr"])

        reference = (f"defaults MAE {default_score['mae']:.2f} "
                     f"({default_score['pct_exact6']:.0f}% exact)"
                     if default_score else "")
        print(f"{batch} {camera} {light:5s}: bg{best['bg_pct']} thr{best['thr']} "
              f"close{best['close_k']} merge{best['merge_dist']} "
              f"minA{best['min_area']} areacount{best['area_count']} -> "
              f"MAE {best['mae']:.2f}  {best['pct_exact6']:.0f}% exact  "
              f"{best['pct_within1']:.0f}% within 1  mean {best['mean_count']:.2f}"
              f"   | {reference}")

    write_csv(os.path.join(outdir, "tuning_grid.csv"), grid_rows)
    write_csv(os.path.join(outdir, "tuning_best.csv"), best_rows)
    write_csv(os.path.join(outdir, "tuning_clips.csv"), clip_rows)
    with open(sc.TUNED_JSON, "w") as handle:
        json.dump(tuned, handle, indent=2)
    print(f"\nWrote {sc.TUNED_JSON} and {os.path.join(outdir, 'tuning_*.csv')} "
          f"({time.time() - started:.0f} s)")


if __name__ == "__main__":
    main()
