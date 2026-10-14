#!/usr/bin/env python3
"""Shared configuration and image analysis for the home-tank recordings.

Every script in this repository imports this module, so detection behaves
identically in parameter selection, validation and the main run, and there is a
single place where settings live.

Contents
    Experiment layout: cohorts, camera to condition mapping, tank regions of
    interest, number of fish.
    Detection parameters: defaults, overridden per cohort, camera and light
    phase by tuned_params.json.
    Video access: frames are read from ffmpeg as raw grey images.
    Light-phase assignment and background estimation.
    Detection: background subtraction, blob merging, fish counting, and the
    frame-to-frame pixel-change index used as the activity measure.

Inputs read from the folder of this file
    tuned_params.json   detection parameters per cohort, camera and light phase
    fish_alive.csv      number of fish alive per tank from a given date
"""

import csv
import json
import os
import re
import subprocess
from datetime import date

import cv2
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TUNED_JSON = os.path.join(HERE, "tuned_params.json")
FISH_ALIVE_CSV = os.path.join(HERE, "fish_alive.csv")

# --- experiment layout -------------------------------------------------------

EXPECTED_FISH = 6

# Camera 003 recorded the sleep-restricted tanks and camera 004 the control
# tanks, in both cohorts. Condition is therefore confounded with camera and
# rack position; see the README.
GROUPS = {"003": "SleepDeprivation", "004": "Control"}

# Cohort name -> (first date, last date), inclusive. A recording is assigned to
# a cohort from the date in its filename.
BATCHES = {
    "b1": (date(2026, 4, 29), date(2026, 5, 6)),
    "b2": (date(2026, 8, 26), date(2026, 9, 1)),
}

# Tank boxes (x0, y0, x1, y1) per cohort and camera; tank1 is on the left.
# The boxes reach the tank floor, where resting fish sit at y of about 675-715.
# Camera 004 sits lower in the frame in the second cohort, hence separate boxes.
CAMERA_ROIS = {
    "b1": {"003": [(45, 350, 600, 710), (640, 350, 1238, 710)],
           "004": [(42, 300, 610, 690), (645, 300, 1240, 690)]},
    "b2": {"003": [(15, 350, 610, 715), (655, 350, 1240, 715)],
           "004": [(20, 305, 622, 715), (672, 305, 1250, 715)]},
}

# Mean frame brightness below which a frame would be called dark. Kept only as a
# fallback: light phase is assigned from the clock hour (see segment_by_clock),
# because camera auto-exposure makes day and night brightness overlap in the
# second cohort.
NIGHT_BRIGHTNESS_CUTOFF = 95.0
NIGHT_CUTOFF = {("b2", "003"): 95.5, ("b2", "004"): 82.0}

# --- detection parameters ----------------------------------------------------
# These are the fallback values; tuned_params.json overrides them per cohort,
# camera and light phase.

DEFAULTS = {
    "day": dict(thr=16, min_area=70, close_k=5, merge_dist=0, motion_thr=16, bg_pct=50),
    "night": dict(thr=12, min_area=60, close_k=5, merge_dist=0, motion_thr=12, bg_pct=50),
}

PROC_FPS = 2.0            # frames per second decoded and analysed
BG_SAMPLE_FRAMES = 60     # frames used to build one background image
MIN_SEGMENT_SEC = 60.0    # shortest light-phase segment kept inside a recording
ILLUM_NORMALIZE = True    # match frame mean to the background (auto-gain drift)
MEDIAN_BLUR_K = 3
MORPH_OPEN_K = 3
MAX_BLOB_AREA = 6000
SPECK_AREA = 10           # pieces smaller than this are dropped before merging
MULTI_FISH_RATIO = 1.6    # a blob above this times one fish counts as several

# --- filenames, cohorts, parameters ------------------------------------------


def parse_name(path):
    """Split a recording filename into (camera, DDMMYYYY, HHMMSS).

    Both separators used during the experiment are accepted:
    '003_R_01052026_060000.mp4' and '003 R 01092026 060000.mp4'.
    """
    parts = re.split(r"[ _]+", os.path.splitext(os.path.basename(path))[0])
    camera = parts[0] if len(parts) > 0 else ""
    day = parts[2] if len(parts) > 2 else ""
    clock = parts[3] if len(parts) > 3 else ""
    return camera, day, clock


def batch_of(path):
    """Cohort a recording belongs to, or None if its date is outside BATCHES."""
    _, day, _ = parse_name(path)
    try:
        recorded = _date_from_ddmmyyyy(day)
    except ValueError:
        return None
    for name, (first, last) in BATCHES.items():
        if first <= recorded <= last:
            return name
    return None


def load_tuned():
    """Parameters from tuned_params.json, or an empty dict if it is absent."""
    if os.path.exists(TUNED_JSON):
        with open(TUNED_JSON) as handle:
            return json.load(handle)
    return {}


def params_for(tuned, batch, camera, light):
    """Full parameter dict for one cohort, camera and light phase ('day'/'night')."""
    params = dict(DEFAULTS[light])
    params.update(tuned.get(batch, {}).get(camera, {}).get(light, {}))
    params.setdefault("single_area", None)
    return params


def cutoff_for(batch, camera):
    """Fallback brightness cutoff for this cohort and camera."""
    return NIGHT_CUTOFF.get((batch, camera), NIGHT_BRIGHTNESS_CUTOFF)


def rois_for(batch, camera):
    """Tank boxes for this cohort and camera, or None if not configured."""
    return CAMERA_ROIS.get(batch, {}).get(camera)


# --- fish alive per tank -----------------------------------------------------
# fish_alive.csv records the fish lost during the experiment:
#   batch,camera,tank,from_date,n_fish
#   b2,003,tank1,29082026,5
# meaning that from 29 August 2026 onward that tank held five fish. The number
# is carried through to the output as metadata. The activity index is a
# tank-level measure and is not divided by it.


def _date_from_ddmmyyyy(text):
    return date(int(text[4:8]), int(text[2:4]), int(text[0:2]))


def load_fish_alive():
    """Rows of fish_alive.csv as (batch, camera, tank, from_date, n_fish)."""
    rows = []
    if not os.path.exists(FISH_ALIVE_CSV):
        return rows
    with open(FISH_ALIVE_CSV, newline="") as handle:
        for record in csv.DictReader(handle):
            if not (record.get("batch") or "").strip():
                continue
            rows.append((record["batch"].strip(),
                         record["camera"].strip().zfill(3),
                         record["tank"].strip(),
                         _date_from_ddmmyyyy(record["from_date"].strip()),
                         int(record["n_fish"])))
    return rows


def n_fish_for(alive_rows, batch, camera, tank, ddmmyyyy):
    """Fish in that tank on that date, from the latest matching entry.

    Defaults to EXPECTED_FISH when no entry applies.
    """
    count, latest = EXPECTED_FISH, None
    try:
        day = _date_from_ddmmyyyy(ddmmyyyy)
    except ValueError:
        return count
    for row_batch, row_camera, row_tank, from_date, n_fish in alive_rows:
        if (row_batch == batch and row_camera == camera and row_tank == tank
                and from_date <= day and (latest is None or from_date >= latest)):
            count, latest = n_fish, from_date
    return count


# --- video access ------------------------------------------------------------


def probe_size(path):
    """Frame width and height of a video, read with ffprobe."""
    result = subprocess.run(
        ["ffprobe", "-v", "error", "-select_streams", "v:0",
         "-show_entries", "stream=width,height", "-of", "csv=p=0", path],
        capture_output=True, text=True, check=True)
    width, height = result.stdout.strip().split("\n")[0].split(",")[:2]
    return int(width), int(height)


def iter_frames(path, fps, max_sec=0):
    """Yield grey uint8 frames at the requested rate, piped from ffmpeg.

    Decoding straight from ffmpeg avoids writing an intermediate file and keeps
    the sampling rate exact for recordings with a variable frame rate.
    """
    width, height = probe_size(path)
    args = ["ffmpeg", "-v", "error"]
    if max_sec:
        args += ["-t", str(max_sec)]
    args += ["-i", path, "-vf", f"fps={fps}", "-f", "rawvideo", "-pix_fmt", "gray", "-"]
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                            bufsize=width * height * 4)
    frame_bytes = width * height
    try:
        while True:
            buf = proc.stdout.read(frame_bytes)
            if len(buf) < frame_bytes:
                break
            yield np.frombuffer(buf, np.uint8).reshape(height, width)
    finally:
        proc.stdout.close()
        proc.wait()


# --- light phase -------------------------------------------------------------

# Default phase of a recording, from the hour it started.
HOUR_LIGHT = {"06": "night", "14": "day", "21": "night"}
# The only switch that can occur inside a recording at that hour.
LIGHT_SWITCH = {"06": "up", "14": None, "21": "down"}
LIGHT_STEP_MIN = 8.0        # grey levels that count as a lights switch
LIGHT_STEP_WIN_SEC = 20.0   # window compared on each side of the step


def segment_by_light(brights, min_seg_frames, cutoff=None):
    """Split a brightness trace into runs above and below a cutoff.

    Returns [[start, end, is_night], ...]. Kept as a fallback for recordings
    whose start hour is unknown; segment_by_clock is what the pipeline uses.
    """
    cutoff = NIGHT_BRIGHTNESS_CUTOFF if cutoff is None else cutoff
    labels = [b < cutoff for b in brights]
    segments, start = [], 0
    for i in range(1, len(labels)):
        if labels[i] != labels[i - 1]:
            segments.append([start, i, labels[i - 1]])
            start = i
    segments.append([start, len(labels), labels[-1]])

    # Absorb segments shorter than the minimum into their neighbour.
    changed = True
    while changed and len(segments) > 1:
        changed = False
        for j in range(len(segments)):
            a, b, _ = segments[j]
            if b - a < min_seg_frames:
                if j == 0:
                    segments[1][0] = a
                else:
                    segments[j - 1][1] = b
                segments.pop(j)
                changed = True
                break

    merged = [segments[0]]
    for a, b, label in segments[1:]:
        if label == merged[-1][2]:
            merged[-1][1] = b
        else:
            merged.append([a, b, label])
    return merged


def segment_by_clock(brights, dt, hour):
    """Split a recording at a lights switch, using the clock hour as the anchor.

    The phase comes from the hour the recording started. A switch inside the
    recording is accepted only when brightness changes abruptly, measured as the
    median of about 20 s after the step minus the median of about 20 s before
    it, by at least LIGHT_STEP_MIN grey levels, and only in the direction that
    the hour allows: downward at 21:00 (lights going off), upward at 06:00
    (lights coming on), never at 14:00. This rejects the slow ramp that camera
    auto-exposure produces after lights-off, which a simple cutoff reads as
    lights coming back on.

    Returns [[start, end, is_night], ...] with at most one switch.
    """
    n = len(brights)
    default_night = HOUR_LIGHT.get(hour, "night") == "night"
    allowed = LIGHT_SWITCH.get(hour)
    min_seg = max(2, int(MIN_SEGMENT_SEC / dt))
    if allowed is None or n < 2 * min_seg + 1:
        return [[0, n, default_night]]

    trace = np.asarray(brights, float)
    window = max(1, int(round(LIGHT_STEP_WIN_SEC / dt)))
    best_index, best_step = None, 0.0
    for i in range(max(min_seg, window), min(n - min_seg, n - window) + 1):
        step = float(np.median(trace[i:i + window]) - np.median(trace[i - window:i]))
        if allowed == "down":
            step = -step
        if step > best_step:
            best_index, best_step = i, step

    if best_index is None or best_step < LIGHT_STEP_MIN:
        return [[0, n, default_night]]
    if allowed == "down":
        return [[0, best_index, False], [best_index, n, True]]
    return [[0, best_index, True], [best_index, n, False]]


def make_segments(brights, reservoir, dt, pct_fn=None, pct_list=None, cutoff=None,
                  hour=None):
    """Light-phase segments of a recording and one background image per segment.

    Args:
        brights: mean brightness of each processed frame.
        reservoir: (frame_index, frame) pairs sampled across the recording.
        dt: seconds between processed frames.
        pct_fn: called with is_night, returns the background percentile to use.
        pct_list: when given, return a {percentile: background} dict per segment
            instead of a single image (used by parameter selection).
        cutoff: brightness cutoff for the fallback segmentation.
        hour: start hour of the recording; when given, segmentation is by clock.

    Returns:
        (segments, backgrounds) with segments as [[start, end, is_night], ...].

    The background is a per-pixel percentile rather than a mean, so that a fish
    resting in the same spot for most of the recording is not absorbed into it.
    A high percentile does this for dark fish, a low one for fish that appear
    bright under infrared light; the percentile is selected per cohort, camera
    and light phase during parameter selection.
    """
    min_seg = max(2, int(MIN_SEGMENT_SEC / dt))
    if hour is not None:
        segments = segment_by_clock(brights, dt, hour)
    else:
        segments = segment_by_light(brights, min_seg, cutoff)

    backgrounds = []
    for a, b, is_night in segments:
        pool = [frame for index, frame in reservoir if a <= index < b]
        if not pool:
            middle = (a + b) // 2
            pool = [min(reservoir, key=lambda item: abs(item[0] - middle))[1]]
        if len(pool) > BG_SAMPLE_FRAMES:
            keep = np.linspace(0, len(pool) - 1, BG_SAMPLE_FRAMES).astype(int)
            pool = [pool[i] for i in keep]
        stack = np.stack(pool)
        if pct_list is not None:
            backgrounds.append({q: np.percentile(stack, q, axis=0).astype(np.uint8)
                                for q in pct_list})
        else:
            q = pct_fn(is_night) if pct_fn else 50
            backgrounds.append(np.percentile(stack, q, axis=0).astype(np.uint8))
    return segments, backgrounds


# --- detection ---------------------------------------------------------------

_OPEN_KERNEL = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (MORPH_OPEN_K, MORPH_OPEN_K))
_CLOSE_KERNELS = {}


def _close_kernel(size):
    if size not in _CLOSE_KERNELS:
        _CLOSE_KERNELS[size] = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (size, size))
    return _CLOSE_KERNELS[size]


def prep(roi, bg_mean):
    """Median filter a region and match its mean to the background.

    Matching the mean removes the slow gain changes of the infrared cameras,
    which would otherwise appear as whole-frame motion.
    """
    grey = cv2.medianBlur(roi, MEDIAN_BLUR_K)
    if ILLUM_NORMALIZE:
        shift = int(round(bg_mean - float(grey.mean())))
        grey = np.clip(grey.astype(np.int16) + shift, 0, 255).astype(np.uint8)
    return grey


def diff_image(roi, bg_roi_prepped, bg_mean):
    """Absolute difference between a prepared region and its background."""
    return cv2.absdiff(prep(roi, bg_mean), bg_roi_prepped)


def raw_pieces(diff, thr, close_k):
    """Threshold and clean a difference image.

    Returns (fraction of the region above threshold, [(area, cx, cy), ...]) for
    pieces larger than SPECK_AREA.
    """
    _, mask = cv2.threshold(diff, thr, 255, cv2.THRESH_BINARY)
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, _OPEN_KERNEL)
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, _close_kernel(close_k))
    occupancy = float(np.count_nonzero(mask)) / mask.size
    count, _, stats, centroids = cv2.connectedComponentsWithStats(mask, connectivity=8)
    pieces = [(int(stats[i, cv2.CC_STAT_AREA]), float(centroids[i, 0]), float(centroids[i, 1]))
              for i in range(1, count) if stats[i, cv2.CC_STAT_AREA] >= SPECK_AREA]
    return occupancy, pieces


def merge_pieces(pieces, merge_dist):
    """Join pieces whose centroids lie within merge_dist pixels.

    Single-link clustering, so that one fish broken into head and tail by the
    threshold is counted once.
    """
    if merge_dist <= 0 or len(pieces) < 2:
        return list(pieces)

    parent = list(range(len(pieces)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    positions = np.array([(p[1], p[2]) for p in pieces])
    squared = ((positions[:, None, :] - positions[None, :, :]) ** 2).sum(-1)
    rows, cols = np.nonzero(np.triu(squared <= merge_dist * merge_dist, 1))
    for i, j in zip(rows.tolist(), cols.tolist()):
        parent[find(i)] = find(j)

    clusters = {}
    for i in range(len(pieces)):
        clusters.setdefault(find(i), []).append(pieces[i])

    merged = []
    for group in clusters.values():
        area = sum(p[0] for p in group)
        cx = sum(p[0] * p[1] for p in group) / area
        cy = sum(p[0] * p[2] for p in group) / area
        merged.append((area, cx, cy))
    return merged


def fish_from_pieces(pieces, min_area, merge_dist, single_area):
    """Blobs kept as fish, as (area, cx, cy, n_fish).

    n_fish is above one when a blob is much larger than a single fish, which
    happens when fish touch or overlap in the side view.
    """
    blobs = [b for b in merge_pieces(pieces, merge_dist)
             if min_area <= b[0] <= MAX_BLOB_AREA]
    fish = []
    for area, cx, cy in blobs:
        count = 1
        if single_area and area > MULTI_FISH_RATIO * single_area:
            count = max(1, int(round(area / single_area)))
        fish.append((area, cx, cy, count))
    return fish


def detect(roi, bg_prepped, bg_mean, params):
    """Detect fish in one tank region. Returns (occupancy, fish list)."""
    diff = diff_image(roi, bg_prepped, bg_mean)
    occupancy, pieces = raw_pieces(diff, params["thr"], params["close_k"])
    fish = fish_from_pieces(pieces, params["min_area"], params["merge_dist"],
                            params.get("single_area"))
    return occupancy, fish


def motion_fraction(current, previous, thr):
    """Fraction of pixels that changed between two consecutive prepared frames.

    This is the activity index analysed in the paper. It is measured from the
    images alone and does not depend on how many fish the detector finds.
    """
    diff = cv2.absdiff(current, previous)
    _, mask = cv2.threshold(diff, thr, 255, cv2.THRESH_BINARY)
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, _OPEN_KERNEL)
    return float(np.count_nonzero(mask)) / mask.size


def prep_bg(bg_full, roi_box):
    """Background region and its mean, ready for detect()."""
    x0, y0, x1, y1 = roi_box
    region = bg_full[y0:y1, x0:x1]
    return cv2.medianBlur(region, MEDIAN_BLUR_K), float(region.mean())
