#!/usr/bin/env python3
"""Re-digitise the measured markers of Nguyen et al. (2022), Figs. 4 and 5 (core/data/nguyen_2022/README.md).

The figures are 300 dpi raster images in the paper's PDF. Extract them first (poppler):

    pdfimages -f 7 -l 7 -png nguyen_2022.pdf p7      # p7-001.png: Fig. 4
    pdfimages -f 9 -l 9 -png nguyen_2022.pdf p9      # p9-000.png: Fig. 5
    python core/data/nguyen_2022/digitise_figures.py --images DIR --out core/data/nguyen_2022

The pixel calibrations (tick positions), the panel boxes and the legend boxes below were
read from those two images. Needs numpy, scipy and Pillow.
"""
import argparse
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

# ---- axis calibration -------------------------------------------------------------
A = np.pi * 0.0282**2 / 4
def lv(ccm): return 0.0786 / (ccm * 1e-6 / 60 / (0.4 * A))   # L/v [s] at eps = 0.4
S4, S5 = lv(50), lv(100)
def lin(pairs):
    p = sorted(pairs); xs = [a for a, _ in p]; vs = [b for _, b in p]
    k = (vs[-1] - vs[0]) / (xs[-1] - xs[0])
    fwd = lambda x: vs[0] + (np.asarray(x, float) - xs[0]) * k
    inv = lambda v: xs[0] + (np.asarray(v, float) - vs[0]) / k
    return fwd, inv
# Fig 5 broken axis (panels a and c): left segment 0..20, right segment linear through the 200..1000 ticks
def broken(left, right):
    lf, li = lin(left); rf, ri = lin(right)
    def fwd(x):
        x = np.asarray(x, float); return np.where(x <= left[-1][0] + 4, lf(x), rf(x))
    def inv(v):
        v = np.asarray(v, float); return np.where(v <= 20.0001, li(v), ri(v))
    return fwd, inv
PANELS = {
  "fig4a": dict(img="p7-001.png", x=lin([(99.5, 0), (670.5, 300)]), y=lin([(514.5, 0.0), (21.5, 1.0)]), box=(101, 22, 669, 513)),
  "fig4b": dict(img="p7-001.png", x=lin([(846, 0), (1393.5, 300)]), y=lin([(518.5, 20), (19, 60)]), box=(848, 20, 1392, 517)),
  "fig5a": dict(img="p9-000.png", x=broken([(114.5, 0), (218, 20)], [(312.5, 200), (668.5, 1000)]), y=lin([(515.5, 0.0), (27.5, 1.0)]), box=(116, 24, 667, 514)),
  "fig5b": dict(img="p9-000.png", x=lin([(836.5, 0), (1386.5, 1000)]), y=lin([(511, 20), (100.5, 45)]), box=(838, 20, 1385, 510)),
  "fig5c": dict(img="p9-000.png", x=broken([(119, 0), (221.5, 20)], [(315.5, 200), (668.5, 1000)]), y=lin([(1115, 0.0), (652.5, 0.30)]), box=(121, 653, 667, 1114)),
}

# ---- marker colours and legend boxes ------------------------------------------------
FILL = {
    "red": lambda a: (a[..., 0] > 170) & (a[..., 1] < 80) & (a[..., 2] < 80),
    "green": lambda a: (a[..., 1] > 110) & (a[..., 0] < 70) & (a[..., 2] < 70),
    "blue": lambda a: (a[..., 2] > 140) & (a[..., 0] < 70) & (a[..., 1] < 110),
}
LEGEND = {"fig4a": (468, 356, 654, 496), "fig4b": (1188, 50, 1374, 192), "fig5a": (395, 356, 658, 494),
          "fig5b": (1112, 39, 1381, 178), "fig5c": (382, 679, 652, 817)}



# ---- isolated markers ----------------------------------------------------------------
def disk(r):
    y, x = np.mgrid[-r:r + 1, -r:r + 1]; return x * x + y * y <= r * r

def series_points(panel, colour, r_open=3, step=6):
    P = PANELS[panel]
    img = np.asarray(Image.open(P["img"]).convert("RGB")).astype(int)
    x0, y0, x1, y1 = P["box"]; a = img[y0:y1, x0:x1]
    dark = a.max(axis=2) < 110
    lx0, ly0, lx1, ly1 = LEGEND[panel]
    excl = np.zeros(dark.shape, bool); excl[max(0, ly0 - y0):max(0, ly1 - y0), max(0, lx0 - x0):max(0, lx1 - x0)] = True
    fill = FILL[colour](a) & ~excl
    core = ndi.binary_opening(fill, structure=disk(r_open))
    near = ndi.binary_dilation(core, structure=disk(9))
    blob = (fill | dark) & near & ~excl
    lab, n = ndi.label(core)
    pts = []
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        h = sl[0].stop - sl[0].start; w = sl[1].stop - sl[1].start
        comp = lab[sl] == i
        if comp.sum() < 12:
            continue
        if w <= 14 and h <= 14:        # one marker
            cy, cx = ndi.center_of_mass(comp); pts.append((cx + sl[1].start, cy + sl[0].start, "m"))
        else:                           # band: sample its centre line
            cmask = np.zeros(core.shape, bool); cmask[sl] = comp
            region = blob & ndi.binary_dilation(cmask, structure=disk(9))
            for cx in range(sl[1].start, sl[1].stop, step):
                rows = np.flatnonzero(region[:, cx])
                if len(rows) == 0:
                    continue
                # split into runs, take runs that contain core pixels
                edges = np.flatnonzero(np.diff(np.r_[0, region[:, cx].astype(int), 0]))
                for s, e in zip(edges[::2], edges[1::2]):
                    if cmask[s:e, max(0, cx - 2):cx + 3].any():
                        if e - s <= 40:
                            pts.append((cx, (s + e - 1) / 2, "b"))
                        else:   # a steep stretch of the band: one point per marker-sized slice
                            for ss in range(s, e, 12):
                                pts.append((cx, min(ss + 6, e - 1), "b"))
    out = []
    for px, py, kind in pts:
        out.append((float(P["x"][0](px + x0)), float(P["y"][0](py + y0)), px + x0, py + y0, kind))
    out.sort()
    return out


# ---- dense marker bands -----------------------------------------------------------------
SPEC = {"fig4a": {"100": "red", "70": "blue", "50": "green"}, "fig4b": {"100": "red", "70": "blue", "50": "green"},
        "fig5a": {"30": "red", "15": "blue", "5": "green"}, "fig5b": {"30": "red", "15": "blue", "5": "green"},
        "fig5c": {"30": "red", "15": "blue", "5": "green"}}
def band_points(panel, series, step=6, r_thick=5, hmax=45, min_fill=3):
    P = PANELS[panel]
    img = np.asarray(Image.open(P["img"]).convert("RGB")).astype(int)
    x0, y0, x1, y1 = P["box"]; a = img[y0:y1, x0:x1]
    dark = a.max(axis=2) < 110
    lx0, ly0, lx1, ly1 = LEGEND[panel]
    excl = np.zeros(dark.shape, bool); excl[max(0, ly0 - y0):max(0, ly1 - y0), max(0, lx0 - x0):max(0, lx1 - x0)] = True
    fills = {name: FILL[c](a) & ~excl for name, c in series.items()}
    union = dark.copy()
    for f in fills.values(): union |= f
    union &= ~excl
    thick = ndi.binary_opening(union, structure=disk(r_thick))
    out = {name: [] for name in series}
    for cx in range(3, thick.shape[1] - 3, step):
        col = thick[:, cx]
        edges = np.flatnonzero(np.diff(np.r_[0, col.astype(int), 0]))
        for s, e in zip(edges[::2], edges[1::2]):
            if e - s > hmax: continue
            counts = {n: f[max(0, s - 2):e + 2, cx - 4:cx + 5].sum() for n, f in fills.items()}
            best = max(counts, key=counts.get)
            if counts[best] < min_fill: continue
            py = (s + e - 1) / 2
            out[best].append((float(P["x"][0](cx + x0)), float(P["y"][0](py + y0)), cx + x0, py + y0, "b"))
    return out

def all_points(panel):
    series = SPEC[panel]
    bands = band_points(panel, series)
    res = {}
    for name, colour in series.items():
        mk = series_points(panel, colour)
        mk = [p for p in mk if p[4] == "m"]
        mxy = np.array([(p[2], p[3]) for p in mk]) if mk else np.zeros((0, 2))
        bd = [p for p in bands[name] if not (len(mxy) and np.min(np.hypot(mxy[:, 0] - p[2], mxy[:, 1] - p[3])) < 9)]
        res[name] = sorted(mk + bd)
    return res


OUTFILES = {"fig4a": ("Fig4", "Y_CO2_{}%.csv", "Time [s], Y_CO2 [-]", "S4", 0.0),
            "fig4b": ("Fig4", "T_CO2_{}%.csv", "Time [s], T [K]", "S4", 273.15),
            "fig5c": ("Fig5", "Y_CO2_{}%.csv", "Time [s], Y_CO2 [-]", "S5", 0.0),
            "fig5a": ("Fig5", "Y_N2_{}%.csv", "Time [s], Y_N2 [-]", "S5", 0.0),
            "fig5b": ("Fig5", "T_CO2_{}%.csv", "Time [s], T [K]", "S5", 273.15)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--images", type=Path, required=True, help="folder with p7-001.png and p9-000.png")
    ap.add_argument("--out", type=Path, required=True, help="folder to write Fig4/ and Fig5/ into")
    args = ap.parse_args()
    for P in PANELS.values():
        P["img"] = str(args.images / P["img"])
    raw = {}
    for panel in SPEC:
        pts = all_points(panel)
        for name in SPEC[panel]:
            raw[f"{panel}/{name}"] = [(p[0], p[1]) for p in pts[name]]
    # Baseline: before breakthrough the markers of all series overlap at y = 0, so the
    # band cannot be told apart by colour. Points at y = 0 are placed every 10 units of
    # dimensionless time up to 5 units before each series' first resolved marker.
    for panel, ymax in (("fig4a", 0.02), ("fig5c", 0.004)):
        for key in [k for k in raw if k.startswith(panel + "/")]:
            first = min(p[0] for p in raw[key] if p[1] > ymax)
            raw[key] = sorted(raw[key] + [(float(x), 0.0) for x in np.arange(0.0, first - 5, 10.0)], key=lambda p: p[0])
    scales = {"S4": S4, "S5": S5}
    for key, pts in raw.items():
        panel, ser = key.split("/")
        fig, name, header, sc, dT = OUTFILES[panel]
        arr = sorted({(round(max(p[0], 0.0), 3), round(p[1], 4)) for p in pts if p[0] >= -0.5})
        (args.out / fig).mkdir(parents=True, exist_ok=True)
        with open(args.out / fig / name.format(ser), "w") as fh:
            fh.write(header + "\n")
            for tb, v in arr:
                fh.write(f"{tb * scales[sc]:.6g},{v + dT:.6g}\n")
        print(fig, name.format(ser), len(arr), "points")


if __name__ == "__main__":
    main()
