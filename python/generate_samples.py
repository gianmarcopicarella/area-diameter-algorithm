import os

import numpy as np
import json
import shutil

PATH_TO_DATA = os.path.join("..", "data")
PATH_TO_EXPERIMENTS = os.path.join(PATH_TO_DATA, "samples", "experiments")
PATH_TO_RAW = os.path.join(PATH_TO_DATA, "raw")
RAW_DATA_FILENAME = "detections_subset.json"


def rand_points_in_square(count, s=20):
    pts = set()
    while len(pts) < count:
        gen_pts = np.random.uniform(low=[0, 0], high=[s, s], size=(1000, 2))
        i = 0
        while len(pts) < count and i < len(gen_pts):
            p = (float(gen_pts[i][0]), float(gen_pts[i][1]))
            if 0 <= gen_pts[i][0] <= s and 0 <= gen_pts[i][1] <= s and p not in pts:
                pts.add(p)
            i += 1
    assert (len(pts) == count)
    return list(pts)


def norm_points_in_square(count, std=1, s=20, c=10):
    pts = set()
    while len(pts) < count:
        gen_pts = np.random.normal(loc=c, scale=std, size=(1000, 2))
        i = 0
        while len(pts) < count and i < len(gen_pts):
            p = (float(gen_pts[i][0]), float(gen_pts[i][1]))
            if 0 <= gen_pts[i][0] <= s and 0 <= gen_pts[i][1] <= s and p not in pts:
                pts.add(p)
            i += 1
    assert (len(pts) == count)
    return list(pts)


def generate_synthetic_data(uniform_x_values, gaussian_x_values, iterations_count):

    PATH_TO_UNIFORM = os.path.join(PATH_TO_EXPERIMENTS, "uniform")
    if os.path.isdir(PATH_TO_UNIFORM):
        shutil.rmtree(PATH_TO_UNIFORM)
    os.makedirs(PATH_TO_UNIFORM)

    for i, n in enumerate(uniform_x_values):
        os.makedirs(os.path.join(PATH_TO_UNIFORM, f"{i}"))
        for j in range(iterations_count):
            data = {"count": int(n),
                    "points": [{"x": x, "y": y} for x, y in rand_points_in_square(n)]}
            json_path = os.path.join(PATH_TO_UNIFORM, f"{i}", f"points_{j}.json")
            with open(json_path, 'w') as file:
                json.dump(data, file)

    PATH_TO_GAUSSIAN = os.path.join(PATH_TO_EXPERIMENTS, "gaussian")
    if os.path.isdir(PATH_TO_GAUSSIAN):
        shutil.rmtree(PATH_TO_GAUSSIAN)
    os.makedirs(PATH_TO_GAUSSIAN)

    gaussian_n = uniform_x_values[0]
    for i, std in enumerate(gaussian_x_values):
        os.makedirs(os.path.join(PATH_TO_GAUSSIAN, f"{i}"))
        for j in range(iterations_count):
            data = {"count": int(gaussian_n),
                    "points": [{"x": x, "y": y} for x, y in norm_points_in_square(gaussian_n, std)]}
            json_path = os.path.join(PATH_TO_GAUSSIAN, f"{i}", f"points_{j}.json")
            with open(json_path, 'w') as file:
                json.dump(data, file)


def generate_real_data(ALGORITHM="midog21_1st_stage", THRESHOLD=0.64):
    def to_mm(v):
        micron_per_pixel = 0.23
        factor = micron_per_pixel / 1000.0
        return v * factor

    with open(os.path.join(PATH_TO_RAW, RAW_DATA_FILENAME), "rb") as file:
        raw_data = json.load(file)
    data = []
    for entry in raw_data:
        if ALGORITHM not in entry:
            continue
        sample = entry[ALGORITHM]
        points = {(to_mm(float(s["x"])), to_mm(float(s["y"]))) for s in sample if s["prob"] >= THRESHOLD}
        if len(points) == 0:
            continue
        data.append(list(points))
    data.sort(key=lambda x: len(x), reverse=True)

    PATH_TO_REAL = os.path.join(PATH_TO_EXPERIMENTS, "real")
    if os.path.isdir(PATH_TO_REAL):
        shutil.rmtree(PATH_TO_REAL)
    os.makedirs(PATH_TO_REAL)

    for i in range(10):
        os.makedirs(os.path.join(PATH_TO_REAL, f"{i}"))
        out_data = {"count": len(data[i]),
                    "points": [{"x": p[0], "y": p[1]} for p in data[i]]}
        json_path = os.path.join(PATH_TO_REAL, f"{i}", "points_0.json")
        with open(json_path, 'w') as file:
            json.dump(out_data, file)


np.random.seed(0)
generate_synthetic_data(np.arange(100, 201, 10), np.arange(0.5, 5.51, 0.5), 100)
generate_real_data()
