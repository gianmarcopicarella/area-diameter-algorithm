import os

import numpy as np
from matplotlib import pyplot as plt

import constants
import utils


def scatter_points(P, side=30, numbers=False):
    pts = np.array(P)
    plt.scatter(pts[:, 0], pts[:, 1], c="black", s=1)
    if numbers:
        for i, p in enumerate(pts):
            plt.text(p[0] - 0.02, p[1], f'{i}', fontsize=3, ha='right')
    plt.gcf().set_dpi(400)
    plt.gca().set_xlim([0, side])
    plt.gca().set_ylim([0, side])
    plt.gca().set_aspect("equal")
    plt.show()


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


def generate_synthetic_data():
    path_to_uniform = os.path.join(constants.PATH_TO_EXPERIMENTS, "uniform")
    path_to_gaussian = os.path.join(constants.PATH_TO_EXPERIMENTS, "gaussian")

    utils.prepare_path(path_to_uniform)
    utils.prepare_path(path_to_gaussian)

    for i, count in enumerate(constants.DENSITY_VALUES):
        path_to_uniform_samples = os.path.join(path_to_uniform, str(i))
        utils.prepare_path(path_to_uniform_samples)
        for j in range(constants.SYNTHETIC_BENCHMARK_UNIFORM_ITERATIONS):
            uniform_data = {"count": int(count),
                            "points": [{"x": x, "y": y} for x, y in rand_points_in_square(count)]}
            utils.write_json(os.path.join(path_to_uniform_samples, f"points_{j}.json"), uniform_data)

    gaussian_point_set_size = 200
    for i, std in enumerate(constants.STD_VALUES):
        path_to_gaussian_samples = os.path.join(path_to_gaussian, str(i))
        utils.prepare_path(path_to_gaussian_samples)
        for j in range(constants.SYNTHETIC_BENCHMARK_GAUSSIAN_ITERATIONS):
            gaussian_data = {"count": int(gaussian_point_set_size),
                             "points": [{"x": x, "y": y} for x, y in
                                        norm_points_in_square(gaussian_point_set_size, std)]}
            utils.write_json(os.path.join(path_to_gaussian_samples, f"points_{j}.json"), gaussian_data)


def generate_real_data():
    def to_mm(v):
        micron_per_pixel = 0.23
        factor = micron_per_pixel / 1000.0
        return v * factor

    raw_dataset = utils.read_json(constants.PATH_TO_RAW_DATASET, "rb")
    dataset = []
    for entry in raw_dataset:
        if constants.DETECTION_METHOD not in entry:
            continue
        sample = entry[constants.DETECTION_METHOD]
        points = {(to_mm(float(s["x"])), to_mm(float(s["y"]))) for s in sample
                  if s["prob"] >= constants.DETECTION_CONFIDENCE_THRESHOLD}
        if len(points) == 0:
            continue
        dataset.append(list(points))
    dataset.sort(key=lambda x: len(x), reverse=True)

    path_to_real = os.path.join(constants.PATH_TO_EXPERIMENTS, "real")
    utils.prepare_path(path_to_real)

    assert (constants.REAL_BENCHMARKS_COUNT < len(dataset))

    for i in range(constants.REAL_BENCHMARKS_COUNT):
        path_to_real_samples = os.path.join(path_to_real, str(i))
        utils.prepare_path(path_to_real_samples)
        ann_info = utils.ann(dataset[i])
        real_data = {"count": len(dataset[i]),
                     "points": [{"x": p[0], "y": p[1]} for p in dataset[i]],
                     "ann_avg": ann_info[0], "ann_std": ann_info[1]}
        utils.write_json(os.path.join(path_to_real_samples, "points_0.json"), real_data)
        # scatter_points(dataset[i])
        # print(real_data["count"], real_data["ann_avg"], real_data["ann_std"])


np.random.seed(0)
generate_synthetic_data()
generate_real_data()
