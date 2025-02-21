import os

import matplotlib.pyplot as plt
import numpy as np

import json
import constants
import utils

from matplotlib.patches import Polygon

data_runs = utils.read_json(constants.PATH_TO_BENCHMARK_BASE_REPORT)
data_results = utils.read_json(constants.PATH_TO_BENCHMARK_CUSTOM_REPORT)


def scatter_points(P, poly):
    plt.gca().add_patch(Polygon(poly, closed=False, edgecolor="red", facecolor='none', linewidth=1))
    pts = np.array(P)
    plt.scatter(pts[:, 0], pts[:, 1], c="black", s=8)
    for i, p in enumerate(pts):
        plt.text(p[0] - 0.2, p[1], f'{i}', fontsize=9, ha='right')
    plt.gcf().set_dpi(200)
    plt.gca().set_aspect("equal")
    plt.show()


def to_kb(bytes):
    return bytes / 1000


def extract_run_data_for(data, distribution, x_values, diameters):
    result = {"Eppstein": {"time": [], "memory": [], "entries": [], "required_entries": []}}
    for d in diameters:
        result["Antipodal_" + str(d)] = {"time": [], "memory": [], "entries": []}

    for r in data["benchmarks"]:
        split_id = r["name"].split("/")
        if distribution == split_id[1]:
            x_step = int(split_id[2])
            tt = (float(x_values[x_step]), float(r["cpu_time"]), 0)
            mt = (float(x_values[x_step]), to_kb(float(r["mem_avg"])), to_kb(float(r["mem_std"])))
            et = (float(x_values[x_step]), float(r["entries_avg"]), float(r["entries_std"]))
            diameter = "_" + split_id[3] if len(split_id) == 5 else ""
            key = split_id[0] + diameter
            result[key]["time"].append(tt)
            result[key]["memory"].append(mt)
            result[key]["entries"].append(et)
            if "required_entries" in result[key]:
                met = (float(x_values[x_step]), float(r["min_entries_avg"]), float(r["min_entries_std"]))
                result[key]["required_entries"].append(met)

    for k in result:
        for ik in result[k]:
            result[k][ik].sort()

    return result


def extract_solutions_data_for(data, distribution, x_values, diameters):
    result = {"Eppstein": {"count": [], "area": [], "diameter": []}}
    for d in diameters:
        result["Antipodal_" + str(d)] = {"count": [], "area": [], "diameter": []}

    for i, r in enumerate(data["results"]):
        split_id = r[0]["id"].split("/")
        if distribution == split_id[1]:
            max_diam = split_id[3]
            algorithm = split_id[0] + "_" + max_diam if split_id[0] == "Antipodal" else split_id[0]
            x_step = split_id[2]
            counts, areas, diams = [], [], []
            for iteration, e in enumerate(r):
                assert ("convex_area" in e)
                sol = e["convex_area"]
                counts.append(sol["count"])
                areas.append(sol["area"])
                path_to_sample = os.path.join(constants.PATH_TO_EXPERIMENTS, distribution.lower(), x_step,
                                              f"points_{iteration}.json")
                with open(path_to_sample, 'r') as file:
                    points = json.load(file)

                fi = int(sol["diameter_indices"][0])
                si = int(sol["diameter_indices"][1])
                fdp = np.array((points["points"][fi]["x"], points["points"][fi]["y"]), dtype=np.float32)
                sdp = np.array((points["points"][si]["x"], points["points"][si]["y"]), dtype=np.float32)

                diams.append(float(np.linalg.norm(fdp - sdp)))

            x_idx = i % len(x_values)
            ct = (float(x_values[x_idx]), float(np.mean(counts)), float(np.std(counts)))
            at = (float(x_values[x_idx]), float(np.mean(areas)), float(np.std(areas)))
            dt = (float(x_values[x_idx]), float(np.mean(diams)), float(np.std(diams)))

            result[algorithm]["count"].append(ct)
            result[algorithm]["area"].append(at)
            result[algorithm]["diameter"].append(dt)

    return result


input_size = np.arange(100, 201, 10)
std_dev = np.arange(0.5, 5.51, 0.5)
diameters = [2, 3, 4, 5]

uniform_run_data = extract_run_data_for(data_runs, "Uniform", input_size, diameters)
gaussian_run_data = extract_run_data_for(data_runs, "Gaussian", std_dev, diameters)

uniform_sol_data = extract_solutions_data_for(data_results, "Uniform", input_size, diameters)
gaussian_sol_data = extract_solutions_data_for(data_results, "Gaussian", std_dev, diameters)


# print(uniform_run_data)
# print(gaussian_run_data)


def base_unroll(data, algorithm, key):
    s = ""
    for e in data[algorithm][key]:
        s += "(" + str(round(e[0], 2)) + ", " + str(round(e[1], 2)) + ") "
    return s


def std_unroll(data, algorithm, key):
    s = ""
    for e in data[algorithm][key]:
        s += "(" + str(round(e[0], 2)) + ", " + str(round(e[1], 2)) + ") +- (0, " + str(round(e[2], 2)) + ") "
    return s


colours = ["Cyan", "Magenta", "Green", "Peach", "Fuchsia"]
rotate = ["0", "0", "0", "180", "0"]
marks = ["diamond", "o", "asterisk", "triangle", "triangle", "triangle"]

print("[Uniform data]")

print("1) Time")

for i, d in enumerate(diameters):
    print("\\addplot+[color=" + colours[i] + ", mark=" + marks[i] + ",mark options={rotate=" + rotate[
        i] + "}] coordinates {" + base_unroll(uniform_run_data, "Antipodal_" + str(d),
                                              "time") + "}; \\addlegendentry{$\\text{Antipodal}_{" + str(d) + "}$}")

print("\\addplot+[color=" + colours[-1] + ", mark=" + marks[-1] + ",mark options={rotate=" + rotate[
    -1] + "}] coordinates {" + base_unroll(uniform_run_data, "Eppstein",
                                           "time") + "}; \\addlegendentry{$\\text{Eppstein}_{\\infty}$}")

print("2) Memory")

for i, d in enumerate(diameters):
    print("\\addplot+[color=" + colours[i] + ", mark=" + marks[i] + ",mark options={rotate=" + rotate[
        i] + "}] coordinates {" + base_unroll(uniform_run_data, "Antipodal_" + str(d),
                                              "memory") + "}; \\addlegendentry{$\\text{Antipodal}_{" + str(d) + "}$}")

print("\\addplot+[color=" + colours[-1] + ", mark=" + marks[-1] + ",mark options={rotate=" + rotate[
    -1] + "}] coordinates {" + base_unroll(uniform_run_data, "Eppstein",
                                           "memory") + "}; \\addlegendentry{$\\text{Eppstein}_{\\infty}$}")

print("3) Entries")

for i, d in enumerate(diameters):
    print("\\addplot+[color=" + colours[i] + ", mark=" + marks[i] + ",mark options={rotate=" + rotate[
        i] + "}, error bars/.cd, y dir=both, y explicit] coordinates {" + std_unroll(uniform_run_data,
                                                                                     "Antipodal_" + str(d),
                                                                                     "entries") + "}; \\addlegendentry{$\\text{Antipodal}_{" + str(
        d) + "}$}")

print("\\addplot+[color=" + colours[-1] + ", mark=" + marks[-1] + ",mark options={rotate=" + rotate[
    -1] + "}] coordinates {" + std_unroll(uniform_run_data, "Eppstein",
                                          "entries") + "}; \\addlegendentry{$\\text{Eppstein}_{\\infty}$}")
print(
    "\\addplot+[color=" + colours[-1] + ",dashed, mark=" + marks[-2] + "] coordinates {" + std_unroll(uniform_run_data,
                                                                                                      "Eppstein",
                                                                                                      "required_entries") + "}; \\addlegendentry{$\\text{Eppstein}_{\\infty}^{\\star}$}")

# -------

print("4) Cardinality")

for i, d in enumerate(diameters):
    print("\\addplot+[color=" + colours[i] + ", mark=" + marks[i] + ",mark options={rotate=" + rotate[
        i] + "}, error bars/.cd, y dir=both, y explicit] coordinates {" + std_unroll(uniform_sol_data,
                                                                                     "Antipodal_" + str(d),
                                                                                     "count") + "}; \\addlegendentry{$\\text{Antipodal}_{" + str(
        d) + "}$}")

print("\\addplot+[color=" + colours[-1] + ", mark=" + marks[-1] + ",mark options={rotate=" + rotate[
    -1] + "}, error bars/.cd, y dir=both, y explicit] coordinates {" + std_unroll(uniform_sol_data, "Eppstein",
                                                                                  "count") + "}; \\addlegendentry{$\\text{Eppstein}_{\\infty}$}")

print("5) Area")

for i, d in enumerate(diameters):
    print("\\addplot+[color=" + colours[i] + ", mark=" + marks[i] + ",mark options={rotate=" + rotate[
        i] + "}, error bars/.cd, y dir=both, y explicit] coordinates {" + std_unroll(uniform_sol_data,
                                                                                     "Antipodal_" + str(d),
                                                                                     "area") + "}; \\addlegendentry{$\\text{Antipodal}_{" + str(
        d) + "}$}")

print("\\addplot+[color=" + colours[-1] + ", mark=" + marks[-1] + ",mark options={rotate=" + rotate[
    -1] + "}, error bars/.cd, y dir=both, y explicit] coordinates {" + std_unroll(uniform_sol_data, "Eppstein",
                                                                                  "area") + "}; \\addlegendentry{$\\text{Eppstein}_{\\infty}$}")

print("6) Diameter")

for i, d in enumerate(diameters):
    print("\\addplot+[color=" + colours[i] + ", mark=" + marks[i] + ",mark options={rotate=" + rotate[
        i] + "}, error bars/.cd, y dir=both, y explicit] coordinates {" + std_unroll(uniform_sol_data,
                                                                                     "Antipodal_" + str(d),
                                                                                     "diameter") + "}; \\addlegendentry{$\\text{Antipodal}_{" + str(
        d) + "}$}")

print("\\addplot+[color=" + colours[-1] + ", mark=" + marks[-1] + ",mark options={rotate=" + rotate[
    -1] + "}, error bars/.cd, y dir=both, y explicit] coordinates {" + std_unroll(uniform_sol_data, "Eppstein",
                                                                                  "diameter") + "}; \\addlegendentry{$\\text{Eppstein}_{\\infty}$}")

print("[Gaussian data]")

print("1) Time")

for i, d in enumerate(diameters):
    print("\\addplot+[color=" + colours[i] + ", mark=" + marks[i] + ",mark options={rotate=" + rotate[
        i] + "}] coordinates {" + base_unroll(gaussian_run_data, "Antipodal_" + str(d),
                                              "time") + "}; \\addlegendentry{$\\text{Antipodal}_{" + str(d) + "}$}")

print("\\addplot+[color=" + colours[-1] + ", mark=" + marks[-1] + ",mark options={rotate=" + rotate[
    -1] + "}] coordinates {" + base_unroll(gaussian_run_data, "Eppstein",
                                           "time") + "}; \\addlegendentry{$\\text{Eppstein}_{\\infty}$}")

print("2) Memory")

for i, d in enumerate(diameters):
    print("\\addplot+[color=" + colours[i] + ", mark=" + marks[i] + ",mark options={rotate=" + rotate[
        i] + "}] coordinates {" + base_unroll(gaussian_run_data, "Antipodal_" + str(d),
                                              "memory") + "}; \\addlegendentry{$\\text{Antipodal}_{" + str(d) + "}$}")

print("\\addplot+[color=" + colours[-1] + ", mark=" + marks[-1] + ",mark options={rotate=" + rotate[
    -1] + "}] coordinates {" + base_unroll(gaussian_run_data, "Eppstein",
                                           "memory") + "}; \\addlegendentry{$\\text{Eppstein}_{\\infty}$}")

print("3) Entries")

for i, d in enumerate(diameters):
    print("\\addplot+[color=" + colours[i] + ", mark=" + marks[i] + ",mark options={rotate=" + rotate[
        i] + "}, error bars/.cd, y dir=both, y explicit] coordinates {" + std_unroll(gaussian_run_data,
                                                                                     "Antipodal_" + str(d),
                                                                                     "entries") + "}; \\addlegendentry{$\\text{Antipodal}_{" + str(
        d) + "}$}")

print("\\addplot+[color=" + colours[-1] + ", mark=" + marks[-1] + ",mark options={rotate=" + rotate[
    -1] + "}] coordinates {" + std_unroll(gaussian_run_data, "Eppstein",
                                          "entries") + "}; \\addlegendentry{$\\text{Eppstein}_{\\infty}$}")
print(
    "\\addplot+[color=" + colours[-1] + ",dashed, mark=" + marks[-2] + "] coordinates {" + std_unroll(
        gaussian_run_data,
        "Eppstein",
        "required_entries") + "}; \\addlegendentry{$\\text{Eppstein}_{\\infty}^{\\star}$}")

print("4) Cardinality")
for i, d in enumerate(diameters):
    print("\\addplot+[color=" + colours[i] + ", mark=" + marks[i] + ",mark options={rotate=" + rotate[
        i] + "}, error bars/.cd, y dir=both, y explicit] coordinates {" + std_unroll(gaussian_sol_data,
                                                                                     "Antipodal_" + str(d),
                                                                                     "count") + "}; \\addlegendentry{$\\text{Antipodal}_{" + str(
        d) + "}$}")

print("\\addplot+[color=" + colours[-1] + ", mark=" + marks[-1] + ",mark options={rotate=" + rotate[
    -1] + "}, error bars/.cd, y dir=both, y explicit] coordinates {" + std_unroll(gaussian_sol_data, "Eppstein",
                                                                                  "count") + "}; \\addlegendentry{$\\text{Eppstein}_{\\infty}$}")

print("5) Area")
for i, d in enumerate(diameters):
    print("\\addplot+[color=" + colours[i] + ", mark=" + marks[i] + ",mark options={rotate=" + rotate[
        i] + "}, error bars/.cd, y dir=both, y explicit] coordinates {" + std_unroll(gaussian_sol_data,
                                                                                     "Antipodal_" + str(d),
                                                                                     "area") + "}; \\addlegendentry{$\\text{Antipodal}_{" + str(
        d) + "}$}")

print("\\addplot+[color=" + colours[-1] + ", mark=" + marks[-1] + ",mark options={rotate=" + rotate[
    -1] + "}, error bars/.cd, y dir=both, y explicit] coordinates {" + std_unroll(gaussian_sol_data, "Eppstein",
                                                                                  "area") + "}; \\addlegendentry{$\\text{Eppstein}_{\\infty}$}")

print("6) Diameter")

for i, d in enumerate(diameters):
    print("\\addplot+[color=" + colours[i] + ", mark=" + marks[i] + ",mark options={rotate=" + rotate[
        i] + "}, error bars/.cd, y dir=both, y explicit] coordinates {" + std_unroll(gaussian_sol_data,
                                                                                     "Antipodal_" + str(d),
                                                                                     "diameter") + "}; \\addlegendentry{$\\text{Antipodal}_{" + str(
        d) + "}$}")

print("\\addplot+[color=" + colours[-1] + ", mark=" + marks[-1] + ",mark options={rotate=" + rotate[
    -1] + "}, error bars/.cd, y dir=both, y explicit] coordinates {" + std_unroll(gaussian_sol_data, "Eppstein",
                                                                                  "diameter") + "}; \\addlegendentry{$\\text{Eppstein}_{\\infty}$}")
