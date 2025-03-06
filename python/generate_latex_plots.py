import os

import numpy as np

import json
import constants
import utils

style_settings = {
    "Antipodal": {
        "0": "color=Cyan,mark=diamond,error bars/.cd, y dir=both, y explicit",
        "1": "color=Magenta,mark=o,error bars/.cd, y dir=both, y explicit",
        "2": "color=Green,mark=asterisk,error bars/.cd, y dir=both, y explicit",
        "3": "color=Peach,mark=triangle,mark options={rotate=180},error bars/.cd, y dir=both, y explicit",
        "4": "color=Blue,mark=+,error bars/.cd, y dir=both, y explicit"
    },
    "Eppstein": {
        "None": "color=Fuchsia,mark=triangle,error bars/.cd, y dir=both, y explicit",
        "*": "color=Fuchsia,mark=triangle,dashed,error bars/.cd, y dir=both, y explicit"
    }
}

data_runs = utils.read_json(constants.PATH_TO_BENCHMARK_RUNS)
data_results = utils.read_json(constants.PATH_TO_BENCHMARK_RESULTS)

input_data = {"run_data": data_runs, "solutions_data": data_results}
data = {"run_data": dict(), "solutions_data": dict()}


def process_benchmark_data():
    for b in input_data["run_data"]["benchmarks"]:
        run = b["name"].split("iterations")[0]
        components = run.split("/")

        distribution = components[1]
        algorithm = components[0]
        diameter = "None" if components[3] == "" else components[3]
        index = int(components[2])

        if distribution not in data["run_data"]:
            data["run_data"][distribution] = dict()
        if algorithm not in data["run_data"][distribution]:
            data["run_data"][distribution][algorithm] = dict()
        if diameter not in data["run_data"][distribution][algorithm]:
            data["run_data"][distribution][algorithm][diameter] = \
                {"time": [], "memory": [], "entries": [], "required_entries": []}

        entry = data["run_data"][distribution][algorithm][diameter]
        entry["time"].append((index, float(b["cpu_time"]), 0))
        entry["memory"].append((index, utils.to_kb(float(b["mem_avg"])), utils.to_kb(float(b["mem_std"]))))
        entry["entries"].append((index, float(b["entries_avg"]), float(b["entries_std"])))
        entry["required_entries"].append((index, float(b["min_entries_avg"]), float(b["min_entries_std"])))

    for results in input_data["solutions_data"]["results"]:
        first_run = results[0]["name"].split("iterations")[0]
        components = first_run.split("/")
        distribution = components[1]
        algorithm = components[0]
        diameter = "None" if components[3] == "" else components[3]
        index = int(components[2])

        if distribution not in data["solutions_data"]:
            data["solutions_data"][distribution] = dict()
        if algorithm not in data["solutions_data"][distribution]:
            data["solutions_data"][distribution][algorithm] = dict()
        if diameter not in data["solutions_data"][distribution][algorithm]:
            data["solutions_data"][distribution][algorithm][diameter] = \
                {"count": [], "area": [], "diameter": [], "indices": []}

        entry = data["solutions_data"][distribution][algorithm][diameter]
        counts, areas, diameters = [], [], []
        for r in results:
            iteration = int(r["name"].split(":")[1])

            assert ("convex_area" in r)
            convex_area = r["convex_area"]

            counts.append(convex_area["count"])
            areas.append(convex_area["area"])

            path_to_sample = os.path.join(constants.PATH_TO_EXPERIMENTS, distribution.lower(),
                                          str(index), f"points_{iteration}.json")
            with open(path_to_sample, 'r') as file:
                points = json.load(file)

            fi = int(convex_area["diameter_indices"][0])
            si = int(convex_area["diameter_indices"][1])
            fdp = np.array((points["points"][fi]["x"], points["points"][fi]["y"]), dtype=np.float32)
            sdp = np.array((points["points"][si]["x"], points["points"][si]["y"]), dtype=np.float32)

            diameters.append(float(np.linalg.norm(fdp - sdp)))

            entry["indices"].append((index, convex_area["hull_indices"]))

        entry["count"].append((index, float(np.mean(counts)), float(np.std(counts))))
        entry["area"].append((index, float(np.mean(areas)), float(np.std(areas))))
        entry["diameter"].append((index, float(np.mean(diameters)), float(np.std(diameters))))


def generate_latex_plot(title, x_label, y_label, y_max, y_mode, x_values, rows, legend_style_key="north west"):
    legend_styles = {"north west": r"legend style={font=\footnotesize,at={(-0.0006,1.0009)},anchor=north west}",
                     "north east": r"legend style={font=\footnotesize,at={(1.0006,1.0009)},anchor=north east}"}
    plot_template = r"""
\begin{center}
\begin{tikzpicture}
\begin{axis}[
    <LEGEND_STYLE>,
    <TITLE>,
    width=12cm, height=9cm,
    <X_LABEL>,
    <Y_LABEL>,
    grid=both,
    mark options={solid},
    <Y_MAX>,
    <Y_MODE>,
    xtick={<X_VALUES>},
    xticklabels={<X_VALUES>},
    xtick pos=bottom,
    ytick pos=left,
    xtick align=outside,
    tick style={thick},
    ytick align=outside
]

<ADD_ROWS_HERE>

\end{axis}
\end{tikzpicture}
\end{center}"""
    plot = plot_template.replace("<LEGEND_STYLE>", legend_styles[legend_style_key])
    plot = plot.replace("<TITLE>", title)
    plot = plot.replace("<X_LABEL>", x_label)
    plot = plot.replace("<Y_LABEL>", y_label)
    plot = plot.replace("<Y_MAX>", y_max)
    plot = plot.replace("<Y_MODE>", y_mode)
    plot = plot.replace("<X_VALUES>", x_values)
    return plot.replace("<ADD_ROWS_HERE>", rows)


def unroll(sequence, dist_key):
    if dist_key == "Uniform":
        return " ".join(
            f"({round(constants.DENSITY_VALUES[t[0]], 2)}, {round(t[1], 2)}) +- (0, {round(t[2], 2)})" for t in
            sequence)
    elif dist_key == "Gaussian":
        return " ".join(
            f"({round(constants.STD_VALUES[t[0]], 2)}, {round(t[1], 2)}) +- (0, {round(t[2], 2)})" for t in sequence)
    else:
        assert False


def legend(algo_key, diam_key):
    if diam_key == "None" or diam_key == "*":
        assert (algo_key == "Eppstein")
        legend_txt = r"\addlegendentry{$\text{" + algo_key + r"}_{\infty}"
        if diam_key == "*": legend_txt += r"^{\star}"
        return legend_txt + "$};\n"
    else:
        return r"\addlegendentry{$\text{" + algo_key + "}_{" + str(
            constants.SYNTHETIC_BENCHMARK_DIAMETERS[int(diam_key)]) + "}$};\n"


def generate_rows(data_type, dist_key, metric_key):
    rows = ""
    entry = data[data_type][dist_key]
    for algo_key in entry:
        for diam_key in entry[algo_key]:
            assert (metric_key in entry[algo_key][diam_key])
            rows += r"\addplot[" + style_settings[algo_key][diam_key] + "] coordinates {" + unroll(
                entry[algo_key][diam_key][metric_key], dist_key) + "};\n"
            rows += legend(algo_key, diam_key)

    if metric_key == "entries" and "Eppstein" in entry:
        rows += r"\addplot[" + style_settings["Eppstein"]["*"] + "] coordinates {" + unroll(
            entry["Eppstein"]["None"]["required_entries"], dist_key) + "};\n"
        rows += legend("Eppstein", "*")

    return rows


process_benchmark_data()

plots = \
    r"""
\usepackage{pgfplots}
\usepackage{subcaption}
\usepackage{pgfplotstable}
\usepackage{amsmath}
\usepackage{comment}
\usepackage{graphicx}

\usepackage[dvipsnames]{xcolor}
\usepackage[a4paper, margin=1in]{geometry}

\pgfplotsset{compat=1.18}
\usepgfplotslibrary{statistics}
\usepgflibrary{plotmarks}

\usepgfplotslibrary{external} 
\tikzexternalize

\begin{document}""" + "\n"

plots += \
    r"""\section{Comparison using synthetic data}
\subsection{Uniform distribution}
\subsubsection{Time, memory and table entries}""" + "\n"

plots += generate_latex_plot("title={Uniform distribution, $100$ repetitions}",
                             "xlabel={Input size (points)}",
                             "ylabel={Time (ms)}", "ymax=10^6", "ymode=log",
                             ",".join(str(x) for x in constants.DENSITY_VALUES),
                             generate_rows("run_data", "Uniform", "time")) + "\n"

plots += generate_latex_plot("title={Uniform distribution, $100$ repetitions}",
                             "xlabel={Input size (points)}",
                             "ylabel={Memory allocation (kb)}", "ymax=10^10", "ymode=log",
                             ",".join(str(x) for x in constants.DENSITY_VALUES),
                             generate_rows("run_data", "Uniform", "memory")) + "\n"

plots += generate_latex_plot("title={Uniform distribution, $100$ repetitions}",
                             "xlabel={Input size (points)}",
                             "ylabel={Table entries}", "ymax=10^10", "ymode=log",
                             ",".join(str(x) for x in constants.DENSITY_VALUES),
                             generate_rows("run_data", "Uniform", "entries")) + "\n"

plots += r"\subsubsection{Enclosed points, area and diameter}" + "\n"

plots += generate_latex_plot("title={Uniform distribution, $100$ repetitions}",
                             "xlabel={Input size (points)}",
                             "ylabel={Cardinality (enclosed points)}", "", "ymode=linear",
                             ",".join(str(x) for x in constants.DENSITY_VALUES),
                             generate_rows("solutions_data", "Uniform", "count")) + "\n"

plots += generate_latex_plot("title={Uniform distribution, $100$ repetitions}",
                             "xlabel={Input size (points)}",
                             r"ylabel={Area ($\text{mm}^2$)}", "ymax=6", "ymode=linear",
                             ",".join(str(x) for x in constants.DENSITY_VALUES),
                             generate_rows("solutions_data", "Uniform", "area")) + "\n"

plots += generate_latex_plot("title={Uniform distribution, $100$ repetitions}",
                             "xlabel={Input size (points)}",
                             "ylabel={Diameter (mm)}", "ymax=30", "ymode=linear",
                             ",".join(str(x) for x in constants.DENSITY_VALUES),
                             generate_rows("solutions_data", "Uniform", "diameter")) + "\n"

plots += r"\subsection{Gaussian distribution}" + "\n" + r"\subsubsection{Time, memory and table entries}" + "\n"

plots += generate_latex_plot("title={Gaussian distribution, $100$ repetitions}",
                             "xlabel={Standard deviation}",
                             "ylabel={Time (ms)}", "ymax=10^6", "ymode=log",
                             ",".join(str(x) for x in constants.STD_VALUES),
                             generate_rows("run_data", "Gaussian", "time"), "north east") + "\n"

plots += generate_latex_plot("title={Gaussian distribution, $100$ repetitions}",
                             "xlabel={Standard deviation}",
                             "ylabel={Memory allocation (kb)}", "ymax=10^9", "ymode=log",
                             ",".join(str(x) for x in constants.STD_VALUES),
                             generate_rows("run_data", "Gaussian", "memory"), "north east") + "\n"

plots += generate_latex_plot("title={Gaussian distribution, $100$ repetitions}",
                             "xlabel={Standard deviation}",
                             "ylabel={Table entries}", "ymax=10^10", "ymode=log",
                             ",".join(str(x) for x in constants.STD_VALUES),
                             generate_rows("run_data", "Gaussian", "entries"), "north east") + "\n"

plots += r"\subsubsection{Enclosed points, area and diameter}" + "\n"

plots += generate_latex_plot("title={Gaussian distribution, $100$ repetitions}",
                             "xlabel={Standard deviation}",
                             "ylabel={Cardinality (enclosed points)}", "", "ymode=linear",
                             ",".join(str(x) for x in constants.STD_VALUES),
                             generate_rows("solutions_data", "Gaussian", "count"), "north east") + "\n"

plots += generate_latex_plot("title={Gaussian distribution, $100$ repetitions}",
                             "xlabel={Standard deviation}",
                             r"ylabel={Area ($\text{mm}^2$)}", "ymax=6", "ymode=linear",
                             ",".join(str(x) for x in constants.STD_VALUES),
                             generate_rows("solutions_data", "Gaussian", "area")) + "\n"

plots += generate_latex_plot("title={Gaussian distribution, $100$ repetitions}",
                             "xlabel={Standard deviation}",
                             "ylabel={Diameter (mm)}", "", "ymode=linear",
                             ",".join(str(x) for x in constants.STD_VALUES),
                             generate_rows("solutions_data", "Gaussian", "diameter")) + "\n"
plots += r"\end{document}"

print(plots)
