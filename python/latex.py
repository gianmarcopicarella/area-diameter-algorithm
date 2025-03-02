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
area_selector_data_results = utils.read_json(constants.PATH_TO_AREA_SELECTOR_RESULTS)

input_data = {
    "run_data": {"ours": data_runs, "area_selector": None},
    "solutions_data": {"ours": data_results, "area_selector": area_selector_data_results}
}

data = {"run_data": {"area_selector": dict(), "ours": dict()},
        "solutions_data": {"area_selector": dict(), "ours": dict()}}


def process_benchmark_data():
    for key in input_data["run_data"]:
        if input_data["run_data"][key] is None:
            continue
        for b in input_data["run_data"][key]["benchmarks"]:
            run = b["name"].split("iterations")[0]
            components = run.split("/")

            distribution = components[1]
            algorithm = components[0]
            diameter = "None" if components[3] == "" else components[3]
            index = int(components[2])

            if distribution not in data["run_data"][key]:
                data["run_data"][key][distribution] = dict()
            if algorithm not in data["run_data"][key][distribution]:
                data["run_data"][key][distribution][algorithm] = dict()
            if diameter not in data["run_data"][key][distribution][algorithm]:
                data["run_data"][key][distribution][algorithm][diameter] = \
                    {"time": [], "memory": [], "entries": [], "required_entries": []}

            entry = data["run_data"][key][distribution][algorithm][diameter]
            entry["time"].append((index, float(b["cpu_time"]), 0))
            entry["memory"].append((index, utils.to_kb(float(b["mem_avg"])), utils.to_kb(float(b["mem_std"]))))
            entry["entries"].append((index, float(b["entries_avg"]), float(b["entries_std"])))
            entry["required_entries"].append((index, float(b["min_entries_avg"]), float(b["min_entries_std"])))

    for key in input_data["solutions_data"]:
        if input_data["solutions_data"][key] is None:
            continue
        for results in input_data["solutions_data"][key]["results"]:
            first_run = results[0]["name"].split("iterations")[0]
            components = first_run.split("/")
            distribution = components[1]
            algorithm = components[0]
            diameter = "None" if components[3] == "" else components[3]
            index = int(components[2])

            if distribution not in data["solutions_data"][key]:
                data["solutions_data"][key][distribution] = dict()
            if algorithm not in data["solutions_data"][key][distribution]:
                data["solutions_data"][key][distribution][algorithm] = dict()
            if diameter not in data["solutions_data"][key][distribution][algorithm]:
                data["solutions_data"][key][distribution][algorithm][diameter] = \
                    {"count": [], "area": [], "diameter": [], "indices": []}

            entry = data["solutions_data"][key][distribution][algorithm][diameter]
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


def generate_latex_comparison_table():
    table = \
        r"""
\begin{table}
    \centering
    \begin{tabular}{|c|c|c|c|c|c|c|c|c|}
        \hline
        \multicolumn{3}{|c|}{\textbf{Dataset}} & \multicolumn{3}{c|}{\textbf{Antipodal approach}} & \multicolumn{3}{c|}{\textbf{Approach used in \cite{umcmodel}}} \\
        \cline{1-9}
        Index & C & ANN & C & A ($\text{mm}^2$) & D (mm) & C & A ($\text{mm}^2$) & D (mm) \\
        \hline

<ADD_ROWS_HERE>

        \hline
    \end{tabular}
    \caption{Comparison of the results from our antipodal approach with those from the method presented in \cite{umcmodel}.}
    \label{tab:comparison}
\end{table}
    """

    ours_solutions_data = data["solutions_data"]["ours"]["Real"]["Antipodal"]["0"]
    area_sel_solutions_data = data["solutions_data"]["area_selector"]["Real"]["AreaSelector"]["0"]

    assert (len(ours_solutions_data["count"]) == constants.REAL_BENCHMARKS_COUNT)
    assert (len(area_sel_solutions_data["count"]) == constants.REAL_BENCHMARKS_COUNT)

    rows = ""
    for i in range(constants.REAL_BENCHMARKS_COUNT):
        path_to_input_data = os.path.join(constants.PATH_TO_EXPERIMENTS, "real", str(i), "points_0.json")
        points = [(p['x'], p['y']) for p in utils.read_json(path_to_input_data)['points']]
        ann = utils.ann(points)
        apr_count, asr_count = int(ours_solutions_data["count"][i][1]), int(area_sel_solutions_data["count"][i][1])
        apr_area, asr_area = ours_solutions_data["area"][i][1], area_sel_solutions_data["area"][i][1]
        apr_diam, asr_diam = ours_solutions_data["diameter"][i][1], area_sel_solutions_data["diameter"][i][1]

        if apr_count > asr_count:
            apr_count = r"\textbf{" + str(apr_count) + "}"
        elif apr_count < asr_count:
            asr_count = r"\textbf{" + str(asr_count) + "}"
        else:
            apr_count = r"\underline{" + str(apr_count) + "}"
            asr_count = r"\underline{" + str(asr_count) + "}"

        rows += f"{i} & {len(points)} & ${round(ann[0], 2)} \\pm {round(ann[1], 2)}$ & " \
                f"{apr_count} & {round(apr_area, 2)} & {round(apr_diam, 2)} & " \
                f"{asr_count} & {round(asr_area, 2)} & {round(asr_diam, 2)}" + r" \\" + "\n"

    return table.replace("<ADD_ROWS_HERE>", rows)


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


def generate_real_data_latex_plots():
    base_plot = \
        r"""
\begin{figure}
    \centering
    \begin{subfigure}{0.45\textwidth}
        \centering
        \begin{tikzpicture}
            \begin{axis}[
                grid=major,
                width=8cm,
                height=8cm,
                axis equal
            ]
            \addplot[thick, <OUR_COLOR>, fill=<OUR_COLOR>!20, opacity=1]
            coordinates {<OUR_POLYGON_VERTICES>};
            \addplot[only marks, mark=*, color=<OUR_COLOR>, mark size=0.6pt]
            coordinates {<OUR_POLYGON_VERTICES>};
            \addplot[only marks, mark=*, color=black, mark size=0.6pt] table {<PATH_TO_DAT_FILE>};
            \end{axis}
        \end{tikzpicture}
        \caption{}
    \end{subfigure}
    \hfill
    \begin{subfigure}{0.45\textwidth}
        \centering
        \begin{tikzpicture}
            \begin{axis}[
                grid=major,
                width=8cm,
                height=8cm,
                axis equal
            ]
            \addplot[thick, <AREA_SELECTOR_COLOR>, fill=<AREA_SELECTOR_COLOR>!20, opacity=1]
            coordinates {<AREA_SELECTOR_POLYGON_VERTICES>};
            \addplot[only marks, mark=*, color=<AREA_SELECTOR_COLOR>, mark size=0.6pt]
            coordinates {<AREA_SELECTOR_POLYGON_VERTICES>};
            \addplot[only marks, mark=*, color=black, mark size=0.6pt] table {<PATH_TO_DAT_FILE>};
            \end{axis}
        \end{tikzpicture}
        \caption{}
    \end{subfigure}
    \caption{Convex area found by the antipodal algorithm (a) and the area selector (b) for point set $<POINT_SET_INDEX>$.}
\end{figure}
    """

    result = ""

    our_data = data["solutions_data"]["ours"]["Real"]["Antipodal"]["0"]
    area_selector_data = data["solutions_data"]["area_selector"]["Real"]["AreaSelector"]["0"]

    for i in range(constants.REAL_BENCHMARKS_COUNT):
        path_to_points = os.path.join(constants.PATH_TO_EXPERIMENTS, "real", str(i), "points_0.json")
        points = utils.read_json(path_to_points)
        dat_format_text = "\n".join(f"{p['x']} {p['y']}"
                                    for j, p in enumerate(points["points"]))
                                    # if j not in our_data["indices"][i][1] and j not in area_selector_data["indices"][i][1])
        path_to_dat = os.path.join(constants.PATH_TO_LATEX, "dat", "real", str(i))
        utils.prepare_path(path_to_dat)
        with open(os.path.join(path_to_dat, "points_0.dat"), "w") as file:
            file.write(dat_format_text)

        plot = base_plot.replace("<PATH_TO_DAT_FILE>", "dat/real/" + str(i) + "/points_0.dat")

        our_polygon_string = " ".join(f"({points['points'][j]['x']}, {points['points'][j]['y']})"
                                      for j in our_data["indices"][i][1] + [our_data["indices"][i][1][0]])
        plot = plot.replace("<OUR_POLYGON_VERTICES>", our_polygon_string)
        area_selector_polygon_string = " ".join(f"({points['points'][j]['x']}, {points['points'][j]['y']})"
                                                for j in area_selector_data["indices"][i][1] + [area_selector_data["indices"][i][1][0]])
        plot = plot.replace("<AREA_SELECTOR_POLYGON_VERTICES>", area_selector_polygon_string)
        plot = plot.replace("<OUR_COLOR>", "red")
        plot = plot.replace("<AREA_SELECTOR_COLOR>", "cyan")
        plot = plot.replace("<POINT_SET_INDEX>", str(i))

        result += plot + "\n"

    return result


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


def generate_rows(data_type, dist_key, metric_key, method="ours"):
    rows = ""
    entry = data[data_type][method][dist_key]
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

latex = \
    r"""
\documentclass{article}
\usepackage{pgfplots}
\usepackage{subcaption}
\usepackage{pgfplotstable}
\usepackage{amsmath}
\usepackage{comment}
\usepackage{graphicx}
\usepackage{array}
\usepackage{booktabs}
\usepackage{pdflscape}
\usepackage{tabularx}
\usepackage{multirow}

\usepackage[dvipsnames]{xcolor}
\usepackage[a4paper, margin=1in]{geometry}

\pgfplotsset{compat=1.18}
\usepgfplotslibrary{statistics}
\usepgflibrary{plotmarks}

\usepgfplotslibrary{external} 
\tikzexternalize

\begin{document}""" + "\n"

latex += r"\section{Uniform distribution}" + "\n" + r"\subsection{Time, memory and table entries}" + "\n"

latex += generate_latex_plot("title={Uniform distribution, $100$ repetitions}",
                             "xlabel={Input size (points)}",
                             "ylabel={Time (ms)}", "ymax=10^6", "ymode=log",
                             ",".join(str(x) for x in constants.DENSITY_VALUES),
                             generate_rows("run_data", "Uniform", "time")) + "\n"

latex += generate_latex_plot("title={Uniform distribution, $100$ repetitions}",
                             "xlabel={Input size (points)}",
                             "ylabel={Memory allocation (kb)}", "ymax=10^10", "ymode=log",
                             ",".join(str(x) for x in constants.DENSITY_VALUES),
                             generate_rows("run_data", "Uniform", "memory")) + "\n"

latex += generate_latex_plot("title={Uniform distribution, $100$ repetitions}",
                             "xlabel={Input size (points)}",
                             "ylabel={Table entries}", "ymax=10^10", "ymode=log",
                             ",".join(str(x) for x in constants.DENSITY_VALUES),
                             generate_rows("run_data", "Uniform", "entries")) + "\n"

latex += r"\subsection{Enclosed points, area and diameter}" + "\n"

latex += generate_latex_plot("title={Uniform distribution, $100$ repetitions}",
                             "xlabel={Input size (points)}",
                             "ylabel={Cardinality (enclosed points)}", "", "ymode=linear",
                             ",".join(str(x) for x in constants.DENSITY_VALUES),
                             generate_rows("solutions_data", "Uniform", "count")) + "\n"

latex += generate_latex_plot("title={Uniform distribution, $100$ repetitions}",
                             "xlabel={Input size (points)}",
                             r"ylabel={Area ($\text{mm}^2$)}", "ymax=6", "ymode=linear",
                             ",".join(str(x) for x in constants.DENSITY_VALUES),
                             generate_rows("solutions_data", "Uniform", "area")) + "\n"

latex += generate_latex_plot("title={Uniform distribution, $100$ repetitions}",
                             "xlabel={Input size (points)}",
                             "ylabel={Diameter (mm)}", "ymax=30", "ymode=linear",
                             ",".join(str(x) for x in constants.DENSITY_VALUES),
                             generate_rows("solutions_data", "Uniform", "diameter")) + "\n"

latex += r"\section{Gaussian distribution}" + "\n" + r"\subsection{Time, memory and table entries}" + "\n"

latex += generate_latex_plot("title={Gaussian distribution, $100$ repetitions}",
                             "xlabel={Standard deviation}",
                             "ylabel={Time (ms)}", "ymax=10^6", "ymode=log",
                             ",".join(str(x) for x in constants.STD_VALUES),
                             generate_rows("run_data", "Gaussian", "time"), "north east") + "\n"

latex += generate_latex_plot("title={Gaussian distribution, $100$ repetitions}",
                             "xlabel={Standard deviation}",
                             "ylabel={Memory allocation (kb)}", "ymax=10^9", "ymode=log",
                             ",".join(str(x) for x in constants.STD_VALUES),
                             generate_rows("run_data", "Gaussian", "memory"), "north east") + "\n"

latex += generate_latex_plot("title={Gaussian distribution, $100$ repetitions}",
                             "xlabel={Standard deviation}",
                             "ylabel={Table entries}", "ymax=10^10", "ymode=log",
                             ",".join(str(x) for x in constants.STD_VALUES),
                             generate_rows("run_data", "Gaussian", "entries"), "north east") + "\n"

latex += r"\subsection{Enclosed points, area and diameter}" + "\n"

latex += generate_latex_plot("title={Gaussian distribution, $100$ repetitions}",
                             "xlabel={Standard deviation}",
                             "ylabel={Cardinality (enclosed points)}", "", "ymode=linear",
                             ",".join(str(x) for x in constants.STD_VALUES),
                             generate_rows("solutions_data", "Gaussian", "count"), "north east") + "\n"

latex += generate_latex_plot("title={Gaussian distribution, $100$ repetitions}",
                             "xlabel={Standard deviation}",
                             r"ylabel={Area ($\text{mm}^2$)}", "ymax=6", "ymode=linear",
                             ",".join(str(x) for x in constants.STD_VALUES),
                             generate_rows("solutions_data", "Gaussian", "area")) + "\n"

latex += generate_latex_plot("title={Gaussian distribution, $100$ repetitions}",
                             "xlabel={Standard deviation}",
                             "ylabel={Diameter (mm)}", "", "ymode=linear",
                             ",".join(str(x) for x in constants.STD_VALUES),
                             generate_rows("solutions_data", "Gaussian", "diameter")) + "\n"

latex += r"\section{Comparison of the antipodal algorithm and the area selector in \cite{umcmodel} when applied to real data}" + "\n"
latex += generate_latex_comparison_table()
latex += generate_real_data_latex_plots()
latex += \
    r"""
\clearpage
\bibliography{sources}
\bibliographystyle{IEEEtran}
\end{document}
"""

with open(os.path.join(constants.PATH_TO_LATEX, "main.tex"), "w") as f:
    f.write(latex)
