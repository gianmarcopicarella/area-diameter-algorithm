import math

import constants
import utils

table = \
"""
\\begin{table}[h]
    \centering
    \\begin{tabular}{|c|c|c|c|c|c|c|c|c|}
        \hline
        \multicolumn{3}{|c|}{\textbf{Dataset}} & \multicolumn{3}{c|}{\textbf{Antipodal approach}} & \multicolumn{3}{c|}{\textbf{\cite{} approach}} \\
        \cline{1-9}
        Index & C & ANN & C & A ($\text{mm}^2$) & D (mm) & C & A ($\text{mm}^2$) & D (mm) \\
        \hline
        
        
        
        \hline
    \end{tabular}
    \caption{Comparison of Methods A and B on Different Datasets}
    \label{tab:results}
\end{table}
"""


data_runs = utils.read_json(constants.PATH_TO_BENCHMARK_BASE_REPORT)
data_results = utils.read_json(constants.PATH_TO_BENCHMARK_CUSTOM_REPORT)


def fill_table():
    assert(len(data_runs["benchmarks"]) == data_results["count"])
    rows = ""

    for i in range(len(data_runs["benchmarks"])):
        prog_info = data_runs["benchmarks"][i]
        sol_info = data_results["results"][i][0]

        # int(math.ceil(utils.to_seconds(prog_info["cpu_time"])))

        rows += ""
        print(prog_info, sol_info)

fill_table()