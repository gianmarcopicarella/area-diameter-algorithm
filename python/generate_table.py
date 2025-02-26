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

antipodal_data_results = utils.read_json(constants.PATH_TO_BENCHMARK_CUSTOM_REPORT)
area_selector_data_results = utils.read_json(constants.PATH_TO_AREA_SELECTOR_REPORT)


def fill_table():
    assert (area_selector_data_results["count"] == antipodal_data_results["count"])
    rows = ""

    for i in range(antipodal_data_results["count"]):

        apr = antipodal_data_results["results"][i][0]
        asr = area_selector_data_results["results"][i][0]

        print(apr["id"][apr["id"].index("/"):], asr["id"][asr["id"].index("/"):])


fill_table()
