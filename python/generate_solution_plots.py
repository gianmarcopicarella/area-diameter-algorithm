import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import constants
import utils


def extract_solutions_plots():
    data = utils.read_json(constants.PATH_TO_BENCHMARK_RESULTS)
    utils.prepare_path(constants.PATH_TO_PLOTS)

    for result in data["results"]:
        assert (len(result) > 0)

        split_id = result[0]["name"].split("/")
        algorithm = split_id[0].lower()
        distribution = split_id[1].lower()
        x_step = split_id[2]
        max_diameter = split_id[3]

        destination_path = os.path.join(constants.PATH_TO_PLOTS, algorithm, distribution, x_step)
        if len(split_id) == 5: destination_path = os.path.join(destination_path, max_diameter)
        if not os.path.isdir(destination_path): os.makedirs(destination_path)

        for iteration in range(len(result)):
            points_path = os.path.join(constants.PATH_TO_EXPERIMENTS, distribution, x_step, f"points_{iteration}.json")
            temp = utils.read_json(points_path)
            points = np.empty((temp["count"], 2), dtype=np.float64)
            for i, p in enumerate(temp["points"]):
                points[i, 0] = float(temp["points"][i]["x"])
                points[i, 1] = float(temp["points"][i]["y"])

            if "convex_area" in result[iteration]:
                convex_area = result[iteration]["convex_area"]
                indices = convex_area["hull_indices"]
                hull_points = [points[i] for i in indices]
                diameter = max([np.linalg.norm(points[i] - points[j])
                                for i in indices for j in indices])
                plt.title("count: " + str(round(convex_area["count"], 2)) +
                          ", area: " + str(round(convex_area["area"], 2)) +
                          ", diameter: " + str(round(diameter, 2)))
                plt.gca().add_patch(Polygon(hull_points, closed=True, edgecolor="red", facecolor='none', linewidth=1))
            else:
                plt.title("No solution found")
            plt.scatter(points[:, 0], points[:, 1], c="black", s=1)
            max_side = 0
            for i, p in enumerate(points):
                plt.text(p[0] - 0.03, p[1], f'{i}', fontsize=3, ha='right')
                max_side = max(max_side, p[0])
                max_side = max(max_side, p[1])
            plt.gcf().set_dpi(400)
            plt.gca().set_aspect("equal")
            plt.gca().set_xlim([0, max_side + 0.2])
            plt.gca().set_ylim([0, max_side + 0.2])
            plt.savefig(os.path.join(destination_path, f"result_{iteration}.png"), dpi=400)
            plt.clf()


extract_solutions_plots()
