# Master Thesis: Finding dense and well-shaped convex clusters in 2d point sets

<div style="display: flex; justify-content: space-between;">

  <img src="https://github.com/gianmarcopicarella/master-thesis/blob/e62b43d1df0f0b7de4b1550c9c84873b0940b9a3/data/ap_0.pdf" alt="Image 1" width="23%" />
  <img src="https://github.com/gianmarcopicarella/master-thesis/blob/e62b43d1df0f0b7de4b1550c9c84873b0940b9a3/data/ap_1.pdf" alt="Image 2" width="23%" />
  <img src="https://github.com/gianmarcopicarella/master-thesis/blob/e62b43d1df0f0b7de4b1550c9c84873b0940b9a3/data/ap_1.pdf" alt="Image 3" width="23%" />
  <img src="https://github.com/gianmarcopicarella/master-thesis/blob/e62b43d1df0f0b7de4b1550c9c84873b0940b9a3/data/ap_1.pdf" alt="Image 4" width="23%" />

</div>

## Introduction

This repository contains the C++ and Python code developed for my Master's thesis at Utrecht University. The research conducted in this thesis explores the usage of several geometric algorithms for the identification of convex regions in Whole-Slide Images (WSI) containing the highest number of mitotic cells while being constrained on area and diameter. We implemented Eppstein et al.'s algorithm which can find such convex regions if the constraint on the diameter is relaxed and a novel antipodal algorithm which is able to find such regions in polynomial time. Our algorithm runs in \$O(kn^7)\$ time and uses \$O(kn^5)\$ space, where \$n\$ is the size of the input set and \$k\$ is the number of points defining the convex hull of the optimal solution. Our experiments show that our algorithm's performance is indeed affected by the point set distribution and the maximum allowed diameter but still is able to process in practice real-world point sets of more than 900 points in less than 30 minutes. There is a lot of room for improvements, including code parallelization (for which our algorithm is a very good candidate) and heuristic search.

For detailed insights into the goals of this project, you can refer to my [research proposal](https://github.com/gianmarcopicarella/master-thesis/blob/cd705a7bf150f72d711a044cddcb1203e70f860c/data/research_proposal_gianmarcopicarella.pdf). A link to the final Master's thesis will be added here upon completion.

## Repository Structure

```
├── data/               # Input samples, reports and LaTeX data
├── external/           # Submodules used in the project
├── python/             # Python scripts used for data generation and reporting
├── src/                # C++ codebase 
└── README.md           # This file
```

## Requirements

Clang == 19.1.6, CMake >= 3.18 and Python >= 3.10.

## Getting Started

Clone the repository with all its submodules

```sh
$ git clone --recurse-submodules https://github.com/gianmarcopicarella/master-thesis.git
```

Build the project in release mode using CMake

```sh
$ cd master-thesis
$ mkdir cmake-build-release
$ cd cmake-build-release
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make
```

Run the test suite

```sh
$ ./src/test
```

Run the benchmarks using synthetic data
```bash
$ ./src/bench
```

Post-process the benchmark data and plot the convex areas found using synthetic data

```
$ cd ../python
$ pip install -r requirements.txt
$ python generate_latex_plots.py
$ python generate_solution_plots.py
```

Run the benchmark using real-world data
```bash
$ cd ../external/Area-Selector
$ pip install -r requirements.txt
$ python run/quality_of_results_comparison.py
```

Most of the provided Python scripts generate LaTeX code to facilitate the visualization of experimental results. The experiments include multiple variants of the antipodal algorithm (tested with maximum diameter values 2mm, 3mm, 4mm, 5mm and 6mm), Eppstein et al.'s algorithm (applied to the entire input set and individual image patches), and the area selector method (tested with patch step sizes of 1.98mm and 0.5mm).

## Generating Data

The dataset used for experiments and real-world runs is available in `data/samples`. However, if you wish to generate the dataset from scratch, follow these steps:

1. Download [detections_subset.json](https://drive.google.com/file/d/1aHM7tw1oLBKeqv6VaCwpLoY8x4KPVu5i/view?usp=drive_link) and place it in the `data/raw` directory.
2. Run the following command to generate the samples:

   ```bash
   $ python generate_input_samples.py
   ```

If any settings affecting data generation are modified, the script should be rerun to update the dataset accordingly.

## License

This repository is licensed under the [MIT License](LICENSE).

## Contact

For any questions, feel free to reach out to me at [[g.picarella@students.uu.nl](mailto\:g.picarella@students.uu.nl)] or via GitHub Issues.
