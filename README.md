# Master Thesis: Finding dense and well-shaped convex clusters in 2d point sets
Supervisors: Marc J. van Kreveld, Frank Staals and Sjoerd de Vries.

Department of Information and Computing Sciences, Utrecht University, The Netherlands.

## Abstract

We address the problem of computing a convex region with bounded area and diameter, enclosing the maximum number of points from a planar point set $P$.

This problem originated from previous work on an automatic pipeline for computing the mitotic count (MC) from histological images. The MC counts mitotic cells within a specific region and serves as a key indicator for determining cancer patient risk levels and treatment regimes. In this context, mitotic cells appear as points in a plane, and we need to identify the optimal mitotic hotspot region---a convex region that is not too elongated (bounded diameter), has a bounded area and contains the maximum number of points. 

We developed a novel dynamic programming algorithm, the Area-Diameter (AD) algorithm, which solves this problem in $O(n^6k)$ time and $O(n^3k)$ space, where $n$ is the size of $P$ and $k$ is the maximum number of enclosed points. To the best of our knowledge, this is the first polynomial time and exact algorithm to solve this problem.

We implemented the AD algorithm and an existing Area-only (A) algorithm that performs the same task without the diameter constraint in $O(n^3k)$ time and $O(n^2k)$ space. We experimentally compared their performance and solutions for three types of point sets: uniformly distributed sets with increasing size, normally distributed sets with increasing standard deviation and real-world medical data. For the AD algorithm, we test different diameters. 

Our results show that despite the higher worst-case complexity of the AD algorithm, the bound on the diameter enables significant pruning that often makes our algorithm practically faster than the A algorithm for moderately dense point sets. The performance of the A algorithm does not change with different point set distributions, but the solutions tend to be too elongated and beyond acceptable limits for medical applications. Finally, we show that our algorithm can process real-world medical datasets in a reasonable time, delivering exact solutions that are better---in terms of the number of enclosed points---than those found by existing methods.

For more information, please refer to my thesis [here](https://github.com/gianmarcopicarella/master-thesis/blob/main/data/picarella-master_thesis_v4.pdf).

### Example output
![Alt Text](https://github.com/gianmarcopicarella/master-thesis/blob/f3332a024c38767f42cb01aae80b2cbe93c10f60/data/example-areas.png)

(a) Mitotic hotspots found by the $\text{AD}_4$ algorithm. (b) Mitotic hotspots found by the AS algorithm using $s=1.98$. (c) Mitotic hotspots found by the AS algorithm using $s=0.5$. (d) Mitotic hotspots found by the A algorithm using $s=0.5$. We report the results for the real-world point set with index $\text{I}=0$. The patch size is $3$ $\times$ $3$ mm. The line spacing is set to $5$ $\text{mm}$.

## Repository Structure

```
├── data/               # Input samples, reports and LaTeX data
├── external/           # Submodules used in the project
├── python/             # Python scripts used for data generation and reporting
└── src/                # C++ codebase
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

Most of these Python scripts generate LaTeX code to facilitate the visualization the experimental results.

## Generating Data

The datasets used in our experiments are available in `data/samples`. However, you can generate the datasets from scratch by following these steps:

1. Download [detections_subset.json](https://drive.google.com/file/d/1aHM7tw1oLBKeqv6VaCwpLoY8x4KPVu5i/view?usp=drive_link) and place it in the `data/raw` directory.
2. Run the data generation and postprocessing script:

   ```bash
   $ python generate_input_samples.py
   ```

## License

This repository is licensed under the [MIT License](LICENSE).

## Contact

For any questions, feel free to reach out to me at [[g.picarella@students.uu.nl](mailto\:g.picarella@students.uu.nl)] or via GitHub Issues.
