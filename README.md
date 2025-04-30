# Master Thesis: Finding dense and well-shaped convex clusters in 2d point sets
Supervisors: Marc J. van Kreveld, Frank Staals and Sjoerd de Vries.
Department of Information and Computing Sciences, Utrecht University, The Netherlands.

## Abstract

We address the problem of computing a convex region with bounded area and diameter, enclosing the maximum number of points from a planar point set $P$.

This problem originated from previous work on an automatic pipeline for computing the mitotic count (MC) from histological images. The MC counts mitotic cells within a specific region and serves as a key indicator for determining cancer patient risk levels and treatment regimes. In this context, mitotic cells appear as points in a plane, and we need to identify the optimal mitotic hotspot region---a convex region that is not too elongated (bounded diameter), has a bounded area and contains the maximum number of points. 

We developed a novel dynamic programming algorithm, the Area-Diameter (AD) algorithm, which solves this problem in $O(n^6k)$ time and $O(n^3k)$ space, where $n$ is the size of $P$ and $k$ is the maximum number of enclosed points. To the best of our knowledge, this is the first polynomial time and exact algorithm to solve this problem.

We implemented the AD algorithm and an existing Area-only (A) algorithm that performs the same task without the diameter constraint in $O(n^3k)$ time and $O(n^2k)$ space. We experimentally compared their performance and solutions for three types of point sets: uniformly distributed sets with increasing size, normally distributed sets with increasing standard deviation and real-world medical data. For the AD algorithm, we test different diameters. 

Our results show that despite the higher worst-case complexity of the AD algorithm, the bound on the diameter enables significant pruning that often makes our algorithm practically faster than the A algorithm for moderately dense point sets. The performance of the A algorithm does not change with different point set distributions, but the solutions tend to be too elongated and beyond acceptable limits for medical applications. Finally, we show that our algorithm can process real-world medical datasets in a reasonable time, delivering exact solutions that are better---in terms of the number of enclosed points---than those found by existing methods.

For more information, please refer to my thesis [here]().

### Example output
![Alt Text](https://github.com/gianmarcopicarella/master-thesis/blob/59369825b71c7e77b649cc473bb48df3fcedce0f/data/example_areas.png)

The optimal convex regions found by our algorithm run with constraints $a_{\text{max}}=4\text{mm}, d_{\text{max}}=4.243\text{mm}$ and the $4$ most populated real-world point sets filtered with probability threshold $t=0.86$. The line spacing is set to $5$ $\text{mm}$.

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
