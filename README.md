# Computing Largest Subsets of Points Whose Convex Hulls have Bounded Area and Diameter

## Authors

### Gianmarco Picarella
Department of Information and Computing Sciences, Utrecht University, The Netherlands
### Marc J. van Kreveld
Department of Information and Computing Sciences, Utrecht University, The Netherlands
### Frank Staals
Department of Information and Computing Sciences, Utrecht University, The Netherlands
### Sjoerd de Vries
Department of Information and Computing Sciences, Utrecht University, The Netherlands

Department of Digital Health, University Medical Center Utrecht, The Netherlands


## Abstract

We study the problem of computing a convex region with bounded area and diameter that contains the maximum number of points from a given point set $P$. We show that this problem can be solved in $O(n^6k)$ time and $O(n^3k)$ space, where $n$ is the size of $P$ and $k$ is the maximum number of points in the found region. 

We experimentally compare this new algorithm with an existing algorithm that does the same but without the diameter constraint, which runs in $O(n^3k)$ time. For the new algorithm, we use different diameters. 

We use both synthetic data and data from an application in cancer detection, which motivated our research.

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
$ git clone --recurse-submodules [REPO_URL]
```

Build the project in release mode using CMake

```sh
$ cd [REPO_NAME]
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

## Data Generation

The dataset used for experiments and real-world runs is available in `data/samples`. The raw version of the real-world data is available [here](https://drive.google.com/file/d/1aHM7tw1oLBKeqv6VaCwpLoY8x4KPVu5i/view?usp=drive_link). Move the file "detections_subset.json" in the `data/raw` directory. Finally, run the data generation and postprocessing script:

   ```bash
   $ python generate_input_samples.py
   ```

We will provide the raw version of the real-world data after receiving official approval from UMC Utrecht.
