# Computing Largest Subsets of Points Whose Convex Hulls have Bounded Area and Diameter
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

The dataset used for experiments and real-world runs is available in `data/samples`. We will provide the raw version of the real-world data after receiving official approval from UMC Utrecht.
