# Thesis Title: Finding dense and well-shaped convex clusters in 2d point sets

## Introduction

Given a point set $P$ in $d$-dimensional space, a discrete convex polygon (DCP) in $P$ is a convex polygon with vertices $V\subseteq P$. The convex hull of a point set is always a DCP. Interestingly, many real-world problems arising in different disciplines (e.g., medicine and the oil industry) can be easily modeled as optimization problems, where a DCP maximizing (or minimizing) a certain objective function must be found. These problems may require taking into consideration additional constraints, usually defined as scalar values limiting specific properties of the polygon (e.g. the maximum allowed DCP area is bounded by $a_{\text{max}}$). The geometric knapsack problem is a generalization of these problems, which has been extensively studied for a variety of different objective functions and constraints.

We now introduce the notion of diameter. Let $P$ be a set of $2$-dimensional points. The diameter $d$ of $P$ is the maximum Euclidean distance between two points $p, q \in P$. The diameter of any convex polygon is always equivalent to its major axis length. 

In this master's thesis project, we focus on a variation of the geometric knapsack problem that, to the best of our knowledge, has not been studied in existing research. The problem is defined as follows.

### Problem Formulation
Given a set of $2$-dimensional points $P$ find a subset $S\subseteq P$ with $|S|>2$ such that the convex hull of $S$ satisfies the following properties:

1. The diameter of $S$ does not exceed $d_{\text{max}}$.
2. The area defined by the convex hull of $S$ does not exceed $a_{\text{max}}$.
3. $S$ has maximum cardinality.

The formulation of this problem comes from the generalization of a real-world need in the field of oncology, specifically in breast cancer prognosis.

### Algorithms in this repository
This repository provides an efficient implementation of [Eppstein's algorithm](https://link.springer.com/article/10.1007/BF02187823) and a novel algorithm designed to tackle the above-mentioned problem in polynomial time. Our new algorithm allows to find a convex region having bounded area and diameter while maximizing the number of enclosed points in O(kn^7) time and O(kn^5) space.

For more information about my project, please refer to the code in this repository and my [research proposal]().

## Getting Started

If you haven't done it yet, install clang == 19.1.6, CMake >= 3.18 and Python >= 3.10. After that, you can follow the following steps.

Clone the project
```bash
$ git clone --recurse-submodules https://github.com/gianmarcopicarella/master-thesis.git
```

Go to the project directory
```bash
$ cd master-thesis
```

Build the project in Release Mode using CMake
```bash
$ mkdir cmake-build-release
$ cd cmake-build-release
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make
```

Run the test suite
```bash
$ ./src/test
```

Run the benchmarks using synthetic and real data
```bash
$ ./src/bench
```

Plot and save the solutions found during the benchmarking phase
```bash
$ cd ../python
$ pip install -r requirements.txt
$ python generate_plots.py
```

Run the [Area-Selector](https://github.com/gianmarcopicarella/Area-Selector)'s algorithms with real data
```bash
$ cd ../external/Area-Selector
$ pip install -r requirements.txt
$ cd run
$ python run_area_selector.py
```

Generate the LaTeX code presenting the results obtained (benchmarking + comparison with [Area-Selector](https://github.com/gianmarcopicarella/Area-Selector)'s results with real data)
```bash
  $ python latex.py
```

An extensive set of settings which can be adjusted to fit your needs can be found in src/CMakeLists.txt. The default settings are the ones used during our experiments.

## Generate the data yourself

The data used for the experiments and real world runs are already available in data/samples. If you want to generate the data from scratch, then you can download [detections_subset.json](https://drive.google.com/file/d/1aHM7tw1oLBKeqv6VaCwpLoY8x4KPVu5i/view?usp=drive_link) to data/raw and run

```bash
  $ python generate_samples.py
```

The same script should be run if any settings affecting data generation is modified.


## Tech Stack
CMake, Python, C++
