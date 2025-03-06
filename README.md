# Master Thesis: Finding dense and well-shaped convex clusters in 2d point sets

This repository contains the C++ and Python code developed for my Master's thesis at Utrecht University. The research conducted in this thesis explores the usage of several geometric algorithms for the identification of convex regions in Whole-Slide Images (WSI) containing the highest number of mitotic cells while being constrained on area and diameter. We implemented Eppstein et al.'s algorithm which can find such convex regions if the constraint on the diameter is relaxed and a novel antipodal algorithm which is able to find such regions in polynomial time. Our algorithm runs in \$O(kn^7)\$ time and uses \$O(kn^5)\$ space, where \$n\$ is the size of the input set and \$k\$ is the number of points defining the convex hull of the optimal solution. Our experiments show that our algorithm's performance is indeed affected by the point set distribution and the maximum allowed diameter but still is able to process in practice real-world point sets of more than 900 points in less than 30 minutes. There is a lot of room for improvements, including code parallelization (for which our algorithm is a very good candidate) and heuristic search.

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

## Example Output





## License

This repository is licensed under the [MIT License](LICENSE).

## Contact

For any questions, feel free to reach out to me at [[g.picarella@students.uu.nl](mailto\:g.picarella@students.uu.nl)] or via GitHub Issues.
