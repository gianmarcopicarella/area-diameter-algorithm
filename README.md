# Thesis Title: Finding dense and well-shaped convex clusters in 2d point sets
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

Run the test suite and benchmark executables.
```bash
$ ./src/test
$ ./src/bench
```

Plot the solutions found during the benchmarks and save them to the data/plots folder (it takes ≈ 10 minutes on my machine)
```bash
$ cd ../python
$ pip install -r requirements.txt
$ python generate_plots.py
```

Generate LateX plots about the data collected during the benchmarking phase
```bash
  $ python latex.py
```

Run the program with the 10 biggest real-world data points (2hrs maximum per file)
```bash
$ python run_with_real_data.py
```

## Generate the data yourself

The data used for the experiments and real world runs are already available in data/samples. If you want to generate the data from scratch, then you should download [detections_subset.json](https://drive.google.com/file/d/1aHM7tw1oLBKeqv6VaCwpLoY8x4KPVu5i/view?usp=drive_link) to the folder data/raw and then run

```bash
  $ python generate_samples.py
```


## Tech Stack
CMake, Python, C++
