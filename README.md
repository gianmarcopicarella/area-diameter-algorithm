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

Run the [Area-Selector](https://github.com/sjoerd-de-vries/Area-Selector)'s algorithms on real world samples
```bash
$ cd ../external/Area-Selector/run
$ python run_area_selector.py
```

Generate the LaTeX code presenting the results obtained (benchmarking + comparison with [Area-Selector](https://github.com/sjoerd-de-vries/Area-Selector)'s results for real-world samples)
```bash
  $ python latex.py
```

An extensive set of settings which can be adjusted to fit your needs can be found in src/CMakeLists.txt. The default settings are the ones used during our experiments.

## Generate the data yourself

The data used for the experiments and real world runs are already available in data/samples. If you want to generate the data from scratch, then you can download [detections_subset.json](https://drive.google.com/file/d/1aHM7tw1oLBKeqv6VaCwpLoY8x4KPVu5i/view?usp=drive_link) to data/raw and run

```bash
  $ python generate_samples.py
```


## Tech Stack
CMake, Python, C++
