# Accelerated Reeds-Shepp Planning

Solves the Reeds-Shepp path planning problem 15x faster on average than [OMPL]'s implementation.

Paper attachment TBD.

Reduces set of cases to consider down to 20 and introduces a new partition that speeds up the computation of the shortest path.

Requires SFML & OMPL. It has a demo in main.cpp that computes and draws paths between randomly generated start/end configurations.

Also contains an implemntation of ["An efficient algorithm to find a shortest path for a car-like robot"].

## Benchmarking Results

This table presents sample benchmarking results of our proposed planner against OMPL's implementation of the original Reeds-Shepp algorithm and against our implementation of [Desaulniers1995].

Performed on idle Linux Intel i9-13980HX.

| Method            | Average Time per State (µs) | Time Ratio to OMPL | Maximum Path Length Error | Average Path Length Error |
| ----------------- | --------------------------- | ------------------ | ------------------------- | ------------------------- |
| Proposed          | 0.0827                      | 15.12              | 9.82e-15                  | 3.28e-16                  |
| [Desaulniers1995] | 0.216                       | 5.79               | 5.82e-15                  | 2.77e-16                  |
| [OMPL]            | 1.25                        | 1                  | -                         | -                         |

Table: Benchmarking results of our proposed planner against OMPL's implementation of the original Reeds-Shepp algorithm and against our implementation of [Desaulniers1995].

[Desaulniers1995]: https://ieeexplore.ieee.org/document/478429
[OMPL]: https://ompl.kavrakilab.org/ReedsSheppStateSpace_8cpp_source.html
["An efficient algorithm to find a shortest path for a car-like robot"]: https://ieeexplore.ieee.org/document/478429

## Sample partition in 3D

3D partition of regions in the configuration space. Each color refers to a unique region (out of 20) where one identified path type is certainly the shortest. 1 million configurations $p_{f}$ spanning $\big[ [-100,100,], [-100,100], [-\pi,\pi)\big]$ for $p_{0} = (0,0,0)$ and $r = 20$. Red axis = $x$ axis, blue axis = $y$ axis, yellow axis = $\theta$ axis. <br>
![alt text](https://github.com/IbrahimSquared/accelerated-RS-planner/blob/master/samples/3D_cases_cropped.png) <br>

## Variability Based on Types

1e8 final configuration samples spanning $\big[ [-1000,1000,], [-1000,1000], [-\pi,\pi)\big]$ for $p_{0} = (0,0,0)$ and $r = 400$.
Computing those speedups for each path type incurs an overhead that makes the speedups slightly slower. This is only to show that the speedup is not uniform across all path types.
The weighted average speedup is around 11.87.

| Type | Occurrences | Speedup |
| ---- | ----------- | ------- |
| 1    | 6829417     | 11.5399 |
| 2    | 8684623     | 13.1266 |
| 3    | 9879506     | 12.3335 |
| 4    | 7630623     | 13.4698 |
| 5    | 6732552     | 12.8536 |
| 6    | 6683699     | 11.8223 |
| 7    | 791193      | 11.4474 |
| 8    | 2287070     | 11.8638 |
| 9    | 554434      | 11.7778 |
| 10   | 2529318     | 11.9787 |
| 11   | 1985832     | 10.639  |
| 12   | 11555       | 10.8744 |
| 13   | 7494149     | 10.3889 |
| 14   | 8013981     | 11.3507 |
| 15   | 7714423     | 9.52923 |
| 16   | 9776168     | 11.3254 |
| 17   | 1118135     | 10.2761 |
| 18   | 9897533     | 11.0233 |
| 19   | 138503      | 9.18825 |
| 20   | 1247286     | 8.42923 |

## How to Use (to be completed)

Set compiler path if needed. Make sure SFML: <br>
`sudo apt install libsfml-dev` <br>
and OMPL libraries: <br>
https://ompl.kavrakilab.org/installation.html <br>
are installed, then: <br>
`sudo apt install cmake` <br>
`sudo apt install g++` <br>
`mkdir build && cd build` <br>
`cmake ..` <br>
`make`

Troubleshooting build issues with OMPL: <br>
By default, OMPL installs in /usr/local/include/ompl-1.6/ompl. Create a symbolic link so that it can be properly located: <br>
`sudo ln -s /usr/local/include/ompl-1.6/ompl /usr/local/include/ompl` <br>
Similarly, create Eigen symbolic links so that required files can be found by OMPL:
`cd /usr/include` <br>
`sudo ln -sf eigen3/Eigen Eigen` <br>
`sudo ln -sf eigen3/unsupported unsupported` <br>
After making any changes to the library paths or creating symbolic links, run the ldconfig command to update the system's cache of shared libraries: <br>
`sudo ldconfig` <br>
After building, you can verify that the required libraries have been properly found using ldd, for example: <br>
`ldd ./acceleratedRSPlanner` <br>

To use this project, follow these simple steps: <br>
