# Accelerated Reeds-Shepp Planning

Solves the Reeds-Shepp path planning problem 11.5x faster on average than OMPL's implementation. 

Paper attachment TBD.

## Benchmarking Results

This table presents sample benchmarking results of our proposed planner against OMPL's implementation of the original Reeds-Shepp algorithm and against our implementation of [Desaulniers1995].

| Method                        | Average Time per State (µs) | Time Ratio to OMPL | Maximum Path Length Error | Average Path Length Error |
|-------------------------------|-----------------------------------------------|--------------------|---------------------------|---------------------------|
| Proposed                      | 0.106                                         | 11.7               | 1.7e-14             | 3.31e-16            |
| [Desaulniers1995]             | 0.298                                         | 4.16               | 1.76e-14            | 2.76e-16            |
| [OMPL]                        | 1.241                                         | 1                  | -                         | -                         |

Table: Benchmarking results of our proposed planner against OMPL's implementation of the original Reeds-Shepp algorithm and against our implementation of [Desaulniers1995].

[Desaulniers1995]: https://ieeexplore.ieee.org/document/478429
[OMPL]: https://ompl.kavrakilab.org/ReedsSheppStateSpace_8cpp_source.html

## Sample partition in 3D
3D partition of regions in the configuration space. Each color refers to a unique region where one identified path type is certainly the shortest. 1 million configurations $p_{f}$ spanning $\big[ [-100,100,], [-100,100], [-\pi,\pi)\big]$ for $p_{0} = (0,0,0)$.
![alt text](https://github.com/IbrahimSquared/accelerated-RS-planner/blob/master/samples/3D_cases_cropped.png) <br>

## How to Use

To use this project, follow these simple steps:

1. Clone the repository:
