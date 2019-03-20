# **EAGO_Differential.jl**
## A Julia Package for Solving Deterministic Global Optimization Problems with Differential Constraints


## Using EAGO_Differential.jl
EAGO_Differential.jl is a work in progress and is currently not a registered Julia package. As such, it should be added or developed using the package url as detailed in the [Pkg documentation](https://julialang.github.io/Pkg.jl/v1/). It has been tested using the [this commit](https://github.com/PSORLab/EAGO.jl/commit/d1744d85991677eb3a9f215416c65c0e16bc55c0) of the EAGO master branch.

The package uses the EAGO architecture to formulate parametric ODE constrained optimization problems  as a series of block-sequential implicit function solution routines and constructs the relaxation using the theory presented in [1,2,3] and a forthcoming paper. This package exports 'solve_ode' function used to solve parametric ODEs of the above form as well as 'ImplicitODELowerEvaluator' and 'ImplicitODEUpperEvaluator' structures which are used to calculated the convex/concave relaxations and associated interval bounds used in the optimization algorithm.

## Forthcoming work

- Variable step-size algorithms.
- Solving parametric ODE constrained semi-infinite programs.
- Additional implicit approaches and mixed explicit-implicit approaches.

## Citing EAGO_Differential.jl

Work contained in this package and associated theory has been submitted for peer review. For now, please cite this package as software.

## Citations

- Stuber, M.D., Scott, J.K., and P.I. Barton. Convex and Concave Relaxations of Implicit Functions. Optimization Methods and Software. 30(3), 424-460, 2014.
- Scott, J.K., Stuber, M.D., and Barton, P.I. Generalized McCormick Relaxations. J Global Optim, 51:569-606, 2011
