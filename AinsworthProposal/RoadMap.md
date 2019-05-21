# Mixed-FEM error estimator
This folder contains the code used to implement Mark Ainsworth proposal for
_a posteriori_ Mixed-FEM estimations.

The reference used is:
>  M. Ainsworth, X. Ma, Non-uniform order mixed FEM approximation: Implementation, post-processing,
> computable error bound and adaptivity. J. Comput. Physc. 231(2) (2012) 436-453.

As stated in the article, the estimation consists of reconstructing the potential which is done through the three following steps
(to be better explained in the next sections):

1. Solution of a local Neumann problem, using the flux given by the FEM as the boundary condition
2. Inter-element smoothing matching the solution obtained by _1_ along the mesh
3. Solution of a local Dirichlet problem, using the solution obtained in _2_ as the boundary condition

The solution obtained at the end of the process is compared to the pressure given by the Mixed-FEM. The bigger the deviation, the worse is the approximation at that element, suggesting it should be refined to achieve a better solution.

## Notes on the Implementation
I'm using this file as a sort of notebook to guide the work.

For now I can think of the following things I need to do:

- Create a class to manage the estimation
- Implement the solution of a Poisson model problem
- Store the multiphysics solution mesh into a member of the created class 
