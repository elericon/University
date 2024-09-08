The files in this folder are used to solve the 2D Navier-Stokes equations for the lid driven cavity with fixed boundary conditions.
The code is a first order solver with explicit treatment of convection, thus giving a poisson equation to be solved at each time step.
It does not support the correct solution for Re > 8500 as it needs a second order solution.
To function it needs the files in the general purpose folder as it uses sparse matrices and sparse linear solvers for the solution 
of the poisson equation.
