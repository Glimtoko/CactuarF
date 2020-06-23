# CactuarF
1D First-order Godunov code - Fortran version

## The Euler Equations
The aim of CactaurF is to provide a simple numerical solver for the Euler equations in 1 dimension. In conservative form, the Euler equations are:

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{U}_t %2B \mathbf{F}(\mathbf{U})_x = 0">

Where:

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{U} = \begin{bmatrix}\rho\\\rho u\\E\end{bmatrix}">

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{F} = \begin{bmatrix}\rho u\\\rho u^2 %2B p\\u(E%2Bp)\end{bmatrix}">

## The Godunov Method
The Godunov method allows us to discretise the conservative equations above to produce a mechanism to step through time to reach a solution. For a given cell at index i, the update to go from timestep n to timestep n+1 is:

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{U}_i^{n%2B1} =  \frac{\Delta t}{\Delta x}  \left ( \mathbf{F}_{i-\frac{1}{2}}^n - \mathbf{F}_{{i%2B\frac{1}{2}}^n \right )">

