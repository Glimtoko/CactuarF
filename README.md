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

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{U}_i^{n%2B1} =  \frac{\Delta t}{\Delta x}  \left ( \mathbf{F}_{i-\frac{1}{2}} - \mathbf{F}_{i%2B\frac{1}{2}} \right )">

Where the intercell numerical flux (F) is given by:

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{F}_{i%2B\frac{1}{2}} = \mathbf{F}(\mathbf{U}_{i%2B\frac{1}{2}}(0))">

To determine these fluxes, we construct Riemann problems on each cell boundary, and solve them at the point:

<img src="https://render.githubusercontent.com/render/math?math=S = \frac{x}{t} = 0">

Note that Riemann solvers generally work in terms of the vector of primitive variables:

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{W} =  \begin{bmatrix}\rho\\u\\p\end{bmatrix}">

We need to convert this vector into the vector of conserved quantities, F. For density and momemtum, this is trivial. For energy, we use:

<img src="https://render.githubusercontent.com/render/math?math=E = \rho(\frac{1}{2}u^2 %2B e)">

Where, for an **ideal gas**, e is given by:

<img src="https://render.githubusercontent.com/render/math?math=e = \frac{p}{\rho(\gamma - 1)}">

Where gamma is the ratio of specific heats. Note that this scheme is stable only where the timestep obeys the relation:

<img src="https://render.githubusercontent.com/render/math?math=\Delta T = \frac{C_{cfl}\Delta x}{S^n_{max}}">

Where:

<img src="https://render.githubusercontent.com/render/math?math=0 < C_{cfl} \le 1">

And S is the maximum wavespeed in the problem. For a 1D code like this, the maximum wavespeed is easy determined by averaging the interface wave speeds from the Riemann solutions. However, extended to multiple dimensions, such an averaging is unstable, so it is more common to use an estimate for S. In this code, the following estimate is used:

<img src="https://render.githubusercontent.com/render/math?math=S^n_{\mathrm{max}} = \mathrm{max} \{|u^n_i| %2B a^n_i \}">
