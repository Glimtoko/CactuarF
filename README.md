# CactuarF
1D First-order Godunov code - Fortran version

# 1. Introduction
## The Euler Equations
The aim of CactaurF is to provide a simple numerical solver for the Euler equations in 1 dimension. In conservative form, the Euler equations are:

$\mathbf{U}_t + \mathbf{F}(\mathbf{U})_x = 0$

Where:

$$\mathbf{U} = \begin{bmatrix}\rho \\ \rho u\\ \E\end{bmatrix}$$

$$\mathbf{F} = \begin{bmatrix}\rho u \\ \rho u^2 + p\\ \u(E+p)\end{bmatrix}$$

## The Godunov Method
The Godunov method allows us to discretise the conservative equations above to produce a mechanism to step through time to reach a solution. For a given cell at index i, the update to go from timestep n to timestep n+1 is:

$\mathbf{U}_i^{n+1} =  \frac{\Delta t}{\Delta x}  \left ( \mathbf{F}_{i-\frac{1}{2}} - \mathbf{F}_{i+\frac{1}{2}} \right )$

Where the intercell numerical flux (F) is given by:

$\mathbf{F}_{i+\frac{1}{2}} = \mathbf{F}(\mathbf{U}_{i+\frac{1}{2}}(0))$

To determine these fluxes, we construct Riemann problems on each cell boundary, and solve them at the point:

$S = \frac{x}{t} = 0$

Note that Riemann solvers generally work in terms of the vector of primitive variables:

$\mathbf{W} =  \begin{bmatrix}\rho\\u\\p\end{bmatrix}$

We need to convert this vector into the vector of conserved quantities, F. For density and momemtum, this is trivial. For energy, we use:

$E = \rho(\frac{1}{2}u^2 + e)$

Where, for an **ideal gas**, e is given by:

$e = \frac{p}{\rho(\gamma - 1)}$

Where gamma is the ratio of specific heats. Note that this scheme is stable only where the timestep obeys the relation:

$\Delta T = \frac{C_{cfl}\Delta x}{S^n_{max}}$

Where:

$0 < C_{cfl} \le 1$

And S is the maximum wavespeed in the problem. For a 1D code like this, the maximum wavespeed is easy determined by averaging the interface wave speeds from the Riemann solutions. However, extended to multiple dimensions, such an averaging is unstable, so it is more common to use an estimate for S. In this code, the following estimate is used:

$S^n_{\mathrm{max}} = \mathrm{max} \{|u^n_i| + a^n_i \}$

Note that estimating S in this way can underestimate the wave speed, and thus create an unstabel timestep. This means that the range of values for C needs to be constrained below one. Testing has suggested that C in the following range is generally stable:

$0 < C_{cfl} \le 0.7$

## Boundary Conditions
Boundary conditions in a 1D Godunov code are easy to deal with. Consider the following mesh:

![Full 1D Mesh](/images/mesh1.png)

Here the cells L and R represent the ends of the regular mesh, so for example in the standard Sod problem setup the left hand edge of cell L would be at x=0.0cm, and the right hand edge of cell R would be at x=1.0cm. We then simply extend the mesh by one cell in each direction, to Lb and Rb. The state in these cells is kept constant, and is initialised at t=0 depending on the nature (reflective or transmissive) of the boundary condition. Currently, CactuarF only supports transmissive boundary conditions, defined as:

$\mathbf{U}(Lb) = \mathbf{U}(L)$

and,

$\mathbf{U}(Rb) = \mathbf{U}(R)$

These boundary cells are used to determine the left-hand and right-hand fluxes on cells L and R respectively.

# 2. The CactuarF code
## Parallelism
Parallelising a 1D code is extremely easy. First of all, we divide the mesh into N sections, where the first section (numbered 0 by convention) corresponds to the left-hand end of the mesh, and the Nth section (numbered N-1) corresponds to the right-hand end of the mesh. For example for N=3 we would have:

![Parallel 1D Mesh](/images/mesh2.png)

Here Lb and Rb correspond to boundary condition cells as before. We have gained new boundary cells labelled Lg and Rg, however, and these are known as "ghost" cells. They will work in a similar way to the boundary condition cells, so that, for example, on mesh section 1 the left-hand flux on L will be determined using Lg. Where things differ is how the values in cells Lg and Rg are determined. These cannot simply be set at time zero, as they correspond to information which changes as time progresses. Instead, they are set via communication from neighbouring processors. This is illustrated in the following image:

![Parallel 1D Comms](/images/mesh5.png)

I.e. on a given processor n, the value of Lg comes from cell R on processor (n-1), and the value of Rg comes from cell L on processor (n+1). This is a very simple scheme to encode in MPI, and requires a total of 6 sends and 6 receives on processors 1 < n < (N-2), and 3 sends and 3 receives on processors 0 and (N-1). These values come from the need to communicate the entire vector **W** of primitive values.

In CactaurF, parallel communication is handled in the *comms.f95* source file.
