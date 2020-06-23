# CactuarF
1D First-order Godunov code - Fortran version

## The Euler Equations
The aim of CactaurF is to provide a simple numerical solver for the Euler equations in 1 dimension. In conservative form, the Euler equations are:

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{U}_t %2B \mathbf{F}(\mathbf{U})_x = 0">

Where:

<img src="https://render.githubusercontent.com/render/math?\mathbf{U} = \begin{bmatrix}\rho\\\rho u\\E\end{bmatrix} , \mathbf{F} = \begin{bmatrix}\rho u\\\rho u^2 %2B p\\u(E%2Bp)\end{bmatrix}">

# Test an equation
<img src="https://render.githubusercontent.com/render/math?math=e^{i %2B\pi} =x%2B1">

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{U}_i^{n%2B1} =  \frac{\mathrm{d} t}{\mathrm{d} x}  \left ( \mathbf{F}_i^n - \mathbf{F}_{i%2B1}^n \right )">
