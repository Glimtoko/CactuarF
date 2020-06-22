program cactaurF
use iso_fortran_env, only: int32, real64
use riemann_solvers
use flux_functions
implicit none

integer(kind=int32) :: ncells       ! Number of cells in problem
real(kind=real64) :: L = 1.0        ! Length of domain
real(kind=real64) :: x0 = 0.5       ! Position of interface/membrane
real(kind=real64) :: gamma = 1.4    ! EoS parameter, ratio of specific heats
real(kind=real64) :: dtmax = 0.1

! Mesh data
real(kind=real64), dimension(:), allocatable :: x   ! Mesh coordinates

! Conserved quantities
real(kind=real64), dimension(:), allocatable :: density
real(kind=real64), dimension(:), allocatable :: momentum
real(kind=real64), dimension(:), allocatable :: energy

! Primitive quantities
real(kind=real64), dimension(:), allocatable :: pressure
real(kind=real64), dimension(:), allocatable :: velocity

! Fluxes
real(kind=real64), dimension(:), allocatable :: density_f
real(kind=real64), dimension(:), allocatable :: momentum_f
real(kind=real64), dimension(:), allocatable :: energy_f

! Other lengths etc.
real(kind=real64) :: dx     ! Cell size
real(kind=real64) :: ein    ! Internal energy for energy calculation
integer(kind=int32) :: i    ! General loop parameters

! Used in main loop
real(kind=real64) :: t      ! Current time
real(kind=real64) :: a      ! Sound speed
real(kind=real64) :: S      ! Characteristic wave speed
real(kind=real64) :: CFL    ! CFL parameter (0 < CFL <= 1.0)
real(kind=real64) :: dt     ! Timestep from CFL condition
real(kind=real64) :: dtdx   ! dt/dx, used in Godunov update
integer(kind=int32) :: step ! Current step

! Output
integer(kind=int32) :: fnum     ! File number

! Problem setup
real(kind=real64) :: rhoL, uL, PL       ! Left hand state
real(kind=real64) :: rhoR, uR, PR       ! Right hand state
real(kind=real64) :: t_end              ! End time
real(kind=real64) :: xupper             ! Cell boundary. Used in setup

procedure(riemann_API), pointer :: model

model => riemann_TSRS
ncells = 250
CFL = 0.6

! Allocate arrays from 0 to ncells+1 to allow for ghosts/boundaries
allocate(x(0:ncells+1)); x = 0.0

allocate(density(0:ncells+1)); density = 0.0
allocate(momentum(0:ncells+1)); momentum = 0.0
allocate(energy(0:ncells+1)); energy = 0.0

allocate(density_f(ncells)); density_f = 0.0
allocate(momentum_f(ncells)); momentum_f = 0.0
allocate(energy_f(ncells)); energy_f = 0.0

allocate(pressure(0:ncells+1)); pressure = 0.0
allocate(velocity(0:ncells+1)); velocity = 0.0

! Initial test - set Sod problem
uL = 0.0
uR = 0.0
rhoL = 1.0
rhoR = 0.125
PL = 1.0
PR = 0.1
t_end = 0.25

! Set cell size (uniform mesh)
dx = L/ncells
write(*,'("Number of cells = ",i4," => cell size = ",f7.3,"cm")') ncells, dx

! Set coordinates
do i = 1, ncells+1
    x(i) = (i-0.5) * dx
end do

! Set initial density and pressure fields
do i = 1, ncells+1
    xupper = x(i) + dx/2.0
    if (xupper <= x0) then
        density(i) = rhoL
        momentum(i) = rhoL*uL
        pressure(i) = PL
        velocity(i) = uL
    else
        density(i) = rhoR
        momentum(i) = rhoR*uR
        pressure(i) = PR
        velocity(i) = uR
    end if
    ein = pressure(i)/((gamma - 1.0)*density(i))
    energy(i) = density(i)*(0.5*velocity(i)*velocity(i) + ein)
end do

! Set boundary cells - will need to change a little for parallel setup
density(0) = density(1)
momentum(0) = momentum(1)
pressure(0) = pressure(1)
velocity(0) = velocity(1)
energy(0) = energy(1)

density(ncells+1) = density(ncells)
momentum(ncells+1) = momentum(ncells)
pressure(ncells+1) = pressure(ncells)
velocity(ncells+1) = velocity(ncells)
energy(ncells+1) = energy(ncells)

t = 0.0
step = 0
do while (.true.)
    step = step + 1
    S = 0.0
    do i = 1, ncells+1
        a = sqrt((gamma*pressure(i))/density(i))
        S = max(S, a + velocity(i))
    end do

    dt = min(dtmax, cfl*dx/S)
    if (dt <= 0.0) then
        write(*,'("ERROR: dt = ",es9.3)') dt
        stop
    end if
    write(*,'("Step: ",i5,", t:",f7.3", dt:",es13.6)') step, t, dt

    ! Get fluxes - need to include right-hand boundary cell
    do i = 1, ncells+2
        call get_flux_from_sample( &
            velocity(i-1), density(i-1), pressure(i-1), &
            velocity(i), density(i), pressure(i), gamma, model, &
            density_f(i), momentum_f(i), energy_f(i) &
        )
    end do

    dtdx = dt/dx

    ! Perform Godunov update
    do i = 1, ncells
        ! Conservative update
        density(i) = density(i) + dtdx*(density_f(i) - density_f(i+1))
        momentum(i) = momentum(i) + dtdx*(momentum_f(i) - momentum_f(i+1))
        energy(i) = energy(i) + dtdx*(energy_f(i) - energy_f(i+1))

        ! Primitive update
        pressure(i) = (gamma - 1.0)*(energy(i) - 0.5*density(i)*velocity(i)*velocity(i))
        velocity(i) = momentum(i)/density(i)

        if (pressure(i) < 0.0) then
            write(*,'("Negative pressure in cell ",i5)') i
            write(*,'("Pressure = ",f7.3)') pressure(i)
            write(*,'("Energy = ",f7.3)') energy(i)
            write(*,'("Density = ",f7.3)') density(i)
            write(*,'("Velocity = ",f7.3)') velocity(i)
            stop
        end if
    end do


    if (t > t_end) exit
    t = t + dt
end do

! Output - Will need to formalise
open(file="output.dat", newunit=fnum, status="replace")
do i = 1, ncells+1
    ein  = energy(i)/density(i) - 0.5*velocity(i)*velocity(i)
    write(fnum,*) i, x(i), density(i), pressure(i), velocity(i), ein
end do


end program cactaurF
