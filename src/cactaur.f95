program cactuarF
use iso_fortran_env, only: int32, real64
use riemann_solvers
use flux_functions
use parallel_comms
use text_output

use mpi
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

! Parallelism
integer(kind=int32) :: idL, idR
integer(kind=int32) :: status   ! Status return from MPI
integer(kind=int32) :: nprocs   ! Number of processors
integer(kind=int32) :: rank     ! MPI rank
integer(kind=int32) :: ncells_per_proc  ! Number of cells per processor

! Runtime
real(kind=real64) :: time1, time2

! Output
integer(kind=int32) :: fnum     ! File number

! Problem setup
real(kind=real64) :: rhoL, uL, PL       ! Left hand state
real(kind=real64) :: rhoR, uR, PR       ! Right hand state
real(kind=real64) :: t_end              ! End time
real(kind=real64) :: xupper             ! Cell boundary. Used in setup

procedure(riemann_API), pointer :: model

model => riemann_TSRS
ncells = 30000
CFL = 0.6

! Initialise MPI
call MPI_Init(status)

! Get number of processors, and this processor's rank
call MPI_Comm_Size(MPI_COMM_WORLD, nprocs, status)
call MPI_Comm_Rank(MPI_COMM_WORLD, rank, status)
if (rank == 0) write(*,'("This is CactuarF running on ",i2," processors")') nprocs

! Get cells per processor, and ensure this is an integer
ncells_per_proc = ncells / nprocs
if (nprocs * ncells_per_proc /= ncells) then
    if (rank == 0) then
        write(*,'("ERROR: Number of cells does not divide by number of processors!")')
        write(*,'("NCELLS =,"i4)') ncells
        write(*,'("NPROCS =,"i4)') nprocs
        write(*,'("NCELLS_PER_PROC =,"i4)') ncells_per_proc
    end if
    call MPI_ABORT(MPI_COMM_WORLD, 0, status)
end if

! For now, serial only
idL = ncells_per_proc*rank + 1
idR = ncells_per_proc*(rank+1)

write(*,'("Processor ",i2,":: idL = ",i4,", idR = ",i4)') rank, idL, idR

! Allocate arrays from 0 to ncells+1 to allow for ghosts/boundaries
allocate(x(idL-1:idR+1)); x = 0.0

allocate(density(idL-1:idR+1)); density = 0.0
allocate(momentum(idL-1:idR+1)); momentum = 0.0
allocate(energy(idL-1:idR+1)); energy = 0.0

allocate(density_f(idL:idR)); density_f = 0.0
allocate(momentum_f(idL:idR)); momentum_f = 0.0
allocate(energy_f(idL:idR)); energy_f = 0.0

allocate(pressure(idL-1:idR+1)); pressure = 0.0
allocate(velocity(idL-1:idR+1)); velocity = 0.0

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
if (rank == 0) then
    write(*,'(/"Total number of cells = ",i4," => cell size = ",f7.3,"cm")') ncells, dx
end if

! Set coordinates
do i = idL, idR+1
    x(i) = (i-0.5) * dx
end do

! Set initial density and pressure fields
do i = idL, idR
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


! Set boundary cells - transmissive
if (idL == 1) then
    print *, rank, "LEFT BOUNDARY"
    density(idL-1) = density(idL)
    momentum(idL-1) = momentum(idL)
    pressure(idL-1) = pressure(idL)
    velocity(idL-1) = velocity(idL)
    energy(idL-1) = energy(idL)
end if

if (idR == ncells) then
    print *, rank, "RIGHT BOUNDARY"
    density(idR+1) = density(idR)
    momentum(idR+1) = momentum(idR)
    pressure(idR+1) = pressure(idR)
    velocity(idR+1) = velocity(idR)
    energy(idR+1) = energy(idR)
end if

t = 0.0
step = 0
time1 = MPI_WTIME()
do while (.true.)
    ! Update ghosts at start of timestep
    call parallel_update(density, pressure, velocity, rank, nprocs)

    step = step + 1
    S = 0.0
    do i = idL, idR+1
        a = sqrt((gamma*pressure(i))/density(i))
        S = max(S, a + velocity(i))
    end do

    dt = min(dtmax, cfl*dx/S)

    ! Communicate smallest timestep to all processors
    call MPI_ALLREDUCE(dt, dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, status)

    if (dt <= 0.0) then
        write(*,'("ERROR: dt = ",es9.3)') dt
        stop
    end if
    if (rank == 0) then
        write(*,'("Step: ",i5,", t:",f7.3", dt:",es13.6)') step, t, dt
    end if

    ! Get fluxes - need to include right-hand boundary cell
    do i = idL, idR+1
        call get_flux_from_sample( &
            velocity(i-1), density(i-1), pressure(i-1), &
            velocity(i), density(i), pressure(i), gamma, model, &
            density_f(i), momentum_f(i), energy_f(i) &
        )
    end do

    dtdx = dt/dx

    ! Perform Godunov update
    do i = idL, idR
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
time2 = MPI_WTIME()

if (rank == 0) then
    write(*,'("Completed ",i5, " steps in ",es13.6,"s")') step, time2 - time1
end if

! Final output
call do_text_output(x, density, pressure, velocity, energy, rank)

call MPI_FINALIZE(status)

end program cactuarF
