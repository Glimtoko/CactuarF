module riemann_exact
! This module contains functions required for an iterative Riemann solver.
! This works by solving the equation:
!
! p[k] = p[k-1] - f(p[k-1])/f'(p[k-1])
!
! where p[k] is the kth iterate. This process is repeated until the relative
! pressure change:
!
! CHA = |p[k] - p[k-1]|/(0.5*|p[k] + p[k-1]|
!
! is less than a specified tolerance.

! The functions f and f' differ for shocks and rarefactions. Note that the forms
! of f and f' are too complex to type in text from here, but can be seen in the
! following functions.
use iso_fortran_env, only: int32, real64
use riemann_sampler
implicit none

! Tolerance for iterative solver
real(kind=real64), parameter :: TOL = 1.0e-6

! Maximum number of iterations
integer(kind=int32), parameter :: MAXITER = 1000

contains
! N.b. in all function and argument names, "k" refers to either the left or
! right side of the Riemann problem.
real(kind=real64) function fk(P, Pk, rhok, gamma)
! Function to calculate f. Calls the appropriate sub-function from fk_shock
! and fk_rarefaction (below).

! Input data
real(kind=real64), intent(in) :: P          ! Current pressure guess
real(kind=real64), intent(in) :: Pk, rhok   ! State on side k of problem
real(kind=real64), intent(in) :: gamma      ! EoS parameter

if (P >= Pk) then
    fk = fk_shock(P, Pk, rhok, gamma)
else
    fk = fk_rarefaction(P, Pk, rhok, gamma)
end if

end function fk


real(kind=real64) function fk_shock(P, Pk, rhok, gamma)
! Function to calculate f for a shock

! Input data
real(kind=real64), intent(in) :: P          ! Current pressure guess
real(kind=real64), intent(in) :: Pk, rhok   ! State on side k of problem
real(kind=real64), intent(in) :: gamma      ! EoS parameter

! Local variables
real(kind=real64) :: Ak, Bk     ! Parameters in equation for f

Ak = 2.0/((gamma + 1.0)*rhok)
Bk = ((gamma - 1.0)/(gamma + 1.0)) * Pk

fK_shock = (P - Pk) * sqrt(Ak/(P + Bk))

end function fk_shock


real(kind=real64) function fk_rarefaction(P, Pk, rhok, gamma)
! Function to calculate f for a rarefaction

! Input data
real(kind=real64), intent(in) :: P          ! Current pressure guess
real(kind=real64), intent(in) :: Pk, rhok   ! State on side k of problem
real(kind=real64), intent(in) :: gamma      ! EoS parameter

! Local variables
real(kind=real64) :: Ak, f     ! Parameters in equation for f

f = (gamma - 1.0)/(2.0*gamma)
ak = sqrt((gamma*Pk)/rhok)

fk_rarefaction = (2.0*ak)/(gamma - 1.0) * ((P/Pk)**f - 1.0)

end function fk_rarefaction


real(kind=real64) function fdashk(P, Pk, rhok, gamma)
! Function to calculate f'. Calls the appropriate sub-function from fdashk_shock
! and fdashk_rarefaction (below).

! Input data
real(kind=real64), intent(in) :: P          ! Current pressure guess
real(kind=real64), intent(in) :: Pk, rhok   ! State on side k of problem
real(kind=real64), intent(in) :: gamma      ! EoS parameter

if (P >= Pk) then
    fdashk = fdashk_shock(P, Pk, rhok, gamma)
else
    fdashk = fdashk_rarefaction(P, Pk, rhok, gamma)
end if

end function fdashk


real(kind=real64) function fdashk_shock(P, Pk, rhok, gamma)
! Function to calculate f' for a shock

! Input data
real(kind=real64), intent(in) :: P          ! Current pressure guess
real(kind=real64), intent(in) :: Pk, rhok   ! State on side k of problem
real(kind=real64), intent(in) :: gamma      ! EoS parameter

! Local variables
real(kind=real64) :: Ak, Bk     ! Parameters in equation for f'

Ak = 2.0/((gamma + 1.0)*rhok)
Bk = ((gamma - 1.0)/(gamma + 1.0)) * Pk

fdashk_shock = sqrt(Ak/(P + Bk)) * (1.0 - (P - Pk)/(2.0*(Bk + P)))

end function fdashk_shock


real(kind=real64) function fdashk_rarefaction(P, Pk, rhok, gamma)
! Function to calculate f' for a rarefaction

! Input data
real(kind=real64), intent(in) :: P          ! Current pressure guess
real(kind=real64), intent(in) :: Pk, rhok   ! State on side k of problem
real(kind=real64), intent(in) :: gamma      ! EoS parameter

! Local variables
real(kind=real64) :: Ak, f     ! Parameters in equation for f'

f = (gamma - 1.0)/(2.0*gamma)
ak = sqrt((gamma*Pk)/rhok)

fdashk_rarefaction = 1.0/(rhok*ak) * (P/Pk)**f

end function fdashk_rarefaction



real(kind=real64) function get_P_estimate(P0, PL, rhoL, uL, PR, rhoR, uR, gamma)
! Calculates an updated estimate for pressure in the star region

! Input data
real(kind=real64), intent(in) :: P0             ! Current pressure guess
real(kind=real64), intent(in) :: PL, rhoL, uL   ! State on left side of problem
real(kind=real64), intent(in) :: PR, rhoR, uR   ! State on right side of problem
real(kind=real64), intent(in) :: gamma      ! EoS parameter

! Local variables
real(kind=real64) :: du     ! Difference in velocities across problem
real(kind=real64) :: f      ! Value of f evaluated on both sides
real(kind=real64) :: fdash  ! Value of f' evaluated on both sides

! N.b.
! f = f(pL, ...) + f(pR, ...) + du
! f' = f'(pL, ...) + f'(pR, ...)

f = fk(P0, PL, rhoL, gamma) + fk(P0, PR, rhoR, gamma) + du
fdash = fdashk(P0, PL, rhoL, gamma) + fdashk(P0, PR, rhoR, gamma)

get_P_estimate = P0 - (f/fdash)

end function get_P_estimate


real(kind=real64) function get_pressure_change(Pnew, Pold)
! Calculates the relative pressure change for this iteration
real(kind=real64), intent(in) :: Pnew   ! Current pressure guess
real(kind=real64), intent(in) :: Pold   ! Previous pressure guess

get_pressure_change = abs(Pnew - Pold)/(0.5*(Pnew + Pold))

end function get_pressure_change


subroutine riemann_iterative( &
    uL, rhoL, PL, uR, rhoR, PR, gamma, &
    Pstar, ustar, rhoLstar, rhoRstar, exact &
)

! Inputs
real(kind=real64), intent(in) :: uL, rhoL, PL   ! Left-hand state
real(kind=real64), intent(in) :: uR, rhoR, PR   ! Right-hand state
real(kind=real64), intent(in) :: gamma   ! EoS parameter

! Outputs
real(kind=real64), intent(out) :: Pstar, ustar          ! Pressure and velocity
real(kind=real64), intent(out) :: rhoLstar, rhoRstar    ! Needed for API, will
                                                        ! always be zero
logical, intent(out) :: exact

! Local
real(kind=real64) :: aL, aR         ! Sound speeds, for initial guess
real(kind=real64) :: Pold, Pnew     ! Previous and current pressure guess
real(kind=real64) :: dP             ! Pressure difference for this iteration
integer(kind=int32) :: iteration    ! Current iteration number

exact = .true.

! Initial pressure guess
aL = sqrt((gamma*PL)/rhoL)
aR = sqrt((gamma*PR)/rhoR)
Pold = 0.5*(PL + PR) - 0.125*(uR - uL)*(rhoR + rhoL)*(aL + aR)

! Iterate to a solution (or until we get bored)
do
    iteration = iteration + 1
    Pnew = get_P_estimate(Pold, PL, rhoL, uL, PR, rhoR, uR, gamma)

    ! Insert test for negative pressure guess

    dP = get_pressure_change(Pnew, Pold)

    if (dP < TOL) exit
    if (iteration > MAXITER) then
        ! Abort here
        exit
    end if

    Pold = Pnew
end do

Pstar = Pnew

! Calculate u
ustar = 0.5*(uL + uR) + 0.5*(fk(Pnew, PR, rhoR, gamma) &
                      - fk(Pnew, PL, rhoL, gamma))

rhoLstar = 0.0
rhoRstar = 0.0

end subroutine riemann_iterative

end module riemann_exact
