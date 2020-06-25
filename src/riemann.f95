module riemann
! "Top-level" module containing routine to call and sample a Riemann solver's
! solution
use iso_fortran_env, only: int32, real64
use riemann_api_mod
use riemann_solvers
use riemann_exact
use riemann_sampler
implicit none

! ! API for all riemann solver functions
! interface
!     subroutine riemann_API( &
!         uL, rhoL, PL, uR, rhoR, PR, gamma, &
!         Pstar, ustar, rhoLstar, rhoRstar, exact &
!     )
!     use iso_fortran_env, only: real64
!     ! Inputs
!     real(kind=real64), intent(in) :: uL, rhoL, PL   ! Left-hand state
!     real(kind=real64), intent(in) :: uR, rhoR, PR   ! Right-hand state
!     real(kind=real64), intent(in) :: gamma   ! EoS parameter
!
!     ! Outputs
!     real(kind=real64), intent(out) :: Pstar, ustar, rhoLstar, rhoRstar
!     logical, intent(out) :: exact
!     end subroutine
! end interface

contains

subroutine solve( &
    uL, rhoL, PL, uR, rhoR, PR, gamma, model, &
    rho, P, u, E &
)
! Generate a solution in terms of primative variables for the Riemann problem
! represented by (uL, rhoL, PL) and (uR, rhoR, PR), using a provided Riemann
! solver (model).

! Inputs
real(kind=real64), intent(in) :: uL, rhoL, PL   ! Left-hand state
real(kind=real64), intent(in) :: uR, rhoR, PR   ! Right-hand state
real(kind=real64), intent(in) :: gamma   ! EoS parameter
procedure(riemann_API), pointer, intent(in) :: model

! Outputs
real(kind=real64), intent(out) :: rho, P, u, E

! Local data
real(kind=real64) :: Pstar, ustar       ! Output from Riemann solver
real(kind=real64) :: rhoLstar, rhoRstar ! Output from approximate Rieman solver
logical :: exact                        ! Is the Riemann solver exact?

real(kind=real64) :: ein                ! Internal energy for E calculation

! Get P, u and maybe rho from Riemann solver
call model(uL, rhoL, PL, uR, rhoR, PR, gamma, Pstar, ustar, rhoLstar, rhoRstar, exact)


! Sample the solution at S = x/t = 0.0
call sample( &
    Pstar, ustar, rhoLstar, rhoRstar, &
    uL, rhoL, PL, &
    uR, rhoR, PR, &
    gamma, exact, &
    rho, P, u &
)

! Calculate energy
ein = P/((gamma - 1.0)*rho)
E = rho * (0.5*u*u + ein)

end subroutine solve

end module
