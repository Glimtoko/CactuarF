module flux_functions
use iso_fortran_env, only: int32, real64
use riemann
implicit none

contains

subroutine get_flux_from_sample( &
    uL, rhoL, PL, &
    uR, rhoR, PR, &
    gamma, model, &
    rho_f, mom_f, E_f &
)
! Returns fluxes calculated via sampling a Riemann solution (exact or
! approximate) at S = x/t = 0.0

! Inputs
real(kind=real64), intent(in) :: uL, rhoL, PL   ! Left-hand state
real(kind=real64), intent(in) :: uR, rhoR, PR   ! Right-hand state
real(kind=real64), intent(in) :: gamma   ! EoS parameter
procedure(riemann_API), pointer, intent(in) :: model

! Outputs
real(kind=real64), intent(out) :: rho_f, mom_f, E_f


! Local data
real(kind=real64) :: rho, P, u, E

! Call Riemann solver
call solve( &
    uL, rhoL, PL, uR, rhoR, PR, gamma, model, &
    rho, P, u, E &
)

! Convert Riemann solution in terms of primitive variables into conserved fluxes
rho_f = rho*u
mom_f = rho*u*u + P
E_f = u*(E + P)

end subroutine

end module flux_functions
