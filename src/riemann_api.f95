module riemann_api_mod
! "Top-level" module containing routine to call and sample a Riemann solver's
! solution
use iso_fortran_env, only: int32, real64
implicit none

! API for all riemann solver functions
interface
    subroutine riemann_API( &
        uL, rhoL, PL, uR, rhoR, PR, gamma, &
        Pstar, ustar, rhoLstar, rhoRstar, exact &
    )
    use iso_fortran_env, only: real64
    ! Inputs
    real(kind=real64), intent(in) :: uL, rhoL, PL   ! Left-hand state
    real(kind=real64), intent(in) :: uR, rhoR, PR   ! Right-hand state
    real(kind=real64), intent(in) :: gamma   ! EoS parameter

    ! Outputs
    real(kind=real64), intent(out) :: Pstar, ustar, rhoLstar, rhoRstar
    logical, intent(out) :: exact
    end subroutine
end interface

end module riemann_api_mod
