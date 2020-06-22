module riemann_solvers
use iso_fortran_env, only: int32, real64
use riemann_sampler
implicit none
! Since Fortran cannot elegantly return more than one value from a function, all
! of these routines will be subroutines, not functions

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

contains
subroutine riemann_PVRS1( &
    uL, rhoL, PL, uR, rhoR, PR, gamma, &
    Pstar, ustar, rhoLstar, rhoRstar, exact &
)

! Inputs
real(kind=real64), intent(in) :: uL, rhoL, PL   ! Left-hand state
real(kind=real64), intent(in) :: uR, rhoR, PR   ! Right-hand state
real(kind=real64), intent(in) :: gamma   ! EoS parameter

! Outputs
real(kind=real64), intent(out) :: Pstar, ustar, rhoLstar, rhoRstar
logical, intent(out) :: exact

! Local
real(kind=real64) :: aL, aR             ! Sound speeds
real(kind=real64) :: rho_bar, a_bar     ! Average density and sound speed

exact = .false.

aL = sqrt((gamma*PL)/rhoL)
aR = sqrt((gamma*PR)/rhoR)

rho_bar = 0.5*(rhoL + rhoR)
a_bar = 0.5*(aL + aR)

Pstar = 0.5*(PL + PR) + 0.5*(uL - uR)*(rho_bar*a_bar)
ustar = 0.5*(uL + uR) + 0.5*(PL - PR)/(rho_bar*a_bar)
rhoLstar = rhoL + (uL - ustar)*(rho_bar/a_bar)
rhoRstar = rhoR + (ustar - uR)*(rho_bar/a_bar)


end subroutine riemann_PVRS1


subroutine riemann_PVRS2( &
    uL, rhoL, PL, uR, rhoR, PR, gamma, &
    Pstar, ustar, rhoLstar, rhoRstar, exact &
)

! Inputs
real(kind=real64), intent(in) :: uL, rhoL, PL   ! Left-hand state
real(kind=real64), intent(in) :: uR, rhoR, PR   ! Right-hand state
real(kind=real64), intent(in) :: gamma   ! EoS parameter

! Outputs
real(kind=real64), intent(out) :: Pstar, ustar, rhoLstar, rhoRstar
logical, intent(out) :: exact

! Local
real(kind=real64) :: aL, aR             ! Sound speeds
real(kind=real64) :: CL, CR             ! Sound speed * density
real(kind=real64) :: rho_bar, a_bar     ! Average density and sound speed
real(kind=real64) :: f                  ! Common factor in equations

exact = .false.

aL = sqrt((gamma*PL)/rhoL)
aR = sqrt((gamma*PR)/rhoR)

CL = rhoL * aL
CR = rhoR * aR

f = 1.0/(CL + CR)

Pstar = f*(CR*PL + CL*PR + CL*CR*(uL - uR))
ustar = f*(CL*uL + CR*uR + (PL - PR))
rhoLstar = rhoL + (Pstar - PL)/(aL**2)
rhoRstar = rhoR + (Pstar - PR)/(aR**2)

end subroutine riemann_PVRS2


subroutine riemann_TRRS( &
    uL, rhoL, PL, uR, rhoR, PR, gamma, &
    Pstar, ustar, rhoLstar, rhoRstar, exact &
)

! Inputs
real(kind=real64), intent(in) :: uL, rhoL, PL   ! Left-hand state
real(kind=real64), intent(in) :: uR, rhoR, PR   ! Right-hand state
real(kind=real64), intent(in) :: gamma   ! EoS parameter

! Outputs
real(kind=real64), intent(out) :: Pstar, ustar, rhoLstar, rhoRstar
logical, intent(out) :: exact

! Local
real(kind=real64) :: aL, aR             ! Sound speeds
real(kind=real64) :: rho_bar, a_bar     ! Average density and sound speed
real(kind=real64) :: PLR                ! Pressure ratio in cell
real(kind=real64) :: z                  ! Common power in equations

exact = .false.

aL = sqrt((gamma*PL)/rhoL)
aR = sqrt((gamma*PR)/rhoR)

z = (gamma - 1.0)/(2.0*gamma)
PLR = (PL/PR)**z

ustar = PLR*uL/aL + uR/aR + 2.0*(PLR - 1.0)/(gamma - 1.0)
ustar = ustar / (PLR/aL + 1.0/aR)

Pstar = 0.5*( &
    PL*(1 + (gamma - 1.0)/(2*aL)*(uL - ustar))**(1.0/z) + &
    PR*(1 + (gamma - 1.0)/(2*aR)*(ustar - uR))**(1.0/z) &
)

rhoLstar = rhoL*(Pstar/PL)**(1/gamma)
rhoRstar = rhoR*(Pstar/PR)**(1/gamma)

end subroutine riemann_TRRS


subroutine riemann_TSRS( &
    uL, rhoL, PL, uR, rhoR, PR, gamma, &
    Pstar, ustar, rhoLstar, rhoRstar, exact &
)

! Inputs
real(kind=real64), intent(in) :: uL, rhoL, PL   ! Left-hand state
real(kind=real64), intent(in) :: uR, rhoR, PR   ! Right-hand state
real(kind=real64), intent(in) :: gamma   ! EoS parameter

! Outputs
real(kind=real64), intent(out) :: Pstar, ustar, rhoLstar, rhoRstar
logical, intent(out) :: exact

! Local
real(kind=real64) :: aaL, aaR               ! Sound speeds
real(kind=real64) :: f                      ! Common factor in equations
real(kind=real64) :: AL, AR, BL, BR, gL, gR ! Terms in equation
real(kind=real64) :: rho_bar, a_bar         ! Average density and sound speed
real(kind=real64) :: Pguess                 ! Initial pressure guess from PVRS

exact = .false.

f = (gamma - 1.0)/(gamma + 1.0)

AL = 2.0/((gamma + 1.0)*rhoL)
AR = 2.0/((gamma + 1.0)*rhoR)
BL = f*PL
BR = f*PR

! Initial guess at pressure from PVRS model
aaL = sqrt((gamma*PL)/rhoL)
aaR = sqrt((gamma*PR)/rhoR)
rho_bar = 0.5*(rhoL + rhoR)
a_bar = 0.5*(aaL + aaR)

Pguess = 0.5*(PL + PR) + 0.5*(uL - uR)*(rho_bar*a_bar)
Pguess = max(Pguess, 0.0)

gL = sqrt(AL/(Pguess + BL))
gR = sqrt(AR/(Pguess + BR))

Pstar = gL*PL + gR*PR - (uR - uL)
Pstar = Pstar / (gL + gR)

ustar = 0.5*(uL + uR) + 0.5*((Pstar - PR)*gR - (Pstar - PL)*gL)

rhoLstar = rhoL * ((Pstar/PL + f)/(f*Pstar/PL + 1))
rhoRstar = rhoR * ((Pstar/PR + f)/(f*Pstar/PR + 1))

end subroutine riemann_TSRS

! subroutine solve( &
!     uL, rhoL, PL, uR, rhoR, PR, gamma, model, &
!     rho, P, u, E &
! )
! ! Generate a solution in terms of primative variables for the Riemann problem
! ! represented by (uL, rhoL, PL) and (uR, rhoR, PR), using a provided Riemann
! ! solver (model).
!
! ! Inputs
! real(kind=real64), intent(in) :: uL, rhoL, PL   ! Left-hand state
! real(kind=real64), intent(in) :: uR, rhoR, PR   ! Right-hand state
! real(kind=real64), intent(in) :: gamma   ! EoS parameter
! procedure(riemann_API), pointer, intent(in) :: model
!
! ! Outputs
! real(kind=real64), intent(out) :: rho, P, u, E
!
! ! Local data
! real(kind=real64) :: Pstar, ustar       ! Output from Riemann solver
! real(kind=real64) :: rhoLstar, rhoRstar ! Output from approximate Rieman solver
! logical :: exact                        ! Is the Riemann solver exact?
!
! real(kind=real64) :: ein                ! Internal energy for E calculation
!
! ! Get P, u and maybe rho from Riemann solver
! call model(uL, rhoL, PL, uR, rhoR, PR, gamma, Pstar, ustar, rhoLstar, rhoRstar, exact)
!
!
! ! Sample the solution at S = x/t = 0.0
! call sample( &
!     Pstar, ustar, rhoLstar, rhoRstar, &
!     uL, rhoL, PL, &
!     uR, rhoR, PR, &
!     gamma, exact, &
!     rho, P, u &
! )
!
! ! Calculate energy
! ein = P/((gamma - 1.0)*rho)
! E = rho * (0.5*u*u + e)
!
! end subroutine solve


end module riemann_solvers
