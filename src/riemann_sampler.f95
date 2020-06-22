module riemann_sampler
use iso_fortran_env, only: int32, real64
implicit none
public :: sample
private

contains

real(kind=real64) function SL_shock(P, PL, rhoL, uL, gamma)
! Calculates wave speed of a left-travelling shock. Here the notation "k"
! can stand for either L (left) or R (right).

! Input
real(kind=real64), intent(in) :: P              ! Pressure in star region
real(kind=real64), intent(in) :: PL, rhoL, uL   ! State on boundary
real(kind=real64), intent(in) :: gamma          ! EoS parameter

! Local data
real(kind=real64) :: aL     ! Sound speed

SL_shock = uL - aL*((gamma + 1.0)/(2.0*gamma)*(P/PL) + &
                (gamma - 1.0)/(2.0*gamma))**0.5

end function SL_shock



real(kind=real64) function SR_shock(P, PR, rhoR, uR, gamma)
! Calculates wave speed of a right-travelling shock.

! Input
real(kind=real64), intent(in) :: P              ! Pressure in star region
real(kind=real64), intent(in) :: PR, rhoR, uR   ! State on right boundary
real(kind=real64), intent(in) :: gamma          ! EoS parameter

! Local data
real(kind=real64) :: aR     ! Sound speed

SR_shock = uR + aR*((gamma + 1.0)/(2.0*gamma)*(P/PR) + &
                (gamma - 1.0)/(2.0*gamma))**0.5

end function SR_shock



real(kind=real64) function rhostar_shock(P, Pk, rhok, gamma)
! Returns density inside the star region for a shock. Here the notation "k"
! can stand for either L (left) or R (right).

! Input
real(kind=real64), intent(in) :: P          ! Pressure in star region
real(kind=real64), intent(in) :: Pk, rhok   ! State on boundary
real(kind=real64), intent(in) :: gamma      ! EoS parameter

! Local data
real(kind=real64) :: f  ! Gamma ratio

f = (gamma - 1.0)/(gamma + 1.0)

rhostar_shock = rhok * (P/Pk + f) / (f*P/Pk + 1.0)

end function rhostar_shock



real(kind=real64) function rhostar_rarefaction(P, Pk, rhok, gamma)
! Returns density inside the star region for a rarefaction. Here the notation
! "k" can stand for either L (left) or R (right).

! Input
real(kind=real64), intent(in) :: P          ! Pressure in star region
real(kind=real64), intent(in) :: Pk, rhok   ! State on boundary
real(kind=real64), intent(in) :: gamma      ! EoS parameter

! Local data
real(kind=real64) :: f  ! Gamma ratio

f = (gamma - 1.0)/(gamma + 1.0)

rhostar_rarefaction = rhok * (P/Pk)**(1.0/gamma)

end function rhostar_rarefaction



subroutine rhoLfan_rarefaction(P, PL, rhoL, uL, gamma, S, rho_out, P_out, u_out)
! Returns the state inside the fan region for a left rarefaction.

! Input
real(kind=real64), intent(in) :: P              ! Pressure in star region
real(kind=real64), intent(in) :: PL, rhoL, uL   ! State on boundary
real(kind=real64), intent(in) :: gamma          ! EoS parameter
real(kind=real64), intent(in) :: S              ! Sample point. Usually 0

! Outputs
real(kind=real64), intent(out) :: rho_out, P_out, u_out  ! Calculated state

! Local data
real(kind=real64) :: aL     ! Sound speed
real(kind=real64) :: f      ! Power used in equations

aL = sqrt((gamma*PL)/rhoL)
f = (2*gamma)/(gamma - 1.0)

rho_out = rhoL*(2.0/(gamma + 1.0) + &
               (gamma - 1.0)/((gamma + 1.0)*aL)*(uL - S))**(2/(gamma - 1.0))

u_out = 2.0/(gamma+1.0)*(aL + (gamma - 1.0)/2.0 * uL + S)

P_out = PL * (2.0/(gamma + 1.0) + &
             (gamma - 1.0)/((gamma + 1.0)*aL) * (uL - S))**f

end subroutine rhoLfan_rarefaction



subroutine rhoRfan_rarefaction(P, PR, rhoR, uR, gamma, S, rho_out, P_out, u_out)
! Returns the state inside the fan region for a right rarefaction.

! Input
real(kind=real64), intent(in) :: P              ! Pressure in star region
real(kind=real64), intent(in) :: PR, rhoR, uR   ! State on boundary
real(kind=real64), intent(in) :: gamma          ! EoS parameter
real(kind=real64), intent(in) :: S              ! Sample point. Usually 0

! Outputs
real(kind=real64), intent(out) :: rho_out, P_out, u_out  ! Calculated state

! Local data
real(kind=real64) :: aR     ! Sound speed
real(kind=real64) :: f      ! Power used in equations

aR = sqrt((gamma*PR)/rhoR)
f = 2*gamma/(gamma - 1.0)

rho_out = rhoR*(2.0/(gamma + 1.0) - &
           (gamma - 1.0)/((gamma + 1.0)*aR)*(uR - S))**(2/(gamma - 1.0))

u_out = 2.0/(gamma+1.0)*(-aR + (gamma - 1.0)/2.0 * uR + S)

P_out = PR * (2.0/(gamma + 1.0) - &
         (gamma - 1.0)/((gamma + 1.0)*aR) * (uR - S))**f


end subroutine rhoRfan_rarefaction


subroutine sample( &
    Pstar, ustar, rhoLstar, rhoRstar, &
    uL, rhoL, PL, &
    uR, rhoR, PR, &
    gamma, exact, &
    rho_out, P_out, u_out &
)

! Inputs
real(kind=real64), intent(in) :: Pstar, ustar        ! Conditions in star region
                                                     ! from Riemann solve
real(kind=real64), intent(in) ::  rhoLstar, rhoRstar ! Conditions in star region
                                                     ! from approx Riemann
real(kind=real64), intent(in) :: uL, rhoL, PL        ! Left-hand boundary state
real(kind=real64), intent(in) :: uR, rhoR, PR        ! Right-hand boundary state
real(kind=real64), intent(in) :: gamma               ! EoS parameter
logical, intent(in) :: exact               ! Was solver exact?

! Outputs
real(kind=real64), intent(out) :: rho_out, P_out, u_out ! State at sample point

! Local data
real(kind=real64) :: Sk     ! Shock speed
real(kind=real64) :: ak     ! Sound speed on boundary
real(kind=real64) :: astar  ! Sound speed in star region
real(kind=real64) :: SHk    ! Speed of wave head
real(kind=real64) :: STk    ! Speed of wave tail

! Note that the sample position, S (=x/t) is always 0.0 for a Godunov solve, and
! so is stored as a parameter
real(kind=real64), parameter :: S = 0.0

! Which side of the solution are we on?
if (S < ustar) then
    ! Left side of solution
    if (Pstar > PL) then
        ! Shock
        Sk = SL_shock(Pstar, PL, rhoL, uL, gamma)
        if (S < Sk) then
            ! Outside star region, return boundary state
            rho_out = rhoL
            P_out = PL
            u_out = uL
        else
            ! Inside star region
            if (exact) then
                ! Exact Riemann solver requires calculation of density
                rho_out = rhostar_shock(Pstar, PL, rhoL, gamma)
                P_out = Pstar
                u_out = ustar
            else
                rho_out = rhoLstar
                P_out = Pstar
                u_out = ustar
            end if
        end if

    else
        ! Rarefaction
        ak = sqrt((gamma*PL)/rhoL)
        astar = ak*(Pstar/PL)**((gamma - 1.0)/(2.0*gamma))

        SHk = uL - ak
        STk = ustar - astar

        if (S <= SHk) then
            ! Outside of rarefaction, return boundary state
            rho_out = rhoL
            P_out = PL
            u_out = uL
        elseif (S <= STk) then
            ! Inside fan region
            call rhoLfan_rarefaction( &
                Pstar, PL, rhoL, uL, gamma, S, &
                rho_out, P_out, u_out &
            )
        else
            ! Inside star region
            if (exact) then
                rho_out = rhostar_rarefaction(Pstar, PL, rhoL, gamma)
                P_out = Pstar
                u_out = ustar
            else
                rho_out = rhoLstar
                P_out = Pstar
                u_out = ustar
            end if
        end if
    end if


else
    ! Right side of solution
    if (Pstar > PR) then
        ! Shock
        Sk = SR_shock(Pstar, PR, rhoR, uR, gamma)
        if (S < Sk) then
            ! Outside star region, return boundary state
            rho_out = rhoR
            P_out = PR
            u_out = uR
        else
            ! Inside star region
            if (exact) then
                rho_out = rhostar_shock(Pstar, PR, rhoR, gamma)
                P_out = Pstar
                u_out = ustar
            else
                rho_out = rhoRstar
                P_out = Pstar
                u_out = ustar
            end if
        end if
    else
        ! Rarefaction
        ak = sqrt((gamma*PR)/rhoR)
        astar = ak*(Pstar/PR)**((gamma - 1.0)/(2.0*gamma))

        SHk = uR + ak
        STk = ustar + astar

        if (S >= SHk) then
            ! Outside of rarefaction, return boundary state
            rho_out = rhoR
            P_out = PR
            u_out = uR
        elseif (S >= STk) then
            ! Inside fan region
            call rhoRfan_rarefaction( &
                Pstar, PR, rhoR, uR, gamma, S, &
                rho_out, P_out, u_out &
            )
        else
            ! Inside star region
            if (exact) then
                rho_out = rhostar_rarefaction(Pstar, PR, rhoR, gamma)
                P_out = Pstar
                u_out = ustar
            else
                rho_out = rhoRstar
                P_out = Pstar
                u_out = ustar
            end if
        end if
    end if
end if

end subroutine sample

end module riemann_sampler
