module text_output
use iso_fortran_env, only: int32, real64
implicit none (type, external)
private
public :: do_text_output
contains

subroutine do_text_output(x, density, pressure, velocity, energy, rank)
! Inputs
! Note, we include the allocatable specifier here as it ensures Fortran
! passes the actual array bounds, rather than assuming 1:size
real(kind=real64), dimension(:), allocatable, intent(in) :: x
real(kind=real64), dimension(:), allocatable, intent(in) :: density
real(kind=real64), dimension(:), allocatable, intent(in) :: pressure
real(kind=real64), dimension(:), allocatable, intent(in) :: velocity
real(kind=real64), dimension(:), allocatable, intent(in) :: energy
integer(kind=int32), intent(in) :: rank

! Local variables
integer(kind=int32) :: L, R                     ! Bounds of main data
integer(kind=int32) :: i                        ! Loop
real(kind=real64) :: ein                        ! Internal energy

integer(kind=int32), save :: filenumber = 1     ! File number
character(len=40) :: filename_stub              ! Base of file name
character(len=4) :: filenumber_str              ! File number as string (fXXX)
character(len=4) :: rank_str                    ! Rank as string (rXXX)
character(len=4), parameter :: ext = ".dat"     ! Output file extension
integer(kind=int32) :: unitnum                  ! Output file unit

! Array bounds
L = lbound(density, 1) + 1
R = ubound(density, 1) - 1

! File name stub
filename_stub = "output."

! Generate the file number string
write(filenumber_str,'("f",I3.3)') filenumber

! Generate the rank number string
write(rank_str,'("r",I3.3)') rank

open(newunit=unitnum, file=trim(filename_stub)//filenumber_str//"."//rank_str//ext)
do i = L, R
    ein  = energy(i)/density(i) - 0.5*velocity(i)*velocity(i)
    write(unitnum,*) i, x(i), density(i), pressure(i), velocity(i), ein
end do
close(unit=unitnum)

end subroutine do_text_output


end module text_output
