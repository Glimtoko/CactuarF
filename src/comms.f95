module parallel_comms
use iso_fortran_env, only: int32, real64
use mpi
implicit none

public :: parallel_update
private

contains

subroutine parallel_update(density, pressure, velocity, rank, nprocs)
! Inputs
! Note, we include the allocatable specifier here as it ensures Fortran
! passes the actual array bounds, rather than assuming 1:size
real(kind=real64), dimension(:), allocatable, intent(in) :: density
real(kind=real64), dimension(:), allocatable, intent(in) :: pressure
real(kind=real64), dimension(:), allocatable, intent(in) :: velocity
integer(kind=int32), intent(in) :: rank
integer(kind=int32), intent(in) :: nprocs

! Local variables
integer(kind=int32) :: error
! integer(kind=int32), dimension(3) :: recvL, recvR
! integer(kind=int32), dimension(3) :: sendL, sendR
integer(kind=int32) :: neighbourL, neighbourR
integer(kind=int32) :: L, R                 ! Bounds of main data
integer(kind=int32) :: Lg, Rg               ! Ghost cells

integer(kind=int32), dimension(:), allocatable :: recv, send
integer(kind=int32), dimension(:,:), allocatable :: status
integer(kind=int32) :: nrequests

if (nprocs == 1) return

Lg = lbound(density, 1)
L = Lg + 1

Rg = ubound(density, 1)
R = Rg - 1

! print *, size(density), lbound(density), ubound(density)
! print *, density(0)

neighbourL = rank - 1
neighbourR = rank + 1

! Post receives
if (rank == 0) then
    ! Left-most processor. Receieve from right-hand neighbour only
    nrequests = 3
    allocate(recv(nrequests), send(nrequests), status(MPI_STATUS_SIZE, nrequests))
    call MPI_Irecv(density(Rg), 1, MPI_DOUBLE, neighbourR, rank*100+1, MPI_COMM_WORLD, recv(1), error)
    call MPI_Irecv(pressure(Rg), 1, MPI_DOUBLE, neighbourR, rank*100+2, MPI_COMM_WORLD, recv(2), error)
    call MPI_Irecv(velocity(Rg), 1, MPI_DOUBLE, neighbourR, rank*100+3, MPI_COMM_WORLD, recv(3), error)
elseif (rank == nprocs-1) then
    ! Right-most processor. Receive from left-hand neighbour only
    nrequests = 3
    allocate(recv(nrequests), send(nrequests), status(MPI_STATUS_SIZE, nrequests))

    call MPI_Irecv(density(Lg), 1, MPI_DOUBLE, neighbourL, rank*100+11, MPI_COMM_WORLD, recv(1), error)
    call MPI_Irecv(pressure(Lg), 1, MPI_DOUBLE, neighbourL, rank*100+12, MPI_COMM_WORLD, recv(2), error)
    call MPI_Irecv(velocity(Lg), 1, MPI_DOUBLE, neighbourL, rank*100+13, MPI_COMM_WORLD, recv(3), error)
else
    ! Non-boundary processor. Receive from left and right
    nrequests = 6
    allocate(recv(nrequests), send(nrequests), status(MPI_STATUS_SIZE, nrequests))

    call MPI_Irecv(density(Rg), 1, MPI_DOUBLE, neighbourR, rank*100+1, MPI_COMM_WORLD, recv(1), error)
    call MPI_Irecv(pressure(Rg), 1, MPI_DOUBLE, neighbourR, rank*100+2, MPI_COMM_WORLD, recv(2), error)
    call MPI_Irecv(velocity(Rg), 1, MPI_DOUBLE, neighbourR, rank*100+3, MPI_COMM_WORLD, recv(3), error)

    call MPI_Irecv(density(Lg), 1, MPI_DOUBLE, neighbourL, rank*100+11, MPI_COMM_WORLD, recv(4), error)
    call MPI_Irecv(pressure(Lg), 1, MPI_DOUBLE, neighbourL, rank*100+12, MPI_COMM_WORLD, recv(5), error)
    call MPI_Irecv(velocity(Lg), 1, MPI_DOUBLE, neighbourL, rank*100+13, MPI_COMM_WORLD, recv(6), error)
end if


! Post sends
if (rank == 0) then
    ! Left-most processor. Send to right-hand neighbour only
    call MPI_Isend(density(R), 1, MPI_DOUBLE, neighbourR, neighbourR*100+11, MPI_COMM_WORLD, send(1), error)
    call MPI_Isend(pressure(R), 1, MPI_DOUBLE, neighbourR, neighbourR*100+12, MPI_COMM_WORLD, send(2), error)
    call MPI_Isend(velocity(R), 1, MPI_DOUBLE, neighbourR, neighbourR*100+13, MPI_COMM_WORLD, send(3), error)
elseif (rank == nprocs-1) then
    ! Right-most processor. Send to left-hand neighbour only
    call MPI_Isend(density(L), 1, MPI_DOUBLE, neighbourL, neighbourL*100+1, MPI_COMM_WORLD, send(1), error)
    call MPI_Isend(pressure(L), 1, MPI_DOUBLE, neighbourL, neighbourL*100+2, MPI_COMM_WORLD, send(2), error)
    call MPI_Isend(velocity(L), 1, MPI_DOUBLE, neighbourL, neighbourL*100+3, MPI_COMM_WORLD, send(3), error)
else
    ! Non-boundary processor. Send to left and right
    call MPI_Isend(density(R), 1, MPI_DOUBLE, neighbourR, neighbourR*100+11, MPI_COMM_WORLD, send(1), error)
    call MPI_Isend(pressure(R), 1, MPI_DOUBLE, neighbourR, neighbourR*100+12, MPI_COMM_WORLD, send(2), error)
    call MPI_Isend(velocity(R), 1, MPI_DOUBLE, neighbourR, neighbourR*100+13, MPI_COMM_WORLD, send(3), error)

    call MPI_Isend(density(L), 1, MPI_DOUBLE, neighbourL, neighbourL*100+1, MPI_COMM_WORLD, send(4), error)
    call MPI_Isend(pressure(L), 1, MPI_DOUBLE, neighbourL, neighbourL*100+2, MPI_COMM_WORLD, send(5), error)
    call MPI_Isend(velocity(L), 1, MPI_DOUBLE, neighbourL, neighbourL*100+3, MPI_COMM_WORLD, send(6), error)
end if

! Wait for receives to complete, since that's what we care about really
call MPI_Waitall(nrequests, recv, status, error)

end subroutine

end module parallel_comms
