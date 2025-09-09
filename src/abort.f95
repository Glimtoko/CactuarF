module abort_mod
use iso_fortran_env, only: int32

use mpi, only: MPI_ABORT, MPI_COMM_WORLD

implicit none (type, external)

private
public :: abort

contains

subroutine abort(reason)
character(len=*), intent(in) :: reason
integer(kind=int32) :: status

write(*,'("CactuarF has encountered an abort condition:")')
write(*,'(a)') reason

call MPI_ABORT(MPI_COMM_WORLD, 0, status)

end subroutine abort

end module abort_mod
