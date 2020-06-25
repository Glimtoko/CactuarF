module abort_mod
use iso_fortran_env, only: int32

use mpi
implicit none

contains

subroutine abort(reason)
character(len=*), intent(in) :: reason
integer(kind=int32) :: status

write(*,'("CactuarF has encountered an abort condition:")')
write(*,'(a)') reason

call MPI_ABORT(MPI_COMM_WORLD, 0, status)

end subroutine

end module abort_mod
