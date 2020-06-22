program cactaurF
use iso_fortran_env, only: int32, real64
use riemann_solvers

integer(kind=int32) :: ncells



procedure(riemann_API), pointer :: model

model => riemann_PVRS1

end program cactaurF
