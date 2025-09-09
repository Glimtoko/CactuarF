set(MY_SRC   abort.f95
             comms.f95
             riemann_api.f95
             riemann_sampler.f95
             riemann.f95
             cactuar.f95
             flux.f95
             riemann_exact.f95
             riemann_solvers.f95
             text_output.f95
             CACHE INTERNAL "")
add_executable(cactuarf "${MY_SRC}")
target_link_libraries(cactuarf MPI::MPI_Fortran)