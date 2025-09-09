cmake_minimum_required(VERSION 3.1)

project(cactuarf VERSION 0.1
        DESCRIPTION "Cactuar Godonov Solver"
        LANGUAGES Fortran)
enable_language(Fortran)

# Currently setting the Fortran compiler to use -std=gnu, change this if you
# want a specific standard
set(FVERSION "-std=f2018")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION}")

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Source code
add_subdirectory(src)
