cmake_minimum_required(VERSION 3.1)

project(cactuarf VERSION 1.0
        DESCRIPTION "Cactuar Godonov Solver"
        LANGUAGES Fortran)
enable_language(Fortran)

# Set Fortran Standard
set(FVERSION "-std=f2018")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FVERSION}")

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Source code
add_subdirectory(src)
