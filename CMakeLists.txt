cmake_minimum_required(VERSION 3.24)
project(SLAE_4term)

set(CMAKE_CXX_STANDARD 20)

enable_testing()
add_subdirectory(googletest)

add_executable(SLAE_4term main.cpp
        "Tridiagonal_matrix/Tridiagonal_matrix.hpp" "Tridiagonal_matrix/Tridiagonal_matrix.h"
        "Tridiagonal_matrix_solver/Tridiag_matrix_solver.hpp" "Tridiagonal_matrix_solver/Tridiag_matrix_solver.h"
        )

target_link_libraries(SLAE_4term gtest gmock)

