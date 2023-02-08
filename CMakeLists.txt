cmake_minimum_required(VERSION 3.24)
project(SLAE_4term)

set(CMAKE_CXX_STANDARD 20)

add_executable(SLAE_4term main.cpp
        "Tridiagonal_matrix/Tridiagonal_matrix.cpp" "Tridiagonal_matrix/Tridiagonal_matrix.h"
        "Tridiagonal_matrix_solver/Tridiag_matrix_solver.cpp" "Tridiagonal_matrix_solver/Tridiag_matrix_solver.h"
        )
