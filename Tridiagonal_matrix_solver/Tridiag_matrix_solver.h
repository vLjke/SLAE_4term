#pragma once

#ifndef SLAE_4TERM_TRIDIAG_MATRIX_SOLVER_H
#define SLAE_4TERM_TRIDIAG_MATRIX_SOLVER_H

#include <vector>
#include "../Tridiagonal_matrix/Tridiagonal_matrix.h"

class Tridiag_matrix_solver {
protected:
    std::vector<double> p;
    std::vector<double> q;
public:
    // Constructors
    Tridiag_matrix_solver() = default;
    // Function to solve SLAE
    std::vector<double> Solution(const Tridiagonal_matrix& m, const std::vector<double>& d);
    // Destructor
    ~Tridiag_matrix_solver() = default;
};


#endif //SLAE_4TERM_TRIDIAG_MATRIX_SOLVER_H
