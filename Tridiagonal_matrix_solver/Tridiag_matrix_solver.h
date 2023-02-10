#ifndef SLAE_4TERM_TRIDIAG_MATRIX_SOLVER_H
#define SLAE_4TERM_TRIDIAG_MATRIX_SOLVER_H

#include "../Tridiagonal_matrix/Tridiagonal_matrix.h"


template<typename T>
class Tridiag_matrix_solver {
protected:
    std::vector<T> p;
    std::vector<T> q;
public:
    // Constructors
    Tridiag_matrix_solver() = default;
    // Function to solve SLAE
    std::vector<T> Solution(const Tridiagonal_matrix<T>& m, const std::vector<T>& d);
    // Destructor
    ~Tridiag_matrix_solver() = default;
};


#endif //SLAE_4TERM_TRIDIAG_MATRIX_SOLVER_H
