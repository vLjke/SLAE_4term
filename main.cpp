#include <iostream>
#include "Tridiagonal_matrix/Tridiagonal_matrix.h"
#include "Tridiagonal_matrix_solver/Tridiag_matrix_solver.h"

int main() {
    Tridiag_matrix_solver solver {};
    std::vector<double> d {1, 2, 3, 4, 5};
    Tridiagonal_matrix m {5, {3, 4, 1, 11}, {6, 8, 9, 10, 13}, {3, 3, 4, 5}};
    std::cout << m << std::endl;
    for (auto v: d)
        std::cout << v << " ";
    std::cout << std::endl;
    for (auto v: solver.Solution(m, d))
        std::cout << v << " ";
    std::cout << std::endl;
    return 0;
}
