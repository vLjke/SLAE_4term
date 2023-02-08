#pragma once

#ifndef SLAE_4TERM_TRIDIAGONAL_MATRIX_H
#define SLAE_4TERM_TRIDIAGONAL_MATRIX_H

#include <iostream>
#include <vector>

// Triplet struct
struct Triplet {
    double a = 0;
    double b = 0;
    double c = 0;
};
// Display triplet
std::ostream& operator<<(std::ostream& os, const Triplet& t);


class Tridiagonal_matrix {
    friend class Tridiag_matrix_solver;
protected:
    std::vector<Triplet> data;
    // Matrix order
    size_t order = 0;
public:
    // Constructors
    Tridiagonal_matrix() = default;
    Tridiagonal_matrix(int N, std::vector<double>& a, std::vector<double>& b, std::vector<double>& c);
    Tridiagonal_matrix(int N, std::vector<double>&& a, std::vector<double>&& b, std::vector<double>&& c);
    Tridiagonal_matrix(std::vector<Triplet>& vec);
    Tridiagonal_matrix(std::vector<Triplet>&& vec);
    // Operators
    Triplet& operator[](unsigned int i);
    friend std::ostream& operator<<(std::ostream& os, const Tridiagonal_matrix& m);
    // Destructor
    ~Tridiagonal_matrix() = default;
};


#endif //SLAE_4TERM_TRIDIAGONAL_MATRIX_H
