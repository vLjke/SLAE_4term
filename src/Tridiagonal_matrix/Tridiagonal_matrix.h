#ifndef SLAE_4TERM_TRIDIAGONAL_MATRIX_H
#define SLAE_4TERM_TRIDIAGONAL_MATRIX_H

#include <iostream>
#include <vector>

// Matrix class
template<typename T>
class Tridiagonal_matrix {
protected:
    struct Triplet {
        T a = 0;
        T b = 0;
        T c = 0;
    };
    std::vector<Triplet> data;
public:
    // Constructors
    Tridiagonal_matrix() = default;
    Tridiagonal_matrix(int N, const std::vector<T>& a, const std::vector<T>& b, const std::vector<T>& c);
    // Operators
    T getElem(size_t i, size_t j) const;
    const Triplet& operator[](size_t i) const;
    Triplet& operator[](size_t i);
    template<typename Q>
    friend std::ostream& operator<<(std::ostream& os, const Tridiagonal_matrix<Q>& m);
    size_t getOrder() const;
    // Method to solve SLAE
    std::vector<T> Solution(const std::vector<T>& d) const;
    // Destructor
    ~Tridiagonal_matrix() = default;
};


#endif //SLAE_4TERM_TRIDIAGONAL_MATRIX_H
