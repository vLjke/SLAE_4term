#ifndef SLAE_4TERM_TRIDIAGONAL_MATRIX_H
#define SLAE_4TERM_TRIDIAGONAL_MATRIX_H

#include <vector>

// Triplet struct
template<typename T>
struct Triplet {
    T a = 0;
    T b = 0;
    T c = 0;
};
// Display triplet
template<typename T>
std::ostream& operator<<(std::ostream& os, const Triplet<T>& t);

template<typename T>
class Tridiagonal_matrix {
protected:
    std::vector<Triplet<T>> data;
public:
    // Constructors
    Tridiagonal_matrix() = default;
    Tridiagonal_matrix(int N, std::vector<T>& a, std::vector<T>& b, std::vector<T>& c);
    Tridiagonal_matrix(int N, std::vector<T>&& a, std::vector<T>&& b, std::vector<T>&& c);
    Tridiagonal_matrix(std::vector<Triplet<T>>& vec);
    Tridiagonal_matrix(std::vector<Triplet<T>>&& vec);
    // Operators
    Triplet<T> operator[](unsigned int i) const;
    Triplet<T>& operator[](unsigned int i);
    template<typename Q>
    friend std::ostream& operator<<(std::ostream& os, const Tridiagonal_matrix<Q>& m);
    // Method to get matrix order
    unsigned int getOrder() const;
    // Destructor
    ~Tridiagonal_matrix() = default;
};


#endif //SLAE_4TERM_TRIDIAGONAL_MATRIX_H
