#ifndef SLAE_4TERM_TRIDIAGONAL_MATRIX_H
#define SLAE_4TERM_TRIDIAGONAL_MATRIX_H

#include <iostream>
#include <vector>

template<typename T>
class Tridiagonal_matrix {
private:
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
    T operator()(size_t i, size_t j) const;
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

// ================ Function implementation ================
// Constructors
template<typename T>
Tridiagonal_matrix<T>::Tridiagonal_matrix(int N, const std::vector<T>& a, const std::vector<T>& b, const std::vector<T>& c) {
    // Vector's size
    this->data.resize(N);
    // First line
    this->data[0].b = b[0];
    this->data[0].c = c[0];
    // Lines between first and last
    for (int i = 1; i < N - 1; ++i) {
        this->data[i].a = a[i - 1];
        this->data[i].b = b[i];
        this->data[i].c = c[i];
    }
    // Last line
    this->data[N - 1].a = a[N - 2];
    this->data[N - 1].b = b[N - 1];
}

// Tridiagonal matrix operators
template<typename T>
T  Tridiagonal_matrix<T>::operator()(size_t i, size_t j) const{
    if (i - j == 0)
        return this->data[i].b;
    else if (i - j == 1)
        return this->data[i].a;
    else if (i - j == -1)
        return this->data[i].c;
    else
        return static_cast<T>(0);
}

template<typename T>
const typename Tridiagonal_matrix<T>::Triplet& Tridiagonal_matrix<T>::operator[](size_t i) const {
    return this->data[i];
}

template<typename T>
typename Tridiagonal_matrix<T>::Triplet& Tridiagonal_matrix<T>::operator[](size_t i) {
    return this->data[i];
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Tridiagonal_matrix<T>& m) {
    if (m.getOrder() >= 3) {
        os << m.data[0].b << " " << m.data[0].c;
        for (int i = 0; i < m.getOrder() - 2; ++i)
            os << " 0";
        os << std::endl;
        for (int i = 1; i < m.getOrder() - 1; ++i) {
            for (int j = 0; j < i - 1; ++j)
                os << "0 ";
            os << m.data[i].a << " " << m.data[i].b << " " << m.data[i].c;;
            for (int j = 0; j < m.getOrder() - i - 2; ++j)
                os << " 0";
            os << std::endl;
        }
        for (int i = 0; i < m.getOrder() - 2; ++i)
            os << "0 ";
        os << m.data[m.getOrder() - 1].a << " " << m.data[m.getOrder() - 1].b;
    }
    else if (m.getOrder() == 2) {
        os << m.data[0].b << " " << m.data[0].c << std::endl;
        os << m.data[1].a << " " << m.data[1].b;
    }
    else
        os << m.data[0].b;
    return os;
}

template<typename T>
size_t Tridiagonal_matrix<T>::getOrder() const {
    return this->data.size();
}

// Tridiagonal matrix Solver
template<typename T>
std::vector<T> Tridiagonal_matrix<T>::Solution(const std::vector<T> &d) const {
    std::vector<T> p(this->getOrder());
    std::vector<T> q(this->getOrder());
    // Setting value of p_0 and q_0
    p[0] = -this->data[0].c / this->data[0].b;
    q[0] = d[0] / this->data[0].b;
    // First stage: finding p and q
    for (int i = 1; i < this->getOrder() - 1; ++i) {
        p[i] = -this->data[i].c / (this->data[i].a * p[i - 1] + this->data[i].b);
        q[i] = (d[i] - this->data[i].a * q[i - 1]) / (this->data[i].a * p[i - 1] + this->data[i].b);
    }
    // Setting value of X_n
    q[this->getOrder() - 1] = (d[this->getOrder() - 1] - this->data[this->getOrder() - 1].a * q[this->getOrder() - 2]) /
                              (this->data[this->getOrder() - 1].a * p[this->getOrder() - 2] + this->data[this->getOrder() - 1].b);
    // Second stage: finding X
    for (int i = this->getOrder() - 2; i >= 0; --i)
        q[i] += p[i] * q[i + 1];
    return q;
}


#endif //SLAE_4TERM_TRIDIAGONAL_MATRIX_H
