#include "Tridiagonal_matrix.h"

// Operator for triplet
template<typename T>
std::ostream& operator<<(std::ostream& os, const Triplet<T>& t){
    os << t.a << " " << t.b << " " << t.c;
    return os;
}

// Constructors
template<typename T>
Tridiagonal_matrix<T>::Tridiagonal_matrix(int N, std::vector<T>& a, std::vector<T>& b, std::vector<T>& c) {
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

template<typename T>
Tridiagonal_matrix<T>::Tridiagonal_matrix(int N, std::vector<T>&& a, std::vector<T>&& b, std::vector<T>&& c) {
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

template<typename T>
Tridiagonal_matrix<T>::Tridiagonal_matrix(std::vector<Triplet<T>>& vec): data(vec) {}

template<typename T>
Tridiagonal_matrix<T>::Tridiagonal_matrix(std::vector<Triplet<T>>&& vec) {
    std::move(vec.begin(), vec.end(), this->data.begin());
}
// Operators
template<typename T>
Triplet<T> Tridiagonal_matrix<T>::operator[](unsigned int i) const {
    if (i < this->getOrder())
        return this->data[i];
}
template<typename T>
Triplet<T>& Tridiagonal_matrix<T>::operator[](unsigned int i) {
    if (i < this->getOrder())
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
            os << m.data[i];
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
unsigned int Tridiagonal_matrix<T>::getOrder() const {
    return this->data.size();
}


