#include "Tridiagonal_matrix.h"

// Operator for triplet
std::ostream& operator<<(std::ostream& os, const Triplet& t){
    os << t.a << " " << t.b << " " << t.c;
    return os;
}

// Constructors
Tridiagonal_matrix::Tridiagonal_matrix(int N, std::vector<double>& a, std::vector<double>& b, std::vector<double>& c) {
    // Vector's size
    this->data.resize(N);
    // Matrix order
    this->order = this->data.size();
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

Tridiagonal_matrix::Tridiagonal_matrix(int N, std::vector<double>&& a, std::vector<double>&& b, std::vector<double>&& c) {
    // Vector's size
    this->data.resize(N);
    // Matrix order
    this->order = this->data.size();
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

Tridiagonal_matrix::Tridiagonal_matrix(std::vector<Triplet>& vec): data(vec) { this->order = this->data.size();}

Tridiagonal_matrix::Tridiagonal_matrix(std::vector<Triplet>&& vec) {
    std::move(vec.begin(), vec.end(), this->data.begin());
    // Matrix order
    this->order = this->data.size();
}
// Operators
Triplet& Tridiagonal_matrix::operator[](unsigned int i) {
    if (i < this->order)
        return this->data[i];
}

std::ostream& operator<<(std::ostream& os, const Tridiagonal_matrix& m) {
    if (m.order >= 3) {
        os << m.data[0].b << " " << m.data[0].c;
        for (int i = 0; i < m.order - 2; ++i)
            os << " 0";
        os << std::endl;
        for (int i = 1; i < m.order - 1; ++i) {
            for (int j = 0; j < i - 1; ++j)
                os << "0 ";
            os << m.data[i];
            for (int j = 0; j < m.order - i - 2; ++j)
                os << " 0";
            os << std::endl;
        }
        for (int i = 0; i < m.order - 2; ++i)
            os << "0 ";
        os << m.data[m.order - 1].a << " " << m.data[m.order - 1].b;
    }
    else if (m.order == 2) {
        os << m.data[0].b << " " << m.data[0].c << std::endl;
        os << m.data[1].a << " " << m.data[1].b;
    }
    else
        os << m.data[0].b;
    return os;
}


