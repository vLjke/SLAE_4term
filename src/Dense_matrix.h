#ifndef SLAE_4TERM_DENSE_MATRIX_H
#define SLAE_4TERM_DENSE_MATRIX_H

#include <iostream>
#include <vector>
#include <ranges>
#include <cmath>
#include <utility>
#include "Vector_operations.h"

template<typename T>
class Dense_matrix{
protected:
    size_t M{}, N{};
    std::vector<T> data;
public:
    // Constructors
    Dense_matrix() = default;
    Dense_matrix(const std::vector<std::vector<T>>& matrix);
    Dense_matrix(const std::vector<T>& n); // Projection constructor
    // Operators
    const T& operator()(size_t i, size_t j) const;
    T& operator()(size_t i, size_t j);
    template<typename Q>
    friend std::ostream& operator<<(std::ostream& os, const Dense_matrix<Q>& m);
    std::vector<T> operator*(const std::vector<T>& vec) const;
    std::pair<size_t, size_t> getOrder() const;
    // QR decomposition
    std::pair<Dense_matrix<T>, Dense_matrix<T>> QR_decomp_HH() const;
    // Descructor
    ~Dense_matrix() = default;
};

// ================ Function implementation ================
// Constructors
template<typename T>
Dense_matrix<T>::Dense_matrix(const std::vector<std::vector<T>>& matrix): M(matrix.size()), N(matrix[0].size()),
    data(std::ranges::join_view(matrix).begin(), std::ranges::join_view(matrix).end()) {}
template<typename T>
Dense_matrix<T>::Dense_matrix(const std::vector<T>& n) {
    this->M = n.size();
    this->N = n.size();
    for (int i = 0; i < this->M; ++i) {
        std::vector<T> e(this->M);
        e[i] = 1;
        e = e - 2 * n * n[i] / (n * n);
        for (int j = 0; j < this->M; ++j)
            this->data[j * this->M + i] = e[j];
    }
}
// Dense matrix operators
template<typename T>
const T& Dense_matrix<T>::operator()(size_t i, size_t j) const {
    return this->data[i * this->N + j];
}

template<typename T>
T& Dense_matrix<T>::operator()(size_t i, size_t j) {
    return this->data[i * this->N + j];
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Dense_matrix<T>& m) {
    for (int i = 0; i < m.getOrder().first; ++i) {
        for (int j = 0; j < m.getOrder().second; ++j)
            os << m(i, j) << " ";
        os << std::endl;
    }
    return os;
}

template<typename T>
std::vector<T> Dense_matrix<T>::operator*(const std::vector<T>& vec) const {
    std::vector<T> res(this->M, 0);
    for (int i = 0; i < this->M; ++i)
        for (int j = 0; j < this->N; ++j)
            res[i] += this->data[i * this->N + j] * vec[j];
    return res;
}

template<typename T>
std::pair<size_t, size_t> Dense_matrix<T>::getOrder() const {
    return std::make_pair(this->M, this->N);
}
// Dense matrix QR decomposition
template<typename T>
std::pair<Dense_matrix<T>, Dense_matrix<T>> Dense_matrix<T>::QR_decomp_HH() const {
    // R matrix
    Dense_matrix<T> R = *this;
    // Dimension of vector space on the first step
    size_t dimension = R.getOrder().first;
    // Getting first column
    std::vector<T> x(dimension);
    for (int i = 0; i < dimension; ++i)
        x[i] = R(i, 0);
    // Basis' vector e
    std::vector<T> e(dimension);
    e[0] = 1;
    // Finding n for the first projection
    std::vector<T> n(dimension);
    if (x[0] > 0)
        n = x + sqrt(x * x) * e;
    else
        n = x - sqrt(x * x) * e;
    // Creating Q_0
    Dense_matrix<T> Q {n};
    // First step, changing R
    R(0, 0) -= 2 * n[0] * (x * n) / (n * n);
    for (int i = 1; i < R.getOrder().first; ++i)
        R(i, 0) = 0;
    // Changing other columns after first one
    for (int j = 1; j < R.getOrder().second; ++j) {
        std::vector<T> x1(dimension);
        for (int i = 0; i < x1.size(); ++i)
            x1[i] = R(i, j);
        for (int i = 0; i < dimension; ++i)
            R(i, j) -= 2 * n[i] * (x1 * n) / (n * n);
    }
    // Cycle through all columns
    for (int k = 1; k < this->getOrder().second; ++k) {
        // Decreasing dimension of projection
        dimension--;
        x.resize(dimension);
        n.resize(dimension);
        e.resize(dimension);
        // Finding n for projection number k
        for (int i = 0; i < dimension; ++i)
            x[i] = R(i + k , k);
        if (x[0] > 0)
            n = x + sqrt(x * x) * e;
        else
            n = x - sqrt(x * x) * e;
        // Step number k, changing R
        R(k, k) -= 2 * n[0] * (x * n) / (n * n);
        for (int i = k + 1; i < R.getOrder().first; ++i)
            R(i, k) = 0;
        // Changing other columns after column number k
        for (int j = k + 1; j < R.getOrder().second; ++j) {
            for (int i = 0; i < dimension; ++i)
                x[i] = R(k + i, j);
            for (int i = 0; i < dimension; ++i)
                R(k + i, j) -= 2 * n[i] * (x * n) / (n * n);
        }
        // Step k, finding Q_k
        for (int i = 0; i < Q.getOrder().first; ++i) {
            for (int j = 0; j < dimension; ++j)
                x[j] = R(i, k + j);
            for (int j = 0; j < dimension; ++j)
                Q(i, j + k) -= 2 * n[i] * (x * n) / (n * n);
        }
    }
    // Transposing Q
    for (int i = 0; i < Q.getOrder().first; ++i)
        for (int j = i + 1; j < Q.getOrder().second; ++j)
            std::swap(Q(i, j), Q(j, i));
    return std::make_pair(Q, R);
}

#endif //SLAE_4TERM_DENSE_MATRIX_H

