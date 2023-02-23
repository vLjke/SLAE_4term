#ifndef SLAE_4TERM_DENSE_MATRIX_HPP
#define SLAE_4TERM_DENSE_MATRIX_HPP

#include "Dense_matrix.h"

// Constructors
template<typename T>
Dense_matrix<T>::Dense_matrix(const std::vector<std::vector<T>>& matrix) {
    this->M = matrix.size();
    this->N = matrix[0].size();
    this->data.resize(this->M * this->N);
    for (int i = 0; i < this->M; ++i)
        for (int j = 0; j < this->N; ++j)
            this->getElem(i, j) = matrix[i][j];
}
// Operators
template<typename T>
const T& Dense_matrix<T>::getElem(size_t i, size_t j) const {
    return this->data[i * this->N + j];
}

template<typename T>
T& Dense_matrix<T>::getElem(size_t i, size_t j) {
    return this->data[i * this->N + j];
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Dense_matrix<T>& m) {
    for (int i = 0; i < m.getOrder().first; ++i) {
        for (int j = 0; j < m.getOrder().second; ++j)
            os << m.getElem(i, j) << " ";
        os << std::endl;
    }
    return os;
}

template<typename T>
std::vector<T> Dense_matrix<T>::operator*(const std::vector<T>& vec) const {
    std::vector<T> res(this->M, 0);
    for (int i = 0; i < this->M; ++i)
        for (int j = 0; j < this->N; ++j)
            res[i] += this->getElem(i, j) * vec[j];
    return res;
}

template<typename T>
std::pair<size_t, size_t> Dense_matrix<T>::getOrder() const {
    return std::make_pair(this->M, this->N);
}

#endif //SLAE_4TERM_DENSE_MATRIX_HPP
