#ifndef SLAE_4TERM_DENSE_MATRIX_H
#define SLAE_4TERM_DENSE_MATRIX_H

#include <iostream>
#include <vector>
#include <ranges>
#include <utility>

template<typename T>
class Dense_matrix{
protected:
    size_t M{}, N{};
    std::vector<T> data;
public:
    // Constructors
    Dense_matrix() = default;
    Dense_matrix(const std::vector<std::vector<T>>& matrix);
    // Operators
    const T& operator()(size_t i, size_t j) const;
    T& operator()(size_t i, size_t j);
    template<typename Q>
    friend std::ostream& operator<<(std::ostream& os, const Dense_matrix<Q>& m);
    std::vector<T> operator*(const std::vector<T>& vec) const;
    std::pair<size_t, size_t> getOrder() const;
    // Descructor
    ~Dense_matrix() = default;
};

// ================ Function implementation ================
// Constructors
template<typename T>
Dense_matrix<T>::Dense_matrix(const std::vector<std::vector<T>>& matrix): M(matrix.size()), N(matrix[0].size()),
    data(std::ranges::join_view(matrix).begin(), std::ranges::join_view(matrix).end()) {}
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

#endif //SLAE_4TERM_DENSE_MATRIX_H

