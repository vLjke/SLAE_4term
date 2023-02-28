#ifndef SLAE_4TERM_CSR_MATRIX_H
#define SLAE_4TERM_CSR_MATRIX_H

#include <iostream>
#include <vector>
#include <ranges>
#include <algorithm>
#include <utility>

namespace DOK_cell_space
{
    template<typename T>
    struct cell {
        size_t i;
        size_t j;
        T value;
        // Comparison operator
        bool operator<(const cell<T>& right);
    };
}

namespace CSR_matrix_space
{
    template<typename T>
    class CSR_matrix {
    private:
        size_t M{}, N{};
        std::vector<T> values;
        std::vector<size_t> cols;
        std::vector<size_t> rows;
    public:
        // Constructors
        CSR_matrix() = default;
        CSR_matrix(size_t M, size_t N, const std::vector<DOK_cell_space::cell<T>>& vec);
        // Operators
        const T& operator()(size_t i, size_t j) const;
        T operator()(size_t i, size_t j);
        std::vector<T> operator*(const std::vector<T>& vec) const;
        std::pair<size_t, size_t> getOrder() const;
        // Destructor
        ~CSR_matrix() = default;
    };
}


// ================ Function implementation ================
// DOK_cell operator
template<typename T>
bool DOK_cell_space::cell<T>::operator<(const DOK_cell_space::cell<T>& right) {
    if (this->i > right.i || (this->i == right.i && this->j > right.j))
        return false;
    return true;
}

// Constructors
template<typename T>
CSR_matrix_space::CSR_matrix<T>::CSR_matrix(size_t M, size_t N, const std::vector<DOK_cell_space::cell<T>>& vec) {
    // Setting matrix dimensions
    this->M = M;
    this->N = N;
    // Sorting input data
    std::vector<DOK_cell_space::cell<T>> d = vec;
    std::sort(d.begin(), d.end());
    // Filling in values and cols
    this->values.resize(d.size());
    this->cols.resize(d.size());
    for (int i = 0; i < d.size(); ++i) {
        this->values[i] = d[i].value;
        this->cols[i] = d[i].j;
    }
    // Filling in rows
    this->rows.resize(M + 1, 0);
    if (!d.empty()) {
        int temp = 0;
        for (int i = 0; i < d[d.size() - 1].i + 1; ++i) {
            if (d[temp].i == i)
                temp += std::count_if(d.begin() + rows[i], d.end(),
                                      [i](const DOK_cell_space::cell<T> &c) { return c.i == i; });
            this->rows[i + 1] = temp;
        }
        // Last lines with only zeros if there are some
        for (int i = d[d.size() - 1].i + 1; i < this->M + 1; ++i)
            this->rows[i] = temp;
    }
}

// CSR_matrix operators
template<typename T>
const T& CSR_matrix_space::CSR_matrix<T>::operator()(size_t i, size_t j) const {
    // Lambda to find correct index j
    auto coincideIndex = [j](size_t k){return k == j;};
    // Finding needed index in correct range
    if (auto it = std::find_if(std::next(this->cols.begin(), this->rows[i]),
                               std::next(this->cols.begin(), this->rows[i + 1]), coincideIndex);
                                it != std::next(this->cols.begin(), this->rows[i + 1]))
        return this->values[std::distance(this->cols.begin(), it)];
    else
        return static_cast<T>(0);
}

template<typename T>
T CSR_matrix_space::CSR_matrix<T>::operator()(size_t i, size_t j) {
    // Lambda to find correct index j
    auto coincideIndex = [j](size_t k){return k == j;};
    // Finding needed index in correct range
    if (auto it = std::find_if(std::next(this->cols.begin(), this->rows[i]),
                               std::next(this->cols.begin(), this->rows[i + 1]), coincideIndex);
            it != std::next(this->cols.begin(), this->rows[i + 1]))
        return this->values[std::distance(this->cols.begin(), it)];
    else
        return static_cast<T>(0);
}

template<typename T>
std::vector<T> CSR_matrix_space::CSR_matrix<T>::operator*(const std::vector<T>& vec) const {
    // Result vector
    std::vector<T> result(this->M);
    // Multiply algorithm
    for (int k = 0; k < M; ++k)
        for (int s = this->rows[k]; s < this->rows[k + 1]; ++s)
            result[k] += vec[this->cols[s]] * this->values[s];
    return result;
}

template<typename T>
std::pair<size_t, size_t> CSR_matrix_space::CSR_matrix<T>::getOrder() const {
    return std::make_pair(this->M, this->N);
}


#endif //SLAE_4TERM_CSR_MATRIX_H