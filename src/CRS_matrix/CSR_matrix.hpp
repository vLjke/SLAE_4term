#ifndef SLAE_4TERM_CSR_MATRIX_HPP
#define SLAE_4TERM_CSR_MATRIX_HPP

#include "CSR_matrix.h"
#include <ranges>
#include <algorithm>

// Comparison operators
template<typename T>
bool DOK_cell_space::cell<T>::operator==(const DOK_cell_space::cell<T>&right) {
    return this->i == right.i && this->j = right.j;
}

template<typename T>
bool DOK_cell_space::cell<T>::operator<(const DOK_cell_space::cell<T>& right) {
    if (this->i > right.i || (this->i == right.i && this->j > right.j))
        return false;
    return true;
}

template<typename T>
bool DOK_cell_space::cell<T>::operator>(const DOK_cell_space::cell<T>& right) {
    return !(*this < right);
}

// Constructors
template<typename T>
CSR_matrix_space::CSR_matrix<T>::CSR_matrix(size_t M, size_t N, const std::vector<DOK_cell_space::cell<T>>& vec) {
    // Setting matrix dimensions
    this->M = M;
    this->N = N;
    // Checking whether input data is correct and sorting it
    auto correctIndex = [this](const DOK_cell_space::cell<T>& c)
            { return c.i < this->getOrder().first && c.j < this->getOrder().second; };
    auto notZeroValue = [](const DOK_cell_space::cell<T>& c)
            { return c.value != 0; };
    auto t = vec
             | std::ranges::views::filter(correctIndex)
             | std::ranges::views::filter(notZeroValue);
    std::vector<DOK_cell_space::cell<T>> d {};
    for (const auto& v: t)
        d.push_back(v);
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
            temp += std::ranges::count_if(d, [i](const DOK_cell_space::cell<T> &c) { return c.i == i; });
            this->rows[i + 1] = temp;
        }
        // Last lines with only zeros if there are some
        for (int i = d[d.size() - 1].i + 1; i < this->M + 1; ++i)
            this->rows[i] = temp;
    }
}

// Operators
template<typename T>
const T& CSR_matrix_space::CSR_matrix<T>::getElem(size_t i, size_t j) const {
    // Lambda to find correct index j
    auto coincideIndex = [j](size_t k){return k == j;};
    // Finding needed index in correct range
    if (auto it = std::ranges::find_if(this->cols | std::ranges::views::drop(this->rows[i])
            | std::ranges::views::take(this->rows[i + 1] - this->rows[i]), coincideIndex); it)
        return this->values[*it];
    else
        return static_cast<T>(0);
}

template<typename T>
T& CSR_matrix_space::CSR_matrix<T>::getElem(size_t i, size_t j) {
    // Lambda to find correct index j
    auto coincideIndex = [j](size_t k){return k == j;};
    // Finding needed index in correct range
    if (auto it = std::ranges::find_if(this->cols | std::ranges::views::drop(this->rows[i])
            | std::ranges::views::take(this->rows[i + 1] - this->rows[i]), coincideIndex); it)
        return this->values[*it];
    else
        return static_cast<T>(0);
}

template<typename T>
std::vector<T> CSR_matrix_space::CSR_matrix<T>::operator*(const std::vector<T>& vec) const {
    // Result vector
    std::vector<T> result(this->M, 0);
    // Multiply algorithm
    for (int k = 0; k < M; ++k) {
        if (this->rows[k + 1] - this->rows[k] > 0) {
            for (int s = this->rows[k]; s < this->rows[k + 1]; ++s)
                result[k] += vec[this->cols[s]] * this->values[s];
        }
    }
    return result;
}

template<typename T>
std::pair<size_t, size_t> CSR_matrix_space::CSR_matrix<T>::getOrder() const {
    return std::make_pair(this->M, this->N);
}

#endif //SLAE_4TERM_CSR_MATRIX_HPP