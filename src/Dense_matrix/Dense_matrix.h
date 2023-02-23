#ifndef SLAE_4TERM_DENSE_MATRIX_H
#define SLAE_4TERM_DENSE_MATRIX_H

#include <iostream>
#include <vector>
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
    const T& getElem(size_t i, size_t j) const;
    T& getElem(size_t i, size_t j);
    template<typename Q>
    friend std::ostream& operator<<(std::ostream& os, const Dense_matrix<Q>& m);
    std::vector<T> operator*(const std::vector<T>& vec) const;
    std::pair<size_t, size_t> getOrder() const;
    // Descructor
    ~Dense_matrix() = default;
};

#endif //SLAE_4TERM_DENSE_MATRIX_H

