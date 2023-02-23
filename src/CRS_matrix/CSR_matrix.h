#ifndef SLAE_4TERM_CSR_MATRIX_H
#define SLAE_4TERM_CSR_MATRIX_H

#include <iostream>
#include <vector>
#include <utility>

namespace DOK_cell_space
{
    template<typename T>
    struct cell {
        size_t i;
        size_t j;
        T value;
        // Comparison operators
        bool operator==(const cell<T>& right);
        bool operator<(const cell<T>& right);
        bool operator>(const cell<T>& right);
    };
}

namespace CSR_matrix_space
{
    template<typename T>
    class CSR_matrix {
    protected:
        size_t M{}, N{};
        std::vector<T> values;
        std::vector<size_t> cols;
        std::vector<size_t> rows;
    public:
        // Constructors
        CSR_matrix() = default;
        CSR_matrix(size_t M, size_t N, const std::vector<DOK_cell_space::cell<T>>& vec);
        // Operators
        const T& getElem(size_t i, size_t j) const;
        T& getElem(size_t i, size_t j);
        std::vector<T> operator*(const std::vector<T>& vec) const;
        std::pair<size_t, size_t> getOrder() const;
        // Destructor
        ~CSR_matrix() = default;
    };
}

#endif //SLAE_4TERM_CSR_MATRIX_H