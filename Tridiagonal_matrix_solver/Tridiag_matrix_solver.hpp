#include "Tridiag_matrix_solver.h"

// Function to solve SLAE
template<typename T>
std::vector<T> Tridiag_matrix_solver<T>::Solution(const Tridiagonal_matrix<T>& m, const std::vector<T>& d) {
    // p[0] and q[0] will not be used
    this->p.resize(m.getOrder());
    this->q.resize(m.getOrder());
    // Setting value of p_1 and q_1
    this->p[1] = -m[0].c / m[0].b;
    this->q[1] = d[0] / m[0].b;
    // First stage: finding p and q
    for (int i = 2; i < m.getOrder(); ++i) {
        this->p[i] = -m[i - 1].c / (m[i - 1].a * this->p[i - 1] + m[i - 1].b);
        this->q[i] = (d[i - 1] - m[i - 1].a * this->q[i - 1]) / (m[i - 1].a * this->p[i - 1] + m[i - 1].b);
    }
    // X vector
    std::vector<T> result(m.getOrder());
    // Setting value of X_n
    result[m.getOrder() - 1] = (d[m.getOrder() - 1] - m[m.getOrder() - 1].a * q[m.getOrder() - 1]) /
                            (m[m.getOrder() - 1].a * p[m.getOrder() - 1] + m[m.getOrder() - 1].b);
    // Second stage: finding X
    for (int i = m.getOrder() - 2; i >= 0; --i)
        result[i] = p[i + 1] * result[i + 1] + q[i + 1];
    return result;
}