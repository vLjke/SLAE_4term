#ifndef SLAE_4TERM_CSR_MATRIX_H
#define SLAE_4TERM_CSR_MATRIX_H

#include <iostream>
#include <vector>
#include <ranges>
#include <algorithm>
#include "Vector_operations.h"

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

        // Additional methods
        void SSOR_make_iteration(std::vector<T>& x, const std::vector<T>& b, double omega) const {
            // Step 0.5
            for (int k = 0; k < x.size(); ++k) {
                // temp var to find Dx later
                double temp = x[k];
                x[k] = omega * b[k];
                double diag;
                for (int s = this->rows[k]; s < this->rows[k + 1]; ++s) {
                    if (this->cols[s] == k)
                        diag = this->values[s];
                    else
                        x[k] -= omega * x[this->cols[s]] * this->values[s];
                }
                x[k] /= diag;
                x[k] += (1 - omega) * temp;
            }
            // Step 1
            for (int k = x.size() - 1; k >= 0; --k) {
                // temp var for to find Dx later
                double temp = x[k];
                x[k] = omega * b[k];
                double diag;
                for (int s = this->rows[k]; s < this->rows[k + 1]; ++s) {
                    if (this->cols[s] == k)
                        diag = this->values[s];
                    else
                        x[k] -= omega * x[this->cols[s]] * this->values[s];
                }
                x[k] /= diag;
                x[k] += (1 - omega) * temp;
            }
        }
        void Sym_GS_make_iteration(std::vector<T>& x, const std::vector<T>& b) const {
            // Step 0.5
            for (int k = 0; k < x.size(); ++k) {
                x[k] = b[k];
                double diag;
                for (int s = this->rows[k]; s < this->rows[k + 1]; ++s) {
                    if (this->cols[s] == k)
                        diag = this->values[s];
                    else
                        x[k] -= x[this->cols[s]] * this->values[s];
                }
                x[k] /= diag;
            }
            // Step 1
            for (int k = x.size() - 1; k >= 0; --k) {
                x[k] = b[k];
                double diag;
                for (int s = this->rows[k]; s < this->rows[k + 1]; ++s) {
                    if (this->cols[s] == k)
                        diag = this->values[s];
                    else
                        x[k] -= x[this->cols[s]] * this->values[s];
                }
                x[k] /= diag;
            }
        }

    public:
        // Constructors
        CSR_matrix() = default;
        CSR_matrix(size_t M, size_t N, const std::vector<DOK_cell_space::cell<T>>& vec);

        // Operators
        const T& operator()(size_t i, size_t j) const;
        T operator()(size_t i, size_t j);
        std::vector<T> operator*(const std::vector<T>& vec) const;
        std::pair<size_t, size_t> getOrder() const;

        // SLAE Solvers
        // Jacobi method
        std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> Jacobi_method
        (const std::vector<T>& x_0, const std::vector<T>& b, double accuracy) const;
        // Gauss-Seidel method
        std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> Gauss_Seidel_method
        (const std::vector<T>& x_0, const std::vector<T>& b, double accuracy) const;
        // Symmetrical Gauss-Seidel w/ Chebyshev acceleration method
        std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> Symmetrical_GS_method
        (const std::vector<T>& x_0, const std::vector<T>& b, double rho, double accuracy) const;
        // SOR method
        std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> SOR_method
        (const std::vector<T>& x_0, const std::vector<T>& b, double omega, double accuracy) const;
        // SSOR w/ Chebyshev acceleration method
        std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> SSOR_method
        (const std::vector<T>& x_0, const std::vector<T>& b, double omega, double rho, double accuracy) const;
        // Simple-iteration method
        std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> Simple_iteration_method
        (const std::vector<T>& x_0, const std::vector<T>& b, double tau, double accuracy) const;
        // Simple-iteration method w/ Chebyshev acceleration
        std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> SIM_Chebyshev_acceleration
        (const std::vector<T>& x_0, const std::vector<T>& b, size_t r, double eig_min, double eig_max, double accuracy) const;
        // Steepest descent method
        std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> Steepest_descent_method
        (const std::vector<T>& x_0, const std::vector<T>& b, double accuracy) const;
        // Heavy ball method
        std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> Heavy_ball_method
        (const std::vector<T>& x_0, const std::vector<T>& b, double accuracy) const;

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
    for (int k = 0; k < this->M; ++k)
        for (int s = this->rows[k]; s < this->rows[k + 1]; ++s)
            result[k] += vec[this->cols[s]] * this->values[s];
    return result;
}

template<typename T>
std::pair<size_t, size_t> CSR_matrix_space::CSR_matrix<T>::getOrder() const {
    return std::make_pair(this->M, this->N);
}

// CSR matrix SLAE solvers
// Jacobi method
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::Jacobi_method
(const std::vector<T>& x_0, const std::vector<T>& b, double accuracy) const
{
    // Initial approximation
    std::vector<T> x0 = x_0;
    std::vector<T> r = b - (*this) * x0;
    // x on next iteration
    std::vector<T> x(x_0.size());
    // Total number of iterations
    size_t count = 0;
    // Vector to collect residual on each iteration
    std::vector<double> nev {sqrt(r * r)};
    // Stop condition
    while (nev[count] > accuracy) {
        // Finding x_(i+1)
        for (int k = 0; k < x.size(); ++k) {
            x[k] = b[k];
            double diag;
            for (int s = this->rows[k]; s < this->rows[k + 1]; ++s) {
                if (this->cols[s] == k)
                    diag = this->values[s];
                else
                    x[k] -= x0[this->cols[s]] * this->values[s];
            }
            x[k] /= diag;
        }
        x0 = x;
        r = b - (*this) * x0;
        nev.push_back(sqrt(r * r));
        count++;
    }
    return std::make_pair(x0, std::make_pair(count, nev));
}

// Gauss-Seidel method
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::Gauss_Seidel_method
(const std::vector<T>& x_0, const std::vector<T>& b, double accuracy) const
{
    // Initial approximation
    std::vector<T> x = x_0;
    std::vector<T> r = b - (*this) * x;
    // Total number of iterations
    size_t count = 0;
    // Vector to collect residual on each iteration
    std::vector<double> nev {sqrt(r * r)};
    // Stop condition
    while (nev[count] > accuracy) {
        // Finding x_(i+1)
        for (int k = 0; k < x.size(); ++k) {
            x[k] = b[k];
            double diag;
            for (int s = this->rows[k]; s < this->rows[k + 1]; ++s) {
                if (this->cols[s] == k)
                    diag = this->values[s];
                else
                    x[k] -= x[this->cols[s]] * this->values[s];
            }
            x[k] /= diag;
        }
        r = b - (*this) * x;
        nev.push_back(sqrt(r * r));
        count++;
    }
    return std::make_pair(x, std::make_pair(count, nev));
}

// Symmetrical Gauss-Seidel w/ Chebyshev acceleration method
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::Symmetrical_GS_method
(const std::vector<T>& x_0, const std::vector<T>& b, double rho, double accuracy) const
{
    // y_0
    std::vector<T> y_even = x_0;
    // y_1
    std::vector<T> y_odd = x_0;
    Sym_GS_make_iteration(y_odd, b);
    // mu_0 && mu_1
    double mu_even = 1;
    double mu_odd = 1 / rho;
    // Vector to keep (Py + c)
    std::vector<T> temp = y_odd;
    // Find (Py_1 + c)
    Sym_GS_make_iteration(temp, b);
    // Vector to collect residual on each iteration
    std::vector<T> r = b - (*this) * y_odd;
    std::vector<double> nev {sqrt(r * r)};
    // Total number of iterations
    size_t count = 0;
    // Stop condition
    while (nev[count] > accuracy) {
        // On (i + 1) iteration we rewrite y_(i-1) vector and mu_(i-1),
        // on different iterations y_(i-1) is in different variables, same thing with mu_(i-1)
        // If count is even, y_(i-1) is in y_even and mu_(i-1) is in mu_even,
        if (count % 2 == 0) {
            y_even = -mu_even * y_even + 2 * mu_odd / rho * temp;
            mu_even = 2 / rho * mu_odd - mu_even;
            y_even /= mu_even;
            r = b - (*this) * y_even;
            temp = y_even;
            // Calculate (Py_i + c)
            Sym_GS_make_iteration(temp, b);
        }
        else {
            y_odd = -mu_odd * y_odd + 2 * mu_even / rho * temp;
            mu_odd = 2 / rho * mu_even - mu_odd;
            y_odd /= mu_odd;
            r = b - (*this) * y_odd;
            temp = y_odd;
            // Calculate (Py_i + c)
            Sym_GS_make_iteration(temp, b);
        }
        nev.push_back(sqrt(r * r));
        count++;
    }
    if (count % 2 == 0)
        return std::make_pair(y_odd, std::make_pair(count, nev));
    else
        return std::make_pair(y_even, std::make_pair(count, nev));
}

// SOR method
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::SOR_method
(const std::vector<T>& x_0, const std::vector<T>& b, double omega, double accuracy) const
{
    // Initial approximation
    std::vector<T> x = x_0;
    std::vector<T> r = b - (*this) * x;
    // Total number of iterations
    size_t count = 0;
    // Vector to collect residual on each iteration
    std::vector<double> nev {sqrt(r * r)};
    // Stop condition
    while (nev[count] > accuracy) {
        // Finding x_(i+1)
        for (int k = 0; k < x.size(); ++k) {
            // temp var to find Dx later
            double temp = x[k];
            x[k] = omega * b[k];
            double diag;
            for (int s = this->rows[k]; s < this->rows[k + 1]; ++s) {
                if (this->cols[s] == k)
                    diag = this->values[s];
                else
                    x[k] -= omega * x[this->cols[s]] * this->values[s];
            }
            x[k] /= diag;
            x[k] += (1 - omega) * temp;
        }
        r = b - (*this) * x;
        nev.push_back(sqrt(r * r));
        count++;
    }
    return std::make_pair(x, std::make_pair(count, nev));
}

// SSOR w/ Chebyshev acceleration method
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::SSOR_method
(const std::vector<T>& x_0, const std::vector<T>& b, double omega, double rho, double accuracy) const
{
    // y_0
    std::vector<T> y_even = x_0;
    // y_1
    std::vector<T> y_odd = x_0;
    SSOR_make_iteration(y_odd, b, omega);
    // mu_0 && mu_1
    double mu_even = 1;
    double mu_odd = 1 / rho;
    // Vector to keep (Py + c)
    std::vector<T> temp = y_odd;
    // Find (Py_1 + c)
    SSOR_make_iteration(temp, b, omega);
    // Vector to collect residual on each iteration
    std::vector<T> r = b - (*this) * y_odd;
    std::vector<double> nev {sqrt(r * r)};
    // Total number of iterations
    size_t count = 0;
    // Stop condition
    while (nev[count] > accuracy) {
        // On (i + 1) iteration we rewrite y_(i-1) vector and mu_(i-1),
        // on different iterations y_(i-1) is in different variables, same thing with mu_(i-1)
        // If count is even, y_(i-1) is in y_even and mu_(i-1) is in mu_even,
        if (count % 2 == 0) {
            y_even = -mu_even * y_even + 2 * mu_odd / rho * temp;
            mu_even = 2 / rho * mu_odd - mu_even;
            y_even /= mu_even;
            r = b - (*this) * y_even;
            temp = y_even;
            // Calculate (Py_i + c)
            SSOR_make_iteration(temp, b, omega);
        }
        // If count is odd, y_(i-1) is in y_odd and mu_(i-1) is in mu_odd
        else {
            y_odd = -mu_odd * y_odd + 2 * mu_even / rho * temp;
            mu_odd = 2 / rho * mu_even - mu_odd;
            y_odd /= mu_odd;
            r = b - (*this) * y_odd;
            temp = y_odd;
            // Calculate (Py_i + c)
            SSOR_make_iteration(temp, b, omega);
        }
        nev.push_back(sqrt(r * r));
        count++;
    }
    if (count % 2 == 0)
        return std::make_pair(y_odd, std::make_pair(count, nev));
    else
        return std::make_pair(y_even, std::make_pair(count, nev));
}

// Simple-iteration method
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::Simple_iteration_method
(const std::vector<T>& x_0, const std::vector<T>& b, double tau, double accuracy) const
{
    // Initial approximation
    std::vector<T> x = x_0;
    std::vector<T> r = b - (*this) * x;
    // Total number of iterations
    size_t count = 0;
    // Vector to collect residual on each iteration
    std::vector<double> nev{sqrt(r * r)};
    // Stop condition
    while (nev[count] > accuracy) {
        // Finding x_(i+1)
        r = b - (*this) * x;
        x = x + tau * r;
        nev.push_back(sqrt(r * r));
        count++;
    }
    return std::make_pair(x, std::make_pair(count, nev));
}

// SIM w/ Chebyshev acceleration
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::SIM_Chebyshev_acceleration
(const std::vector<T>& x_0, const std::vector<T>& b, size_t R, double eig_min, double eig_max, double accuracy) const
{
    // Amount of polynomial's roots
    double n = pow(2, R);
    // Vector of polynomial's roots
    std::vector<double> z(static_cast<size_t>(n));
    // Constants
    const double sinA = sin(M_PI / n);
    const double cosA = cos(M_PI / n);
    double sinB = sin(M_PI / (2 * n));
    z[0] = cos(M_PI / (2 * n));
    // Algorithm to find roots
    for (int i = 1; i < n; ++i) {
        z[i] = z[i - 1] * cosA - sinB * sinA;
        sinB = z[i - 1] * sinA + sinB * cosA;
        z[i - 1] = (eig_min + eig_max) / 2 + (eig_max - eig_min) / 2 * z[i - 1];
    }
    z[static_cast<size_t>(n) - 1] = (eig_min + eig_max) / 2 + (eig_max - eig_min) / 2 * z[static_cast<size_t>(n) - 1];
    // Vector of indexes
    std::vector<size_t> indexes(static_cast<size_t>(n));
    indexes[static_cast<size_t>(pow(2, R - 1))] = 1;
    // Indexes reshuffle
    for (int i = 2; i <= R; ++i) {
        auto step = static_cast<int>(pow(2, R - i));
        for (int j = 0; j < n; j += 2 * step) {
            indexes[j + step] = static_cast<size_t>(pow(2, i)) - indexes[j] - 1;
        }
    }
    // Initial approximation
    std::vector<T> x = x_0;
    std::vector<T> r = b - (*this) * x;
    // Total number of iterations
    size_t count = 0;
    // Vector to collect residual on each iteration
    std::vector<double> nev{sqrt(r * r)};
    while (nev[count] > accuracy) {
        // Unique tau on each iteration
        for (auto& id: indexes) {
            // Finding x_(i+1)
            r = b - (*this) * x;
            x = x + (1 / z[id]) * r;
            nev.push_back(sqrt(r * r));
            count++;
            if (nev[count] < accuracy)
                break;
        }
    }
    return std::make_pair(x, std::make_pair(count, nev));
}

// Steepest descent method
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::Steepest_descent_method
(const std::vector<T>& x_0, const std::vector<T>& b, double accuracy) const
{
    // Initial approximation
    std::vector<T> x = x_0;
    std::vector<T> r = b - (*this) * x;
    // Total number of iterations
    size_t count = 0;
    // Vector to collect residual on each iteration
    std::vector<double> nev{sqrt(r * r)};
    // Each iteration's temporary variables
    std::vector<T> Ar(x_0.size());
    double alpha;
    // Stop condition
    while (nev[count] > accuracy) {
        // Ar will be used twice, we should calculate it only once
        Ar = (*this) * r;
        // Finding alpha and x_(i+1)
        alpha = r * r / (r * Ar);
        x = x + alpha * r;
        // Finding r_(i+1)
        r = r - alpha * Ar;
        nev.push_back(sqrt(r * r));
        count++;
    }
    return std::make_pair(x, std::make_pair(count, nev));
}

// Heavy ball method
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::Heavy_ball_method
(const std::vector<T>& x_0, const std::vector<T>& b, double accuracy) const {
    // Vector x_i
    std::vector<T> x = x_0;
    // Vector x_(i-1)
    std::vector<T> x0 = x;
    // Vector d = x_i - x_(i-1) -- direction on previous iteration
    std::vector<T> d(x_0.size());
    // Residual
    std::vector<T> r = (*this) * x - b;
    // Total number of iterations
    size_t count = 0;
    // Vector to collect residual on each iteration
    std::vector<double> nev{sqrt(r * r)};
    // Each iteration's temporary variables
    T rAr, rAd, rr;
    // Find x_1 using steepest descent method
    double alpha = (r * r) / (r * ((*this) * r));;
    x = x - alpha * r;
    // Beta variable
    double beta;
    // Stop condition
    while (nev[count] > accuracy) {
        // Direction on previous iteration
        d = x - x0;
        x0 = x;
        // Variables that are used several times
        rAd = r * ((*this) * d);
        rAr = r * ((*this) * r);
        rr = r * r;
        // Finding alpha and beta
        beta = (rr * rAd - r * d * rAr) / (d * ((*this) * d) * rAr - rAd * rAd);
        alpha = (rr + beta * rAd) / rAr;
        // Finding x_(i+1)
        x = x - alpha * r + beta * d;
        // Finding r_(i+1)
        r = (*this) * x - b;
        nev.push_back(sqrt(r * r));
        count++;
    }
    return std::make_pair(x, std::make_pair(count, nev));
}
#endif //SLAE_4TERM_CSR_MATRIX_H