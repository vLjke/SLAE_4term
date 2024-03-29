#ifndef SLAE_4TERM_CSR_MATRIX_H
#define SLAE_4TERM_CSR_MATRIX_H

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <ranges>
#include <algorithm>
#include "Vector_operations.h"
#include "Dense_matrix.h"

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
        std::vector<T> SSOR_make_iteration(const std::vector<T>& x_0, const std::vector<T>& b, double omega) const {
            std::vector<T> x = x_0;
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
            return x;
        }
        std::vector<T> Sym_GS_make_iteration(const std::vector<T>& x_0, const std::vector<T>& b) const {
            std::vector<T> x = x_0;
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
            return x;
        }

        std::vector<T> make_Arnoldi_iteration(std::vector<std::vector<T>>& v, Dense_matrix<T>& H, size_t i) const {
            std::vector<T> res = (*this) * v[i];
            for (int j = 0; j < i+1; ++j) {
                H(j, i) = v[j] * res;
                res = res - H(j, i) * v[j];
            }
            H(i+1, i) = std::sqrt(res * res);
            res = res / H(i+1, i);
            return res;
        }
        void make_GMRES_rotation(Dense_matrix<T>& H, const std::vector<std::pair<double, double>>& SinCos, size_t k, size_t j) const {
            double temp = H(k, j);
            H(k, j) = SinCos[k].second * H(k, j) + SinCos[k].first * H(k+1, j);
            H(k+1, j) = -SinCos[k].first * temp + SinCos[k].second * H(k+1, j);
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
        // Getters
        const std::vector<T>& getValues() const {return this->values;}
        const std::vector<size_t>& getCols() const {return this->cols;}
        const std::vector<size_t>& getRows() const {return this->rows;}

        // Transposed matrix
        CSR_matrix<T> Transposed() const;

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
        // SSOR method
        std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> SSOR_method
        (const std::vector<T>& x_0, const std::vector<T>& b, double omega, double accuracy) const;
        // SSOR w/ Chebyshev acceleration method
        std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> SSOR_Cheb_accel_method
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
        // Conjugate gradient method
        std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CG_method
        (const std::vector<T>& x_0, const std::vector<T>& b, double accuracy) const;
        // BiCG method
        std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> BiCG_method
        (const std::vector<T>& x_0, const std::vector<T>& b, double accuracy) const;
        // CMRES(n) method
        std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> GMRESn_method
        (const std::vector<T>& x_0, const std::vector<T>& b, size_t n, double accuracy) const;

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

template<typename T>
CSR_matrix_space::CSR_matrix<T> CSR_matrix_space::CSR_matrix<T>::Transposed() const {
    // Vector to construct transposed matrix
    std::vector<DOK_cell_space::cell<T>> cells;
    cells.reserve(this->values.size());
    //
    for (int k = 0; k < this->M; ++k)
        for (int s = this->rows[k]; s < this->rows[k + 1]; ++s)
            cells.push_back({this->cols[s], static_cast<size_t>(k), this->values[s]});
    return CSR_matrix<T>{this->N, this->M, cells};
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
    std::vector<double> nev {std::sqrt(r * r)};
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
        nev.push_back(std::sqrt(r * r));
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
    std::vector<double> nev {std::sqrt(r * r)};
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
        nev.push_back(std::sqrt(r * r));
        count++;
    }
    return std::make_pair(x, std::make_pair(count, nev));
}

// Symmetrical Gauss-Seidel w/ Chebyshev acceleration method
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::Symmetrical_GS_method
(const std::vector<T>& x_0, const std::vector<T>& b, double rho, double accuracy) const
{
    // y_(i-1) && y_i
    std::vector<T> y_0 = x_0;
    std::vector<T> y_1 = Sym_GS_make_iteration(y_0, b);
    // mu_(i-1) && mu_i
    double mu_0 = 1;
    double mu_1 = 1 / rho;
    // Vector to keep (Py + c)
    std::vector<T> temp = Sym_GS_make_iteration(y_1, b);
    // Vector to collect residual on each iteration
    std::vector<T> r = b - (*this) * y_1;
    std::vector<double> nev {std::sqrt(r * r)};
    // Total number of iterations
    size_t count = 0;
    // Stop condition
    while (nev[count] > accuracy) {
        // On (i + 1) iteration we rewrite y_(i-1) vector and mu_(i-1),
        y_0 = -mu_0 * y_0 + 2 * mu_1 / rho * temp;
        mu_0 = 2 / rho * mu_1 - mu_0;
        y_0 /= mu_0;
        // Calculate (Py_i + c)
        temp = Sym_GS_make_iteration(y_0, b);
        // Swap vectors and mu
        y_1.swap(y_0);
        std::swap(mu_0, mu_1);
        // Calculate residual
        r = b - (*this) * y_1;
        nev.push_back(std::sqrt(r * r));
        count++;
    }
    return std::make_pair(y_1, std::make_pair(count, nev));
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
    std::vector<double> nev {std::sqrt(r * r)};
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
        nev.push_back(std::sqrt(r * r));
        count++;
    }
    return std::make_pair(x, std::make_pair(count, nev));
}

// SSOR method
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::SSOR_method
(const std::vector<T>& x_0, const std::vector<T>& b, double omega, double accuracy) const {
    // Initial approximation
    std::vector x = x_0;
    std::vector<T> r = b - (*this) * x;
    // Total number of iterations
    size_t count = 0;
    // Vector to collect residual on each iteration
    std::vector<double> nev {std::sqrt(r * r)};
    // Stop condition
    while (nev[count] > accuracy) {
        // Finding x_(i+1)
        x = SSOR_make_iteration(x, b, omega);
        r = b - (*this) * x;
        // Finding residual
        nev.push_back(std::sqrt(r * r));
        count++;
    }
    return std::make_pair(x, std::make_pair(count, nev));
}

// SSOR w/ Chebyshev acceleration method
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::SSOR_Cheb_accel_method
(const std::vector<T>& x_0, const std::vector<T>& b, double omega, double rho, double accuracy) const
{
    // y_(i-1) && y_i
    std::vector<T> y_0 = x_0;
    std::vector<T> y_1 = SSOR_make_iteration(y_0, b, omega);
    // mu_(i-1) && mu_i
    double mu_0 = 1;
    double mu_1 = 1 / rho;
    // Vector to keep (Py + c)
    std::vector<T> temp = SSOR_make_iteration(y_1, b, omega);
    // Vector to collect residual on each iteration
    std::vector<T> r = b - (*this) * y_1;
    std::vector<double> nev {std::sqrt(r * r)};
    // Total number of iterations
    size_t count = 0;
    // Stop condition
    while (nev[count] > accuracy) {
        // On (i + 1) iteration we rewrite y_(i-1) vector and mu_(i-1),
        y_0 = -mu_0 * y_0 + 2 * mu_1 / rho * temp;
        mu_0 = 2 / rho * mu_1 - mu_0;
        y_0 /= mu_0;
        // Calculate (Py_i + c)
        temp = SSOR_make_iteration(y_0, b, omega);
        // Swap vectors and mu
        y_1.swap(y_0);
        std::swap(mu_0, mu_1);
        // Calculate residual
        r = b - (*this) * y_1;
        nev.push_back(std::sqrt(r * r));
        count++;
    }
    return std::make_pair(y_1, std::make_pair(count, nev));
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
    std::vector<double> nev{std::sqrt(r * r)};
    // Stop condition
    while (nev[count] > accuracy) {
        // Finding x_(i+1)
        r = b - (*this) * x;
        x = x + tau * r;
        nev.push_back(std::sqrt(r * r));
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
    double n = std::pow(2, R);
    // Vector of polynomial's roots
    std::vector<double> z(static_cast<size_t>(n));
    // Constants
    const double sinA = std::sin(M_PI / n);
    const double cosA = std::cos(M_PI / n);
    double sinB = std::sin(M_PI / (2 * n));
    z[0] = std::cos(M_PI / (2 * n));
    // Algorithm to find roots
    for (int i = 1; i < n; ++i) {
        z[i] = z[i - 1] * cosA - sinB * sinA;
        sinB = z[i - 1] * sinA + sinB * cosA;
        z[i - 1] = (eig_min + eig_max) / 2 + (eig_max - eig_min) / 2 * z[i - 1];
    }
    z[static_cast<size_t>(n) - 1] = (eig_min + eig_max) / 2 + (eig_max - eig_min) / 2 * z[static_cast<size_t>(n) - 1];
    // Vector of indexes
    std::vector<size_t> indexes(static_cast<size_t>(n));
    indexes[static_cast<size_t>(std::pow(2, R - 1))] = 1;
    // Indexes reshuffle
    for (int i = 2; i <= R; ++i) {
        auto step = static_cast<int>(std::pow(2, R - i));
        for (int j = 0; j < n; j += 2 * step) {
            indexes[j + step] = static_cast<size_t>(std::pow(2, i)) - indexes[j] - 1;
        }
    }
    // Initial approximation
    std::vector<T> x = x_0;
    std::vector<T> r = b - (*this) * x;
    // Total number of iterations
    size_t count = 0;
    // Vector to collect residual on each iteration
    std::vector<double> nev{std::sqrt(r * r)};
    while (nev[count] > accuracy) {
        // Unique tau on each iteration
        for (auto& id: indexes) {
            // Finding x_(i+1)
            r = b - (*this) * x;
            x = x + (1 / z[id]) * r;
            nev.push_back(std::sqrt(r * r));
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
    std::vector<double> nev{std::sqrt(r * r)};
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
        nev.push_back(std::sqrt(r * r));
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
    std::vector<double> nev{std::sqrt(r * r)};
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
        nev.push_back(std::sqrt(r * r));
        count++;
    }
    return std::make_pair(x, std::make_pair(count, nev));
}

// Conjugate gradient method
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::CG_method
(const std::vector<T>& x_0, const std::vector<T>& b, double accuracy) const {
    // Vector x_i
    std::vector<T> x = x_0;
    // Residual
    std::vector<T> r = (*this) * x - b;
    // Vector d, alpha and beta variables
    std::vector<T> d = r;
    double beta;
    double alpha;
    // We use Ad twice, will calculate it only once
    std::vector<T> Ad(x_0.size());
    // Total number of iterations
    size_t count = 0;
    // Vector to collect residual on each iteration
    std::vector<double> nev{std::sqrt(r * r)};
    while (nev[count] > accuracy) {
        // N iterations to find solution (theoretically it's enough)
        for (int i = 0; i < x_0.size(); ++i) {
            // Check whether we reached solution
            if (!isZero(r)) {
                Ad = *this * d;
                // Finding alpha_i and x_(i+1) via d_i and r_i
                alpha = r * r / (d * Ad);
                x = x - alpha * d;
                // Finding beta_(i+1) via r_i and r_(i+1)
                beta = 1 / (r * r);
                // Calculating r_(i+1)
                r = r - alpha * Ad;
                beta *= (r * r);
                // Finding d_(i+1)
                d = r + beta * d;
                // Collecting residual and increasing counter
                nev.push_back(std::sqrt(r * r));
                count++;
            }
            else
                return std::make_pair(x, std::make_pair(count, nev));
        }
    }
    return std::make_pair(x, std::make_pair(count, nev));
}

// BiCG method
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::BiCG_method
(const std::vector<T>& x_0, const std::vector<T>& b, double accuracy) const {
    // Transposed matrix
    CSR_matrix<T> A_T = this->Transposed();
    // Vector x_i
    std::vector<T> x = x_0;
    // Initial residual
    std::vector<T> r_0 = (*this) * x - b;
    // Total number of iterations
    size_t count = 0;
    // Vector to collect residual on each iteration
    std::vector<double> nev{std::sqrt(r_0 * r_0)};
    // Restart if needed
    while (nev[count] > accuracy) {
        // r and p vectors
        std::vector<T> r = (*this) * x - b;
        std::vector<T> r_wave = r, p = r, p_wave = r;
        // Temporary variables
        T q, t, rr_current, rr_prev;
        std::vector<T> Ap;
        // n iterations cycle
        for (int i = 0; i < this->M; ++i) {
            // Updating rr_(i-1)
            rr_current = r * r_wave;
            // Restart if r*r_current == 0
            if (rr_current == 0) break;
            // Except for the first iteration
            if (i != 0) {
                t = rr_current / rr_prev;
                p = r + t * p;
                p_wave = r_wave + t * p_wave;
            }
            // On the first iteration
            else
                rr_prev = rr_current;
            // Finding temporary variables
            Ap = (*this) * p;
            q = rr_current / (r_wave * Ap);
            // Finding x_i and r_i
            x = x - q * p;
            r = r - q * Ap;
            r_wave = r_wave - q * (A_T * p_wave);
            // Updating rr_(i-2)
            rr_prev = rr_current;
            // Collecting residual and increasing counter
            nev.push_back(std::sqrt(r * r));
            count++;
        }
    }
    return std::make_pair(x, std::make_pair(count, nev));
}

// GMRES(n) method
template<typename T>
std::pair<std::vector<T>, std::pair<size_t, std::vector<double>>> CSR_matrix_space::CSR_matrix<T>::GMRESn_method
(const std::vector<T> &x_0, const std::vector<T> &b, size_t n, double accuracy) const {
    std::vector<T> x = x_0;
    // Vector to keep basis vectors
    std::vector<std::vector<T>> v(n+1, std::vector<T>(this->M));
    // Hessenberg matrix
    Dense_matrix<T> H{n+1, n};
    // Vector to keep sinA && cosA
    std::vector<std::pair<double, double>> SinCos(n);
    // Initial residual
    std::vector<T> r_0 = (*this) * x_0 - b;
    std::vector<double> nev{std::sqrt(r_0 * r_0)};
    // Restart if needed
    while (nev[nev.size() - 1] > accuracy) {
        // Right side vector
        std::vector<T> e(n+1);
        // Residual vector
        std::vector<T> r = (*this) * x - b;
        e[0] = std::sqrt(r * r);
        // First basis vector
        v[0] = r / e[0];
        // Make n iterations
        for (int i = 0; i < n; ++i) {
            // Finding next basis vector
            v[i+1] = make_Arnoldi_iteration(v, H, i);
            // Making previous rotations
            for (int k = 0; k < i; ++k) {
                make_GMRES_rotation(H, SinCos, k, i);
            }
            // Making current iteration rotation
            double temp = std::sqrt(H(i + 1, i) * H(i + 1, i) + H(i, i) * H(i, i));
            SinCos[i].first = H(i + 1, i) / temp;
            SinCos[i].second = H(i, i) / temp;
            H(i, i) = SinCos[i].second * H(i, i) + SinCos[i].first * H(i + 1, i);
            H(i + 1, i) = 0;
            // Rotating right side vector
            temp = e[i];
            e[i] = SinCos[i].second * e[i] + SinCos[i].first * e[i + 1];
            e[i + 1] = -SinCos[i].first * temp + SinCos[i].second * e[i + 1];
            // Collecting residual on each iteration
            nev.push_back(std::abs(e[i + 1]));
        }
        // Vector to keep coordinates in Krylov's space
        std::vector<T> y(n);
        // Reverse Gauss technique
        for (int i = n - 1; i >= 0; --i) {
            // Finding y_i
            y[i] = e[i];
            for (int j = n - 1; j > i; --j) {
                y[i] -= H(i, j) * y[j];
            }
            y[i] /= H(i, i);
            // Finding solution
            x = x - y[i] * v[i];
        }
    }
    return std::make_pair(x, std::make_pair(nev.size() - 1, nev));
}


#endif //SLAE_4TERM_CSR_MATRIX_H