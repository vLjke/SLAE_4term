#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "gtest/gtest.h"
#include "../../../src/CSR_matrix.h"

using namespace CSR_matrix_space;

TEST(Tasks, Task_1) {
    double accuracy = 1e-13;

    double a = 10;
    double b = 25;

    size_t N = 289;
    size_t L = 17;

    std::vector<DOK_cell_space::cell<double>> cells;
    std::vector<double> x_0(N);
    std::vector<double> c(N, 2);
    for (size_t i = 1; i < N - 1; ++i) {
        cells.push_back({i, i, 2 * b});
        cells.push_back({i, i - 1, a});
        cells.push_back({i, i + 1, a});
    }
    cells.push_back({0, 0, 2 * b});
    cells.push_back({N - 1, N - 1, 2 * b});
    cells.push_back({0, 1, a});
    cells.push_back({N - 1, N - 2, a});
    for (size_t i = 0; i < N - L; ++i) {
        cells.push_back({i + L, i, a});
        cells.push_back({i, i + L, a});
    }
    CSR_matrix<double> A {N, N, cells};

    double eig_max = 2 * (b + 2 * a * std::cos(M_PI / static_cast<double>(L + 1)));
    double eig_min = 2 * (b - 2 * a * std::cos(M_PI / static_cast<double>(L + 1)));

    std::ofstream out;
    out.open("/home/vljke/Documents/Clion projects/SLAE_4term/tests/Tasks/Test_2/Task1/residuals.txt");

    // 1)
    double tau = 1 / eig_max;
    auto resGrad = A.Simple_iteration_method(x_0, c, tau, accuracy);
    for (auto& i : resGrad.second.second)
        out << i << " ";
    out << std::endl;

    // 2)
    tau = 2 / (eig_max + eig_min);
    auto SIMOpt = A.Simple_iteration_method(x_0, c, tau, accuracy);
    for (auto& i : SIMOpt.second.second)
        out << i << " ";
    out << std::endl;

    // 3)
    size_t r = 3;
    auto SIMChebAcceleration = A.SIM_Chebyshev_acceleration(x_0, c, r, eig_min, eig_max, accuracy);
    for (auto& i : SIMChebAcceleration.second.second)
        out << i << " ";
    out << std::endl;

    // 4)
    double omega = 0.5;
    auto Ssor = A.SSOR_method(x_0, c, omega, accuracy);
    for (auto& i : Ssor.second.second)
        out << i << " ";
    out << std::endl;

    out.close();
    // different eig_min eig_max
    out.open("/home/vljke/Documents/Clion projects/SLAE_4term/tests/Tasks/Test_2/Task1/different_eig_max.txt");
    double eig = eig_max;
    while (eig < eig_max * 50) {
        auto tempRes = A.SIM_Chebyshev_acceleration(x_0, c, r, eig_min, eig, accuracy);
        out << eig - eig_min << " " << tempRes.second.first << std::endl;
        eig += eig_max / 100;
    }

}

TEST(Tasks, Task_2) {
    size_t N = 4;
    std::vector<double> x_0(N);
    double accuracy = 1e-13;
    std::vector<double> b(4, 6);
    std::vector<double> c = b;
    std::vector<DOK_cell_space::cell<double>> cells;
    cells.push_back({0, 0, 12});
    cells.push_back({1, 1, 14});
    cells.push_back({2, 2, 16});
    cells.push_back({3, 3, 18});
    CSR_matrix<double> A{N, N, cells};
    double eig_max = 18;
    double eig_min = 12;
    // 1)
    double tau = 0.9 * 2 / eig_max;
    auto resGrad = A.Simple_iteration_method(x_0, b, tau, accuracy);
    std::cout << resGrad.second.first << std::endl;
    // 2)
    tau = 2 / (eig_max + eig_min);
    auto SIMOpt = A.Simple_iteration_method(x_0, b, tau, accuracy);
    std::cout << SIMOpt.second.first << std::endl;
    // 3)
    auto SteepestGrad = A.Steepest_descent_method(x_0, b, accuracy);
    std::cout << SteepestGrad.second.first << std::endl;
    // 4)
    size_t r = 3;
    auto SIMChebAcceleration = A.SIM_Chebyshev_acceleration(x_0, b, r, eig_min, eig_max, accuracy);
    std::cout << SIMChebAcceleration.second.first << std::endl;
    // 5)
    auto CG = A.CG_method(x_0, b, accuracy);
    std::cout << CG.second.first << std::endl;
}