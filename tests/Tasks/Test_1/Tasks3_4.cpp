#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "gtest/gtest.h"
#include "../../../src/CSR_matrix.h"

using namespace CSR_matrix_space;

TEST(Tasks, Task_3) {
    CSR_matrix<double> m {3, 3, {{0, 0, 10}, {0, 1, 1}, {1, 0, 1}, {1, 1, 7}, {2, 1, 0.1}, {2, 2, 1}}};
    std::vector<double> x0(3);
    std::vector<double> b {20, 30, 1};
    double accuracy = 1e-12;
    double tau = 1e-3;

    std::ofstream out;
    out.open("/home/vljke/Documents/Clion projects/SLAE_4term/tests/Tasks/Task3/data.txt");
    if (out.is_open())
        while (tau < 0.1) {
            auto res = m.Simple_iteration_method(x0, b, tau, accuracy);
            out << res.second.first << " " << tau;
            out << std::endl;
            tau += 1e-4;
    }
}

TEST(Tasks, Task_4) {
    CSR_matrix<double> m {3, 3, {{0, 0, 12}, {0, 1, 17}, {0, 2, 3}, {1, 0, 17}, {1, 1, 15825}, {1, 2, 28},
                                 {2, 0, 3}, {2, 1, 28}, {2, 2, 238}}};
    std::vector<double> x0 {1, 1, 1};
    std::vector<double> b {1, 2, 3};
    std::vector<double> c {1e-10, 1e-12, 1e5};
    double accuracy = 1e-12;

    std::ofstream out1;
    out1.open("/home/vljke/Documents/Clion projects/SLAE_4term/tests/Tasks/Task4/data_Fixed_point.txt");
    if (out1.is_open()) {
        double tau = 1e-4;
        auto res = m.Simple_iteration_method(x0, b, tau, accuracy);
        for (int i = 0; i < res.second.second.size(); ++i) {
            out1 << i << " " << res.second.second[i];
            out1 << std::endl;
        }
        std::cout << res.first << std::endl;
    }

    std::ofstream out2;
    out2.open("/home/vljke/Documents/Clion projects/SLAE_4term/tests/Tasks/Task4/data_Jacobi.txt");
    if (out2.is_open()) {
        auto res = m.Jacobi_method(x0, b, accuracy);
        for (int i = 0; i < res.second.second.size(); ++i) {
            out2 << i << " " << res.second.second[i];
            out2 << std::endl;
        }
        std::cout << res.first << std::endl;
    }

    std::ofstream out3;
    out3.open("/home/vljke/Documents/Clion projects/SLAE_4term/tests/Tasks/Task4/data_Gauss_Seidel.txt");
    if (out3.is_open()) {
        auto res = m.Gauss_Seidel_method(x0, b, accuracy);
        for (int i = 0; i < res.second.second.size(); ++i) {
            out3 << i << " " << res.second.second[i];
            out3 << std::endl;
        }
        std::cout << res.first << std::endl;
    }
}