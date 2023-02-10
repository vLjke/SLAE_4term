#include <iostream>
#include <cmath>
#include "Tridiagonal_matrix/Tridiagonal_matrix.hpp"
#include "Tridiagonal_matrix_solver/Tridiag_matrix_solver.hpp"
#include "googletest/googlemock/include/gmock/gmock.h"
#include "googletest/googletest/include/gtest/gtest.h"

// Solver for tests
Tridiag_matrix_solver<double> solver {};

// For the first test
Tridiagonal_matrix<double> m1 {5, {2, 3, 1, 4}, {6, 8, 9, 10, 13}, {3, 3, 4, 5}};
std::vector<double> d1 {1, 2, 3, 4, 5};
std::vector<double> correctSolution1 {0.08620, 0.16093, 0.18006, 0.22418, 0.31564};
// For the second test
Tridiagonal_matrix<double> m2 {4, {4, 13, 31}, {99, 154, 33, 34}, {44, 67, 13}};
std::vector<double> d2 {1, 5, 6, 7};
std::vector<double> correctSolution2 {0.03240, -0.05017, 0.18800, 0.03447};
// For the third test
Tridiagonal_matrix<double> m3 {10, {455, 200, 34, 455, 12, 34, 4, 3, 100}, {555, 1234, 788, 45, 522, 34, 88, 9, 12, 122},
                               {111, 123, 333, 1, 50, 12, 22, 5, 9}};
std::vector<double> d3 {1, 2, 3, 1, 2, 3, 1, 2, 3, 1};
std::vector<double> correctSolution3 {0.0013932624, 0.0020426968, -0.0093871732, 0.0299956550, -0.0306405875, 0.0869272728,
                                      0.0343466480, -0.2262741953, 0.7798162331, -0.6309969124};


// 5x5 matrix
TEST(Tridiag_matrix_solver_tests, Subtest_1) {
    for (int i = 0; i < d1.size(); ++i) {
        ASSERT_NEAR(correctSolution1[i], solver.Solution(m1, d1)[i], pow(10, -4))
        << "!!! TEST FAILED ON COORDINATE NUMBER " << i << " !!!" << std::endl;
    }
}
// 4x4 matrix
TEST(Tridiag_matrix_solver_tests, Subtest_2) {
    for (int i = 0; i < d2.size(); ++i) {
        ASSERT_NEAR(correctSolution2[i], solver.Solution(m2, d2)[i], pow(10, -4))
        << "!!! TEST FAILED ON COORDINATE NUMBER " << i << " !!!" << std::endl;
    }
}
// 10x10 matrix
TEST(Tridiag_matrix_solver_tests, Subtest_3) {
    for (int i = 0; i < d3.size(); ++i) {
        ASSERT_NEAR(correctSolution3[i], solver.Solution(m3, d3)[i], pow(10, -10))
        << "!!! TEST FAILED ON COORDINATE NUMBER " << i << " !!!" << std::endl;
    }
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
