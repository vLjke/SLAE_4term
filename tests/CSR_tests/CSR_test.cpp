#include <iostream>
#include <cmath>
#include "../../src/Vector_operations.h"
#include "gtest/gtest.h"
#include "../../src/CSR_matrix.h"

using namespace DOK_cell_space;
using namespace CSR_matrix_space;

std::vector<cell<double>> d1 {{1, 1, 1}, {1, 12, 3}, {0, 0, 1}, {18, 0, 1}, {0, 1, 0}};
std::vector<double> r1 {1, 2};
CSR_matrix<double> m1 {2, 2, d1};

TEST(CSR_matrix_tests, DOK_cell_sort) {
    std::sort(d1.begin(), d1.end());
    for (int i = 0; i < d1.size() - 1; ++i) {
        std::cout << d1[i].i << " " << d1[i].j << std::endl;
        ASSERT_TRUE(d1[i] < d1[i + 1]) <<
        "!!! TEST FAILED ON COORDINATE NUMBER " << i << " !!!" << std::endl;
    }
}

TEST(CSR_matrix_tests, matrix_mult_vector) {
    auto res = m1 * r1;
    for (int i = 0; i < r1.size(); ++i) {
        std::cout <<"Initial coordinate: " << r1[i] << "; Result coordinate: " << res[i] << std::endl;
        ASSERT_NEAR(r1[i], res[i], pow(10, -10)) <<
        "!!! TEST FAILED ON COORDINATE NUMBER " << i << " !!!" << std::endl;
    }
}