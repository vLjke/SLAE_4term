#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include "../../src/Vector_operations.h"
#include "gtest/gtest.h"
#include "../../src/CSR_matrix.h"

using namespace DOK_cell_space;
using namespace CSR_matrix_space;

// DOK cells sort test
TEST(CSR_matrix_tests, DOK_cell_sort) {
    std::string filename("/home/vljke/Documents/Clion projects/SLAE_4term/tests/CSR_tests/DOK_cell_sort.txt");
    std::string temp, s;
    std::ifstream file(filename);
    if (file.is_open()) {
        // Getting number of elements to get from file
        int N;
        getline(file, temp);
        std::stringstream ssn(temp);
        ssn >> s;
        N = std::stoi(s);
        // Getting indexes and values from file
        std::vector<DOK_cell_space::cell<double>> cells(N);
        std::vector<double> i_j_val(3);
        for (int i = 0; i < N; ++i) {
            getline(file, temp);
            std::stringstream ssj(temp);
            for (int j = 0; j < 3; ++j) {
                ssj >> s;
                i_j_val[j] = std::stod(s);
            }
            cells[i].i = i_j_val[0];
            cells[i].j = i_j_val[1];
            cells[i].value = i_j_val[2];
        }
        std::sort(cells.begin(), cells.end());
        std::cout << "(i," << " j)" << std::endl;
        std::cout << "(" << cells[0].i << ", " << cells[1].j << ")" << std::endl;
        for (int i = 1; i < N; ++i) {
            std::cout << "(" << cells[i].i << ", " << cells[i].j << ")" << std::endl;
            ASSERT_TRUE(cells[i - 1].i < cells[i].i || (cells[i - 1].i == cells[i].i && cells[i - 1].j<= cells[i].j))
            << "!!! TEST FAILED ON CELL NUMBER " << i << " !!!" << std::endl;
        }
    }
    file.close();
}

// CSR matrix 20x20 multiply test
TEST(CSR_matrix_tests, matrix_mult_vector_1) {
    std::string filename("/home/vljke/Documents/Clion projects/SLAE_4term/tests/CSR_tests/CSR_matrix_mult_vector_1.txt");
    std::string temp, s;
    std::ifstream file(filename);
    if (file.is_open()) {
        // Getting matrix dimensions
        std::vector<size_t> MN(2);
        getline(file, temp);
        std::stringstream ssi(temp);
        for (int i = 0; i < 2; ++i) {
            ssi >> s;
            MN[i] = std::stoi(s);
        }
        // Getting number of non-zero elements
        int N;
        getline(file, temp);
        std::stringstream ssN(temp);
        ssN >> s;
        N = std::stoi(s);
        // Getting DOK cells from file
        std::vector<DOK_cell_space::cell<double>> cells(N);
        std::vector<double> i_j_val(3);
        for (int i = 0; i < N; ++i) {
            getline(file, temp);
            std::stringstream ssj(temp);
            for (int j = 0; j < 3; ++j) {
                ssj >> s;
                i_j_val[j] = std::stod(s);
            }
            cells[i].i = i_j_val[0];
            cells[i].j = i_j_val[1];
            cells[i].value = i_j_val[2];
        }
        // Creating CSR matrix
        CSR_matrix<double> m {MN[0], MN[1], cells};
        // Getting vector to multiply from file
        std::vector<double> v(MN[1]);
        getline(file, temp);
        std::stringstream ssv(temp);
        for (int i = 0; i < MN[1]; ++i) {
            ssv >> s;
            v[i] = std::stod(s);
        }
        // Getting multiply result vector from file
        std::vector<double> r(MN[0]);
        getline(file, temp);
        std::stringstream ssr(temp);
        for (int i = 0; i < MN[1]; ++i) {
            ssr >> s;
            r[i] = std::stod(s);
        }
        // Checking multiply operator
        std::vector<double> res = m * v;
        for (int i = 0; i < MN[0]; ++i)
            ASSERT_NEAR(res[i], r[i], pow(10, -10))
            << "!!! TEST FAILED ON COORDINATE NUMBER " << i << " !!!" << std::endl;
    }
    file.close();
}

// CSR matrix 100x200 multiply test
TEST(CSR_matrix_tests, matrix_mult_vector_2) {
    std::string filename("/home/vljke/Documents/Clion projects/SLAE_4term/tests/CSR_tests/CSR_matrix_mult_vector_2.txt");
    std::string temp, s;
    std::ifstream file(filename);
    if (file.is_open()) {
        // Getting matrix dimensions
        std::vector<size_t> MN(2);
        getline(file, temp);
        std::stringstream ssi(temp);
        for (int i = 0; i < 2; ++i) {
            ssi >> s;
            MN[i] = std::stoi(s);
        }
        // Getting number of non-zero elements
        int N;
        getline(file, temp);
        std::stringstream ssN(temp);
        ssN >> s;
        N = std::stoi(s);
        // Getting DOK cells from file
        std::vector<DOK_cell_space::cell<double>> cells(N);
        std::vector<double> i_j_val(3);
        for (int i = 0; i < N; ++i) {
            getline(file, temp);
            std::stringstream ssj(temp);
            for (int j = 0; j < 3; ++j) {
                ssj >> s;
                i_j_val[j] = std::stod(s);
            }
            cells[i].i = i_j_val[0];
            cells[i].j = i_j_val[1];
            cells[i].value = i_j_val[2];
        }
        // Creating CSR matrix
        CSR_matrix<double> m {MN[0], MN[1], cells};
        // Getting vector to multiply from file
        std::vector<double> v(MN[1]);
        getline(file, temp);
        std::stringstream ssv(temp);
        for (int i = 0; i < MN[1]; ++i) {
            ssv >> s;
            v[i] = std::stod(s);
        }
        // Getting multiply result vector from file
        std::vector<double> r(MN[0]);
        getline(file, temp);
        std::stringstream ssr(temp);
        for (int i = 0; i < MN[0]; ++i) {
            ssr >> s;
            r[i] = std::stod(s);
        }
        // Checking multiply operator
        std::vector<double> res = m * v;
        for (int i = 0; i < MN[0]; ++i)
            ASSERT_NEAR(res[i], r[i], pow(10, -10))
            << "!!! TEST FAILED ON COORDINATE NUMBER " << i << " !!!" << std::endl;
    }
    file.close();
}

// CSR matrix get N elements
TEST(CSR_matrix_tests, matrix_get_element) {
    std::string filename("/home/vljke/Documents/Clion projects/SLAE_4term/tests/CSR_tests/CSR_matrix_get_elems.txt");
    std::string temp, s;
    std::ifstream file(filename);
    if (file.is_open()) {
        // Getting matrix dimensions
        std::vector<size_t> MN(2);
        getline(file, temp);
        std::stringstream ssi(temp);
        for (int i = 0; i < 2; ++i) {
            ssi >> s;
            MN[i] = std::stoi(s);
        }
        // Getting number of non-zero elements
        int N;
        getline(file, temp);
        std::stringstream ssN(temp);
        ssN >> s;
        N = std::stoi(s);
        // Getting DOK cells from file
        std::vector<DOK_cell_space::cell<double>> cells(N);
        std::vector<double> i_j_val(3);
        for (int i = 0; i < N; ++i) {
            getline(file, temp);
            std::stringstream ssi(temp);
            for (int j = 0; j < 3; ++j) {
                ssi >> s;
                i_j_val[j] = std::stod(s);
            }
            cells[i].i = i_j_val[0];
            cells[i].j = i_j_val[1];
            cells[i].value = i_j_val[2];
        }
        // Creating CSR matrix
        CSR_matrix<double> m {MN[0], MN[1], cells};
        // Getting number of elements to get from file
        int n;
        getline(file, temp);
        std::stringstream ssn(temp);
        ssn >> s;
        n = std::stoi(s);
        // Getting indexes and values from file
        for (int i = 0; i < n; ++i) {
            getline(file, temp);
            std::stringstream ssj(temp);
            for (int j = 0; j < 3; ++j) {
                ssj >> s;
                i_j_val[j] = std::stod(s);
            }
            ASSERT_NEAR(i_j_val[2], m(i_j_val[0], i_j_val[1]), pow(10, -10))
            << "!!! TEST FAILED ON INDEXES i = " << i_j_val[0] << ", j = " << i_j_val[1] << " !!!" << std::endl;
        }
    }
    file.close();
}