#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "gtest/gtest.h"
#include "../../src/Dense_matrix.h"

// Dense matrix 100x100 multiply test
TEST(Dense_matrix_tests, matrix_mult_vector_1) {
    std::string filename("/home/vljke/Documents/Clion projects/SLAE_4term/tests/dense_tests/Dense_matrix_mult_vector_1.txt");
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
        // Getting matrix from file
        std::vector<std::vector<double>> construct_matrix(MN[0]);
        for (int i = 0; i < MN[0]; ++i) {
            construct_matrix[i].resize(MN[1]);
            getline(file, temp);
            std::stringstream ssm(temp);
            for (int j = 0; j < MN[1]; ++j) {
                ssm >> s;
                construct_matrix[i][j] = std::stod(s);
            }
        }
        Dense_matrix<double> m{construct_matrix};

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

// Dense matrix 30x50 multiply test
TEST(Dense_matrix_tests, matrix_mult_vector_2) {
    std::string filename("/home/vljke/Documents/Clion projects/SLAE_4term/tests/dense_tests/Dense_matrix_mult_vector_2.txt");
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
        // Getting matrix from file
        std::vector<std::vector<double>> construct_matrix(MN[0]);
        for (int i = 0; i < MN[0]; ++i) {
            construct_matrix[i].resize(MN[1]);
            getline(file, temp);
            std::stringstream ssm(temp);
            for (int j = 0; j < MN[1]; ++j) {
                ssm >> s;
                construct_matrix[i][j] = std::stod(s);
            }
        }
        Dense_matrix<double> m{construct_matrix};

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

// Dense matrix 10x10 get N elements test
TEST(Dense_matrix_tests, matrix_get_elements) {
    std::string filename("/home/vljke/Documents/Clion projects/SLAE_4term/tests/dense_tests/Dense_matrix_get_elems.txt");
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
        // Getting matrix from file
        std::vector<std::vector<double>> construct_matrix(MN[0]);
        for (int i = 0; i < MN[0]; ++i) {
            construct_matrix[i].resize(MN[1]);
            getline(file, temp);
            std::stringstream ssm(temp);
            for (int j = 0; j < MN[1]; ++j) {
                ssm >> s;
                construct_matrix[i][j] = std::stod(s);
            }
        }
        Dense_matrix<double> m{construct_matrix};
        // Getting number of elements to get from file
        int N;
        std::vector<double> i_j_val(3);
        getline(file, temp);
        std::stringstream ssn(temp);
        ssn >> s;
        N = std::stoi(s);
        // Getting indexes and values from file
        for (int i = 0; i < N; ++i) {
            getline(file, temp);
            std::stringstream ssj(temp);
            for (int j = 0; j < 3; ++j) {
                ssj >> s;
                i_j_val[j] = std::stod(s);
            }
            ASSERT_NEAR(i_j_val[2], m(i_j_val[0], i_j_val[1]), pow(10, -10))
            << "!!! TEST FAILED ON INDEXES NUMBER " << i_j_val[0] << " " << i_j_val[0] << " !!!" << std::endl;
        }
    }
    file.close();
}

// QR decomposition HH algorithm test
TEST(Dense_matrix_tests, QR_decomp_HH) {

}