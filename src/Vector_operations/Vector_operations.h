#ifndef SLAE_4TERM_VECTOR_OPERATIONS_H
#define SLAE_4TERM_VECTOR_OPERATIONS_H

#include<iostream>
#include<vector>

template<typename T>
std::vector<T>& operator+(std::vector<T>& right) {
    return right;
}

template<typename T>
std::vector<T>& operator-(std::vector<T>& right) {
    for (auto& v: right)
        v = -v;
    return right;
}

template<typename T>
std::vector<T> operator+(const std::vector<T>& left, const std::vector<T>& right) {
    std::vector<T> res(left.size());
    for (int i = 0; i < left.size(); ++i)
        res[i] = left[i] + right[i];
    return res;
}

template<typename T>
std::vector<T> operator-(const std::vector<T>& left, const std::vector<T>& right) {
    std::vector<T> res(left.size());
    for (int i = 0; i < left.size(); ++i)
        res[i] = left[i] - right[i];
    return res;
}

template<typename T>
std::vector<T>& operator*(T& left, std::vector<T>& right) {
    for (auto& v: right)
        v *= left;
    return right;
}

template<typename T>
std::vector<T>& operator*(std::vector<T>& left, T& right) {
    for (auto& v: left)
        v *= right;
    return left;
}

template<typename T>
std::vector<T>& operator*(T&& left, std::vector<T>& right) {
    for (auto& v: right)
        v *= left;
    return right;
}

template<typename T>
std::vector<T>& operator*(std::vector<T>& left, T&& right) {
    for (auto& v: left)
        v *= right;
    return left;
}

template<typename T>
T operator*(std::vector<T>& left, std::vector<T>& right) {
    T res = 0;
    for (int i = 0; i < left.size(); ++i)
        res += left[i] * right[i];
    return res;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    for (const auto& v: vec)
        os << v << " ";
    return os;
}


#endif //SLAE_4TERM_VECTOR_OPERATIONS_H
