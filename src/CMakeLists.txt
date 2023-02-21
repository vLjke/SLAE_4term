cmake_minimum_required(VERSION 3.24)
project(SLAE_4term)

set(CMAKE_CXX_STANDARD 20)

file(GLOB_RECURSE source *.h *.hpp)

add_library(SLAE_4term INTERFACE ${source})