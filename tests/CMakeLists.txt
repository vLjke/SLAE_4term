include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_subdirectory(Googletest)

add_executable(tridiag_test tridiag_test.cpp) # здесь название исполняемого файла и файл в котором лежит тест
add_test(NAME tridiag_test COMMAND tridiag_test.cpp)
target_link_libraries(tridiag_test SLAE_4term gtest gtest_main pthread)

add_executable(CSR_test CSR_test.cpp)
add_test(NAME tCSR_test COMMAND CSR_test.cpp)
target_link_libraries(CSR_test SLAE_4term gtest gtest_main pthread)  # Здесь линковка, тот же исполняемый файл, название проекта и эти же либы


