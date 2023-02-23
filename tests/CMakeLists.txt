include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_subdirectory(googletest)

file(GLOB_RECURSE tests_cpp *.cpp)
foreach(test ${tests_cpp})
    get_filename_component(name ${test} NAME_WLE)
    add_executable("${name}" ${test})
    add_test(NAME "${name}" COMMAND ${test})
    target_link_libraries("${name}" SLAE_4term gtest gtest_main pthread)
endforeach()


