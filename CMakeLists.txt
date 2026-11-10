cmake_minimum_required(VERSION 3.20)
project(options_pricing_engine CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra -Wpedantic)

add_library(ope INTERFACE)
target_include_directories(ope INTERFACE include)

add_executable(dump_cases tools/dump_cases.cpp)
target_link_libraries(dump_cases PRIVATE ope)

enable_testing()
add_executable(test_options tests/test_options.cpp)
target_link_libraries(test_options PRIVATE ope)
add_test(NAME options_scenarios COMMAND test_options)
