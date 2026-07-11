cmake_minimum_required(VERSION 3.16)
project(NRLinkSimulator VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

option(NR_BUILD_TESTS "Build unit tests" ON)
option(NR_BUILD_EXAMPLES "Build examples" ON)
option(NR_USE_MKL "Use Intel MKL for linear algebra" OFF)
option(NR_USE_OPENMP "Use OpenMP for parallelization" ON)

if(NR_USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        link_libraries(OpenMP::OpenMP_CXX)
        add_compile_definitions(NR_HAS_OPENMP)
    endif()
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(src)

if(NR_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(NR_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
