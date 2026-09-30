cmake_minimum_required(VERSION 3.15)

project(Lab1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../../bin")

add_executable(lab1
    main.cpp
    eigen.cpp
    matrixtools.cpp
    solvemethods.cpp
)

target_include_directories(lab1 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
