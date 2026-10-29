cmake_minimum_required(VERSION 3.15)

project(Lab3 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../../bin")

add_executable(lab3
    main.cpp
    interpolation.cpp
    spline.cpp
    lsm.cpp
    differentiation.cpp
    integration.cpp
)

target_include_directories(lab3 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})