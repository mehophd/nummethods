cmake_minimum_required(VERSION 3.15)

project(Lab4 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../../bin")

add_executable(lab4
    main.cpp
    ode.cpp
    task41.cpp
    task42.cpp
)

target_include_directories(lab4 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Корректный вывод кириллицы (исходники в UTF-8) в MSVC.
if(MSVC)
    target_compile_options(lab4 PRIVATE /utf-8)
endif()
