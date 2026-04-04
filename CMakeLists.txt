cmake_minimum_required(VERSION 3.10)
project(NanoSat_MEKF C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_executable(mekf_test
    source/main.c
    source/mekf_wb.c
)

target_include_directories(mekf_test PRIVATE source)

target_link_libraries(mekf_test m)