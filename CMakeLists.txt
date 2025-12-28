cmake_minimum_required(VERSION 3.24)
project(main)

add_executable(${PROJECT_NAME}
    src/ntt.cpp
    src/mult.cpp
    src/main.cpp
)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)
