cmake_minimum_required(VERSION 3.24)
project(main)

add_executable(${PROJECT_NAME}
    src/fft.cpp
    src/product.cpp
    src/tests.cpp
    src/main.cpp
	src/utils.cpp
)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)

target_compile_options(${PROJECT_NAME} PRIVATE
    -O3
    -march=native
)
