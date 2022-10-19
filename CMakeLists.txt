cmake_minimum_required(VERSION 3.22.1)

project(AOATsetlinMachine)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Compiler-specific flags and definitions
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -O3 -mavx512f")
endif()

add_subdirectory(src/tsetlin)
add_subdirectory(src/AOA)

include_directories(src/tsetlin)
include_directories(src/AOA)

add_executable(xordemo demo/xor.cpp)

target_link_libraries(xordemo tmLib)

