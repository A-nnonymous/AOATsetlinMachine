cmake_minimum_required(VERSION 3.22.1)

project(RSATsetlinMachine)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Compiler-specific flags and definitions
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -O3 -mavx512f")
endif()

add_subdirectory(src/tsetlin)
add_subdirectory(src/RSA)

include_directories(src/tsetlin)
include_directories(src/RSA)

add_executable(xordemo  demo/xor.cpp)
add_executable(rsademo  demo/rsaDemo.cpp)
add_executable(siRNAdemo demo/siRNAdemo.cpp)

target_link_libraries(xordemo tmLib rsaLib)
target_link_libraries(rsademo tmLib rsaLib)
target_link_libraries(siRNAdemo tmLib rsaLib)


