cmake_minimum_required(VERSION 3.22.a)
project(AOATsetlinMachine VERSION 0.1.0)

aux_source_directory(. AOA_LIB_SRCS)
add_library(aoaLib ${AOA_LIB_SRCS})