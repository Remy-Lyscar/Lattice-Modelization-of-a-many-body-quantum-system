cmake_minimum_required(VERSION 3.20) # old CMake version 
project(QuantumProject) # name of the executable of the project 

# cmake_policy(SET CMP0167 NEW)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ======== Dependencies ========
include_directories(include)   # include the header files (including Eigen library)

# find ARPACK and LAPACK 


# ======== Build ========
file(GLOB_RECURSE SRC_FILES src/*.cpp)

add_executable(${PROJECT_NAME} ${SRC_FILES})
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -march=native)


