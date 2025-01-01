cmake_minimum_required(VERSION 3.20) # old CMake version 
project(QuantumProject) # name of the executable of the project 

# cmake_policy(SET CMP0167 NEW)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ======== Dependencies ========
include_directories(include)   # include the header files (including 
                               # Eigen and arpack++ library)

# find ARPACK and LAPACK 
find_package(PkgConfig REQUIRED)
pkg_check_modules(LAPACK REQUIRED lapack)  # if found, create LAPACK_INCLUDE_DIRS and other variables
pkg_check_modules(ARPACK REQUIRED arpack)

include_directories(${LAPACK_INCLUDE_DIRS})
include_directories(${ARPACK_INCLUDE_DIRS})


# Find SuperLU
pkg_check_modules(SuperLU REQUIRED superlu)
include_directories(${SuperLU_INCLUDE_DIRS})



# ======== Build ========
file(GLOB_RECURSE SRC_FILES src/*.cpp)

add_executable(${PROJECT_NAME} ${SRC_FILES})
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -march=native)

# Link against ARPACK and LAPACK 
target_link_libraries(${PROJECT_NAME} ${LAPACK_LIBRARIES} ${ARPACK_LIBRARIES})
# During the linking process, the linker combines object files (machine code, 
# but the references, variables, ... have not been dealt with yet) from 
# our source code with object files from the ARPACK and LAPACK libraries to 
# create the final executable. 