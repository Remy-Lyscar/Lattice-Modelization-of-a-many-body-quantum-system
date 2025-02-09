cmake_minimum_required(VERSION 3.20) # old CMake version 
project(QuantumProject) # name of the executable of the project 

# cmake_policy(SET CMP0167 NEW)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# ======== Dependencies ========
include_directories(${CMAKE_SOURCE_DIR}/include)  # include the header files

include_directories(${CMAKE_SOURCE_DIR}/external/spectra/include)


# ======== Doxygen ========
find_package(Doxygen)

if (DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(doc_doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
else (DOXYGEN_FOUND)
    message(STATUS "Doxygen not found, unable to generate documentation")
endif (DOXYGEN_FOUND)


# ======== Build ========
file(GLOB_RECURSE SRC_FILES src/*.cpp)

add_executable(${PROJECT_NAME} ${SRC_FILES})

# Add compiler warning flags
add_compile_options(-Wall -Wextra -Wunused-variable)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -march=native)


# During the linking process, the linker combines object files (machine code, 
# but the references, variables, ... have not been dealt with yet) from 
# our source code with object files from the ARPACK and LAPACK libraries to 
# create the final executable. 