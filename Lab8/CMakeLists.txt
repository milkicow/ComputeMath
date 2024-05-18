cmake_minimum_required(VERSION 3.15)

project(exp_computing
        LANGUAGES CXX)

if (${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
    message(FATAL_ERROR "In-source build is forbidden")
endif()

find_package(Boost REQUIRED
             COMPONENTS MPI PROGRAM_OPTIONS HEADERS)


set(CMAKE_CXX_STANDARD          20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

set(SRC_DIR ${PROJECT_SOURCE_DIR})
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR})

add_compile_options(-g)

add_executable(lab1
               ${SRC_DIR}/main.cpp)

target_link_libraries(lab1
                      Boost::program_options
                      Boost::headers
                      Boost::mpi)

target_include_directories(lab1
                           PRIVATE ${INCLUDE_DIR})
