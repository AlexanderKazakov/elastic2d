cmake_minimum_required(VERSION 2.8.4)
project(elastic2d)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Profiling (default is not)
option(PROFILE "PROFILE" OFF)
IF(PROFILE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
ENDIF()

# Debug/Release (default is Release)
IF(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
ELSE()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
ENDIF()

# MPI library
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

# The program itself
set(SOURCE_MAIN_FILE src/main.cpp)
file(GLOB_RECURSE LIB_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/lib/*.cpp")
include_directories(src)
add_executable(elastic2d ${SOURCE_MAIN_FILE} ${LIB_SOURCES})
target_link_libraries(elastic2d ${MPI_CXX_LIBRARIES})

# Compile tests or not (default is not)
option(COMPILE_TESTS "COMPILE_TESTS" OFF)
IF(COMPILE_TESTS)
    # Google test
    find_package(GTest REQUIRED)
    include_directories(${GTEST_INCLUDE_DIRS})
    # This is for gtest too
    find_package (Threads REQUIRED)

    set(TEST_MAIN_FILE src/test/TestMain.cpp)
    set(TEST_MPI_FILE src/test/TestMPI.cpp)

    add_executable(elastic2d_tests ${TEST_MAIN_FILE} ${LIB_SOURCES})
    add_executable(elastic2d_mpi_tests ${TEST_MPI_FILE} ${LIB_SOURCES})

    target_link_libraries(elastic2d_tests ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    target_link_libraries(elastic2d_tests ${MPI_CXX_LIBRARIES})
    target_link_libraries(elastic2d_mpi_tests ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    target_link_libraries(elastic2d_mpi_tests ${MPI_CXX_LIBRARIES})
ENDIF()
