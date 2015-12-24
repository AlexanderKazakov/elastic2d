cmake_minimum_required(VERSION 2.8.4)
project(gcm)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Profiling (default is not)
option(PROFILE "PROFILE" OFF)
IF(PROFILE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
ENDIF()

# Debug/Release (default is Release)
IF(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
ELSE()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
ENDIF()

# MPI library
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

# GSL library
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

# logging
find_package(Log4cxx REQUIRED)
include_directories(${LOG4CXX_INCLUDE_DIRS})

# The program itself
set(SOURCE_MAIN_FILE src/main.cpp)
file(GLOB_RECURSE LIB_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/lib/*.cpp")
include_directories(src)
add_executable(gcm ${SOURCE_MAIN_FILE} ${LIB_SOURCES})
target_link_libraries(gcm ${MPI_CXX_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES})

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

    add_executable(gcm_tests ${TEST_MAIN_FILE} ${LIB_SOURCES})
    add_executable(gcm_mpi_tests ${TEST_MPI_FILE} ${LIB_SOURCES})

    target_link_libraries(gcm_tests ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES})
    target_link_libraries(gcm_tests ${MPI_CXX_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES})
    target_link_libraries(gcm_mpi_tests ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES})
    target_link_libraries(gcm_mpi_tests ${MPI_CXX_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES})
ENDIF()
