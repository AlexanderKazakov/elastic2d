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

# Log4cxx library
find_package(Log4cxx REQUIRED)
include_directories(${LOG4CXX_INCLUDE_DIRS})

# Google test library
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
# This is for gtest too
find_package (Threads REQUIRED)

# source files
include_directories(src)
set(SOURCE_MAIN_FILE src/main.cpp)
file(GLOB_RECURSE LIB_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/lib/*.cpp")
set(TEST_MAIN_FILE src/test/TestMain.cpp)
set(TEST_MPI_FILE src/test/TestMPI.cpp)

# libgcm
add_library(gcm SHARED ${LIB_SOURCES})

# The program itself
add_executable(gcm_exe ${SOURCE_MAIN_FILE})
target_link_libraries(gcm_exe gcm ${MPI_CXX_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES})

add_executable(gcm_tests ${TEST_MAIN_FILE})
add_executable(gcm_mpi_tests ${TEST_MPI_FILE})

target_link_libraries(gcm_tests     gcm ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES} ${MPI_CXX_LIBRARIES})
target_link_libraries(gcm_mpi_tests gcm ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES} ${MPI_CXX_LIBRARIES})
