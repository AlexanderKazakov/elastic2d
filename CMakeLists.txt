cmake_minimum_required(VERSION 2.8.4)
project(elastic2d)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
include_directories(src)

# Google test
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
# This is for gtest too
find_package (Threads REQUIRED)

# MPI
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

set(SOURCE_MAIN_FILE src/main.cpp)
set(TEST_MAIN_FILE src/test/TestMain.cpp)
file(GLOB_RECURSE LIB_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/lib/*.cpp")

add_executable(elastic2d ${SOURCE_MAIN_FILE} ${LIB_SOURCES})
add_executable(elastic2d_tests ${TEST_MAIN_FILE} ${LIB_SOURCES})

target_link_libraries(elastic2d ${MPI_CXX_LIBRARIES})

target_link_libraries(elastic2d_tests ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(elastic2d_tests ${MPI_CXX_LIBRARIES})
