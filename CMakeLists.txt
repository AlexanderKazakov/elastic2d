cmake_minimum_required(VERSION 2.8.4)
project(gcm)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")

set(WARNING_FLAGS "-Wall -Wextra -Wsign-conversion -Wconversion -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -Woverloaded-virtual")
set(SAFE_MATH_OPERATIONS "-frounding-math -fsignaling-nans") # it's also necessary for CGAL working
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -DBOOST_PARAMETER_MAX_ARITY=12 ${WARNING_FLAGS} ${SAFE_MATH_OPERATIONS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp -lmpfr -lgmp")

# Show make output (default is not)
option(VERBOSE_MAKE "VERBOSE_MAKE" OFF)
IF(VERBOSE_MAKE)
    set(CMAKE_VERBOSE_MAKEFILE ON)
ENDIF()
# Additional unsafe optimization (default is not)
option(ADDITIONAL_UNSAFE_OPTIMIZE "ADDITIONAL_UNSAFE_OPTIMIZE" OFF)
set(ADDITIONAL_OPTIMIZATION_FLAGS "-march=native -Ofast -msse3 -ffast-math -funroll-loops")
IF(ADDITIONAL_UNSAFE_OPTIMIZE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_OPTIMIZATION_FLAGS}")
ENDIF()

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
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
# CGAL library
set(CGAL_FIND_VERSION "4.8")
find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS})
# VTK library
find_package(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
if(VTK_MAJOR_VERSION EQUAL 5)
    set(CONFIG_VTK_5 1)
elseif (VTK_MAJOR_VERSION EQUAL 6)
    set(CONFIG_VTK_6 1)
elseif (VTK_MAJOR_VERSION EQUAL 7)
    set(CONFIG_VTK_7 1)
else()
    message(FATAL_ERROR "Unsupported VTK version found: ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
endif()

# cmake generated config
set(CMAKE_CONFIG_SHARE_GCM "${CMAKE_INSTAL_PREFIX}/share/gcm")

set(CMAKE_CONFIG_ENABLE_LOGGING       1)
set(CMAKE_CONFIG_ENABLE_LOGGING_FATAL 1)
set(CMAKE_CONFIG_ENABLE_LOGGING_ERROR 1)
set(CMAKE_CONFIG_ENABLE_LOGGING_WARN  1)
set(CMAKE_CONFIG_ENABLE_LOGGING_INFO  1)
set(CMAKE_CONFIG_ENABLE_LOGGING_DEBUG 0)
set(CMAKE_CONFIG_ENABLE_LOGGING_TRACE 0)

IF(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    set(CMAKE_CONFIG_ENABLE_LOGGING_INFO  1)
ELSEIF(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CONFIG_ENABLE_LOGGING_INFO  1)
    set(CMAKE_CONFIG_ENABLE_LOGGING_DEBUG 1)
    set(CMAKE_CONFIG_ENABLE_LOGGING_TRACE 1)
ENDIF()

set(CMAKE_CONFIG_ENABLE_ASSERTIONS    1)
set(LIBGCM_DOUBLE_PRECISION           1)
configure_file("${CMAKE_CURRENT_LIST_DIR}/cmake/config.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/generated_sources/libgcm/config.hpp")

# source files
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/" "${CMAKE_CURRENT_BINARY_DIR}/generated_sources/")
FILE(GLOB_RECURSE HEADERS "${CMAKE_CURRENT_LIST_DIR}/src/*.hpp")
file(GLOB_RECURSE GCM_EXE_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/launcher/*.cpp")
file(GLOB_RECURSE LIBGCM_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/libgcm/*.cpp")
file(GLOB_RECURSE LIBCGALMESHER_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/libcgalmesher/*.cpp")
file(GLOB_RECURSE SEQUENCE_TEST_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/test/sequence/*.cpp")
set(TEST_MPI_FILE src/test/TestMPI.cpp)

#libcgalmesher
add_library(cgalmesher SHARED ${LIBCGALMESHER_SOURCES})

# libgcm
add_library(gcm SHARED ${LIBGCM_SOURCES})

# gcm_exe
add_executable(gcm_exe ${GCM_EXE_SOURCES} ${HEADERS})
target_link_libraries(gcm_exe gcm cgalmesher ${CGAL_LIBRARIES} ${MPI_CXX_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES} ${VTK_LIBRARIES})

# tests
add_executable(gcm_tests ${SEQUENCE_TEST_SOURCES})
target_link_libraries(gcm_tests gcm cgalmesher ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CGAL_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES} ${MPI_CXX_LIBRARIES} ${VTK_LIBRARIES})

# tmp
#add_executable(tmp ${CMAKE_CURRENT_LIST_DIR}/src/test/sequence/TestLineWalkSearch3D.cpp ${CMAKE_CURRENT_LIST_DIR}/src/test/sequence/TestMain.cpp)
add_executable(tmp ${CMAKE_CURRENT_LIST_DIR}/src/test/sequence/TestLinal.cpp ${CMAKE_CURRENT_LIST_DIR}/src/test/sequence/TestMain.cpp)
target_link_libraries(tmp gcm cgalmesher ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CGAL_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES} ${MPI_CXX_LIBRARIES} ${VTK_LIBRARIES})

# tests on mpi
add_executable(gcm_mpi_tests ${TEST_MPI_FILE})
target_link_libraries(gcm_mpi_tests gcm cgalmesher ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CGAL_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES} ${MPI_CXX_LIBRARIES} ${VTK_LIBRARIES})

install(
    TARGETS gcm gcm_exe
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
)
