cmake_minimum_required(VERSION 2.8.4)
project(gcm)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")

set(WARNING_FLAGS "-Wall -Wextra -Wsign-conversion -Wconversion -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -Woverloaded-virtual")
set(SAFE_MATH_OPERATIONS "-frounding-math -fsignaling-nans") # it's also necessary for CGAL working
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${WARNING_FLAGS} ${SAFE_MATH_OPERATIONS}")

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
find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS})
# VTK library
find_package(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
if(VTK_MAJOR_VERSION EQUAL 5)
    set(CONFIG_VTK_5 1)
elseif (VTK_MAJOR_VERSION EQUAL 6)
    set(CONFIG_VTK_6 1)
else()
    message(FATAL_ERROR "Unsupported VTK version found: ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
endif()

# cmake generated config
set(CMAKE_CONFIG_SHARE_GCM "${CMAKE_INSTAL_PREFIX}/share/gcm")

set(CMAKE_CONFIG_ENABLE_LOGGING       1)
set(CMAKE_CONFIG_ENABLE_LOGGING_FATAL 1)
set(CMAKE_CONFIG_ENABLE_LOGGING_ERROR 1)
set(CMAKE_CONFIG_ENABLE_LOGGING_WARN  1)
set(CMAKE_CONFIG_ENABLE_LOGGING_INFO  0)
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
configure_file("${CMAKE_CURRENT_LIST_DIR}/cmake/config.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/generated_sources/lib/config.hpp")

# source files
include_directories("${CMAKE_CURRENT_LIST_DIR}/src/" "${CMAKE_CURRENT_BINARY_DIR}/generated_sources/")
FILE(GLOB_RECURSE HEADERS "${CMAKE_CURRENT_LIST_DIR}/src/*.hpp")
file(GLOB_RECURSE GCM_EXE_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/launcher/gcm_exe/*.cpp")
file(GLOB_RECURSE SEISMO_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/launcher/seismo/*.cpp")
file(GLOB_RECURSE INVERSE_PROBLEM_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/launcher/inverse_problem/*.cpp")
file(GLOB_RECURSE LIB_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/lib/*.cpp")
file(GLOB_RECURSE SEQUENCE_TEST_SOURCES "${CMAKE_CURRENT_LIST_DIR}/src/test/sequence/*.cpp")
set(TEST_MPI_FILE src/test/TestMPI.cpp)

# libgcm
add_library(gcm SHARED ${LIB_SOURCES})

# gcm_exe
add_executable(gcm_exe ${GCM_EXE_SOURCES} ${HEADERS})
target_link_libraries(gcm_exe gcm ${CGAL_LIBRARIES} ${MPI_CXX_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES} ${VTK_LIBRARIES})

# gcm_seismo
add_executable(gcm_seismo ${SEISMO_SOURCES})
target_link_libraries(gcm_seismo gcm ${MPI_CXX_LIBRARIES} ${CGAL_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES} ${VTK_LIBRARIES})

# gcm_inverse_problem
add_executable(gcm_inverse_problem ${INVERSE_PROBLEM_SOURCES})
target_link_libraries(gcm_inverse_problem gcm ${MPI_CXX_LIBRARIES} ${CGAL_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES} ${VTK_LIBRARIES})

# tests
add_executable(gcm_tests ${SEQUENCE_TEST_SOURCES})
target_link_libraries(gcm_tests gcm ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CGAL_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES} ${MPI_CXX_LIBRARIES} ${VTK_LIBRARIES})

## tmp
#add_executable(tmp ${CMAKE_CURRENT_LIST_DIR}/src/test/sequence/TestEngine.cpp ${CMAKE_CURRENT_LIST_DIR}/src/test/sequence/TestMain.cpp)
#target_link_libraries(tmp gcm ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CGAL_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES} ${MPI_CXX_LIBRARIES} ${VTK_LIBRARIES})

# tests on mpi
add_executable(gcm_mpi_tests ${TEST_MPI_FILE})
target_link_libraries(gcm_mpi_tests gcm ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CGAL_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES} ${MPI_CXX_LIBRARIES} ${VTK_LIBRARIES})

install(
    TARGETS gcm gcm_exe
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
)
