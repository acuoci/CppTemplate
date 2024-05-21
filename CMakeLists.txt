# Works with 3.14 and tested through 3.29
cmake_minimum_required(VERSION 3.15)

# Project name and a few useful settings. 
project(
  ProjectTemplate
  DESCRIPTION "Project Template C++ including static anlysis, code coverage, automatic tests, and documentation"
)

# Library name
set (TEMPLATE_LIB "Template")

# Enable C++ and C languages
enable_language(CXX C)

# Version number
set( VERSION_MAJOR 0 )
set( VERSION_MINOR 1 )
set( VERSION_PATCH 0 )
set( VERSION_COMPLETE "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

# Configure header file
configure_file( ${PROJECT_SOURCE_DIR}/cmake/config.h.in ${PROJECT_BINARY_DIR}/config.h )

# User-defined modules
set( CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

# Detect operating system
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    add_definitions(-DSYSTEM_LINUX)
endif()
if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    add_definitions(-DSYSTEM_DARWIN)
endif()
if(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
    add_definitions(-DSYSTEM_AIX)
endif()
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions(-DSYSTEM_WINDOWS)
endif()

# Detect SSE/AVX instructions
include(FindSSE)

# Detect host processor
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}")

# C++ compiler options
message(STATUS "The compiler is ${CMAKE_CXX_COMPILER_ID}")
include(CompilerOptions)

# Add Eigen3 C++ Libraries
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED)

message(STATUS "Boost libraries: ${Boost_LIBRARY_DIRS}")


# Docs only available if this is the main app
find_package(Doxygen)
if(Doxygen_FOUND)
    add_subdirectory(docs)
else()
    message(STATUS "Doxygen not found, not building docs")
endif()


# Google tests
option(ENABLE_GTEST "Build all tests according to the Google Test projects" OFF)
if(ENABLE_GTEST)
	include(CTest)
	add_subdirectory(test/gtest)
endif(ENABLE_GTEST)


# Code Coverage
option(ENABLE_CODECOV "Apply code coverage analysis through gcov (require Debug build type)" OFF)
if(ENABLE_CODECOV)
	add_subdirectory(test/codecov)
endif(ENABLE_CODECOV)


# CppCheck utility
option(ENABLE_CPPCHECK "Apply CppCheck for static analysis of code quality" OFF)
if(ENABLE_CPPCHECK)
	add_subdirectory(test/analysis)
endif(ENABLE_CPPCHECK)


# Include directories
include_directories("${PROJECT_SOURCE_DIR}/include")

# Add the source files
file(GLOB SOURCES src/*.cpp)


# Add the library
add_library(${TEMPLATE_LIB} ${SOURCES})

# Add the library links
target_link_libraries(${TEMPLATE_LIB} Boost::boost Eigen3::Eigen)

# Add targets
install(TARGETS ${TEMPLATE_LIB} DESTINATION lib)