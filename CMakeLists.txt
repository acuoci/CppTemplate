# Minimum required version
cmake_minimum_required(VERSION 3.15)

# Project name and metadata
project(
  WildFireSimulator
  VERSION 0.1.0
  DESCRIPTION "Wildfire Simulator based on the Cellular Automata technique"
  LANGUAGES CXX C
)

# Library name
set(WFS_LIB "wildfiresimulator")

# Set C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # Ensures that no GNU extensions are used

# Configure the config.h file
configure_file(${PROJECT_SOURCE_DIR}/cmake/config.h.in ${PROJECT_BINARY_DIR}/config.h)

# User-defined modules
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Detect operating system and define SYSTEM macros
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(SYSTEM_DEFINITIONS SYSTEM_LINUX)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(SYSTEM_DEFINITIONS SYSTEM_DARWIN)
elseif(CMAKE_SYSTEM_NAME STREQUAL "AIX")
    set(SYSTEM_DEFINITIONS SYSTEM_AIX)
elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(SYSTEM_DEFINITIONS SYSTEM_WINDOWS)
endif()

# Apply system-specific definitions to the target
if(SYSTEM_DEFINITIONS)
    message(STATUS "System definition: ${SYSTEM_DEFINITIONS}")
endif()

# Detect SSE/AVX instructions
include(FindSSE)

# Detect host processor
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}")

# Compiler-specific options
include(CompilerOptions)

# Eigen3 C++ Libraries
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Boost filesystem library
find_package(Boost REQUIRED COMPONENTS filesystem program_options)

message(STATUS "Boost include directory: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost libraries: ${Boost_LIBRARIES}")

# Documentation with Doxygen
find_package(Doxygen)
if(Doxygen_FOUND)
    add_subdirectory(docs)
else()
    message(STATUS "Doxygen not found, skipping documentation")
endif()

# Enable Google Test if requested
option(ENABLE_GTEST "Enable building of Google Tests" OFF)
if(ENABLE_GTEST)
    include(CTest)
    add_subdirectory(test/gtest)
endif()

# Enable Code Coverage if requested
option(ENABLE_CODECOV "Enable code coverage analysis (requires Debug build)" OFF)
if(ENABLE_CODECOV)
    add_subdirectory(test/codecov)
endif()

# Enable CppCheck if requested
option(ENABLE_CPPCHECK "Enable CppCheck static analysis" OFF)
if(ENABLE_CPPCHECK)
    add_subdirectory(test/analysis)
endif()

# Add include directories (prefer target_include_directories for target scope)
set(INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/src")

# Collect source files
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

# Add the library target
add_library(${WFS_LIB} STATIC ${SOURCES})

# Set target properties
target_include_directories(${WFS_LIB}
    PRIVATE ${INCLUDE_DIRS}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
           $<INSTALL_INTERFACE:include>
)

# Link external libraries to the library target
target_link_libraries(${WFS_LIB} PUBLIC Boost::filesystem Boost::program_options Eigen3::Eigen)

# Apply system-specific definitions to the target
if(SYSTEM_DEFINITIONS)
    target_compile_definitions(${WFS_LIB} PRIVATE ${SYSTEM_DEFINITIONS})
endif()

# Install the library and header files
install(TARGETS ${WFS_LIB}
    EXPORT ${WFS_LIB}Targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/
    DESTINATION include
)

# Add standalone executable subdirectory
add_subdirectory(standalone)
