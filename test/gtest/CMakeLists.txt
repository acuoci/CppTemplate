# Set path to installed Google Test library
set(GTEST_ROOT REQUIRED CACHE PATH "Directory to Google Test library")

# Prepare Google Tests
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
enable_testing()

# Search for Google Test library
find_package(GTest REQUIRED)
include_directories( ${GTEST_INCLUDE_DIRS} ${GMOCK_INCLUDE_DIRS} )

# Unit Tests
add_executable (test_template test_template/test.cpp)

# Standard linking to gtest stuff.
target_link_libraries(test_template GTest::gtest GTest::gtest_main)

# Extra linking for the project.
target_link_libraries(test_template ${TEMPLATE_LIB} )

# Add test
add_test(NAME AllTests COMMAND test_template ${PROJECT_SOURCE_DIR})