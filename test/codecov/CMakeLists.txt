# Recognize gcov package
include(CodeCoverage)

# Add instructions to gcc compiler
APPEND_COVERAGE_COMPILER_FLAGS()

# Exclude directories and files from coverage report
set(COVERAGE_EXCLUDES 'Core/*' 'plugins/*' 'usr/*' 'gtest/internal/*')

# Setup for basic code coverage
SETUP_TARGET_FOR_COVERAGE( NAME test_codecov EXECUTABLE test_linear ${PROJECT_SOURCE_DIR}/data/)

# Setup for Cobertura output
SETUP_TARGET_FOR_COVERAGE_COBERTURA(NAME test_codecov_cobertura EXECUTABLE test_linear ${PROJECT_SOURCE_DIR}/data/ --gtest_output=xml:coverage.junit.xml)