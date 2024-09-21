# Add the executable target for standalone
add_executable(${PROJECT_NAME} main.cpp)  # Replace main.cpp with your actual source files
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17 OUTPUT_NAME "WildFireSimulator")

# Link the wildfiresimulator library (WFS_LIB) to the standalone executable
target_link_libraries(${PROJECT_NAME} PRIVATE wildfiresimulator)

# Add include directories for the standalone executable (in case it needs to access public headers from Template)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)
