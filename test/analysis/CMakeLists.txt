# Find the CppCheck binary
find_package(CppCheck)

# Static analysis. Should be set before adding subprojects
set(ALL_ANALYSIS_TARGETS)

# Add sub directories
#add_subdirectory(linear)

# Add the "make analysis" target
if( CPPCHECK_FOUND )
	add_custom_target(analysis)
	ADD_DEPENDENCIES(analysis ${ALL_ANALYSIS_TARGETS})
	set_target_properties(analysis PROPERTIES EXCLUDE_FROM_ALL TRUE)
	message(STATUS "Analysis analysis targets are ${ALL_ANALYSIS_TARGETS}")
endif()