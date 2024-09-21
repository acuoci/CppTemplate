# Determine architecture (32/64 bit)
set(X64 OFF)
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(X64 ON)
endif()
message(STATUS "Architecture X64: " ${X64})


# GNU or CLANG compilers
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")

    set(CMAKE_CXX_STANDARD 11)

    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++20")
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wall")
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wextra")				# enables some extra warning flags that are not enabled by -Wall
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wunused")				# all the -Wunused-xxx options combined
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wswitch-default")			# warn whenever a switch statement does not have a default case
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wuninitialized")			# warn if an automatic variable is used without first being initialized
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Woverloaded-virtual ")		# Warn when a function declaration hides virtual functions from a base class
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wwrite-strings")			# warn about strings
			
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wignored-qualifiers")		# warn if the return type of a function has a type qualifier such as const (-Wextra)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wmissing-braces")			# warn if an aggregate or union initializer is not fully bracketed (-Wall)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wreturn-type")			# warn whenever a function is defined with a return type that defaults to int (-Wall)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wswitch")				# warn about switch statement (-Wall)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wmissing-field-initializers")	# warn if a structure’s initializer has some fields missing (-Wextra)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")		# warn when a #pragma directive is encountered that is not understood (-Wall)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -Wsign-compare ")			# warn comparison between signed and unsigned values (-Wall)

    # Debug options		
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g --coverage -fprofile-arcs -ftest-coverage")

    # Release options
    if(AVX2_FOUND)
      set(CMAKE_CXX_FLAGS_RELEASE      "-O3 -march=core-avx2")
    else(AVX2_FOUND)
      set(CMAKE_CXX_FLAGS_RELEASE      "-O3 -march=native")
    endif(AVX2_FOUND)

endif()


# INTEL C++ Compiler
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")

    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -std=c++20")

    # Release options
    if(AVX2_FOUND)
      set(CMAKE_CXX_FLAGS_RELEASE      "-O3 -qopenmp -mcore-avx2")
    else(AVX2_FOUND)
      set(CMAKE_CXX_FLAGS_RELEASE      "-O3 -qopenmp -mavx")
    endif(AVX2_FOUND)

endif()


# Microsoft Visual Studio Compiler
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")

	set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} /MP")	# -> build with multiple processes
	set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} /W4")	# -> warning level 4        
	set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} /WX")	# -> treat warnings as errors
	set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} /wd4251") # -> disable warning: 'identifier': class 'type' needs dll-interface to be used by clients of class 'type2'
	set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} /wd4592") # -> disable warning: 'identifier': symbol will be dynamically initialized (implementation limitation)
	set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} /wd4201") # -> disable warning: nonstandard extension used: nameless struct/union (caused by GLM)
	set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} /wd4127") # -> disable warning: conditional expression is constant (caused by Qt)

    	set(CMAKE_CXX_FLAGS_RELEASE 	"${CMAKE_CXX_FLAGS_RELEASE} /Gw")  # -> whole program global optimization
	set(CMAKE_CXX_FLAGS_RELEASE 	"${CMAKE_CXX_FLAGS_RELEASE} /GS-")  # -> buffer security check: no 
	set(CMAKE_CXX_FLAGS_RELEASE 	"${CMAKE_CXX_FLAGS_RELEASE} /GL")  # -> whole program optimization: enable link-time code generation (disables Zi)
	set(CMAKE_CXX_FLAGS_RELEASE 	"${CMAKE_CXX_FLAGS_RELEASE} /GF")  # -> enable string pooling)

endif()


# Linker options 
set(DEFAULT_LINKER_OPTIONS)

# Use pthreads on mingw and linux
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" OR "${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
    set(DEFAULT_LINKER_OPTIONS
        -pthread
    )
endif()
