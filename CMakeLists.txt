cmake_minimum_required(VERSION 3.13)

project(ppic Fortran) 

# Unit testing
set(ENABLE_UNIT_TESTS_PPIC False CACHE BOOL "build and run ppic tests")
if (ENABLE_UNIT_TESTS_PPIC)
  add_subdirectory(tests)
endif ()

# Get dependencies
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  # Update submodules as needed
  option(GIT_SUBMODULE "Check submodules during build" ON)
  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()
endif()

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message("-- Build type is set to ${CMAKE_BUILD_TYPE}")
endif()

message("-- Using ${CMAKE_Fortran_COMPILER_ID} Fortran compiler")

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback -check bounds -debug all")
  set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
  set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native")
endif()


# Add directiories
add_subdirectory(examples)
add_subdirectory(external/Optimization_f90)
add_subdirectory(src)