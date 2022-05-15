cmake_minimum_required(VERSION 3.13)

project(ppic C Fortran) 

# Unit testing
set(ENABLE_UNIT_TESTS_PPIC False CACHE BOOL "build and run ppic tests")
if (ENABLE_UNIT_TESTS_PPIC)
  add_subdirectory(tests)
endif ()

# Add directiories
add_subdirectory(examples)
add_subdirectory(external/Optimization_f90)
add_subdirectory(external/r3d)
add_subdirectory(src)
