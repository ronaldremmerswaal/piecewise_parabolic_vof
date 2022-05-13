cmake_minimum_required(VERSION 3.13)

project(ppic C Fortran) 

# Unit testing
if (ENABLE_UNIT_TESTS)
  add_subdirectory(tests)
endif ()
add_library(ppic)

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(ppic::ppic ALIAS ppic)
set(ppic_LIBRARIES ppic::ppic CACHE STRING "ppic library target")

# find_package(Git QUIET)
# if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# # Update submodules as needed
#     option(GIT_SUBMODULE "Check submodules during build" ON)
#     if(GIT_SUBMODULE)
#         message(STATUS "Submodule update")
#         execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
#                         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#                         RESULT_VARIABLE GIT_SUBMOD_RESULT)
#         if(NOT GIT_SUBMOD_RESULT EQUAL "0")
#             message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
#         endif()
#     endif()
# endif()

# Add directiories
add_subdirectory(external/Optimization_f90)
add_subdirectory(external/r3d)
add_subdirectory(src)

# Wonton targets
install(TARGETS ppic
  EXPORT ppic_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include
  INCLUDES DESTINATION include
  )

#-----------------------------------------------------------------------------
# Prepare output for configuration files to be used by projects importing Wonton
#-----------------------------------------------------------------------------

# export targets

install(EXPORT ppic_LIBRARIES
  FILE ppicTargets.cmake
  NAMESPACE ppic::
  EXPORT_LINK_INTERFACE_LIBRARIES
  DESTINATION lib/cmake/ppic)
