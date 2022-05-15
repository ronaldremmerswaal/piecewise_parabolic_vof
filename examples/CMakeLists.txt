add_executable(demo
  reconstruction_demo.f90
  )

set_target_properties (demo PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(demo PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(demo ppic)