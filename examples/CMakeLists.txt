add_executable(demo
  recon_demo.f90
  )
add_executable(benchmark
  recon_benchmark.f90
  )

set_target_properties (demo PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(demo PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(demo ppic)

set_target_properties (benchmark PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(benchmark PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(benchmark ppic)