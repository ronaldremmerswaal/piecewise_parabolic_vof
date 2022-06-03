add_library(ppic
  m_reconstruction.f90
  m_reconstruction_util.f90
  m_polygon.f90
  m_r2d.f90
  m_r2d_parabolic.f90
  r2d_parabolic.h
  r2d_parabolic.c
  m_common.f90
  )

set_target_properties (ppic PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(ppic PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(ppic optim r3d)