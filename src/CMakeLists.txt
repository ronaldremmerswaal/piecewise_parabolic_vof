add_library(ppic
  m_r2d_parabolic.f90
  r2d_parabolic.h
  r2d_parabolic.c
  )

set_target_properties (ppic PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(ppic PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
