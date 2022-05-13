find_package(PFUNIT REQUIRED)
enable_testing()

# tests
add_pfunit_ctest (ppic_tests
  TEST_SOURCES reconstruction.pf intersection.pf
  LINK_LIBRARIES ppic
  )
