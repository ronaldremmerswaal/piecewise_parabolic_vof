project(tests Fortran) 

if(DEFINED PFUNIT_INSTALL)
    message(STATUS "Manual setup of variable PFUNIT_INSTALL: ${PFUNIT_INSTALL}")
    set(PFUNIT_DIR ${PFUNIT_INSTALL})
else()
    include(ExternalProject)

    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external/pfunit
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        )
    ExternalProject_Add(pfunit
        DOWNLOAD_COMMAND git submodule update
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/pfunit
        BINARY_DIR ${PROJECT_BINARY_DIR}/external/pfunit-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/pfunit-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/pfunit-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${ExternalProjectCMakeArgs}
        )
    include_directories(${PROJECT_BINARY_DIR}/external/pfunit/mod)
    set(PFUNIT_DIR ${PROJECT_BINARY_DIR}/external/pfunit)
endif()

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/generated)
file(WRITE ${PROJECT_BINARY_DIR}/generated/testSuites.inc "")

include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_BINARY_DIR}/generated
    ${PFUNIT_DIR}/mod
    )

set(_test_sources)
foreach(_test
    temp
    )
    if(DEFINED PFUNIT_INSTALL)
        set(test_dependency ${PROJECT_SOURCE_DIR}/${_test}.pf)
    else()
        set(test_dependency pfunit ${PROJECT_SOURCE_DIR}/${_test}.pf)
    endif()
    add_custom_command(
        OUTPUT ${PROJECT_BINARY_DIR}/generated/${_test}.pf
        COMMAND python3 ${PFUNIT_DIR}/bin/funit/pFUnitParser.py ${PROJECT_SOURCE_DIR}/${_test}.pf ${PROJECT_BINARY_DIR}/generated/${_test}.pf
        DEPENDS ${test_dependency}
        )
    set(_test_sources ${_test_sources} ${PROJECT_BINARY_DIR}/generated/${_test}.pf)
    file(APPEND ${PROJECT_BINARY_DIR}/generated/testSuites.inc "ADD_TEST_SUITE(${_test}_suite)\n")
endforeach()

set_source_files_properties(${PFUNIT_DIR}/include/driver.F90 PROPERTIES GENERATED 1)

add_executable(
    pftest_alltests
    ${PFUNIT_DIR}/include/driver.F90
    ${_test_sources}
    )

target_link_libraries(
    pftest_alltests
    ${PFUNIT_DIR}/lib/libpfunit.a
    ppic
    )

enable_testing()
add_test(pftest_alltests ${PROJECT_BINARY_DIR}/pftest_alltests)