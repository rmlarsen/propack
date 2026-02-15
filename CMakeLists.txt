cmake_minimum_required(VERSION 3.14)
project(PROPACK VERSION 2.2 LANGUAGES Fortran)

# C++ API wrappers need a C++ compiler for header validation / tests.
# The Fortran wrappers themselves are pure Fortran (iso_c_binding).

include(GNUInstallDirs)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(PROPACK_BUILD_EXAMPLES    "Build example programs" ON)
option(PROPACK_BUILD_BENCHMARKS  "Build benchmark programs" OFF)
if(PROPACK_BUILD_BENCHMARKS AND CMAKE_BUILD_TYPE AND
   NOT CMAKE_BUILD_TYPE STREQUAL "Release" AND
   NOT CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  message(WARNING
    "Benchmarks are most meaningful in Release mode (current: ${CMAKE_BUILD_TYPE}).")
endif()
option(PROPACK_BUILD_TESTING     "Build CTest targets for examples" ON)
option(PROPACK_USE_OPENMP        "Enable OpenMP parallelization" OFF)

option(PROPACK_BUILD_SINGLE    "Build single-precision (s) library"         ON)
option(PROPACK_BUILD_DOUBLE    "Build double-precision (d) library"         ON)
option(PROPACK_BUILD_COMPLEX8  "Build complex single-precision (c) library" ON)
option(PROPACK_BUILD_COMPLEX16 "Build complex double-precision (z) library" ON)
option(PROPACK_BUILD_CXX_API  "Build C/C++ API wrappers (iso_c_binding)"   ON)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
# Prefer OpenBLAS when available; fall back to any BLAS/LAPACK otherwise.
# Users can override with -DBLA_VENDOR=<vendor> on the cmake command line.
if(NOT BLA_VENDOR)
  set(BLA_VENDOR OpenBLAS)
  find_package(LAPACK QUIET)
  if(NOT LAPACK_FOUND)
    unset(BLA_VENDOR)
    find_package(LAPACK REQUIRED)
  endif()
else()
  find_package(LAPACK REQUIRED)
endif()

if(PROPACK_USE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

# ---------------------------------------------------------------------------
# Default Fortran flags (only when not already set by the user)
# ---------------------------------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(_base_flags "-fno-second-underscore -Wall")
  set(_release_flags "-O3 -funroll-all-loops")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  set(_base_flags "-warn general")
  set(_release_flags "-O3")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(_base_flags "-warn general")
  set(_release_flags "-O3")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang|Flang")
  set(_base_flags "-Wall")
  set(_release_flags "-O3")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
  set(_base_flags "-Minform=warn")
  set(_release_flags "-O3 -fast")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  set(_base_flags "-w=all")
  set(_release_flags "-O3")
else()
  set(_base_flags "")
  set(_release_flags "-O3")
endif()

if(_base_flags AND NOT CMAKE_Fortran_FLAGS)
  set(CMAKE_Fortran_FLAGS "${_base_flags}" CACHE STRING "" FORCE)
endif()
if(NOT CMAKE_Fortran_FLAGS_RELEASE)
  set(CMAKE_Fortran_FLAGS_RELEASE "${_release_flags}" CACHE STRING "" FORCE)
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# CTest — must come before add_test() calls inside the function
# ---------------------------------------------------------------------------
if(PROPACK_BUILD_TESTING)
  enable_testing()
endif()

# ---------------------------------------------------------------------------
# Helper: build a PROPACK precision variant
#
#   _add_propack_variant(<prefix> <type_dir> <real_dir> <real_prefix>)
#
#   e.g.  _add_propack_variant(d double double d)
#         _add_propack_variant(z complex16 double d)
# ---------------------------------------------------------------------------
function(_add_propack_variant PREFIX TYPE_DIR REAL_DIR REAL_PREFIX)

  # --- sources that live in every precision directory ---
  set(_own_sources
    ${TYPE_DIR}/${PREFIX}lanbpro.F
    ${TYPE_DIR}/${PREFIX}reorth.F
    ${TYPE_DIR}/${PREFIX}getu0.F
    ${TYPE_DIR}/${PREFIX}safescal.F
    ${TYPE_DIR}/${PREFIX}blasext.F
    ${TYPE_DIR}/${PREFIX}lansvd.F
    ${TYPE_DIR}/${PREFIX}gemm_ovwr.F
    ${TYPE_DIR}/${PREFIX}lansvd_irl.F
    ${TYPE_DIR}/${PREFIX}ritzvec.F
  )

  if(PROPACK_BUILD_CXX_API)
    list(APPEND _own_sources ${TYPE_DIR}/${PREFIX}propack_c.f90)
  endif()

  # --- real-valued helpers needed by every variant ---
  # For real variants (double, single) these partially overlap with
  # _own_sources; only add ortho_rec and bsvd which are not in the list
  # above.  For complex variants the real helpers come from a different
  # directory so all four files are unique.
  if(TYPE_DIR STREQUAL REAL_DIR)
    set(_real_sources
      ${REAL_DIR}/${REAL_PREFIX}ortho_rec.F
      ${REAL_DIR}/${REAL_PREFIX}bsvd.F
    )
  else()
    set(_real_sources
      ${REAL_DIR}/${REAL_PREFIX}blasext.F
      ${REAL_DIR}/${REAL_PREFIX}ortho_rec.F
      ${REAL_DIR}/${REAL_PREFIX}bsvd.F
      ${REAL_DIR}/${REAL_PREFIX}gemm_ovwr.F
    )
  endif()

  set(_lib_name ${PREFIX}propack)
  add_library(${_lib_name} ${_own_sources} ${_real_sources})
  add_library(PROPACK::${_lib_name} ALIAS ${_lib_name})

  target_include_directories(${_lib_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/${TYPE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/common
  )
  if(PROPACK_BUILD_CXX_API)
    target_include_directories(${_lib_name} PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
  endif()
  target_link_libraries(${_lib_name} PUBLIC propack_common LAPACK::LAPACK)

  install(TARGETS ${_lib_name} EXPORT PROPACKTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
  if(PROPACK_USE_OPENMP)
    target_link_libraries(${_lib_name} PUBLIC OpenMP::OpenMP_Fortran)
  endif()

  # --- example programs ---
  if(PROPACK_BUILD_EXAMPLES)
    set(_ex_dir ${TYPE_DIR}/Examples)

    foreach(_prog example example_irl compare)
      set(_tgt ${PREFIX}propack_${_prog})
      add_executable(${_tgt}
        ${_ex_dir}/${_prog}.F
        ${_ex_dir}/matvec.F
      )
      target_include_directories(${_tgt} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/${_ex_dir}
        ${CMAKE_CURRENT_SOURCE_DIR}/${TYPE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/common
      )
      target_link_libraries(${_tgt} PRIVATE ${_lib_name})
    endforeach()

    # ------- CTest targets -------
    if(PROPACK_BUILD_TESTING)
      # Copy data files into a per-variant directory under the build tree
      # so the executables can find them at test time.
      set(_runtime_dir ${CMAKE_CURRENT_BINARY_DIR}/testing/${TYPE_DIR})
      file(MAKE_DIRECTORY ${_runtime_dir})
      file(MAKE_DIRECTORY ${_runtime_dir}/Output)

      file(GLOB _data_files
        ${CMAKE_CURRENT_SOURCE_DIR}/${_ex_dir}/*.in
        ${CMAKE_CURRENT_SOURCE_DIR}/${_ex_dir}/*.rra
        ${CMAKE_CURRENT_SOURCE_DIR}/${_ex_dir}/*.cua
        ${CMAKE_CURRENT_SOURCE_DIR}/${_ex_dir}/*.coord
        ${CMAKE_CURRENT_SOURCE_DIR}/${_ex_dir}/*.diag
      )
      file(GLOB _ref_files
        ${CMAKE_CURRENT_SOURCE_DIR}/${_ex_dir}/Output/*.ascii
      )

      foreach(_f ${_data_files})
        get_filename_component(_fname ${_f} NAME)
        configure_file(${_f} ${_runtime_dir}/${_fname} COPYONLY)
      endforeach()
      foreach(_f ${_ref_files})
        get_filename_component(_fname ${_f} NAME)
        configure_file(${_f} ${_runtime_dir}/Output/${_fname} COPYONLY)
      endforeach()

      # The Fortran programs read from stdin.  Use a CMake script wrapper
      # for portable input redirection (works on Unix and native Windows).
      set(_stdin_script ${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_with_stdin.cmake)

      add_test(
        NAME ${PREFIX}propack_example_run
        COMMAND ${CMAKE_COMMAND}
          -DEXE=$<TARGET_FILE:${PREFIX}propack_example>
          -DINPUT=${_runtime_dir}/example.in
          -P ${_stdin_script}
        WORKING_DIRECTORY ${_runtime_dir}
      )
      set_tests_properties(${PREFIX}propack_example_run PROPERTIES
        FIXTURES_SETUP ${PREFIX}propack_example_fixture
      )

      add_test(
        NAME ${PREFIX}propack_example_irl_run
        COMMAND ${CMAKE_COMMAND}
          -DEXE=$<TARGET_FILE:${PREFIX}propack_example_irl>
          -DINPUT=${_runtime_dir}/example_irl.in
          -P ${_stdin_script}
        WORKING_DIRECTORY ${_runtime_dir}
      )
      set_tests_properties(${PREFIX}propack_example_irl_run PROPERTIES
        FIXTURES_SETUP ${PREFIX}propack_example_irl_fixture
      )

      add_test(
        NAME ${PREFIX}propack_verify
        COMMAND ${CMAKE_COMMAND}
          -DEXE=$<TARGET_FILE:${PREFIX}propack_compare>
          -DINPUT=${_runtime_dir}/compare.in
          -P ${_stdin_script}
        WORKING_DIRECTORY ${_runtime_dir}
      )
      set_tests_properties(${PREFIX}propack_verify PROPERTIES
        FIXTURES_REQUIRED ${PREFIX}propack_example_fixture
      )

      add_test(
        NAME ${PREFIX}propack_verify_irl
        COMMAND ${CMAKE_COMMAND}
          -DEXE=$<TARGET_FILE:${PREFIX}propack_compare>
          -DINPUT=${_runtime_dir}/compare_irl.in
          -P ${_stdin_script}
        WORKING_DIRECTORY ${_runtime_dir}
      )
      set_tests_properties(${PREFIX}propack_verify_irl PROPERTIES
        FIXTURES_REQUIRED ${PREFIX}propack_example_irl_fixture
      )
    endif()

    # --- benchmark program ---
    if(PROPACK_BUILD_BENCHMARKS)
      set(_bench_tgt ${PREFIX}propack_benchmark)
      add_executable(${_bench_tgt}
        ${_ex_dir}/benchmark.F
        ${_ex_dir}/matvec.F
      )
      target_include_directories(${_bench_tgt} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/${_ex_dir}
        ${CMAKE_CURRENT_SOURCE_DIR}/${TYPE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/common
      )
      target_link_libraries(${_bench_tgt} PRIVATE ${_lib_name})

      if(PROPACK_BUILD_TESTING)
        add_test(
          NAME ${PREFIX}propack_benchmark_run
          COMMAND ${_bench_tgt}
          WORKING_DIRECTORY ${_runtime_dir}
        )
        set_tests_properties(${PREFIX}propack_benchmark_run PROPERTIES
          LABELS "benchmark"
          TIMEOUT 600
        )
      endif()

      # --- generated-matrix scaling benchmark ---
      set(_benchgen_tgt ${PREFIX}propack_benchmark_gen)
      add_executable(${_benchgen_tgt} ${_ex_dir}/benchmark_gen.F)
      target_include_directories(${_benchgen_tgt} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/${TYPE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/common
      )
      target_link_libraries(${_benchgen_tgt} PRIVATE ${_lib_name})
      if(PROPACK_BUILD_TESTING)
        add_test(NAME ${PREFIX}propack_benchmark_gen_run
          COMMAND ${_benchgen_tgt})
        set_tests_properties(${PREFIX}propack_benchmark_gen_run PROPERTIES
          LABELS "benchmark" TIMEOUT 1800)
      endif()
    endif()
  endif()
endfunction()

# ---------------------------------------------------------------------------
# Common library (shared integer utilities + statistics printing)
# ---------------------------------------------------------------------------
add_library(propack_common
  common/iblasext.F
  common/printstat.F
)
target_include_directories(propack_common PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/common
)
target_link_libraries(propack_common PUBLIC LAPACK::LAPACK)

install(TARGETS propack_common EXPORT PROPACKTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# ---------------------------------------------------------------------------
# Build requested precision variants
# ---------------------------------------------------------------------------
#                       prefix  type_dir    real_dir  real_prefix
if(PROPACK_BUILD_DOUBLE)
  _add_propack_variant(d       double      double    d)
endif()
if(PROPACK_BUILD_SINGLE)
  _add_propack_variant(s       single      single    s)
endif()
if(PROPACK_BUILD_COMPLEX16)
  _add_propack_variant(z       complex16   double    d)
endif()
if(PROPACK_BUILD_COMPLEX8)
  _add_propack_variant(c       complex8    single    s)
endif()

# ---------------------------------------------------------------------------
# C++ API unit tests (GTest)
# ---------------------------------------------------------------------------
if(PROPACK_BUILD_CXX_API AND PROPACK_BUILD_TESTING)
  enable_language(CXX)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  find_package(GTest QUIET)
  if(GTest_FOUND)
    add_executable(propack_cxx_test tests/test_propack_cxx.cpp)
    target_link_libraries(propack_cxx_test PRIVATE
      GTest::gtest_main
      dpropack spropack zpropack cpropack
    )
    include(GoogleTest)
    gtest_discover_tests(propack_cxx_test)
  else()
    message(STATUS "GTest not found — skipping C++ API tests")
  endif()
endif()

# ---------------------------------------------------------------------------
# C++ API benchmarks (Google Benchmark)
# ---------------------------------------------------------------------------
if(PROPACK_BUILD_CXX_API AND PROPACK_BUILD_BENCHMARKS)
  enable_language(CXX)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  find_package(benchmark QUIET)
  if(benchmark_FOUND)
    add_executable(propack_cxx_benchmark tests/bench_propack_cxx.cpp)
    target_link_libraries(propack_cxx_benchmark PRIVATE
      benchmark::benchmark benchmark::benchmark_main
      dpropack spropack zpropack cpropack)
    if(PROPACK_BUILD_TESTING)
      add_test(NAME propack_cxx_benchmark_run
        COMMAND propack_cxx_benchmark --benchmark_min_time=0.1s)
      set_tests_properties(propack_cxx_benchmark_run PROPERTIES
        LABELS "benchmark" TIMEOUT 600)
    endif()
  else()
    message(STATUS "Google Benchmark not found — skipping C++ API benchmarks")
  endif()
endif()

# ---------------------------------------------------------------------------
# Install C++ header
# ---------------------------------------------------------------------------
if(PROPACK_BUILD_CXX_API)
  install(FILES include/propack.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

# ---------------------------------------------------------------------------
# CMake package config for find_package(PROPACK)
# ---------------------------------------------------------------------------
include(CMakePackageConfigHelpers)

install(EXPORT PROPACKTargets
  FILE PROPACKTargets.cmake
  NAMESPACE PROPACK::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PROPACK
)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/PROPACKConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/PROPACKConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PROPACK
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/PROPACKConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/PROPACKConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/PROPACKConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PROPACK
)
