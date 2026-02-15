cmake_minimum_required(VERSION 3.14)
project(PROPACK VERSION 2.2 LANGUAGES Fortran)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(PROPACK_BUILD_EXAMPLES "Build example programs" ON)
option(PROPACK_BUILD_TESTING  "Build CTest targets for examples" ON)
option(PROPACK_USE_OPENMP     "Enable OpenMP parallelization" OFF)

option(PROPACK_BUILD_SINGLE    "Build single-precision (s) library"         ON)
option(PROPACK_BUILD_DOUBLE    "Build double-precision (d) library"         ON)
option(PROPACK_BUILD_COMPLEX8  "Build complex single-precision (c) library" ON)
option(PROPACK_BUILD_COMPLEX16 "Build complex double-precision (z) library" ON)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
find_package(LAPACK REQUIRED)   # also finds BLAS via transitive dependency

if(PROPACK_USE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

# ---------------------------------------------------------------------------
# Default Fortran flags (only when not already set by the user)
# ---------------------------------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(_base_flags "-fno-second-underscore -Wall")
  if(NOT CMAKE_Fortran_FLAGS)
    set(CMAKE_Fortran_FLAGS "${_base_flags}" CACHE STRING "" FORCE)
  endif()
  if(NOT CMAKE_Fortran_FLAGS_RELEASE)
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops" CACHE STRING "" FORCE)
  endif()
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
  add_library(${_lib_name} STATIC ${_own_sources} ${_real_sources})

  target_include_directories(${_lib_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/${TYPE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/common
  )
  target_link_libraries(${_lib_name} PUBLIC propack_common LAPACK::LAPACK)
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

      # The Fortran programs read from stdin.  CTest does not natively
      # support input redirection, so we use a small shell wrapper.
      # This works on all Unix-like systems and on Windows with Git Bash
      # or MSYS2.  For MSVC/NMake on Windows you would need a .bat
      # wrapper instead, but Fortran + LAPACK users are overwhelmingly
      # on Unix or using MinGW/MSYS2.
      add_test(
        NAME ${PREFIX}propack_example_run
        COMMAND sh -c "$<TARGET_FILE:${PREFIX}propack_example> < example.in"
        WORKING_DIRECTORY ${_runtime_dir}
      )

      add_test(
        NAME ${PREFIX}propack_example_irl_run
        COMMAND sh -c "$<TARGET_FILE:${PREFIX}propack_example_irl> < example_irl.in"
        WORKING_DIRECTORY ${_runtime_dir}
      )

      add_test(
        NAME ${PREFIX}propack_verify
        COMMAND sh -c "$<TARGET_FILE:${PREFIX}propack_compare> < compare.in"
        WORKING_DIRECTORY ${_runtime_dir}
      )
      set_tests_properties(${PREFIX}propack_verify PROPERTIES
        DEPENDS "${PREFIX}propack_example_run"
      )

      add_test(
        NAME ${PREFIX}propack_verify_irl
        COMMAND sh -c "$<TARGET_FILE:${PREFIX}propack_compare> < compare_irl.in"
        WORKING_DIRECTORY ${_runtime_dir}
      )
      set_tests_properties(${PREFIX}propack_verify_irl PROPERTIES
        DEPENDS "${PREFIX}propack_example_irl_run"
      )
    endif()
  endif()
endfunction()

# ---------------------------------------------------------------------------
# Common library (shared integer utilities + statistics printing)
# ---------------------------------------------------------------------------
add_library(propack_common STATIC
  common/iblasext.F
  common/printstat.F
)
target_include_directories(propack_common PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/common
)
target_link_libraries(propack_common PUBLIC LAPACK::LAPACK)

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
