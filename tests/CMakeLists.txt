# (C) Copyright 2022- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

foreach( prec sp dp )
  if( HAVE_${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_aqua_era5_O48 CONFIG aqua_era5_O48.yml PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_aqua_era5_O48 CONFIG aqua_era5_O48.yml MPI 4 OMP 1 PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_aqua_era5_O48 CONFIG aqua_era5_O48.yml MPI 1 OMP 4 PREC ${prec} )

    ecwam_add_test( ${PROJECT_NAME}_test_aqua_oper_an_fc_O48 CONFIG aqua_oper_an_fc_O48.yml PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_aqua_oper_an_fc_O48 CONFIG aqua_oper_an_fc_O48.yml MPI 4 OMP 1 PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_aqua_oper_an_fc_O48 CONFIG aqua_oper_an_fc_O48.yml MPI 1 OMP 4 PREC ${prec} )

    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_era5_O48 CONFIG etopo1_era5_O48.yml PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_era5_O48 CONFIG etopo1_era5_O48.yml MPI 4 OMP 1 PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_era5_O48 CONFIG etopo1_era5_O48.yml MPI 1 OMP 4 PREC ${prec} )

    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48 CONFIG etopo1_oper_an_fc_O48.yml PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48 CONFIG etopo1_oper_an_fc_O48.yml MPI 4 OMP 1 PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48 CONFIG etopo1_oper_an_fc_O48.yml MPI 1 OMP 4 PREC ${prec} )

    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_grib_restart_O96 CONFIG etopo1_oper_an_fc_grib_restart_O96.yml PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_grib_restart_O96 CONFIG etopo1_oper_an_fc_grib_restart_O96.yml MPI 4 OMP 1 PREC ${prec} )
    ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_grib_restart_O96 CONFIG etopo1_oper_an_fc_grib_restart_O96.yml MPI 1 OMP 4 PREC ${prec} )

    # This configuration is retained for the test-base purely for coverage reasons, and is not used for
    # research or operational forecasts. As such, Loki transformations remove the relevant code-paths as
    # this can yield significant performance benefits, especially on GPU. 
    # (for further info see the discussion here: https://github.com/ecmwf-ifs/ecwam/pull/92)
    if( NOT HAVE_LOKI )
      ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48_iphys_0 CONFIG etopo1_oper_an_fc_O48_iphys_0.yml PREC ${prec} )
      ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48_iphys_0 CONFIG etopo1_oper_an_fc_O48_iphys_0.yml MPI 4 OMP 1 PREC ${prec} )
      ecwam_add_test( ${PROJECT_NAME}_test_etopo1_oper_an_fc_O48_iphys_0 CONFIG etopo1_oper_an_fc_O48_iphys_0.yml MPI 1 OMP 4 PREC ${prec} )
    endif()

  endif()
endforeach()

# Check CPU fusion independently of per-source vector-math overrides. Both
# worker objects use identical precise/no-FMA settings, including on a build
# whose production SNONLIN was compiled with different exception assumptions.
if( HAVE_TESTS AND HAVE_sp AND NOT HAVE_ACC AND NOT HAVE_LOKI AND
    CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND CMAKE_BUILD_TYPE MATCHES "^[Bb][Ii][Tt]$" )
  ecbuild_add_test(
    TARGET ecwam_test_snonlin_fusion_sp
    SOURCES snonlin_fusion_check.F90 snonlin_fused_reference.F90 snonlin_split_reference.F90
    LIBS ecwam_sp
    DEFINITIONS WAM_HAVE_SINGLE_PRECISION
    FFLAGS "-fp-model precise -fp-speculation=safe -no-fma -fpe0 -check bounds"
  )
  ecwam_target_compile_definitions_FILENAME( ecwam_test_snonlin_fusion_sp )
endif()

if( HAVE_TESTS AND HAVE_sp AND HAVE_MPI AND NOT HAVE_ACC )
  ecbuild_add_test(
    TARGET ecwam_test_rank_timer_logs_sp
    COMMAND ${ECWAM_PYTHON_INTERP}
    ARGS ${CMAKE_CURRENT_SOURCE_DIR}/rank_timer_log_check.py
         ${CMAKE_CURRENT_BINARY_DIR}/ecwam_test_etopo1_oper_an_fc_O48_sp_mpi4_omp1 4
    TEST_REQUIRES ecwam_test_etopo1_oper_an_fc_O48_sp_mpi4_omp1
  )
endif()


if( HAVE_TESTS AND HAVE_sp AND HAVE_PROPAG_FACTORED_PAIRS )
  set( pair_check_flags "" )
  if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
    set( pair_check_flags "-fpe0 -check bounds" )
  endif()
  ecbuild_add_test(
    TARGET ecwam_test_propag_pairs_sp
    SOURCES propags2_tiling_check.F90
    LIBS ecwam_sp
    DEFINITIONS WAM_HAVE_SINGLE_PRECISION
    FFLAGS "${pair_check_flags}"
    OMP 4
  )
  ecwam_target_compile_definitions_FILENAME( ecwam_test_propag_pairs_sp )
endif()
