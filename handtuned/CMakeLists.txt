# This project owns the handtuned build. It consumes installed dependencies and
# shares only checked-in upstream interfaces, fallback routines, and chief source.

include( ecwam_source_file_compile_flags )

set(ECWAM_HANDTUNED_GPU_MEMORY_MODE "separate" CACHE STRING
    "GPU memory mode inherited from the original ECWAM build")
set_property(CACHE ECWAM_HANDTUNED_GPU_MEMORY_MODE PROPERTY STRINGS separate unified)
if(NOT ECWAM_HANDTUNED_GPU_MEMORY_MODE MATCHES "^(separate|unified)$")
    message(FATAL_ERROR
        "ECWAM_HANDTUNED_GPU_MEMORY_MODE must be separate or unified")
endif()

file(GLOB FORTRAN_SOURCES CONFIGURE_DEPENDS "*.F90")
list(APPEND FORTRAN_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/../../original/source/wam_plume_module.F90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../original/source/timing.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../../original/source/timing_output_mpi.f90"
)
set(GPU_FORTRAN_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/cireduce_loki_gpu.sccs_hoist.F90"
    "${CMAKE_CURRENT_SOURCE_DIR}/outbs_loki_gpu.sccs_hoist.F90"
    "${CMAKE_CURRENT_SOURCE_DIR}/wamintgr_loki_gpu.sccs_hoist.F90"
)

if(HAVE_ACC)
    if(NOT CMAKE_CUDA_COMPILER_LOADED)
        enable_language(CUDA)
    endif()
    file(GLOB CUDA_SOURCES CONFIGURE_DEPENDS "*.cu")
else()
    list(REMOVE_ITEM FORTRAN_SOURCES ${GPU_FORTRAN_SOURCES})
    list(APPEND FORTRAN_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/../../original/source/ciwabr.F90"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../original/source/ciwaf.F90"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../original/source/wamintgr.F90"
    )
    set(CUDA_SOURCES)
endif()

foreach(source_file ${FORTRAN_SOURCES})
    get_filename_component(filename ${source_file} NAME)
    set_source_files_properties(${source_file} PROPERTIES
        COMPILE_DEFINITIONS "__FILENAME__=\"${filename}\""
    )
endforeach()

function(add_handtuned_precision precision)
    if(precision STREQUAL "sp")
        set(precision_definition WAM_HAVE_SINGLE_PRECISION)
    elseif(precision STREQUAL "dp")
        set(precision_definition WAM_HAVE_DOUBLE_PRECISION)
    else()
        message(FATAL_ERROR "Unsupported precision '${precision}'. Use sp or dp.")
    endif()

    set(library_target handtuned_${precision})
    set(accelerator_libraries)
    if(HAVE_ACC)
        list(APPEND accelerator_libraries
            OpenACC::OpenACC_Fortran
        )
        find_library(ECWAM_NVTX_FORTRAN_LIBRARY NAMES nvhpcwrapnvtx nvhpcnvtx
            HINTS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES} REQUIRED)
        list(APPEND accelerator_libraries ${ECWAM_NVTX_FORTRAN_LIBRARY})
    endif()

    ecbuild_add_library(
        TARGET           ${library_target}
        SOURCES          ${FORTRAN_SOURCES} ${CUDA_SOURCES}
        PUBLIC_LIBS      fiat parkind_${precision} ecwam_intfb
                         ${MPI_Fortran_LIBRARIES}
                         field_api_${precision}
        PRIVATE_LIBS     eccodes_f90
                         ${OpenMP_Fortran_LIBRARIES}
                         ${accelerator_libraries}
        PRIVATE_INCLUDES ${CMAKE_CURRENT_LIST_DIR}/../../original/source
        PRIVATE_DEFINITIONS
                         ${ECWAM_PRIVATE_DEFINITIONS}
                         ${precision_definition}
                         $<${HAVE_ACC}:WAM_GPU>
                         $<${HAVE_CUDA}:WAM_HAVE_CUDA>
                         $<${HAVE_CUDA}:WAM_HAVE_HIC>
                         $<${HAVE_GPU_AWARE_MPI}:WITH_GPU_AWARE_MPI>
        PUBLIC_DEFINITIONS ${ECWAM_DEFINITIONS}
    )

    set_property(TARGET ${library_target} APPEND PROPERTY BUILD_RPATH
        "$<TARGET_FILE_DIR:fiat>"
        "$<TARGET_FILE_DIR:parkind_${precision}>"
        "$<TARGET_FILE_DIR:field_api_${precision}>"
        "$<TARGET_FILE_DIR:eccodes_f90>"
    )
    set_target_properties(${library_target} PROPERTIES
        BUILD_WITH_INSTALL_RPATH OFF
        SKIP_BUILD_RPATH OFF
    )

    ecbuild_target_fortran_module_directory(
        TARGET                   ${library_target}
        MODULE_DIRECTORY         ${CMAKE_BINARY_DIR}/module/${library_target}
        INSTALL_MODULE_DIRECTORY module/${library_target}
    )

    if(precision STREQUAL "dp")
        foreach(autopromote_flag IN LISTS autopromote_flags)
            target_compile_options(${library_target} PRIVATE
                "$<$<COMPILE_LANGUAGE:Fortran>:${autopromote_flag}>"
            )
        endforeach()
    endif()

    if(HAVE_ACC AND CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
        if(NOT ECWAM_HANDTUNED_GPU_ARCH MATCHES "^[0-9]+$")
            message(FATAL_ERROR
                "ECWAM_HANDTUNED_GPU_ARCH must name the numeric NVIDIA architecture")
        endif()
        # Match the original implementation's OpenACC compile/link contract.
        target_link_options(${library_target} INTERFACE SHELL:${OpenACC_Fortran_FLAGS})
        target_compile_options(${library_target} PRIVATE
            "$<$<COMPILE_LANGUAGE:Fortran>:-O3>"
            "$<$<COMPILE_LANGUAGE:Fortran>:-fast>"
            "$<$<COMPILE_LANGUAGE:Fortran>:-gpu=cc${ECWAM_HANDTUNED_GPU_ARCH}>"
            "$<$<COMPILE_LANGUAGE:Fortran>:-gpu=fastmath>"
            "$<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math;-O3>"
        )
        if(HAVE_CUDA)
            target_compile_options(${library_target} PRIVATE
                "$<$<COMPILE_LANGUAGE:Fortran>:-cuda>")
            target_link_options(${library_target} PUBLIC
                "-cuda;-gpu=mem:separate:pinnedalloc")
        endif()
    endif()

    ecwam_target_compile_definitions_FILENAME(${library_target})
endfunction()

foreach(precision sp dp)
    if(HAVE_${precision})
        add_handtuned_precision(${precision})
    endif()
endforeach()
