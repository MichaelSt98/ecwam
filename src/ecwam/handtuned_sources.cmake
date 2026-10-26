# Complete checked-in handtuned SP source set (Lukas revision 283d1e6273ed).
# Keep the RAPS project, interface generation, coupling dependencies, UNWAM
# sources, target name and public definitions; replace only the wave sources.
file( GLOB handtuned_sources CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/handtuned/*.F90"
  "${PROJECT_SOURCE_DIR}/handtuned/*.f90" )

# The handtuned project builds these timing files from its original tree.
# This vendored copy includes them in handtuned/ instead.
set( ecwam_srcs ${handtuned_sources} )
# Use the shared obstruction reader and helper so the MPI decoding backport
# also applies to the handtuned SP build.
list( REMOVE_ITEM ecwam_srcs "${PROJECT_SOURCE_DIR}/handtuned/getgrbobstrct.F90" )
list( APPEND ecwam_srcs
  "${CMAKE_CURRENT_SOURCE_DIR}/getgrbobstrct.F90"
  "${CMAKE_CURRENT_SOURCE_DIR}/cpgrbobstrct.F90" )
if( HAVE_UNWAM )
  list( APPEND ecwam_srcs ${unwam_srcs} )
endif()

# RAPS optional plugin source is not part of the upstream handtuned tree.
if( HAVE_PLUME )
  list( APPEND ecwam_srcs ${CMAKE_CURRENT_SOURCE_DIR}/wam_plume_module.F90 )
endif()
