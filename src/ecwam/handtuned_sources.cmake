# Complete checked-in handtuned SP source set (Lukas revision 283d1e6273ed).
# Keep the RAPS project, interface generation, coupling dependencies, UNWAM
# sources, target name and public definitions; replace only the wave sources.
file( GLOB handtuned_sources CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/handtuned/*.F90"
  "${PROJECT_SOURCE_DIR}/handtuned/*.f90" )

# The handtuned project builds these timing files from its original tree.
# This vendored copy includes them in handtuned/ instead.
set( ecwam_srcs ${handtuned_sources} )
if( HAVE_UNWAM )
  list( APPEND ecwam_srcs ${unwam_srcs} )
endif()

# RAPS optional plugin source is not part of the upstream handtuned tree.
if( HAVE_PLUME )
  list( APPEND ecwam_srcs ${CMAKE_CURRENT_SOURCE_DIR}/wam_plume_module.F90 )
endif()
