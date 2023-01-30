ecbuild_add_executable (
  TARGET RadApp.x
  SOURCES RadApp.F90 GEOS_RadiationEnvGridComp.F90 RootRad_GridComp.F90
  LIBS GEOSradiation_GridComp GEOSchem_GridComp esmf OpenMP::OpenMP_Fortran
  )
target_link_libraries (RadApp.x ${MKL_LIBRARIES})

set (CMAKE_Fortran_FLAGS_RELEASE "${GEOS_Fortran_FLAGS_VECT}")

file (GLOB templates CONFIGURE_DEPENDS *.tmpl)

set (programs
   rad_run.j
   construct_extdata_yaml_list.py
   )

install (
   FILES ${templates}
   DESTINATION etc
   )

install (
   PROGRAMS ${programs}
   DESTINATION bin
   )

install (
   FILES logging.yaml
   DESTINATION etc
   )

if(INSTALL_SOURCE_TARFILE)
  set(CFG_INSTALL_SOURCE_TARFILE TRUE)
else()
  set(CFG_INSTALL_SOURCE_TARFILE FALSE)
endif()

set (setup_scripts
   rad_setup
   )

foreach (file ${setup_scripts})
   configure_file(${file} ${file} @ONLY)
   install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${file} DESTINATION bin)
endforeach ()

configure_file(.RADAPP_VERSION .RADAPP_VERSION @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/.RADAPP_VERSION DESTINATION etc)

