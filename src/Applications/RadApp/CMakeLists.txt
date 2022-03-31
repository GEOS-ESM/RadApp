ecbuild_add_executable (
  TARGET RadApp.x
  SOURCES RadApp.F90 GEOS_RadiationEnvGridComp.F90 RootRad_GridComp.F90
  LIBS GEOSradiation_GridComp GEOSchem_GridComp esmf OpenMP::OpenMP_Fortran
  )
target_link_libraries (RadApp.x ${MKL_LIBRARIES})
