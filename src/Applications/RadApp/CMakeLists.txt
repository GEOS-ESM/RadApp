ecbuild_add_executable (
  TARGET RadApp.x
  SOURCES RadApp.F90 GEOS_RadiationEnvGridComp.F90 RootRad_GridComp.F90
  LIBS GEOSradiation_GridComp GEOSchem_GridComp
  )

set (CMAKE_Fortran_FLAGS_RELEASE "${GEOS_Fortran_FLAGS_VECT}")

#link_directories (${MKL_LIBRARIES})
target_link_libraries (RadApp.x ${OpenMP_Fortran_LIBRARIES} ${MKL_LIBRARIES})
set_target_properties(RadApp.x PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
target_include_directories (RadApp.x PUBLIC ${INC_ESMF})

