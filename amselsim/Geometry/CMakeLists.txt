art_make_library( LIBRARY_NAME amselsim_Geometry
                  SOURCE AmSelGeometry.cxx
                  LIBRARIES larcorealg_Geometry
                            ${MF_MESSAGELOGGER}
                            ${MF_UTILITIES}
                            ${FHICLCPP}
                            ${ROOT_BASIC_LIB_LIST} 
                            ${ROOT_GEOM}
                            cetlib cetlib_except
)

simple_plugin( AmSelGeometryService service
                amselsim_Geometry
		${ART_FRAMEWORK_CORE}
		${ART_FRAMEWORK_IO_SOURCES}
		${ART_FRAMEWORK_PRINCIPAL}
		art_Persistency_Common canvas
		art_Persistency_Provenance canvas
		art_Utilities canvas
		${ART_FRAMEWORK_SERVICES_REGISTRY}
		${MF_MESSAGELOGGER}
		${MF_UTILITIES}
		${ROOT_BASIC_LIB_LIST}                
		cetlib cetlib_except
)

install_headers()
install_fhicl()
install_source()

# testing
add_subdirectory(gdml)
