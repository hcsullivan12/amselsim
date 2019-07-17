cet_enable_asserts()

cet_make(         LIBRARIES lardata_Utilities
                            lardata_Utilities_LArFFT_service
                            ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            art_Persistency_Common canvas
                            art_Persistency_Provenance canvas
                            ${ART_FRAMEWORK_SERVICES_BASIC}
                            ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
		            ${ART_ROOT_IO_TFILESERVICE_SERVICE}
		            ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
	 	            art_Utilities canvas
		            ${MF_MESSAGELOGGER}
		            cetlib cetlib_except
		            ${ROOT_BASIC_LIB_LIST}
		            ${ROOT_GEOM}
                            ${FHICLCPP}
                            nutools_IFDatabase
)

simple_plugin( DetectorPropertiesServiceAmSel  "service"
               amselsim_Services
               ${SQLITE3}
               ${ART_ROOT_IO_ROOTDB}
               wda
               ${MF_MESSAGELOGGER}       
               cetlib cetlib_except
               ${ROOT_BASIC_LIB_LIST}
               ${Boost_SYSTEM_LIBRARY}
)

simple_plugin( AmSelSignalShapingService  "service"
               amselsim_Services
               lardata_Utilities_LArFFT_service
)


install_headers()
install_fhicl()
install_source()
