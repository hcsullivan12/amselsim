include_directories( $ENV{LARDATA_INC}    )
include_directories( $ENV{LARCORE_INC}    )
include_directories( $ENV{LARSIM_INC}     )
include_directories( $ENV{POSTGRESQL_INC} )
art_make(LIB_LIBRARIES amselsim_Geometry
         MODULE_LIBRARIES 
                          larevt_Filters
                          lardataobj_Simulation
                          amselsim_Geometry
                          amselsim_Services
                          amselsim_LArG4
                          amselsim_Services_AmSelSignalShapingService_service
                          lardataobj_RecoBase
			  lardataobj_RawData
                          lardata_ArtDataHelper
                          lardata_Utilities
                          nurandom_RandomUtils_NuRandomService_service
			  lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
                          nusimdata_SimulationBase
                          ${ART_FRAMEWORK_CORE}
                          ${ART_FRAMEWORK_PRINCIPAL}
                          ${ART_FRAMEWORK_SERVICES_REGISTRY}
                          ${ART_FRAMEWORK_SERVICES_BASIC}
                          ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                          ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                          ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR}
			  ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                          art_Persistency_Common canvas
                          art_Persistency_Provenance canvas
                          art_Utilities canvas
                          ${MF_MESSAGELOGGER}
                          
                          ${FHICLCPP}
                          cetlib cetlib_except
                          ${CLHEP}
                          ${ROOT_BASIC_LIB_LIST}
)
install_headers()
install_fhicl()
install_source()
