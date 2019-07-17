art_make(
  LIBRARY_NAME amselsim_PropagationAlgs
  LIB_LIBRARIES nurandom_RandomUtils_NuRandomService_service
                amselsim_LArG4
                ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
                ${ART_ROOT_IO_TFILESERVICE_SERVICE}
        
)

install_headers()
install_fhicl()
install_source()
