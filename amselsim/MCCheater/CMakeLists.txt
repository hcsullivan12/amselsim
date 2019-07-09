cet_make_library(
  LIBRARY_NAME amselsim_MCCheater_ParticleInventory
  SOURCE ParticleInventory.cc ParticleInventory.tcc
  LIBRARIES
  canvas
  lardata_Utilities
  larsim_Simulation
  nutools_ParticleNavigation
  lardataobj_Simulation
  nusimdata_SimulationBase
  ${MF_MESSAGELOGGER}
  )

cet_make_library(
  LIBRARY_NAME amselsim_MCCheater_PhotonBackTracker
  SOURCE PhotonBackTracker.cc PhotonBackTracker.tcc
  LIBRARIES
  canvas
  lardata_Utilities
  larsim_Simulation
  amselsim_MCCheater_ParticleInventory
  larcorealg_Geometry
  nutools_ParticleNavigation
  lardataobj_Simulation
  lardataobj_RawData
  nusimdata_SimulationBase
  )

cet_make_library(
  LIBRARY_NAME amselsim_MCCheater_BackTracker
  SOURCE BackTracker.cc BackTracker.tcc
  LIBRARIES
  canvas
  lardata_Utilities
  larsim_Simulation
  amselsim_MCCheater_ParticleInventory
  larcorealg_Geometry
  nutools_ParticleNavigation
  lardataobj_Simulation
  nusimdata_SimulationBase
  )

simple_plugin(ParticleInventoryService "service"
  amselsim_MCCheater_ParticleInventory
  larsim_Simulation
  lardata_Utilities
  nutools_ParticleNavigation
  lardataobj_Simulation
  nusimdata_SimulationBase
  ${MF_MESSAGELOGGER}
  )

simple_plugin(PhotonBackTrackerService "service"
  amselsim_MCCheater_PhotonBackTracker
  amselsim_MCCheater_ParticleInventory
  amselsim_MCCheater_ParticleInventoryService_service
  lardata_Utilities
  larsim_Simulation
  nutools_ParticleNavigation
  lardataobj_Simulation
  lardataobj_RecoBase
  lardataobj_RawData
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  nusimdata_SimulationBase
  ${MF_MESSAGELOGGER}
  )

simple_plugin(BackTrackerService "service"
  amselsim_MCCheater_BackTracker
  amselsim_MCCheater_ParticleInventory
  amselsim_MCCheater_ParticleInventoryService_service
  lardata_Utilities
  larsim_Simulation
  nutools_ParticleNavigation
  lardataobj_Simulation
  lardataobj_RecoBase
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  nusimdata_SimulationBase
  ${MF_MESSAGELOGGER}
  )

simple_plugin(CheckBackTracking "module"
  amselsim_MCCheater_BackTrackerService_service
  amselsim_MCCheater_ParticleInventoryService_service
  lardataobj_RecoBase
  larsim_Simulation
  nutools_ParticleNavigation
  lardataobj_Simulation
  ${MF_MESSAGELOGGER}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  )

simple_plugin(RecoCheckAna "module"
  amselsim_MCCheater_BackTrackerService_service
  amselsim_MCCheater_ParticleInventoryService_service
  lardataobj_RecoBase
  larsim_Simulation
  nutools_ParticleNavigation
  lardataobj_Simulation
  nusimdata_SimulationBase
  ${MF_MESSAGELOGGER}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  art_Persistency_Provenance
  ${ROOT_TREE}
  ${ROOT_HIST}
  )

simple_plugin(BackTrackerLoader "module"
  amselsim_MCCheater_BackTrackerService_service
  amselsim_MCCheater_ParticleInventoryService_service
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  )
simple_plugin(PhotonBackTrackerLoader "module"
  amselsim_MCCheater_ParticleInventoryService_service
  amselsim_MCCheater_PhotonBackTrackerService_service
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  )


install_headers()
install_fhicl()
install_source()
