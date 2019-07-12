////////////////////////////////////////////////////////////////////////
// Class:       SimQPix
// Plugin Type: analyzer (art v3_02_06)
// File:        SimQPix_module.cc
//
// Generated at Thu Jul 11 17:08:14 2019 by Hunter Sullivan using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "amselsim/Geometry/DetectorGeometryService.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/Simulation/SimChannel.h"

namespace amselsim
{

class SimQPix;


class SimQPix : public art::EDAnalyzer {
public:
  explicit SimQPix(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimQPix(SimQPix const&) = delete;
  SimQPix(SimQPix&&) = delete;
  SimQPix& operator=(SimQPix const&) = delete;
  SimQPix& operator=(SimQPix&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;


private:

  // Declare member data here.
  std::string fDriftEModuleLabel = "largeant";

};

//--------------------------------------------------------------------
SimQPix::SimQPix(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} 
{

}

//--------------------------------------------------------------------
void SimQPix::analyze(art::Event const &e)
{
   auto const *detprop = art::ServiceHandle<detinfo::DetectorPropertiesService const>{}->provider();

  // In case trigger simulation is run in the same job...
  //FIXME: you should never call preProcessEvent
  auto const *ts = lar::providerFrom<detinfo::DetectorClocksService>();

  // get the geometry to be able to figure out signal types and chan -> plane mappings
  auto const *geom = art::ServiceHandle<geo::DetectorGeometryService>()->provider();

  //unsigned int signalSize = fNTicks;
  art::Handle< std::vector<sim::SimChannel> > SimListHandle;
  std::vector<art::Ptr<sim::SimChannel> > Simlist;
  if(e.getByLabel("largeant", SimListHandle))
     { art::fill_ptr_vector(Simlist, SimListHandle); }
  
  const auto NChannels = geom->Nchannels();

  // Loop over channels
  for (int iCh = 0; iCh < (int)Simlist.size(); iCh++)
  {
    float totalQ(0);
    const auto& TDCIDEs = Simlist.at(iCh)->TDCIDEMap(); 
    for (const auto& TDCinfo : TDCIDEs)
    {
      //std::cout << TDCinfo.first << " " << detprop->ConvertTDCToTicks(TDCinfo.first) << std::endl;
      for (const auto& ide : TDCinfo.second)
      {
        totalQ += ide.numElectrons;
      }
    }
  } // end loop over channels
}

//--------------------------------------------------------------------
void SimQPix::beginJob()
{
  // Implementation of optional member function here.
}

//--------------------------------------------------------------------
void SimQPix::endJob()
{
  // Implementation of optional member function here.
}


DEFINE_ART_MODULE(SimQPix)
}
