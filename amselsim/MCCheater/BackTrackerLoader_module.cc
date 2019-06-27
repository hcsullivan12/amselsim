////////////////////////////////////////////////////////////////////////
// Class:       BackTrackerLoader
// Module Type: producer
// File:        BackTrackerLoader.h
//
// Generated at Thu Jun 14 06:49:31 2012 by Brian Rebel using artmod
// from art v1_00_11.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

namespace cheat {
  class BackTrackerLoader;
}

class cheat::BackTrackerLoader : public art::EDProducer {
public:
  explicit BackTrackerLoader(fhicl::ParameterSet const & p);

  virtual void produce(art::Event & e);


private:

  // Declare member data here.

};

//------------------------------------------------------------------------------
cheat::BackTrackerLoader::BackTrackerLoader(fhicl::ParameterSet const & p)
  : EDProducer{p}
{
  // Call appropriate Produces<>() functions here.
}

//------------------------------------------------------------------------------
// the sole purpose of this module is to issue the Rebuild command to the
// BackTracker service.  It should be put after all simulation data producing
// modules have run in the job, and only in jobs that create the simulation and
// then make use of the BackTracker in either cheating reconstruction modules
// or analyzers
void cheat::BackTrackerLoader::produce(art::Event & e)
{
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  pi_serv->Rebuild(e);
  bt_serv->Rebuild(e);
}


DEFINE_ART_MODULE(cheat::BackTrackerLoader)
