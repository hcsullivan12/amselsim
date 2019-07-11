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
  void beginRun(art::Run const& r) override;
  void beginSubRun(art::SubRun const& sr) override;
  void endJob() override;
  void endRun(art::Run const& r) override;
  void endSubRun(art::SubRun const& sr) override;
  void respondToCloseInputFile(art::FileBlock const& fb) override;
  void respondToCloseOutputFiles(art::FileBlock const& fb) override;
  void respondToOpenInputFile(art::FileBlock const& fb) override;
  void respondToOpenOutputFiles(art::FileBlock const& fb) override;

private:

  // Declare member data here.

};


SimQPix::SimQPix(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void SimQPix::analyze(art::Event const& e)
{
  // Implementation of required member function here.
}

void SimQPix::beginJob()
{
  // Implementation of optional member function here.
}

void SimQPix::beginRun(art::Run const& r)
{
  // Implementation of optional member function here.
}

void SimQPix::beginSubRun(art::SubRun const& sr)
{
  // Implementation of optional member function here.
}

void SimQPix::endJob()
{
  // Implementation of optional member function here.
}

void SimQPix::endRun(art::Run const& r)
{
  // Implementation of optional member function here.
}

void SimQPix::endSubRun(art::SubRun const& sr)
{
  // Implementation of optional member function here.
}

void SimQPix::respondToCloseInputFile(art::FileBlock const& fb)
{
  // Implementation of optional member function here.
}

void SimQPix::respondToCloseOutputFiles(art::FileBlock const& fb)
{
  // Implementation of optional member function here.
}

void SimQPix::respondToOpenInputFile(art::FileBlock const& fb)
{
  // Implementation of optional member function here.
}

void SimQPix::respondToOpenOutputFiles(art::FileBlock const& fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(SimQPix)
