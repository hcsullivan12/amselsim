/**
 * @file SimQPix_module.cc
 * @brief Module to simulate signal for QPix electronics
 * 
 * @author H. Sullivan (hsulliva@fnal.gov)
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/Utilities/LArFFT.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include "amselsim/Geometry/DetectorGeometryService.h"
#include "amselsim/Services/AmSelSignalShapingService.h"

#include "CLHEP/Random/RandGaussQ.h"

#include <memory>

namespace qpixsim
{

class SimQPix;


class SimQPix : public art::EDProducer {
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
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const& p);

  void GenNoiseInTime(std::vector<float> &noise);

private:

  CLHEP::HepRandomEngine& fEngine; ///< reference to art-managed random-number engine
  float                   fNoiseMean;
  float                   fNoiseSigma;
};

//--------------------------------------------------------------------
SimQPix::SimQPix(fhicl::ParameterSet const& p)
  : EDProducer{p}
   , fEngine(art::ServiceHandle<rndm::NuRandomService>{}
                ->createEngine(*this, "HepJamesRandom", "propagation", p, "PropagationSeed"))  
{
  this->reconfigure(p);

  produces< std::vector<raw::RawDigit>   >();
  //produces< std::vector<recob::Wire> >(fSpillName);
  //produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
}

//--------------------------------------------------------------------
void SimQPix::reconfigure(fhicl::ParameterSet const& p)
{
  fNoiseMean = p.get< double >("NoiseMean");
  fNoiseSigma = p.get< double >("NoiseSigma");
}

//--------------------------------------------------------------------
void SimQPix::produce(art::Event& e)
{
  auto const *detprop = art::ServiceHandle<detinfo::DetectorPropertiesService const>{}->provider();
  auto const *ts      = lar::providerFrom<detinfo::DetectorClocksService>();
  auto const *geom    = art::ServiceHandle<geo::DetectorGeometryService>()->provider();
  art::ServiceHandle<util::AmSelSignalShapingService> sss;
  art::ServiceHandle<util::LArFFT> fft;
  auto nTicks = fft->FFTSize();

  art::Handle<std::vector<sim::SimChannel>> SimListHandle;
  std::vector<art::Ptr<sim::SimChannel>> Simlist;
  if (e.getByLabel("largeant", SimListHandle))
  {
    art::fill_ptr_vector(Simlist, SimListHandle);
  }

  const auto NChannels = geom->Nchannels();

  // Containers for storing signals
  std::vector<short> adcvec(nTicks, 0);
  std::vector<float> chargeWork(nTicks,0.);

  // Make our collections
  std::unique_ptr< std::vector<raw::RawDigit>> digcol(new std::vector<raw::RawDigit>);
  digcol->reserve(NChannels);

  // Loop over channels
  for (int iCh = 0; iCh < (int)Simlist.size(); iCh++)
  {
    auto channel = Simlist.at(iCh)->Channel();

    // Clear out our containers
    std::fill(chargeWork.begin(), chargeWork.end(), 0.);

    const auto& TDCIDEs = Simlist.at(iCh)->TDCIDEMap(); 
    for (const auto& TDCinfo : TDCIDEs)
    {
      auto tick = detprop->ConvertTDCToTicks(TDCinfo.first);
      float totalQ(0);
      //std::cout << TDCinfo.first << " " << detprop->ConvertTDCToTicks(TDCinfo.first) << std::endl;
      for (const auto& ide : TDCinfo.second) totalQ += ide.numElectrons;

      chargeWork[tick] = totalQ;
    }

    // Convolve charge with appropriate response function
    sss->Convolute(chargeWork);

    // Generate noise
    std::vector<float> noisetmp(nTicks,0.);
    GenNoiseInTime(noisetmp);

    for (int i = 0; i < nTicks; ++i)
    {
      float adcval = noisetmp[i] + chargeWork[i];
      adcvec[i] = (unsigned short)(adcval);
    } // end loop over signal size

    raw::Compress(adcvec, raw::kNone); 
      
    // add this digit to the collection
    raw::RawDigit rd(channel,nTicks,adcvec,raw::kNone);// fNTimeSamples, adcvec, fCompression);
    rd.SetPedestal(0);
    digcol->push_back(rd);
  } // end loop over channels

  e.put(std::move(digcol));
  return;

}

//--------------------------------------------------------------------                                                                                                 
void SimQPix::GenNoiseInTime(std::vector<float> &noise)
{
  CLHEP::RandGaussQ rGauss(fEngine,fNoiseMean,fNoiseSigma);
  for (unsigned int i=0; i<noise.size(); i++) noise[i] = rGauss.fire();
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
