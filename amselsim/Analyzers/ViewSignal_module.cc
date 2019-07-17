////////////////////////////////////////////////////////////////////////
// Class:       ViewSignal
// Plugin Type: analyzer (art v3_02_06)
// File:        ViewSignal_module.cc
//
// Generated at Fri Jul 12 13:56:27 2019 by Hunter Sullivan using cetskelgen
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
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "lardata/Utilities/LArFFT.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/Simulation/SimChannel.h"

class ViewSignal;


class ViewSignal : public art::EDAnalyzer {
public:
  explicit ViewSignal(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ViewSignal(ViewSignal const&) = delete;
  ViewSignal(ViewSignal&&) = delete;
  ViewSignal& operator=(ViewSignal const&) = delete;
  ViewSignal& operator=(ViewSignal&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  TH1F* hRawWaveform = nullptr;
  TH1F* hRawCharge   = nullptr;
};


ViewSignal::ViewSignal(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void ViewSignal::analyze(art::Event const& e)
{
  auto const *detprop = art::ServiceHandle<detinfo::DetectorPropertiesService const>{}->provider();

  // Read in the digit List object(s). 
  art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
  e.getByLabel("simpixels", "", digitVecHandle);

  //unsigned int signalSize = fNTicks;
  art::Handle< std::vector<sim::SimChannel> > SimListHandle;
  std::vector<art::Ptr<sim::SimChannel> > Simlist;
  if(e.getByLabel("largeant", SimListHandle))
     { art::fill_ptr_vector(Simlist, SimListHandle); }

  mf::LogInfo("ViewSignal") << "DigitVecHandle size is " << digitVecHandle->size();

  art::ServiceHandle<util::LArFFT> fFFT;
  int transformSize = fFFT->FFTSize();
  std::vector<short> rawadc(transformSize);  // vector holding uncompressed adc values

  float max(-10000);
  size_t maxCh(0);
  for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter)
  { 
    art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
    auto channel = digitVec->Channel();
    auto dataSize = digitVec->Samples();
    raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());

    float sum(0);
    for (int iBin = 0; iBin < dataSize; iBin++) 
    {
      if (rawadc[iBin] > 1) {std::cout << rawadc[iBin] << std::endl;sum += rawadc[iBin];}
      if (sum > max) {max = sum; maxCh = rdIter;}
    }
  }

  art::Ptr<raw::RawDigit> digitVec(digitVecHandle, maxCh);
  auto dataSize = digitVec->Samples();

  std::cout << dataSize << std::endl;
  raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
  for (int i = 1; i <= dataSize; i++) {std::cout << rawadc[i-1] << std::endl; hRawWaveform->SetBinContent(i, rawadc[i-1]);}

  maxCh = 0;
  max = 0;
  for ( int iCh = 0; iCh < (int)Simlist.size(); iCh++)
  {
    const auto& TDCIDEs = Simlist.at(iCh)->TDCIDEMap();
    if (TDCIDEs.size() > max) {max = TDCIDEs.size(); maxCh=iCh;}
  }

  const auto& TDCIDEs = Simlist.at(maxCh)->TDCIDEMap();
  for (const auto& TDCinfo : TDCIDEs)
  {
    auto tick = detprop->ConvertTDCToTicks(TDCinfo.first);
    for (const auto& ide : TDCinfo.second) 
    {
      auto content = hRawCharge->GetBinContent(tick+1);
      hRawCharge->SetBinContent(tick+1, content+ide.numElectrons);
    }
  }  
}

void ViewSignal::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  hRawWaveform = tfs->make<TH1F>("Raw", "Raw", 110000, 0, 110000);
  hRawCharge   = tfs->make<TH1F>("RawCharge", "Raw Charge", 110000, 0, 110000);
}

void ViewSignal::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ViewSignal)
