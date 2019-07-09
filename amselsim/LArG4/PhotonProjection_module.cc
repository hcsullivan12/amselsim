///////////////////////////////////////////////////////////////
/// \file PhotonProjection_module.cc
/// \brief Analyzer module that projects scintillation light 
///        onto the readout plane.
/// 
/// \author H. Sullivan (hsulliva@fnal.gov)
///////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "art_root_io/TFileService.h"

#include "nurandom/RandomUtils/NuRandomService.h"

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include "amselsim/LArG4/IonizationAndScintillation.h"
#include "amselsim/Geometry/AmSelGeometryService.h"

#include "TTree.h"

namespace amselg4 
{

class PhotonProjection;

class PhotonProjection : public art::EDAnalyzer {
public:
  explicit PhotonProjection(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhotonProjection(PhotonProjection const&) = delete;
  PhotonProjection(PhotonProjection&&) = delete;
  PhotonProjection& operator=(PhotonProjection const&) = delete;
  PhotonProjection& operator=(PhotonProjection&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const& p);
  void resetVars();

private:

  CLHEP::HepRandomEngine& fEngine;
  TTree* fTree;

  ULong64_t fNScint;
  ULong64_t fNHits;
  int       fNPixels;
  float     fPixelSpacing;
  float     fPrimX0;
  std::vector<int> fPixelIdVec;
  std::vector<int> fPixelCountVec;
};


PhotonProjection::PhotonProjection(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fEngine(art::ServiceHandle<rndm::NuRandomService>{}
                ->createEngine(*this, "HepJamesRandom", "propagation", p, "PropagationSeed"))

{
  //this->reconfigure(p);
}

void PhotonProjection::reconfigure(fhicl::ParameterSet const& p)
{
}

void PhotonProjection::analyze(art::Event const& e)
{
  resetVars();

  std::cout << "////////////////////////////////////////////////\n"
            << "Starting photon projection...\n";

  // Get the data from IS action
  auto stepPoints = IonizationAndScintillation::Instance()->StepPoints();
  auto stepScint  = IonizationAndScintillation::Instance()->StepScint();

  amselgeo::AmSelGeometry const* geom = art::ServiceHandle<amselgeo::AmSelGeometryService>()->provider();
  double detHalfHeight = geom->DetHalfHeight();
  double detHalfWidth  = geom->DetHalfWidth();
  double detLength     = geom->DetLength();
 
  CLHEP::RandFlat flat(fEngine);

  if (stepPoints.size() != stepScint.size()) throw cet::exception("PhotonProjection") << "StepPoints and StepScint sizes are different!\n";

  if (!stepPoints.size()) std::cout << "No steps recorded!\n";
  
  // For each step, isotropically fire stepScint, keeping only those that are incident on the readout plane
  
  std::cout << "Number of steps = " << stepPoints.size() << "\n";
  size_t nHit(0);
  std::map<int, int> pixelMap;
  for (size_t iStep = 0; iStep < stepPoints.size(); iStep++)
  {
    // convert to cm
    G4ThreeVector posTemp = stepPoints[iStep]/CLHEP::cm;
    TVector3      pos(posTemp.x(), posTemp.y(), posTemp.z());
    size_t        nScint = stepScint[iStep];
    for (size_t iPh = 0; iPh < nScint; iPh++)
    {
      double cosTheta = 2*flat.fire() - 1;
      double sinTheta = std::pow(1-std::pow(cosTheta,2),0.5);
      double phi      = 2*M_PI*flat.fire();
      TVector3 pHat(sinTheta*std::cos(phi),  sinTheta*std::sin(phi), cosTheta);

      // No chance...
      if (pHat.X() >= 0) continue;

      // Project this onto the x = 0 plane
      double t  = -1*pos.X()/pHat.X();
      TVector3 projPos = pos + t * pHat;

      // Assuming coordinate system is centered on front face of argon slab (it should be)
      if (projPos.Z() <= 0               || projPos.Z() >= detLength ||
          projPos.Y() <= -1*detHalfWidth || projPos.Y() >= detHalfWidth) continue;

      // This hit the readout plane
      nHit++;
    
      /// \todo Handle the x coordinate better here 
      TVector3 point(-0.01, projPos.Y(), projPos.Z());
      auto nearestPixelId = geom->NearestPixelID(point); 
      if (nearestPixelId < 0) continue;
 
      // Add this photon to this pixel
      auto pixelIter = pixelMap.find(nearestPixelId);
      if (pixelIter != pixelMap.end()) pixelIter->second++;
      else pixelMap.emplace(nearestPixelId, 1);
    }
  }
  fNHits = nHit;
  std::cout << "Number of photons incident on plane: " << nHit << std::endl;
  std::cout << "////////////////////////////////////////////////\n";

  // Fill our containers
  fPixelIdVec.reserve(pixelMap.size());
  fPixelCountVec.reserve(pixelMap.size());
  for (const auto& p : pixelMap)
  {
    fPixelIdVec.push_back(p.first);
    fPixelCountVec.push_back(p.second);
  }


  // Now fill some info on primary
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();
  art::Handle< std::vector<simb::MCParticle> > plistHandle;
  e.getByLabel("largeant", plistHandle);


  for (size_t p = 0; p < plistHandle->size(); p++)
  {
    simb::MCParticle part = (*plistHandle)[p];
    if (part.Process() != "primary") continue;
    fPrimX0 = part.Vx();
  }

  fTree->Fill();
}

void PhotonProjection::resetVars()
{
  fNScint = 0;
  fNHits = 0;
  fPrimX0 = -99999;
  fPixelIdVec.clear();
  fPixelCountVec.clear();
}

void PhotonProjection::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("anatree","analysis tree");

  fTree->Branch("nScint",       &fNScint,       "nScint/l");
  fTree->Branch("nHits",        &fNHits,        "nHits/l");
  fTree->Branch("nPixels",      &fNPixels,      "nPixels/l");
  fTree->Branch("pixelSpacing", &fPixelSpacing, "pixelSpacing/D");
  fTree->Branch("primX0", &fPrimX0, "primX0/D");
  fTree->Branch("pixelIDs" ,    &fPixelIdVec    );
  fTree->Branch("pixelHits" ,   &fPixelCountVec );

  // How many pixels are we dealing with here?
  amselgeo::AmSelGeometry const* geom = art::ServiceHandle<amselgeo::AmSelGeometryService>()->provider();
  fPixelSpacing = geom->PixelSpacing();
  fNPixels      = geom->NPixels();
}

void PhotonProjection::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(PhotonProjection)
}
