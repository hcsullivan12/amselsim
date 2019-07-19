/**
 * @file PhotonProjection_module.cc
 * @brief Module to project scintillation photon onto readout plane
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
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include "amselsim/Geometry/DetectorGeometryService.h"
#include "amselsim/LArG4/IonizationAndScintillation.h"

#include "TH1.h"

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include <memory>

namespace amselsim
{

class PhotonProjection;


class PhotonProjection : public art::EDProducer {
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
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const& p);

private:

  CLHEP::HepRandomEngine& fEngine; ///< reference to art-managed random-number engine

  TH1I* hTotScint = nullptr;
  TH1I* hTotHits  = nullptr;

};

//--------------------------------------------------------------------
PhotonProjection::PhotonProjection(fhicl::ParameterSet const& p)
  : EDProducer{p}
   , fEngine(art::ServiceHandle<rndm::NuRandomService>{}
                ->createEngine(*this, "HepJamesRandom", "propagation", p, "PropagationSeed"))  
{
  this->reconfigure(p);

  produces< std::map<int, int> >();
}

//--------------------------------------------------------------------
void PhotonProjection::reconfigure(fhicl::ParameterSet const& p)
{
}

//--------------------------------------------------------------------
void PhotonProjection::produce(art::Event& e)
{
  // Get our services 
  auto const *detprop = art::ServiceHandle<detinfo::DetectorPropertiesService const>{}->provider();
  auto const *ts      = lar::providerFrom<detinfo::DetectorClocksService>();
  auto const *geom    = art::ServiceHandle<geo::DetectorGeometryService>()->provider();

  // Make our new collection
  std::unique_ptr< std::map<int, int>> photcol(new std::map<int, int>);

 
  std::cout << "////////////////////////////////////////////////\n"
            << "Starting photon projection...\n";

  // Get the data from IS action
  auto stepPoints = amselg4::IonizationAndScintillation::Instance()->StepPoints();
  auto stepScint  = amselg4::IonizationAndScintillation::Instance()->StepScint();

  double detHalfHeight = geom->DetHalfHeight();
  double detLength     = geom->DetLength();
 
  CLHEP::RandFlat flat(fEngine);

  if (stepPoints.size() != stepScint.size()) throw cet::exception("PhotonProjection") << "StepPoints and StepScint sizes are different!\n";
  if (!stepPoints.size())                    mf::LogWarning("PhotonProjection") << "No steps recorded!\n";
  
  std::cout << "Number of steps = " << stepPoints.size() << "\n";

  int totHits(0);
  int totScint(0);
  photcol->clear();
  for (size_t iStep = 0; iStep < stepPoints.size(); iStep++)
  {
    // convert to cm
    G4ThreeVector posTemp = stepPoints[iStep]/CLHEP::cm;
    TVector3      pos(posTemp.x(), posTemp.y(), posTemp.z());
    // How much scintillation in this step?
    int           nScint = stepScint[iStep];
    totScint += nScint;

    // Only do half since there's no chance for the other half
    // In this case, limit the range of phi to be [90,270]
    for (size_t iPh = 0; iPh < 0.5*nScint; iPh++)
    {
      double cosTheta = 2*flat.fire() - 1;
      double sinTheta = std::pow(1-std::pow(cosTheta,2),0.5);
      double phi      = M_PI*flat.fire() + 0.5*M_PI;
      TVector3 pHat(sinTheta*std::cos(phi),  sinTheta*std::sin(phi), cosTheta);

      // Still check...
      if (pHat.X() >= 0) continue;

      // Project this onto the x = 0 plane
      double t  = -1*pos.X()/pHat.X();
      TVector3 projPos = pos + t * pHat;

      // Assuming coordinate system is centered in y
      if (projPos.Z() <= 0                || projPos.Z() >= detLength ||
          projPos.Y() <= -1*detHalfHeight || projPos.Y() >= detHalfHeight) continue;

      // This hit the readout plane
      totHits++;
  
      /// \todo Handle the x coordinate better here 
      TVector3 point(-0.01, projPos.Y(), projPos.Z());
      auto nearestPixelId = geom->NearestReadoutNodeID(point); 
      if (nearestPixelId < 0) continue;
 
      // Add this photon to this pixel
      (*photcol)[nearestPixelId]++;
      //auto pixelIter = fPixelMap.find(nearestPixelId);
      //if (pixelIter != fPixelMap.end()) pixelIter->second++;
      //else fPixelMap.emplace(nearestPixelId, 1);
    }
  }
  std::cout << "Number of scintillation photons:     " << totScint << std::endl;
  std::cout << "Number of photons incident on plane: " << totHits  << std::endl;
  std::cout << "////////////////////////////////////////////////\n";

  hTotScint->Fill(totScint);
  hTotHits->Fill(totHits);

  e.put(std::move(photcol));
  return;
}

//--------------------------------------------------------------------
void PhotonProjection::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  // In 10s of 1000s
  hTotScint = tfs->make<TH1I>("hTotScint", "Total amount of scintillation", 1000, 0, 10000);
  hTotScint->GetXaxis()->SetTitle("nScint x 10^{4}");
  hTotHits  = tfs->make<TH1I>("hTotHits", "Total amount of scintillation incident on readout plane", 1000, 0, 10000);
  hTotHits->GetXaxis()->SetTitle("nHits x 10^{4}");
}

//--------------------------------------------------------------------
void PhotonProjection::endJob()
{
  // Implementation of optional member function here.
}


DEFINE_ART_MODULE(PhotonProjection)
}
