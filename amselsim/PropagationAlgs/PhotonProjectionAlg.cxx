/**
 * @file PhotonProjectionAlg.cxx
 * @brief Algorithm for emission and projection of 
 *        scintillation light onto the readout plane.
 * 
 * @todo Switch largeant to use scintillation by particle type
 * 
 * @author H. Sullivan (hsulliva@fnal.gov)
 */

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"

#include "amselsim/LArG4/IonizationAndScintillation.h"
#inclued "amselsim/PropagationAlgs/PhotonProjectionAlg.h"
#include "amselsim/Geometry/AmSelGeometryService.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace amselsim
{

//--------------------------------------------------------------------
PhotonProjectionAlg::PhotonProjectionAlg()
  : fEngine(art::ServiceHandle<rndm::NuRandomService>{}
                ->createEngine(*this, "HepJamesRandom", "propagation", p, "PropagationSeed")),
    fNScint(0),
    fNHits(0)
{}

//--------------------------------------------------------------------
void PhotonProjectionAlg::reconfigure(fhicl::ParameterSet const& p)
{}

/**
 * @brief Do projection of scintillation photons. 
 * 
 * For each G4 step, isotropically fire stepScint photons from the 
 * scintillation vertex, keeping only those that are incident on the 
 * readout plane. Presumably, the incident distribution is sparse, 
 * so filling a container of length = Npixels is unecessary and 
 * impractical for large pixelated detectors. Instead, a map 
 * container keeps track of only those pixels which have non zero 
 * incident pixels. 
 * 
 */
void PhotonProjectionAlg::doProjection()
{
  std::cout << "////////////////////////////////////////////////\n"
            << "Starting photon projection...\n";

  // Get the data from IS action
  auto stepPoints = IonizationAndScintillation::Instance()->StepPoints();
  auto stepScint  = IonizationAndScintillation::Instance()->StepScint();

  amselgeo::AmSelGeometry const* geom = art::ServiceHandle<amselgeo::AmSelGeometryService>()->provider();
  double detHalfHeight = geom->DetHalfHeight();
  double driftLength   = geom->DetDriftLength();
  double detLength     = geom->DetLength();
 
  CLHEP::RandFlat flat(fEngine);

  if (stepPoints.size() != stepScint.size()) throw cet::exception("PhotonProjectionAlg") << "StepPoints and StepScint sizes are different!\n";
  if (!stepPoints.size())                    mf::LogWarning("PhotonProjectionAlg") << "No steps recorded!\n";
  
  std::cout << "Number of steps = " << stepPoints.size() << "\n";
  size_t nHit(0);
  std::map<int, int> pixelMap;
  for (size_t iStep = 0; iStep < stepPoints.size(); iStep++)
  {
    // convert to cm
    G4ThreeVector posTemp = stepPoints[iStep]/CLHEP::cm;
    TVector3      pos(posTemp.x(), posTemp.y(), posTemp.z());
    size_t        nScint = stepScint[iStep];
    fNScint += nScint;
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
          projPos.Y() <= -1*driftLength || projPos.Y() >= driftLength) continue;

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
  std::cout << "Number of scintillation photons:     " << fNScint << std::endl;
  std::cout << "Number of photons incident on plane: " << fNHits  << std::endl;
  std::cout << "////////////////////////////////////////////////\n";
}

}