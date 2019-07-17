/**
 * @file PhotonProjectionAlgAlg.h
 * @brief Algorithm for emission and projection of 
 *        scintillation light onto the readout plane.
 * 
 * @author H. Sullivan (hsulliva@fnal.gov)
 */

#ifndef AMSELSIM_PHOTONPROJECTIONALG_H
#define AMSELSIM_PHOTONPROJECTIONALG_H

#include "CLHEP/Random/JamesRandom.h"
#include "fhiclcpp/ParameterSet.h"

namespace amselsim
{

class PhotonProjectionAlg 
{
public:
  PhotonProjectionAlg(fhicl::ParameterSet const& p);

  void doProjection(CLHEP::HepRandomEngine& engine);
  void reconfigure(fhicl::ParameterSet const& p);
  void reset() { fPixelMap.clear(); fNScint=0; fNHits=0; }

  std::map<int, int> const& GetMap() const { return fPixelMap; }; 
  ULong64_t NumIncidentPhotons() const { return fNHits; };
  ULong64_t NumScintPhotons() const { return fNScint; };

private:
  ULong64_t          fNScint;
  ULong64_t          fNHits;
  std::map<int, int> fPixelMap;
};
}
#endif
