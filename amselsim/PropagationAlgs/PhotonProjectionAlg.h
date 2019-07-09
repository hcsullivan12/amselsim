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

  std::map<int, int> const& GetMap() const { return fPixelMap; }; 

private:
  ULong64_t          fNScint;
  ULong64_t          fNHits;
  std::map<int, int> fPixelMap;
};
}
#endif
