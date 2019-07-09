/**
 * @file PhotonProjectionAlgAlg.h
 * @brief Algorithm for emission and projection of 
 *        scintillation light onto the readout plane.
 * 
 * @author H. Sullivan (hsulliva@fnal.gov)
 */

#ifndef AMSELSIM_PHOTONPROJECTIONALG_H
#define AMSELSIM_PHOTONPROJECTIONALG_H

#include "fhiclcpp/ParameterSet.h"
#include "nurandom/RandomUtils/NuRandomService.h"

namespace amselg4 
{

class PhotonProjectionAlg 
{
public:
  PhotonProjectionAlg();

  void doProjection();
  void reconfigure(fhicl::ParameterSet const& p);

private:
  CLHEP::HepRandomEngine& fEngine;
  ULong64_t        fNScint;
  ULong64_t        fNHits;
};
}
#endif