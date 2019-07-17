////////////////////////////////////////////////////////////////////////
/// \file  IonizationAndScintillationAction.cxx
/// \brief Use Geant4's user "hooks" to determine the number of ionization
///        electrons and scintillation photons for each step
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "amselsim/LArG4/IonizationAndScintillationAction.h"
#include "amselsim/LArG4/IonizationAndScintillation.h"


#include "Geant4/G4Step.hh"

#include <algorithm>

namespace amselg4 {

  //----------------------------------------------------------------------------
  // Constructor.
  IonizationAndScintillationAction::IonizationAndScintillationAction()
  {
  }

  //----------------------------------------------------------------------------
  // Destructor.
  IonizationAndScintillationAction::~IonizationAndScintillationAction()
  {
  }

  //----------------------------------------------------------------------------
  // With every step, calculate the number of ionization electrons and
  // scintillation photons using the IonizationAndScintillation singleton.
  void IonizationAndScintillationAction::SteppingAction(const G4Step* step)
  {
    amselg4::IonizationAndScintillation::Instance()->Reset(step);

    return;
  }

} // namespace AmSelG4
