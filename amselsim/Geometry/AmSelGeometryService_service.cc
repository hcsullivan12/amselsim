//////////////////////////////////////////////////////////////////////
/// \file  AmSelGeometryService_service.cxx
/// \brief Service for AmSel geometry information.
///
/// \author  hsulliva@fnal.gov
//////////////////////////////////////////////////////////////////////

#include "AmSelGeometryService.h"

namespace amselgeo 
{

AmSelGeometryService::AmSelGeometryService
  (fhicl::ParameterSet const& pset, art::ActivityRegistry&)
  : fProvider(new AmSelGeometry(pset))
  {}

}
DEFINE_ART_SERVICE(amselgeo::AmSelGeometryService)
