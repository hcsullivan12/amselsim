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
  {
    fProvider = std::make_unique<amselgeo::AmSelGeometry>(pset);
  }


void AmSelGeometryService::reconfigure(fhicl::ParameterSet const& pset)
{

}
}
DEFINE_ART_SERVICE_INTERFACE_IMPL(amselgeo::AmSelGeometryService, geo::DetectorGeometryService)

