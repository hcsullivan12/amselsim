
#include "AmSelGeometryService.h"

namespace amselgeo 
{

AmSelGeometryService::AmSelGeometryService
  (fhicl::ParameterSet const& pset, art::ActivityRegistry&)
  : fProvider(new AmSelGeometry(pset))
  {}

}
DEFINE_ART_SERVICE(amselgeo::AmSelGeometryService)
