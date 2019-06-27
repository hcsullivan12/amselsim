

#ifndef AMSELGEO_AMSELGEOMETRYSERVICE_H
#define AMSELGEO_AMSELGEOMETRYSERVICE_H

#include "amselsim/Geometry/AmSelGeometry.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace amselgeo 
{

class AmSelGeometryService 
{
  public:
    using provider_type = AmSelGeometry; ///< type of the service provider

    // Standard art service constructor
    AmSelGeometryService(fhicl::ParameterSet const&, art::ActivityRegistry&);

    /// Return a pointer to a (constant) detector properties provider
    provider_type const* provider() const { return fProvider.get(); }
 
  private:
    std::unique_ptr<AmSelGeometry> fProvider; ///< owned provider
    std::string fGDML;

}; // class AmSelGeometryService

}
DECLARE_ART_SERVICE(amselgeo::AmSelGeometryService, LEGACY)

#endif
