//////////////////////////////////////////////////////////////////////
/// \file  AmSelGeometryService.h
/// \brief Service for AmSel geometry information.
///
/// \author  hsulliva@fnal.gov
//////////////////////////////////////////////////////////////////////

#ifndef AMSELGEO_AMSELGEOMETRYSERVICE_H
#define AMSELGEO_AMSELGEOMETRYSERVICE_H

#include "amselsim/Geometry/AmSelGeometry.h"
#include "amselsim/Geometry/DetectorGeometryService.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace amselgeo 
{

class AmSelGeometryService : public DetectorGeometryService
{
  public:
    using provider_type = AmSelGeometry; ///< type of the service provider

    // Standard art service constructor
    AmSelGeometryService(fhicl::ParameterSet const&, art::ActivityRegistry&);

    /// Return a pointer to a (constant) detector properties provider
   // provider_type const* provider() const { return fProvider.get(); }
    virtual const provider_type* provider() const override { return fProvider.get();}
    
    virtual void   reconfigure(fhicl::ParameterSet const& pset) override;

 
  private:
    std::unique_ptr<AmSelGeometry> fProvider; ///< owned provider
    std::string fGDMLPath;

}; // class AmSelGeometryService

}
DECLARE_ART_SERVICE_INTERFACE_IMPL(amselgeo::AmSelGeometryService, amselgeo::DetectorGeometryService, LEGACY)

#endif
