////////////////////////////////////////////////////////////////////////
// \file DetectorGeometryService.h
// \brief Pure virtual service interface for Geometry functions
//
// \author hsulliva@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef GEO_DETECTORGEOMETRYSERVICE_H
#define GEO_DETECTORGEOMETRYSERVICE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "amselsim/Geometry/DetectorGeometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace geo{
  class DetectorGeometryService {

    public:
    typedef geo::DetectorGeometry provider_type;

    public:
      virtual ~DetectorGeometryService() = default;

      virtual void   reconfigure(fhicl::ParameterSet const& pset) = 0;
      virtual const  geo::DetectorGeometry* provider() const = 0;

    }; 
} 
DECLARE_ART_SERVICE_INTERFACE(geo::DetectorGeometryService, LEGACY)
#endif 
