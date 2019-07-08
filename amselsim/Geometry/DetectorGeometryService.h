////////////////////////////////////////////////////////////////////////
// \file DetectorGeometryService.h
// \brief Pure virtual service interface for Geometry functions
//
// \author hsulliva@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETECTORGEOMETRYSERVICE_H
#define DETECTORGEOMETRYSERVICE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "amselsim/Geometry/DetectorGeometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace amselgeo{
  class DetectorGeometryService {

    public:
    typedef amselgeo::DetectorGeometry provider_type;

    public:
      virtual ~DetectorGeometryService() = default;

      virtual void   reconfigure(fhicl::ParameterSet const& pset) = 0;
      virtual const  amselgeo::DetectorGeometryService* provider() const = 0;

    }; 
} 
DECLARE_ART_SERVICE_INTERFACE(amselgeo::DetectorGeometryService, LEGACY)
#endif 
