////////////////////////////////////////////////////////////////////////
// \file DetectorGeometry.h
// \brief Pure virtual base interface for detector geometry
//
// \author hsulliva@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#ifndef AMSELGEO_DETECTORGEOMETRY_H
#define AMSELGEO_DETECTORGEOMETRY_H

#include "fhiclcpp/ParameterSet.h"

#include "TVector3.h"

#include <stdexcept> // std::runtime_error()


///General LArSoft Utilities
namespace amselgeo{

  class DetectorGeometry {
    public:

      DetectorGeometry(const DetectorGeometry &) = delete;
      DetectorGeometry(DetectorGeometry &&) = delete;
      DetectorGeometry& operator = (const DetectorGeometry &) = delete;
      DetectorGeometry& operator = (DetectorGeometry &&) = delete;
      virtual ~DetectorGeometry() = default;

      virtual double      DetHalfWidth() const = 0;
      virtual double      DetDriftLength() const = 0;
      virtual double      DetHalfHeight() const = 0;
      virtual double      DetLength()    const = 0;
      virtual std::string GetLArTPCVolumeName() const = 0;
      virtual std::string VolumeName(TVector3 const& point) const = 0;
      virtual std::string GDMLFile() const = 0;
      virtual std::string OpDetGeoName() const = 0;
      virtual size_t      Ncryostats() const = 0;
      virtual size_t      NTPC() const = 0;
      virtual size_t      NOpDets() const = 0;
      virtual size_t      NAuxDets() const = 0;
      virtual size_t      NSensitiveVolume() const = 0;
      virtual int         NReadoutNodes() const = 0;
      virtual int         NearestReadoutNodeID(TVector3 const& point) const = 0;

    protected:
      DetectorGeometry() = default;

    }; // class DetectorGeometry
} //namespace detinfo

#endif // AMSELGEO_DETECTORGEOMETRY_H
