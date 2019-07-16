////////////////////////////////////////////////////////////////////////
// \file DetectorGeometry.h
// \brief Pure virtual base interface for detector geometry
//
// \author hsulliva@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#ifndef GEO_DETECTORGEOMETRY_H
#define GEO_DETECTORGEOMETRY_H

#include "fhiclcpp/ParameterSet.h"

#include "TVector3.h"
#include "TGeoManager.h"

#include <stdexcept> // std::runtime_error()

///General LArSoft Utilities
namespace geo{

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
      virtual std::string ROOTFile() const = 0;
      virtual std::string OpDetGeoName() const = 0;
      virtual std::string DetectorName() const = 0;
      virtual size_t      Ncryostats() const = 0;
      virtual size_t      NTPC() const = 0;
      virtual size_t      NOpDets() const = 0;
      virtual size_t      NAuxDets() const = 0;
      virtual size_t      NSensitiveVolume() const = 0;
      virtual int         NReadoutNodes() const = 0;
      virtual size_t      Nchannels() const = 0;
      virtual int         NearestReadoutNodeID(TVector3 const& point) const = 0;
      virtual int         NearestChannel(double* xyz, 
                                         size_t const& p, 
                                         size_t const& tpc,
                                         size_t const& cryo) const = 0;
      virtual double      TotalMass(std::string const& vol) const = 0;
      virtual TGeoManager* ROOTGeoManager() const = 0; 

      virtual const std::vector<double> PlaneLocation(size_t const& p) const = 0;
      virtual int           DriftDirection() const = 0;
      virtual size_t        Nplanes() const = 0;
      virtual double        PlanePitch(size_t const& p1=0, size_t const& p2=1) const = 0;
      virtual void          GetOpDetCenter(double* xyz) const = 0;

    protected:
      DetectorGeometry() = default;

    }; // class DetectorGeometry
} //namespace detinfo

#endif // GEO_DETECTORGEOMETRY_H
