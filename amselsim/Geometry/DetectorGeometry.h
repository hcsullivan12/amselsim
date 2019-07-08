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

      // If we wanted to completely refactor how geometry is handled,
      // the functions to define here would high level, e.g. NCryos,
      // NTPCs

      //virtual double Efield(unsigned int planegap=0) const = 0;

    protected:
      DetectorGeometry() = default;

    }; // class DetectorGeometry
} //namespace detinfo

#endif // AMSELGEO_DETECTORGEOMETRY_H
