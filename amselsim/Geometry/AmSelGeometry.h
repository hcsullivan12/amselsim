
#ifndef AMSELGEO_AMSELGEOMETRY_H
#define AMSELGEO_AMSELGEOMETRY_H

#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"

namespace amselgeo
{

/// Configuration parameter documentation goes here
class AmSelGeometry 
{
  public:
    /// Constructor: reads the configuration from a parameter set
    AmSelGeometry(fhicl::ParameterSet const& pset);
    ~AmSelGeometry();

    std::string GDMLFile() const { return fGDML; }
    std::string OpDetGeoName() const { return "opDetector"; }
    size_t NOpDets() const { return 1; }
    size_t Ncryostats() const { return 1; }
    size_t NTPC() const { return 1; }
    size_t NAuxDets() const { return 0; }
    size_t NSensitiveVolume() const { return 1; }
    double Density(double temperature) const;
    double Efield() const { return fEfield; }
    double Temperature() const { return fTemperature; }
    double ElectronLifetime() const { return fElectronLifetime; }
    double DriftVelocity(double efield, double temperature) const;
    double DetHalfHeight() const { return fDetHalfHeight; }
    double DetHalfWidth() const  { return fDetHalfWidth; }
    double DetLength() const     { return fDetLength; }
    std::string GetLArTPCVolumeName() const { return fLArTPCVolName; }
    std::string VolumeName(geo::Point_t const& point) const;
    std::string VolumeName(TVector3 const& point) const
      { return VolumeName(geo::vect::toPoint(point)); }




  private:
    void Initialize();

    std::string fGDML;
    std::string fLArTPCVolName;
    double fTemperature;
    double fEfield;
    double fElectronLifetime;
    double fDetHalfHeight;
    double fDetHalfWidth;
    double fDetLength;


}; // class AmSelGeometry

}

#endif
