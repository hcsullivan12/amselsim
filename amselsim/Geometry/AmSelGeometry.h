//////////////////////////////////////////////////////////////////////
/// \file  AmSelGeometry.h
/// \brief Interface to AmSel geometry information.
///
/// \author  hsulliva@fnal.gov
//////////////////////////////////////////////////////////////////////

#ifndef AMSELGEO_AMSELGEOMETRY_H
#define AMSELGEO_AMSELGEOMETRY_H

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h" // Point_t

// amselsim includes
#include "amselsim/Geometry/DetectorGeometry.h"

// C++ includes
#include <set>

namespace amselgeo
{
  using UShort2_t = unsigned short;
  using ULong4_t = unsigned long;
  using ULong8_t = unsigned long long;

/// Configuration parameter documentation goes here
class AmSelGeometry : public DetectorGeometry
{
  public:

    /// Structure for configuration parameters
    struct Configuration_t 
    {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<std::string> GDML         
      {
        Name("GDML"),
        Comment("Name of GDML file for AmSel geometry")
      };

      fhicl::Atom<double> PixelSpacing
      {
        Name("PixelSpacing"),
        Comment("Pixel spacing")
      };
    };

    /// Constructor: reads the configuration from a parameter set
    AmSelGeometry(fhicl::ParameterSet const& pset,
                  std::set<std::string> const& ignore_params = {});
    ~AmSelGeometry();

    /// Method to validate parameter set
    void ValidateAndConfigure(
        fhicl::ParameterSet const& p,
        std::set<std::string> const& ignore_params = {});

    Configuration_t ValidateConfiguration(
        fhicl::ParameterSet const& p,
        std::set<std::string> const& ignore_params = {});

    /// Extracts the relevant configuration from the specified object
    void Configure(Configuration_t const& config);        

    std::string GDMLFile() const { return fGDMLPath; }
    std::string OpDetGeoName() const { return "opDetector"; }
    size_t NOpDets() const { return 0; }
    size_t Ncryostats() const { return 1; }
    size_t NTPC() const { return 1; }
    size_t NAuxDets() const { return 0; }
    size_t NSensitiveVolume() const { return 1; }
    double DetHalfHeight() const { return fDetHalfHeight; }
    double DetHalfWidth() const  { return fDetHalfWidth; }
    double DetLength() const     { return fDetLength; }
    float  PixelSpacing() const  { return fPixelSpacing; }
    ULong8_t NPixels() const { return fNPixels; }
    ULong8_t NearestPixelID(const std::vector<double>& point) const;
    std::string GetLArTPCVolumeName() const { return fLArTPCVolName; }
    std::string VolumeName(geo::Point_t const& point) const;
    std::string VolumeName(TVector3 const& point) const
      { return VolumeName(geo::vect::toPoint(point)); }




  private:
    void Initialize();

    std::string fGDMLPath;
    std::string fLArTPCVolName;
    double fDetHalfHeight;
    double fDetHalfWidth;
    double fDetLength;
    float  fPixelSpacing;
    std::vector<float> fPixelLimitsY;
    std::vector<float> fPixelLimitsZ;
    ULong8_t fNPixels;
    TGeoVolume* fPixelPlane;
}; // class AmSelGeometry

}

#endif
