//////////////////////////////////////////////////////////////////////
/// \file  AmSelGeometry.cxx
/// \brief Interface to AmSel geometry information.
///
/// \author  hsulliva@fnal.gov
//////////////////////////////////////////////////////////////////////

#include "amselsim/Geometry/AmSelGeometry.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcorealg/CoreUtils/ProviderUtil.h" // lar::IgnorableProviderConfigKeys()

#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoVolume.h"

namespace amselgeo
{

//--------------------------------------------------------------------
AmSelGeometry::AmSelGeometry(fhicl::ParameterSet const& pset,
                             std::set<std::string> const& ignore_params)
 : fNPixels(0),
   fPixelPlane(0)
{
  ValidateAndConfigure(pset, ignore_params);
  Initialize();
}

//--------------------------------------------------------------------
AmSelGeometry::~AmSelGeometry()
{}

//--------------------------------------------------------------------
void AmSelGeometry::ValidateAndConfigure(
    fhicl::ParameterSet const& p,
    std::set<std::string> const& ignore_params /* = {} */) 
{
  Configure(ValidateConfiguration(p, ignore_params));
}

//--------------------------------------------------------------------
AmSelGeometry::Configuration_t
AmSelGeometry::ValidateConfiguration(
    fhicl::ParameterSet const& p,
    std::set<std::string> const& ignore_params /* = {} */) 
{
  std::set<std::string> ignorable_keys = lar::IgnorableProviderConfigKeys();
  ignorable_keys.insert(ignore_params.begin(), ignore_params.end());
  // parses and validates the parameter set:
  fhicl::Table<Configuration_t> config_table { p, ignorable_keys };
  return std::move(config_table());
}

//--------------------------------------------------------------------
void AmSelGeometry::Configure(Configuration_t const& config) 
{
  fGDMLPath = config.GDML();
  fPixelSpacing = config.PixelSpacing();
}

//--------------------------------------------------------------------
void AmSelGeometry::Initialize()
{
  // We first need to validate the GDML file path
  cet::search_path sp("FW_SEARCH_PATH");
  std::string GDMLFilePath;
  if( !sp.find_file(fGDMLPath, GDMLFilePath) ) 
  {
    throw cet::exception("AmSelGeometry")
      << "Can't find geometry file '" << fGDMLPath << "'!\n";
  }

  // Reset the gdml path which now contains the full path
  fGDMLPath = GDMLFilePath;

  // Load the geometry from the gdml file
  if (gGeoManager) TGeoManager::UnlockGeometry();
  TGeoManager::Import(fGDMLPath.c_str());
  gGeoManager->LockGeometry();

  // Initialize high level information
  // Let's loop over our volumes, since we could potentially have 
  // a large number of pixels
  TObjArray* volumes = gGeoManager->GetListOfVolumes();
  size_t nVols = volumes->GetEntries();
  ULong4_t maxId = 0;
  for (size_t iVol = 0; iVol < nVols; iVol++)
  {
    TGeoVolume* vol = (TGeoVolume*)volumes->At(iVol);

    if (std::string(vol->GetName()) == "volPixelPlane")
    {
      fPixelPlane = vol;
      TObjArray* nodes = fPixelPlane->GetNodes();
      fNPixels = nodes->GetEntries();
    }
    if (std::string(vol->GetName()) == "volLArActive")
    { 
      fDetHalfHeight = ((TGeoBBox*)vol->GetShape())->GetDX(); 
      fDetHalfWidth  = ((TGeoBBox*)vol->GetShape())->GetDY(); 
      fDetLength     = 2*((TGeoBBox*)vol->GetShape())->GetDZ(); 

      fLArTPCVolName = "volLArActive";
    }
  }
  if (fLArTPCVolName.find("volLArActive") == std::string::npos) throw cet::exception("AmSelGeometry") << "Couldn't find LAr active volume!\n";
  if (!fPixelPlane)                                             throw cet::exception("AmSelGeometry") << "Couldn't find pixel plane volume!\n";


  // Now let's load our pixels
  float zMin = std::numeric_limits<float>::max();
  float zMax = std::numeric_limits<float>::min();
  float yMin = std::numeric_limits<float>::max();
  float yMax = std::numeric_limits<float>::min();
  for (ULong8_t i = 0; i < fNPixels; i++)
  {
    TGeoNode*   pixelNode = fPixelPlane->GetNode(i);
    TGeoVolume  pixelVol  = pixelNode->GetVolume();
    std::string pixelName = std::string(pixelVol->GetName()); 

    // Get the pixel ID
    size_t iD(0);
    for (; iD < pixelName.size(); iD++) {if(std::isdigit(pixelName[iD])) break;}
    size_t pixelID = std::stoi(pixelName.substr(iD));
    std::cout << pixelName << " " << pixelID << std::endl;


    auto o = ((TGeoBBox*)pixelNode->GetVolume()->GetShape())->GetOrigin();
    Double_t m[3];
    pixelNode->LocalToMaster(o,m);

    yMin = m[1] < yMin ? m[1] : yMin;
    yMax = m[1] > yMax ? m[1] : yMax;

    zMin = m[2] < zMin ? m[2] : zMin;
    zMax = m[2] > zMax ? m[2] : zMax;
  }


  // For now, add the half length to these to convert to world coordinates
  fPixelLimitsY.push_back(yMin); fPixelLimitsY.push_back(yMax);
  fPixelLimitsZ.push_back(zMin+fDetLength*0.5); fPixelLimitsZ.push_back(zMax+fDetLength*0.5);

  mf::LogInfo("AmSelGeometry")<<"Initialized geometry:"
                              <<"\nNpixels = "<<fNPixels
                              <<"\nz limits = ["<<fPixelLimitsZ[0]<<", "<<fPixelLimitsZ[1]<<"]"
                              <<"\ny limits = ["<<fPixelLimitsY[0]<<", "<<fPixelLimitsY[1]<<"]";
}

//--------------------------------------------------------------------
ULong8_t AmSelGeometry::NearestPixelID(const std::vector<double>& point) const
{
  return 1;
}

//--------------------------------------------------------------------
std::string AmSelGeometry::VolumeName(geo::Point_t const& point) const
{
  // check that the given point is in the World volume at least
  TGeoVolume const*volWorld = gGeoManager->FindVolumeFast("volWorld");
  double halflength = ((TGeoBBox*)volWorld->GetShape())->GetDZ();
  double halfheight = ((TGeoBBox*)volWorld->GetShape())->GetDY();
  double halfwidth  = ((TGeoBBox*)volWorld->GetShape())->GetDX();
  if(std::abs(point.x()) > halfwidth  ||
     std::abs(point.y()) > halfheight ||
     std::abs(point.z()) > halflength
     ){
    mf::LogWarning("GeometryCoreBadInputPoint") << "point (" << point.x() << ","
                                              << point.y() << "," << point.z() << ") "
                                              << "is not inside the world volume "
                                              << " half width = " << halfwidth
                                              << " half height = " << halfheight
                                              << " half length = " << halflength
                                              << " returning unknown volume name";
      const std::string unknown("unknownVolume");
      return unknown;
  }

  return gGeoManager->FindNode(point.X(), point.Y(), point.Z())->GetName();
}

}
