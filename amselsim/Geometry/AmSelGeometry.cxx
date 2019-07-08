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

  TGeoNode* topNode = gGeoManager->GetTopNode();
  LookAtNode(topNode);

  if (fLArTPCVolName.find("volLArActive") == std::string::npos) throw cet::exception("AmSelGeometry") << "Couldn't find LAr active volume!\n";
  if (!fPixelPlane)                                             throw cet::exception("AmSelGeometry") << "Couldn't find pixel plane volume!\n";

  mf::LogInfo("AmSelGeometry")<<"Initialized geometry with " << fNPixels << " pixels\n";
}

//--------------------------------------------------------------------
void AmSelGeometry::LookAtNode(const TGeoNode* currentNode) 
{
  // Get the volume of this node
  TGeoVolume* nodeVol = currentNode->GetVolume();
  std::string volName = std::string(nodeVol->GetName());

  // Leave if this is a pixel
  if (volName.find("volPixelPad") != std::string::npos) return;
  if (volName == "volPixelPlane")
  {
    fPixelPlane = nodeVol;
    TObjArray* pixelNodes = fPixelPlane->GetNodes();
    fNPixels = pixelNodes->GetEntries();

    // We have everything we need
    return;
  }
  if (volName == "volLArActive")
  { 
    fDetHalfHeight = ((TGeoBBox*)nodeVol->GetShape())->GetDX(); 
    fDetHalfWidth  = ((TGeoBBox*)nodeVol->GetShape())->GetDY(); 
    fDetLength     = 2*((TGeoBBox*)nodeVol->GetShape())->GetDZ(); 
  
    fLArTPCVolName = "volLArActive";
  }  

  // Check the nodes
  TObjArray* nodes = currentNode->GetNodes();
  if (!nodes) return;
  for (int iN = 0; iN < nodes->GetEntries(); iN++) LookAtNode(nodeVol->GetNode(iN));
}


//--------------------------------------------------------------------
ULong8_t AmSelGeometry::NearestPixelID(geo::Point_t const& point) const
{
  std::string pixelName = std::string(VolumeName(point));
  if (pixelName.find("PixelPad") == std::string::npos) 
  {
    mf::LogWarning("AmSelGeometry") << "point (" << point.x() << ","
                                    << point.y() << "," << point.z() << ") "
                                    << "is not on a pixel pad. Returning "
                                    << " pixelID 1\n";
    return 1;
  }

  // Get the pixel ID
  size_t iD(0);
  for (; iD < pixelName.size(); iD++) {if(std::isdigit(pixelName[iD])) break;}

  return std::stoi(pixelName.substr(iD));
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
    mf::LogWarning("AmSelGeometry") << "point (" << point.x() << ","
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
