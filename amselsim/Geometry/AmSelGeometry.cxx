//////////////////////////////////////////////////////////////////////
/// \file  AmSelGeometry.cxx
/// \brief Interface to AmSel geometry information.
///
/// \todo Change world_vol to volWorld 
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

  mf::LogInfo("AmSelGeometry")<<"Initialized geometry:"
                              <<"\nNpixels = "<<fNPixels;
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
ULong8_t AmSelGeometry::NearestPixelID(const std::vector<double>& point) const
{
  if (point.size() != 3) throw cet::exception("AmSelGeometry") << "Point must be 3D!\n";

  TGeoNode* pixelNode = gGeoManager->FindNode(point[0], point[1], point[2]);
  if (!pixelNode) {mf::LogWarning("AmSelGeometry") << "Pixel node is null! Return 1.\n"; return 1;}

  TGeoVolume* pixelVol  = pixelNode->GetVolume();
  std::string pixelName = std::string(pixelVol->GetName()); 
  // Get the pixel ID
  size_t iD(0);
  for (; iD < pixelName.size(); iD++) {if(std::isdigit(pixelName[iD])) break;}
  
  return std::stoi(pixelName.substr(iD));
  //Double_t m[3];
  //auto o = ((TGeoBBox*)pixelVol->GetShape())->GetOrigin();
  //pixelNode->LocalToMaster(o,m);
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
