//////////////////////////////////////////////////////////////////////
/// \file  AmSelGeometry.cxx
/// \brief Interface to AmSel geometry information.
///
/// There is a question of how to handle pixels. For even small active
/// volumes, the number of pixels is 10s of 1000s. Therefore, there 
/// are two options allowed. The user has the option to load a GDML
/// file containing all pixel pads or a GDML file containing the minimum
/// information necessary to assume a pixelization scheme. The latter
/// is assumed to be a file '*simple.gdml'. 
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
   fPixelPlane(0),
   fIsSimpleGeometry(false)
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
  fGDMLPath     = config.GDML();

  // Look to see if this is GDML with simplified pixels
  if (fGDMLPath.find("simple") != std::string::npos) 
  {
    fIsSimpleGeometry = true; 
    fPixelSpacing = config.PixelSpacing();
  }
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
  // Create new path
  std::string path = topNode->GetName();
  fNodePaths.push_back(path);
  LookAtNode(topNode, path);
 
  if (fLArTPCVolName.find("volLArActive") == std::string::npos) throw cet::exception("AmSelGeometry") << "Couldn't find LAr active volume!\n";
  if (!fPixelPlane)                                             throw cet::exception("AmSelGeometry") << "Couldn't find pixel plane volume!\n";

  // Load simplified geometry
  if (fIsSimpleGeometry) LoadSimpleGeometry();

  for (const auto& np : fNodePaths) std::cout << np << std::endl;

  mf::LogInfo("AmSelGeometry")<<"Initialized geometry with " << fNPixels << " pixels\n";
}

//--------------------------------------------------------------------
void AmSelGeometry::LookAtNode(TGeoNode const* currentNode, std::string const& currentPath) 
{
  // Get the volume of this node
  TGeoVolume* nodeVol = currentNode->GetVolume();
  std::string volName = std::string(nodeVol->GetName());

  // Leave if this is a pixel
  if (volName.find("volPixelPad") != std::string::npos) return;
  if (volName == "volPixelPlane")
  {
    fPixelPlane = nodeVol; 
    fNPixels = fPixelPlane->GetNodes()->GetEntries();
  }
  if (volName == "volLArActive")
  { 
    fDetHalfHeight = ((TGeoBBox*)nodeVol->GetShape())->GetDX(); 
    fDetHalfWidth  = ((TGeoBBox*)nodeVol->GetShape())->GetDY(); 
    fDetLength     = 2*((TGeoBBox*)nodeVol->GetShape())->GetDZ(); 
  
    fLArTPCVolName = "volLArActive";
  }  

  // Check nodes
  TObjArray* nodes = currentNode->GetNodes();
  if (!nodes) return;
  for (int iN = 0; iN < nodes->GetEntries(); iN++) 
  {
    std::string nextNodePath = currentPath+"/"+nodeVol->GetNode(iN)->GetName();
    fNodePaths.push_back(nextNodePath);   
    LookAtNode(nodeVol->GetNode(iN), nextNodePath);
  }
}


//--------------------------------------------------------------------
void AmSelGeometry::LoadSimpleGeometry()
{
  // It is assumed that the gdml file contains the positions of the 
  // columns and rows for the pixels. The ID will be identified by 
  // the row/column combination.
  std::vector<float> y, z;
  TObjArray* pixelNodes = fPixelPlane->GetNodes();
  for (int iP = 0; iP < pixelNodes->GetEntries(); iP++)
  {
    TGeoNode* currentNode = fPixelPlane->GetNode(iP); 
  }
}

//--------------------------------------------------------------------
int AmSelGeometry::NearestPixelID(geo::Point_t const& point) const
{
  std::string pixelName = std::string(VolumeName(point));
  if (pixelName.find("volPixelPad") == std::string::npos) return -1;

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
