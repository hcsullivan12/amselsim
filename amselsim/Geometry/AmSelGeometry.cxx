/**
 * @file AmSelGeometry.cxx
 * @brief Interface to AmSel geometry information.
 * 
 * There is a question of how to handle pixels. For even small active
 * volumes, the number of pixels is 10s of 1000s. Therefore, there 
 * are two options. The user has the option to load a GDML file 
 * containing all pixel pads or to load a simplified geometry. A
 * few assumptions are made in the simplified geometry construction.
 *  
 * @author H. Sullivan (hsulliva@fnal.gov)
 */

#include "amselsim/Geometry/AmSelGeometry.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcorealg/CoreUtils/ProviderUtil.h" // lar::IgnorableProviderConfigKeys()

#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoVolume.h"

namespace amselgeo
{

AmSelGeometry::AmSelGeometry()
{}

//--------------------------------------------------------------------
AmSelGeometry::AmSelGeometry(fhicl::ParameterSet const& pset,
                             std::set<std::string> const& ignore_params)
 : fDetectorName("none"),
   fLArTPCVolName("none"),
   fOpDetVolName("none"),
   fNPixels(0),
   fPixelPlane(0),
   fNOpDets(0),
   fNCryo(1),
   fNTpc(1),
   fUseSimpleGeometry(false)
{
  ValidateAndConfigure(pset, ignore_params);
  Initialize();
}

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
  fDetectorName      = config.DetectorName();
  fGDMLPath          = config.GDML();
  fUseSimpleGeometry = config.UseSimpleGeometry(); 
  fPixelSpacing      = config.PixelSpacing();

  std::transform(fDetectorName.begin(), fDetectorName.end(),
      fDetectorName.begin(), ::tolower);
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
 
  // Force the gdml to have the optical and active volumes 
  if (fOpDetVolName.find("volOpDetSensitive") == std::string::npos) throw cet::exception("AmSelGeometry") << "Couldn't find optical detector volume!\n";
  if (fLArTPCVolName.find("volLArActive") == std::string::npos)     throw cet::exception("AmSelGeometry") << "Couldn't find LAr active volume!\n";
  if (!fPixelPlane)                                                 throw cet::exception("AmSelGeometry") << "Couldn't find pixel plane volume!\n";

  // Load simplified geometry
  if (fUseSimpleGeometry) 
  {
    std::cout <<"\n";
    mf::LogInfo("AmSelGeometry")<<"Loading simple geometry\n";
    LoadSimpleGeometry();
  }
  std::cout << "\n";
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
  if (volName == "volOpDetSensitive") 
  {
    fNOpDets++;
    fOpDetVolName = "volOpDetSensitive";
  }
  if (volName == "volPixelPlane")
  {
    fPixelPlane = nodeVol;
    // Try looking for pixels 
    if (fPixelPlane->GetNodes()) fNPixels = fPixelPlane->GetNodes()->GetEntries();
    // We're done
    return;
  }
  if (volName == "volLArActive")
  { 
    fDriftLength = 2*((TGeoBBox*)nodeVol->GetShape())->GetDX(); 
    fDetHalfY    =   ((TGeoBBox*)nodeVol->GetShape())->GetDY(); 
    fDetLength   = 2*((TGeoBBox*)nodeVol->GetShape())->GetDZ(); 
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
  fSimpleGeoY.clear();
  fSimpleGeoZ.clear();

  // Build a simplified pixelization scheme based on the pixel spacing
  // The pixel placements preserve the symmetry of the pixel plane
  float edge = 0.1*fPixelSpacing;
  float currentZ = 0.5 * fPixelSpacing;
  while ( (currentZ+0.5*fPixelSpacing) < 0.5 * (fDetLength-edge)) 
  {
    fSimpleGeoZ.push_back(currentZ);
    fSimpleGeoZ.push_back(-1 * currentZ);
    currentZ += fPixelSpacing;
  }
  float currentY = 0.5 * fPixelSpacing;

  while ( (currentY+0.5*fPixelSpacing) < (fDetHalfY-edge)) 
  {
    fSimpleGeoY.push_back(currentY);
    fSimpleGeoY.push_back(-1 * currentY);
    currentY += fPixelSpacing;
  }

  fNPixels = fSimpleGeoY.size() * fSimpleGeoZ.size();

  // Sorting so that top left hand corner is pixel 1
  std::sort(fSimpleGeoZ.begin(), fSimpleGeoZ.end(), [](float const& l, float const& r) {return l<r;});
  std::sort(fSimpleGeoY.begin(), fSimpleGeoY.end(), [](float const& l, float const& r) {return l>r;});

  // Convert to world coordinates
  std::string planeNodePath = *std::find_if(fNodePaths.begin(), fNodePaths.end(), [](std::string const& s) {return s.find("PixelPlane") != std::string::npos;}); 
  gGeoManager->cd(planeNodePath.c_str());
  TGeoNode* node = gGeoManager->GetCurrentNode();
  Double_t transl[3];
  auto o = ((TGeoBBox*)node->GetVolume()->GetShape())->GetOrigin();
  gGeoManager->LocalToMaster(o,transl);

  for (auto& z : fSimpleGeoZ) z+=transl[2];
  for (auto& y : fSimpleGeoY) y+=transl[1];

  gGeoManager->cd(fNodePaths[0].c_str());
}

//--------------------------------------------------------------------
void AmSelGeometry::GetOpDetCenter(double* xyz) const
{
  std::string planeNodePath = *std::find_if(fNodePaths.begin(), fNodePaths.end(), [](std::string const& s) {return s.find("PixelPlane") != std::string::npos;}); 
  gGeoManager->cd(planeNodePath.c_str());
  TGeoNode* node = gGeoManager->GetCurrentNode();
  auto o = ((TGeoBBox*)node->GetVolume()->GetShape())->GetOrigin();
  gGeoManager->LocalToMaster(o,xyz);
  gGeoManager->cd(fNodePaths[0].c_str()); 
}

//--------------------------------------------------------------------
int AmSelGeometry::NearestPixelID(geo::Point_t const& point) const
{
  // Check for simplified geometry
  if (fUseSimpleGeometry) return FindSimpleID(point);

  std::string pixelName = std::string(VolumeName(point));
  if (pixelName.find("volPixelPad") == std::string::npos) return -1;

  // Get the pixel ID
  size_t iD(0);
  for (; iD < pixelName.size(); iD++) {if(std::isdigit(pixelName[iD])) break;}

  return std::stoi(pixelName.substr(iD));
}

//--------------------------------------------------------------------
//
// \warning Assumes containers are already sorted in increasing Z and 
//          decreasing y.
//
int AmSelGeometry::FindSimpleID(geo::Point_t const& point) const
{
  // Find the first column closest to our test point
  float testZ        = point.z();
  auto  highZiter    = std::find_if(fSimpleGeoZ.begin(), fSimpleGeoZ.end(), [testZ](float const& z) {return z > testZ;});
  auto  lowZiter     = (highZiter-1) != fSimpleGeoZ.begin() ? (highZiter-1) : highZiter;
  auto  closestZiter = std::abs(*highZiter-testZ) < std::abs(*lowZiter-testZ) ? highZiter : lowZiter;
  size_t column      = std::distance(fSimpleGeoZ.begin(), closestZiter)+1; // 1,2,3...

  // Find the first row closest to our test point
  float testY        = point.y();
  auto  lowYiter     = std::find_if(fSimpleGeoY.begin(), fSimpleGeoY.end(), [testY](float const& y) {return y < testY;});
  auto  highYiter    = (lowYiter-1) != fSimpleGeoY.begin() ? (lowYiter-1) : lowYiter;
  auto  closestYiter = std::abs(*highYiter-testY) < std::abs(*lowYiter-testY) ? highYiter : lowYiter;
  size_t row         = std::distance(fSimpleGeoY.begin(), closestYiter)+1; // 1,2,3...

  //
  // @todo Handle this better
  //

  // If the distance is greater than some threshold, return -1
  double diffZ = testZ - *closestZiter;
  double diffY = testY - *closestYiter;
  if (std::sqrt(diffZ*diffZ+diffY*diffY) > 2*fPixelSpacing) return -1;

  return ((row-1) * fSimpleGeoZ.size() + column) - 1; // 0,1,2...
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
