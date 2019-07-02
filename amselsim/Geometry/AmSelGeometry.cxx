
#include "AmSelGeometry.h"

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
  fGDML = config.GDML();
  fPixelSpacing = config.PixelSpacing();
}

//--------------------------------------------------------------------
void AmSelGeometry::Initialize()
{
  // We first need to validate the GDML file path
  cet::search_path sp("FW_SEARCH_PATH");
  std::string GDMLFilePath;
  if( !sp.find_file(fGDML, GDMLFilePath) ) 
  {
    throw cet::exception("AmSelGeometry")
      << "Can't find geometry file '" << fGDML << "'!\n";
  }

  // Reset the gdml path which now contains the full path
  fGDML = GDMLFilePath;

  // Load the geometry from the gdml file
  if (gGeoManager) TGeoManager::UnlockGeometry();
  TGeoManager::Import(fGDMLPath.c_str());
  gGeoManager->LockGeometry();

  // Initialize high level information
  // Let's loop over our volumes, since we could potentially have 
  // a large number of pixels
  TObjArray* volumes = gGeoManager->GetListOfVolumes();
  int nVols = volumes->GetEntries();
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

      fLArTPCVolName = volLArActive;
    }
  }
  if (fLArTPCVolName.find("volLArActive") == std::string::npos) throw cet::exception("AmSelGeometry") << "Couldn't find LAr active volume!\n";
  if (!fPixelPlane)                                             throw cet::exception("AmSelGeometry") << "Couldn't find pixel plane volume!\n";


  // Now let's try to load our pixels...
  TObjArray* pixelNodes = fPixelPlane->GetNodes();
  for (int i = 0; i < fNPixels; i++)
  {
    auto o = nodes->At(i)->GetShape()->GetOrigin();
    cout << o[0] << " " << o[1] << " " << o[2] << "\n";
    //TGeoMatrix* matrix = nodes->At(i)->GetMatrix();
    //auto t = matrix->GetTranslation();
    //std::cout << 
  }

  mf::LogInfo("AmSelGeometry") << "Initialized geometry:"
                               << "\nNpixels = " << fNPixels;
}

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
