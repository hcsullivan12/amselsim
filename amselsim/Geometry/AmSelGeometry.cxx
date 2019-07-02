
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
 : fGDML(pset.get< std::string >("GDML")),
   fTemperature(pset.get< double >("Temperature", 90.7)),
   fEfield(pset.get<double>("Efield", 0.5)),
   fElectronLifetime(pset.get<double>("ElectronLifetime")),
   fPixelSpacing(pset.get<float>("PixelSpacing")),
   fNPixels(0),
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
  fGDML = config.GDML;
  fPixelSpacing = config.PixelSpacing;
}

//--------------------------------------------------------------------
void AmSelGeometry::Initialize()
{
  TGeoManager::Import(fGDML.c_str());

  // Initialize detector properties
  TObjArray* volumes = gGeoManager->GetListOfVolumes();
  int nvols = volumes->GetEntries();
  ULong4_t maxId = 0;
  std::string volPixelPad("volPixelPad");
  std::string volPixelPlane("volPixelPlane");
  std::string volLArActive("volLArActive");
  for (int i = 0; i < nvols; i++)
  {
    TGeoVolume* vol = (TGeoVolume*)volumes->At(i);

    if (std::string(vol->GetName()).find(volPixelPlane) != std::string::npos)
    {
      std::cout << "\n\nHEEYYYY\n\n\n";
      auto tempVol = vol->GetNode(0)->GetVolume();
      auto b = ((TGeoBBox*)vol->GetShape())->GetOrigin();
      auto m = vol->GetNode(0)->GetMatrix();
      auto t = m->GetTranslation();
      Double_t d[3];
      m->LocalToMaster(t, d);
      std::cout << t[0] << " " << t[1] << " " << t[2] << std::endl;
      std::cout << d[0] << " " << d[1] << " " << d[2] << std::endl;
      std::cout << b[0] << " " << b[1] << " " << b[2] << std::endl;
    }


    if (std::string(vol->GetName()).find(volPixelPlane) != std::string::npos)
    {
      fPixelPlane = vol;
      TObjArray* nodes = fPixelPlane->GetNodes();
      fNPixels = nodes->GetEntries();
    }
    if (std::string(vol->GetName()).find(volLArActive) != std::string::npos)
    { 
      fDetHalfHeight = ((TGeoBBox*)vol->GetShape())->GetDX(); 
      fDetHalfWidth  = ((TGeoBBox*)vol->GetShape())->GetDY(); 
      fDetLength     = 2*((TGeoBBox*)vol->GetShape())->GetDZ(); 

      fLArTPCVolName = volLArActive;
    }
  }

  if (fLArTPCVolName.find(volLArActive) == std::string::npos) throw cet::exception("AmSelGeometry") << "Couldn't find LAr active volume!\n";

  if (!fPixelPlane) throw cet::exception("AmSelGeometry") << "Couldn't find pixel plane!\n";

  // Not sure about this, but let's do it anyway...
  //TObjArray* nodes = fPixelPlane->GetNodes();
  /*for (int i = 0; i < nodes->GetEntries(); i++)
  {
    TGeoMatrix* matrix = nodes->At(i)->GetMatrix();
    auto t = matrix->GetTranslation();
    std::cout << 
  }*/

  mf::LogInfo("AmSelGeometry") << "Initialized geomertry:"
                               << "\nNpixels = " << fNPixels;
}

//--------------------------------------------------------------------
double AmSelGeometry::Density(double temperature) const
{
  // Default temperature use internal value.
  if(temperature == 0.) temperature = fTemperature; //Temperature();
  double density = -0.00615*temperature + 1.928;
  return density;
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

//--------------------------------------------------------------------
double AmSelGeometry::DriftVelocity(double efield, double temperature) const 
{
  // Drift Velocity as a function of Electric Field and LAr Temperature
  // from : W. Walkowiak, NIM A 449 (2000) 288-294
  //
  // Efield should have units of kV/cm
  // Temperature should have units of Kelvin

  // Default Efield, use internal value.
  if(efield == 0.)
    efield = Efield();
  //
  if(efield > 4.0)
    mf::LogWarning("DetectorPropertiesAmSel") << "DriftVelocity Warning! : E-field value of "
				    << efield
				    << " kV/cm is outside of range covered by drift"
				    << " velocity parameterization. Returned value"
				    << " may not be correct";


  // Default temperature use internal value.
  if(temperature == 0.)
    temperature = Temperature();

  if(temperature < 87.0 || temperature > 94.0)
    mf::LogWarning("DetectorPropertiesAmSel") << "DriftVelocity Warning! : Temperature value of "
				    << temperature
				    << " K is outside of range covered by drift velocity"
				    << " parameterization. Returned value may not be"
				    << " correct";




  double tshift = -87.203+temperature;
  double xFit = 0.0938163-0.0052563*tshift-0.0001470*tshift*tshift;
  double uFit = 5.18406+0.01448*tshift-0.003497*tshift*tshift-0.000516*tshift*tshift*tshift;
  double vd;


// Icarus Parameter Set, use as default
  double  P1 = -0.04640; // K^-1
  double  P2 = 0.01712;  // K^-1
  double  P3 = 1.88125;   // (kV/cm)^-1
  double  P4 =  0.99408;    // kV/cm
  double  P5 =  0.01172;   // (kV/cm)^-P6
  double  P6 =  4.20214;
  double  T0 =  105.749;  // K
      // Walkowiak Parameter Set
  double    P1W = -0.01481; // K^-1
  double  P2W = -0.0075;  // K^-1
  double   P3W =  0.141;   // (kV/cm)^-1
  double   P4W =  12.4;    // kV/cm
  double   P5W =  1.627;   // (kV/cm)^-P6
  double   P6W =  0.317;
  double   T0W =  90.371;  // K

// From Craig Thorne . . . currently not documented
// smooth transition from linear at small fields to 
//     icarus fit at most fields to Walkowiak at very high fields
   if (efield < xFit) vd=efield*uFit;
   else if (efield<0.619) { 
     vd = ((P1*(temperature-T0)+1)
	       *(P3*efield*std::log(1+P4/efield) + P5*std::pow(efield,P6))
	       +P2*(temperature-T0));
   }
   else if (efield<0.699) {
     vd = 12.5*(efield-0.619)*((P1W*(temperature-T0W)+1)
	       *(P3W*efield*std::log(1+P4W/efield) + P5W*std::pow(efield,P6W))
	       +P2W*(temperature-T0W))+
       12.5*(0.699-efield)*((P1*(temperature-T0)+1)
	       *(P3*efield*std::log(1+P4/efield) + P5*std::pow(efield,P6))
	       +P2*(temperature-T0));
   }
   else {
     vd = ((P1W*(temperature-T0W)+1)
	       *(P3W*efield*std::log(1+P4W/efield) + P5W*std::pow(efield,P6W))
	       +P2W*(temperature-T0W));     
   }

  vd /= 10.;

  return vd; // in cm/us
}

}
