////////////////////////////////////////////////////////////////////////
//
//  \file DetectorPropertiesAmSel.cxx
//
// Separation of service from Detector info class:
// hsulliva@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#include "amselsim/Services/DetectorPropertiesAmSel.h"

#include <cassert>

// LArSoft includes
#include "larcorealg/CoreUtils/ProviderUtil.h" // lar::IgnorableProviderConfigKeys()
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"

#include "nutools/IFDatabase/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Art includes
#include "fhiclcpp/make_ParameterSet.h"

namespace {
  
  template <typename T>
  inline T sqr(T v) { return v*v; }
  
} // local namespace

namespace ldp{

  //--------------------------------------------------------------------
  DetectorPropertiesAmSel::DetectorPropertiesAmSel() :
    fLP(0), fClocks(0), fGeo(0)
  {

  }
  
  //--------------------------------------------------------------------
  DetectorPropertiesAmSel::DetectorPropertiesAmSel(fhicl::ParameterSet const& pset,
					 const amselgeo::AmSelGeometry* geo,
					 const detinfo::LArProperties* lp,
					 const detinfo::DetectorClocks* c,
					 std::set<std::string> ignore_params /* = {} */
					 ):
    fLP(lp), fClocks(c), fGeo(geo)
  {
    fTPCClock = fClocks->TPCClock();
    
    ValidateAndConfigure(pset, ignore_params);
    
    // initialize prev run number
    fPrevRunNumber = 0;
    fCachedElectronLifetimes.reserve(10000);
    
  }
    
  //--------------------------------------------------------------------
  DetectorPropertiesAmSel::DetectorPropertiesAmSel(fhicl::ParameterSet const& pset,
					 providers_type providers,
					 std::set<std::string> ignore_params /* = {} */
					 ):
    DetectorPropertiesAmSel(pset,
      providers.get<amselgeo::AmSelGeometry>(),
      providers.get<detinfo::LArProperties>(),
      providers.get<detinfo::DetectorClocks>(),
      ignore_params
      )
    {}
  
  //--------------------------------------------------------------------
  bool DetectorPropertiesAmSel::Update(uint64_t t) 
  {

    bool retVal = true;
   
    // if the run number changed, we will need to update the electron lifetime
    if (fGetElectronlifetimeFromDB && t != fPrevRunNumber) {
    
      // let's avoid DB queries if we can (which can take up to 30 seconds on 
      // bad days!) by checking cached electron lifetimes to see if we've already 
      // found this run's lifetime value
      fElectronlifetime = 0;
      for(size_t i=0; i<fCachedElectronLifetimes.size(); i++){
        if( fCachedElectronLifetimes[i].first == t ) {
          fElectronlifetime = fCachedElectronLifetimes[i].second;
          break;
        }
      }
      if( fElectronlifetime == 0 ){
        retVal = UpdateElectronLifetime(t);
        fCachedElectronLifetimes.push_back(std::make_pair(t,fElectronlifetime));
      }
      
      fPrevRunNumber = t;
    }
    
    std::cout<<"Using electron lifetime "<<fElectronlifetime<<" microseconds.\n";
    
    return retVal;
  }

  //--------------------------------------------------------------------
  bool DetectorPropertiesAmSel::UpdateElectronLifetime(uint64_t t) 
  {

    std::string tableName = "elifetime";
    nutools::dbi::Table tbl;
    
    tbl.SetDetector("AmSel");
    tbl.SetTableName(tableName);
    tbl.SetTableType(nutools::dbi::kConditionsTable);
    tbl.SetDataTypeMask(nutools::dbi::kDataOnly);
    if( fElectronlifetimeTag != "" ) tbl.SetTag(fElectronlifetimeTag);

    int ltIdx = tbl.AddCol("lt","float");
    //int ltErrPlusIdx = tbl.AddCol("ltsigplus","float");
    //int ltErrMinusIdx = tbl.AddCol("ltsigminus","float");
    
    tbl.SetMinTSVld(t);
    tbl.SetMaxTSVld(t);

    tbl.SetVerbosity(100);
    tbl.Load();
    
    if (tbl.NRow() == 0) {
      std::cout << "No lifetime found from database!  Defaulting to nominal fhicl setting.";
      return false;
    }
    if (tbl.NRow() > 1) {
      std::cout << "More than one lifetime found from database!  This should NEVER happen, aborting.";
      abort();
    }

    nutools::dbi::Row* row;
    float lt;
    row = tbl.GetRow(0);
    row->Col(ltIdx).Get(lt);
    //row->Col(ltErrPlusIdx).Get(lterr1);
    //row->Col(ltErrMinusIdx).Get(lterr2);
    std::cout << "Setting electron lifetime to database value of " << lt << std::endl;

    fElectronlifetime = lt;
    
    return true;
  }
  
  //--------------------------------------------------------------------
  bool DetectorPropertiesAmSel::UpdateClocks(const detinfo::DetectorClocks* clks) 
  {
    fClocks = clks;
    
    fTPCClock = fClocks->TPCClock();
    return true;
  }
  
  //------------------------------------------------------------
  double DetectorPropertiesAmSel::ConvertTDCToTicks(double tdc) const
  {
    return fClocks->TPCTDC2Tick(tdc);
  }
  
  //--------------------------------------------------------------
  double DetectorPropertiesAmSel::ConvertTicksToTDC(double ticks) const
  {
    return fClocks->TPCTick2TDC(ticks);
  }
  
  //--------------------------------------------------------------------
  void DetectorPropertiesAmSel::Configure(Configuration_t const& config) {
   
    fEfield                     = config.Efield();
    fElectronlifetime           = config.Electronlifetime();
    fElectronlifetimeTag        = config.ElectronlifetimeTag();
    fGetElectronlifetimeFromDB  = config.GetElectronlifetimeFromDB();
    fTemperature                = config.Temperature();
    fElectronsToADC             = config.ElectronsToADC();
    fNumberTimeSamples          = config.NumberTimeSamples();
    fReadOutWindowSize          = config.ReadOutWindowSize();
    fTimeOffsetU                = config.TimeOffsetU();
    fTimeOffsetV                = config.TimeOffsetV();
    fTimeOffsetZ                = config.TimeOffsetZ();
    
    fSternheimerParameters.a    = config.SternheimerA();
    fSternheimerParameters.k    = config.SternheimerK();
    fSternheimerParameters.x0   = config.SternheimerX0();
    fSternheimerParameters.x1   = config.SternheimerX1();
    fSternheimerParameters.cbar = config.SternheimerCbar();
    
    fSamplingRate               = config.SamplingRate();
    if( fSamplingRate <= 0 )    fSamplingRate = fTPCClock.TickPeriod() * 1.e3;
  } // DetectorPropertiesAmSel::Configure()
  
  //--------------------------------------------------------------------
  DetectorPropertiesAmSel::Configuration_t
  DetectorPropertiesAmSel::ValidateConfiguration(
    fhicl::ParameterSet const& p, std::set<std::string> ignore_params /* = {} */
  ) {
    
    std::set<std::string> ignorable_keys = lar::IgnorableProviderConfigKeys();
    ignorable_keys.insert(ignore_params.begin(), ignore_params.end());
    
    // parses and validates the parameter set:
    fhicl::Table<Configuration_t> config_table { p, ignorable_keys };
    
    return std::move(config_table());
    
  } // DetectorPropertiesAmSel::ValidateConfiguration()
  
  //--------------------------------------------------------------------
  void DetectorPropertiesAmSel::ValidateAndConfigure(
    fhicl::ParameterSet const& p, std::set<std::string> ignore_params /* = {} */
  ) {
    Configure(ValidateConfiguration(p, ignore_params));
  } // ValidateAndConfigure()
  
  
  //------------------------------------------------------------------------------------//
  void DetectorPropertiesAmSel::Setup(providers_type providers) {
    
    SetGeometry(providers.get<amselgeo::AmSelGeometry>());
    SetLArProperties(providers.get<detinfo::LArProperties>());
    SetDetectorClocks(providers.get<detinfo::DetectorClocks>());
    
  } // DetectorPropertiesAmSel::Setup()
  
  
  //------------------------------------------------------------------------------------//
  double DetectorPropertiesAmSel::Efield(unsigned int planegap) const
  {
    if(planegap >= fEfield.size())
      throw cet::exception("DetectorPropertiesAmSel") << "requesting Electric field in a plane gap that is not defined\n";
    
    return fEfield[planegap];
  }
    
  //------------------------------------------------
  double DetectorPropertiesAmSel::Density(double temperature) const
  {
    // Default temperature use internal value.
    if(temperature == 0.)
      temperature = Temperature();
  
    double density = -0.00615*temperature + 1.928;
  
    return density;
  } // DetectorPropertiesAmSel::Density()
  
  
  //----------------------------------------------------------------------------------
  // Restricted mean energy loss (dE/dx) in units of MeV/cm.
  //
  // For unrestricted mean energy loss, set tcut = 0, or tcut large.
  //
  // Arguments:
  //
  // mom  - Momentum of incident particle in GeV/c.
  // mass - Mass of incident particle in GeV/c^2.
  // tcut - Maximum kinetic energy of delta rays (MeV).
  //
  // Returned value is positive.
  //
  // Based on Bethe-Bloch formula as contained in particle data book.
  // Material parameters (stored in larproperties.fcl) are taken from
  // pdg web site http://pdg.lbl.gov/AtomicNuclearProperties/
  //
  double DetectorPropertiesAmSel::Eloss(double mom, double mass, double tcut) const
  {
    // Some constants.
  
    double K = 0.307075;     // 4 pi N_A r_e^2 m_e c^2 (MeV cm^2/mol).
    double me = 0.510998918; // Electron mass (MeV/c^2).
  
    // Calculate kinematic quantities.
  
    double bg = mom / mass;           // beta*gamma.
    double gamma = sqrt(1. + bg*bg);  // gamma.
    double beta = bg / gamma;         // beta (velocity).
    double mer = 0.001 * me / mass;   // electron mass / mass of incident particle.
    double tmax = 2.*me* bg*bg / (1. + 2.*gamma*mer + mer*mer);  // Maximum delta ray energy (MeV).
  
    // Make sure tcut does not exceed tmax.
  
    if(tcut == 0. || tcut > tmax)
      tcut = tmax;
  
    // Calculate density effect correction (delta).
  
    double x = std::log10(bg);
    double delta = 0.;
    if(x >= fSternheimerParameters.x0) {
      delta = 2. * std::log(10.) * x - fSternheimerParameters.cbar;
      if(x < fSternheimerParameters.x1)
        delta += fSternheimerParameters.a * std::pow(fSternheimerParameters.x1 - x, fSternheimerParameters.k);
    }
  
    // Calculate stopping number.
  
    double B = 0.5 * std::log(2.*me*bg*bg*tcut / (1.e-12 * sqr(fLP->ExcitationEnergy())))
      - 0.5 * beta*beta * (1. + tcut / tmax) - 0.5 * delta;
  
    // Don't let the stopping number become negative.
  
    if(B < 1.)
      B = 1.;
  
    // Calculate dE/dx.
  
    double dedx = Density() * K*fLP->AtomicNumber()*B / (fLP->AtomicMass() * beta*beta);
  
    // Done.
  
    return dedx;
  } // DetectorPropertiesAmSel::Eloss()
  
  //----------------------------------------------------------------------------------
  double DetectorPropertiesAmSel::ElossVar(double mom, double mass) const
  {
    // Some constants.
  
    double K = 0.307075;     // 4 pi N_A r_e^2 m_e c^2 (MeV cm^2/mol).
    double me = 0.510998918; // Electron mass (MeV/c^2).
  
    // Calculate kinematic quantities.
  
    double bg = mom / mass;          // beta*gamma.
    double gamma2 = 1. + bg*bg;      // gamma^2.
    double beta2 = bg*bg / gamma2;   // beta^2.
  
    // Calculate final result.
  
    double result = gamma2 * (1. - 0.5 * beta2) * me * (fLP->AtomicNumber() / fLP->AtomicMass()) * K * Density();
    return result;
  } // DetectorPropertiesAmSel::ElossVar()


  //------------------------------------------------------------------------------------//
  double DetectorPropertiesAmSel::DriftVelocity(double efield, double temperature) const 
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

  //----------------------------------------------------------------------------------
  // The below function assumes that the user has applied the lifetime correction and
  // effective pitch between the wires (usually after 3D reconstruction). Using with
  // mean wire pitch will not give correct results.
  // parameters:
  //  dQdX in electrons/cm, charge (amplitude or integral obtained) divided by
  //         effective pitch for a given 3D track.
  // returns dEdX in MeV/cm
  double DetectorPropertiesAmSel::BirksCorrection(double dQdx) const
  {
    // Correction for charge quenching using parameterization from
    // S.Amoruso et al., NIM A 523 (2004) 275
    
    double  A3t    = util::kRecombA;
    double  K3t    = util::kRecombk;                     // in KV/cm*(g/cm^2)/MeV
    double  rho    = Density();                    // LAr density in g/cm^3
    double Wion    = 1000./util::kGeVToElectrons;        // 23.6 eV = 1e, Wion in MeV/e
    double E_field  = Efield();                           // Electric Field in the drift region in KV/cm
    K3t           /= rho;                                // KV/MeV
    double dEdx    = dQdx/(A3t/Wion-K3t/E_field*dQdx);    //MeV/cm
    
    return dEdx;
  }  
  
  //----------------------------------------------------------------------------------
  // Modified Box model correction 
  double DetectorPropertiesAmSel::ModBoxCorrection(double dQdx) const
  {
    // Modified Box model correction has better behavior than the Birks
    // correction at high values of dQ/dx.
    double  rho    = Density();                    // LAr density in g/cm^3
    double Wion    = 1000./util::kGeVToElectrons;        // 23.6 eV = 1e, Wion in MeV/e
    double E_field  = Efield();                           // Electric Field in the drift region in KV/cm
    double Beta    = util::kModBoxB / (rho * E_field);
    double Alpha   = util::kModBoxA;
    double dEdx = (exp(Beta * Wion * dQdx ) - Alpha) / Beta;
    
    return dEdx;
  }
  
  //------------------------------------------------------------------------------------//
  int  DetectorPropertiesAmSel::TriggerOffset()     const 
  {
    return fTPCClock.Ticks(fClocks->TriggerOffsetTPC() * -1.);
  }
  
  
  //--------------------------------------------------------------------
  //  x<--> ticks conversion methods 
  //
  //  Ben Jones April 2012, 
  //  based on code by Herb Greenlee in SpacePointService
  //  
  
  


  //--------------------------------------------------------------------
  // Take an X coordinate, and convert to a number of ticks, the
  // charge deposit occured at t=0
 
  double DetectorPropertiesAmSel::ConvertXToTicks(double X, int p, int t, int c) const
  {
    return (X / (fXTicksCoefficient * fDriftDirection.at(c).at(t)) +  fXTicksOffsets.at(c).at(t).at(p) );
  }
  


  //-------------------------------------------------------------------
  // Take a cooridnate in ticks, and convert to an x position
  // assuming event deposit occured at t=0
 
  double  DetectorPropertiesAmSel::ConvertTicksToX(double ticks, int p, int t, int c) const
  {
    return (ticks - fXTicksOffsets.at(c).at(t).at(p)) * fXTicksCoefficient * fDriftDirection.at(c).at(t);  
  }
  

  //--------------------------------------------------------------------
  void DetectorPropertiesAmSel::CheckIfConfigured() const
  {
    if (!fGeo) throw cet::exception(__FUNCTION__) << "Geometry is uninitialized!";
    if (!fLP) throw cet::exception(__FUNCTION__) << "LArPropertiesAmSel is uninitialized!";
    if (!fClocks) throw cet::exception(__FUNCTION__) << "DetectorClocks is uninitialized!";
  }

  //--------------------------------------------------------------------
  double DetectorPropertiesAmSel::GetXTicksOffset(int p, int t, int c) const
  {
    return fXTicksOffsets.at(c).at(t).at(p);
  }

  //--------------------------------------------------------------------
  double DetectorPropertiesAmSel::GetXTicksCoefficient() const
  {
    return fXTicksCoefficient;
  }

  //--------------------------------------------------------------------
  double DetectorPropertiesAmSel::GetXTicksCoefficient(int t, int c) const
  {
    return fXTicksCoefficient * fDriftDirection.at(c).at(t);
  }

} // namespace
