/**
 * @file AnaTree_module.cc
 * @brief Analyzer module for * 
 * @author H. Sullivan (hsulliva@fnal.gov)
 */

// Framework includes 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes 
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT includes 
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"

// amselsim includes
#inclued "amselsim/PropagationAlgs/PhotonProjectionAlg.h"

namespace namespace
{

class AnaTree : public art::EDAnalyzer 
{
public:
  explicit AnaTree(fhicl::ParameterSet const & p);
  virtual ~AnaTree();

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob();
  void reconfigure(fhicl::ParameterSet const & p);

private:

  // Function used to reset all the variables  
  void ResetVars();
  
  // Storing information into TTree
  TTree* fTree;
   
  // Storing Run Information 
  int    run;			    ///< Run Number
  int    subrun;			///< SubRun Number
  int    event;			  ///< Event Number
  double evttime;		  ///< Event Time Stamp
  double efield[3];		///< Electric Field 

  // Storing Geant4 MC Truth Information 
  int    no_primaries;			///< Number of primary Geant4 particles in the event
  int    geant_list_size;		///< Number of Geant4 particles tracked
  double primary_p;				  ///< Primary particle momentum
  std::vector<int>    PDG;
  std::vector<double> StartPointx; 
  std::vector<double> StartPointy;
  std::vector<double> StartPointz;
  std::vector<double> StartEnergy;
  std::vector<double> StartKE;
  std::vector<double> LastKE;
  std::vector<double> StartPx;
  std::vector<double> StartPy;
  std::vector<double> StartPz;
  std::vector<double> EndPointx; 
  std::vector<double> EndPointy;
  std::vector<double> EndPointz;
  std::vector<double> EndEnergy;
  std::vector<double> EndPx;
  std::vector<double> EndPy;
  std::vector<double> EndPz;
  std::vector<int>    InteractionPoint;         ///< Geant 4 Primary Trj Point Corresponding to the Interaction
  std::vector<int>    InteractionPointType;     ///< Geant 4 Primary Interaction Type

  std::vector<int> Process;
  std::vector<int> NumberDaughters;
  std::vector<int> TrackId;
  std::vector<int> Mother;
  std::vector<int> process_primary;	          ///< Is this particle primary (primary = 1, non-primary = 0)
  std::vector<std::string> G4Process;         ///< The process which created this particle
  std::vector<std::string> G4FinalProcess;    ///< The last process which this particle went under

  // Storing additional Geant4 MC Truth Information for the primary track only 	   
  std::vector<int> NTrTrajPts;
  std::vector< std::vector<double> > MidPosX;
  std::vector< std::vector<double> > MidPosY;
  std::vector< std::vector<double> > MidPosZ;
  std::vector< std::vector<double> > MidPx;
  std::vector< std::vector<double> > MidPy;
  std::vector< std::vector<double> > MidPz;

  // Storing additional Geant4 MC Truth Information for the daughter tracks 
  std::vector<double> NDTrTrajPts;
  std::vector<int>    DTrackId;
  std::vector<int>    DPdgCode;
  std::vector<double> DStartKE;
  std::vector<double> DStartEnergy;
  std::vector<double> DStartP;
  std::vector< std::vector<double> > DMidPosX;
  std::vector< std::vector<double> > DMidPosY;
  std::vector< std::vector<double> > DMidPosZ; 

  // Scintillation information
  // @todo Confirm nPrimScint
  ULong64_t nPrimScint;           ///< Number of scintillation photons produced by primary track 
  ULong64_t nReadoutIncPhotons;   ///< Number of scintillation photons incident on readout plane
  std::vector<int> pixelIdVec;    ///< Container of the pixel IDs which detected > 0 photons
  std::vector<int> pixelCountVec; ///< Container of the counts for each pixel ID

  std::string fTreeName;
  std::string fG4ModuleLabel;
  bool        fDoPhotonProjection;
  PhotonProjectionAlg fProjAlg;
};

//--------------------------------------------------------------------
AnaTree::AnaTree(fhicl::ParameterSet const & pset) 
 : EDAnalyzer{p}
{
  this->reconfigure(pset);
}

//--------------------------------------------------------------------
AnaTree::~AnaTree()
{}

//--------------------------------------------------------------------
void AnaTree::reconfigure(fhicl::ParameterSet const & pset)
{
  fTreeName           = pset.get< std::string >("TreeName", "anatree");
  fG4ModuleLabel      = pset.get< std::string >("G4ModuleLabel", "largeant");
  fDoPhotonProjection = pset.get< std::string >("DoPhotonProjection", true);

  if (fDoPhotonProjection) fProjAlg.reconfigure(pset);
  return;
}

//--------------------------------------------------------------------
void AnaTree::analyze(art::Event const & evt)
{
  // Reset variables
  ResetVars();

  // Geometry Service 
  art::ServiceHandle<geo::Geometry> geom;
  // Detector properties service 
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  // MC particle list
  art::Handle< std::vector<simb::MCParticle> > plistHandle;
  e.getByLabel(fG4ModuleLabel, plistHandle);
  auto plist = *plistHandle;

  run    = evt.run();
  subrun = evt.subRun();
  event  = evt.id().event();
   
  std::cout<<std::endl;
  std::cout<<"========================================="<<std::endl;
  std::cout << "Run = "         << run 
            << ", SubRun = "    << subrun 
            << ", Evt = "       << event    << std::endl;
  std::cout<<"========================================="<<std::endl;
  std::cout<<std::endl;

  // Apply our photon projection alg
  if (fDoPhotonProjection)
  {
    pixelIdVec.reserve(pixelMap.size());
    pixelCountVec.reserve(pixelMap.size());
    for (const auto& p : pixelMap)
    {
      pixelIdVec.push_back(p.first);
      pixelCountVec.push_back(p.second);
    }
  }

  // Electric Field 
  efield[0] = detprop->Efield(0);
  efield[1] = detprop->Efield(1);
  efield[2] = detprop->Efield(2);
   
  // Setting a boolian to only output MC info if this is MC-info 
  bool isdata = false;
  if (evt.isRealData()) {isdata = true;}
  else isdata = false;
   
  //
  // Filling MCTruth information
  //
  if(!isdata)
  {   
    // Setting a string for process type 
    std::string pri("primary");
    std::string hadElastic("hadElastic");
    std::string CoulombScat("CoulombScat");
    std::string ProtonInelastic("protonInelastic");
    
    int primary(0);
    int geant_particle(0);
    
    // Determine the number of primary particles from geant 
    for( unsigned int i = 0; i < plist.size(); ++i )
	  {
	    geant_particle++;
	    if(plist[i].Process()==pri) {primary++;}
	  }//<---End i loop
    no_primaries    = primary;
    geant_list_size = plisticle;
     
    // Looping over all the Geant4 particles 
    int iPrim = 0;
    for(unsigned int i = 0; i < plist.size(); ++i)
    {
      // Check if primary
	    if(plist[i].Process()==pri){process_primary.push_back(1);}
	    else                        {process_primary.push_back(0);}
	   
      // Get the process
      Process.push_back(plist[i].Process());
    
	    // Saving other info 
      PDG.push_back(plist[i].PdgCode());
      Mother.push_back(plist[i].Mother());
      TrackId.push_back(plist[i].TrackId());
      StartEnergy.push_back(plist[i].E());
	    EndEnergy.push_back(plist[i].EndE());

      double startp = sqrt( pow(plist[i].Px(),2) + pow(plist[i].Py(),2) + pow(plist[i].Pz(),2));
      double mass = plist[i].Mass(); 
      double startke = sqrt( pow(startp,2) + pow(mass,2) ) - mass;
      StartKE.push_back(startke);

      double lastke = 0;
	    simb::MCTrajectory truetraj = plist[i].Trajectory();
      int counter_ii = 0;
	    for(auto itTraj = truetraj.begin(); itTraj != truetraj.end(); ++itTraj) 
      {
        if(truetraj.Z(counter_ii) > 0   && truetraj.Z(counter_ii) < 90   &&
           truetraj.X(counter_ii) > 0   && truetraj.X(counter_ii) < 47.5 &&
           truetraj.Y(counter_ii) > -20 && truetraj.Y(counter_ii) < 20) 
        {
          double thisp;
          if(counter_ii) 
          {
            thisp = sqrt( pow(truetraj.Px(counter_ii-1),2) + 
                          pow(truetraj.Py(counter_ii-1),2) + 
                          pow(truetraj.Pz(counter_ii-1),2));
          }
          else 
          {
            thisp = sqrt( pow(truetraj.Px(counter_ii),2) + 
                          pow(truetraj.Py(counter_ii),2) + 
                          pow(truetraj.Pz(counter_ii),2));
          }
          lastke = sqrt( pow(thisp,2) + pow(mass,2) ) - mass;
        }//<--End if in tpc
        counter_ii++;
	    }//<--End loop on true trajectory points
      LastKE.push_back(lastke);
 
	    // Saving the start and end Px, Py, Pz info 
	    StartPx.push_back(plist[i].Px());
	    StartPy.push_back(plist[i].Py());
	    StartPz.push_back(plist[i].Pz());
	    EndPx.push_back(plist[i].EndPx());
	    EndPy.push_back(plist[i].EndPy());
	    EndPz.push_back(plist[i].EndPz());
	   
	    // Saving the Start and End Point for this particle 
	    StartPointx.push_back(plist[i].Vx());
	    StartPointy.push_back(plist[i].Vx());
	    StartPointz.push_back(plist[i].Vx());
	    EndPointx.push_back(plist[i].EndPosition()[0]);
	    EndPointy.push_back(plist[i].EndPosition()[1]);
	    EndPointz.push_back(plist[i].EndPosition()[2]);

	    // Saving the processes for this particle 
	    G4Process.push_back( plist[i].Process() );
	    G4FinalProcess.push_back( plist[i].EndProcess() );
 	   
	    // Saving the number of Daughters for this particle 
	    NumberDaughters.push_back(plist[i].NumberDaughters());

  	  // Save intermediary information for the primary track
  	  if (plist[i].Process() == pri)
      {
        primary_p = plist[i].P();
  	    NTrTrajPts.push_back(plist[i].NumberTrajectoryPoints());
        std::vector<double> midx;  std::vector<double> midy;  std::vector<double> midz;
        std::vector<double> midpx; std::vector<double> midpy; std::vector<double> midpz;
   
  	    int iPrimPt = 0;	   
  	    for(auto itTraj = truetraj.begin(); itTraj != truetraj.end(); ++itTraj)
  	    {
          // Pushing back vars 
          midx.push_back(truetraj.X(iPrimPt));
          midy.push_back(truetraj.Y(iPrimPt));
          midz.push_back(truetraj.Z(iPrimPt));
          midpx.push_back(truetraj.Px(iPrimPt));
          midpy.push_back(truetraj.Py(iPrimPt));
          midpz.push_back(truetraj.Pz(iPrimPt));
  	      iPrimPt++;
  	    }//<--End loop on true trajectory points

        MidPosX.push_back(midx); MidPosY.push_back(midy); MidPosZ.push_back(midz);
        MidPx.push_back(midpx);  MidPy.push_back(midpy);  MidPz.push_back(midpz);
	    
        // Look at interesting points
  	    auto thisTracjectoryProcessMap =  truetraj.TrajectoryProcesses();
  	    // Ok, if the size of the map is 0, all the action might happen at the end of the track
  	    // So we check the daugthers:
  	    //    - Case 1. There are daugthers:
  	    //               * The interesting point is the last one
  	    //               * The interaction type is the one that created the first daugther (this might be improved)
  	    //    - Case 2. There are NO daugthers:
  	    //              * We assign the interaction type to be 0: nothing happens, thought going particle
  	    //              * The interesting point is the last one (might not be in the TPC)
  	    if(!thisTracjectoryProcessMap.size())
  	    {
  	      int interestingPoint = (int) (NTrTrajPts[i] - 1);
  	      InteractionPoint.push_back(interestingPoint);
  	      if(NumberDaughters[i])
  	      {		 
  	        auto thePrimaryDaughterID = plist[i]. Daughter(0); 
  	        for(unsigned int iD = 0; iD < geant_part.size(); ++iD )
  	        {
  	          if(plist[iD].TrackId() == thePrimaryDaughterID) {InteractionPointType.push_back(plist[iD].Process());}
  	        }//<--- End particle loop
  	      }//<--- End if there are daughters
          else {InteractionPointType.push_back("nodaughters");}
  	    }
        else
   	    {
   	      // The map is not zero: somthing interesting might happen in the middle of the track!!
   	      for(auto const& couple: thisTracjectoryProcessMap) 
   	      {
   	        int interestingPoint = (int) couple.first;
   	        InteractionPoint.push_back(interestingPoint);
            InteractionPointType.push_back(truetraj.KeyToProcess(couple.second));
   	      }
   	    }	
   
	      iPrim++;
	    }//<--End if primary

  	  if (plist[i].Process() != pri)
      {
        if (plist[i].Mother() == 1) 
        {
          if (plist[i].PdgCode() < 10000)
          {
  	        NDTrTrajPts.push_back(plist[i].NumberTrajectoryPoints());
            DTrackId.push_back(plist[i].TrackId());
            DPdgCode.push_back(plist[i].PdgCode());
            DStartEnergy.push_back(plist[i].E());
            DStartP.push_back(plist[i].P());
  
  	        simb::MCTrajectory truetraj = plist[i].Trajectory();
  	        int jDaughtPt = 0;	
            std::vector<double> midx; std::vector<double> midy; std::vector<double> midz;
  	        for(auto itTraj = truetraj.begin(); itTraj != truetraj.end(); ++itTraj)
  	        {
              midx.push_back(truetraj.X(jDaughtPt));
              midy.push_back(truetraj.Y(jDaughtPt));
              midz.push_back(truetraj.Z(jDaughtPt));
  	          jDaughtPt++;
  	        }//<--End loop on true trajectory points
            DMidPosX.push_back(midx); DMidPosY.push_back(midy); DMidPosZ.push_back(midz);
          }
        }//<-- End if primary daughter
      }//<-- End if not primary
    }//<--End loop on geant particles   
  }//<---End checking if this is MC 

  fTree->Fill();
}//<---End analyze()

//--------------------------------------------------------------------
void AnaTree::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>(fTreeName.c_str(),  fTreeName.c_str());
  fTree->Branch("run",                         &run,"run/I");
  fTree->Branch("subrun",                      &subrun,"subrun/I");
  fTree->Branch("event",                       &event,"event/I");
  fTree->Branch("evttime",                     &evttime,"evttime/D");
  fTree->Branch("efield",                      efield,"efield[3]/D"); 
  fTree->Branch("no_primaries",                &no_primaries,"no_primaries/I");
  fTree->Branch("geant_list_size",             &geant_list_size,"geant_list_size/I");
  fTree->Branch("primary_p",                   &primary_p,"primary/D");
  fTree->Branch("PDG",                         &PDG);
  fTree->Branch("StartKE",                     &StartKE);
  fTree->Branch("LastKE",                      &LastKE);
  fTree->Branch("StartEnergy",                 &StartEnergy);
  fTree->Branch("StartPx",                     &StartPx);
  fTree->Branch("StartPy",                     &StartPy);
  fTree->Branch("StartPz",                     &StartPz);
  fTree->Branch("EndEnergy",                   &EndEnergy);
  fTree->Branch("EndPx",                       &EndPx);
  fTree->Branch("EndPy",                       &EndPy);
  fTree->Branch("EndPz",                       &EndPz);
  fTree->Branch("StartPointx",                 &StartPointx);
  fTree->Branch("StartPointy",                 &StartPointy);
  fTree->Branch("StartPointz",                 &StartPointz);
  fTree->Branch("EndPointx",                   &EndPointx);
  fTree->Branch("EndPointy",                   &EndPointy);
  fTree->Branch("EndPointz",                   &EndPointz);
  fTree->Branch("Process",                     &Process);
  fTree->Branch("NumberDaughters",             &NumberDaughters);
  fTree->Branch("Mother",                      &Mother);
  fTree->Branch("TrackId",                     &TrackId);
  fTree->Branch("process_primary",             &process_primary);
  fTree->Branch("G4Process",                   &G4Process);
  fTree->Branch("G4FinalProcess",              &G4FinalProcess);  
  fTree->Branch("NTrTrajPts",                  &NTrTrajPts);
  fTree->Branch("NProtonDaughters",            &NProtonDaughters,"NProtonDaughters/I");
  fTree->Branch("NNeutronDaughters",           &NNeutronDaughters,"NNeutronDaughters/I");
  fTree->Branch("NDTrTrajPts",                 &NDTrTrajPts);
  fTree->Branch("DTrackId",                    &DTrackId);
  fTree->Branch("DPdgCode",                    &DPdgCode);
  fTree->Branch("DStartEnergy",                &DStartEnergy);
  fTree->Branch("DStartP",                     &DStartP);
  fTree->Branch("MidPosX",                     &MidPosX);
  fTree->Branch("MidPosY",                     &MidPosY);
  fTree->Branch("MidPosZ",                     &MidPosZ);
  fTree->Branch("MidPx",                       &MidPx);
  fTree->Branch("MidPy",                       &MidPy);
  fTree->Branch("MidPz",                       &MidPz);
  fTree->Branch("DMidPosX",                   &DMidPosX);
  fTree->Branch("DMidPosY",                   &DMidPosY);
  fTree->Branch("DMidPosZ",                   &DMidPosZ);
  fTree->Branch("InteractionPoint",            &InteractionPoint);
  fTree->Branch("InteractionPointType",        &InteractionPointType);
  fTree->Branch("nPrimScint",         &nPrimScint,         "nPrimScint/l");
  fTree->Branch("nReadoutIncPhotons", &nReadoutIncPhotons, "nReadoutIncPhotons/l");
  fTree->Branch("pixelIDs" ,    &pixelIdVec    );
  fTree->Branch("pixelHits" ,   &pixelCountVec );
}

//--------------------------------------------------------------------
void AnaTree::ResetVars()
{
  G4Process.clear();
  G4FinalProcess.clear();
  InteractionPoint.clear();
  InteractionPointType.clear();
  PDG.clear();
  Mother.clear();
  TrackId.clear();
  NumberDaughters.clear();
  process_primary.clear();
  StartPointx.clear();  
  StartPointy.clear();
  StartPointz.clear();
  StartKE.clear();
  LastKE.clear();
  StartEnergy.clear();
  StartPx.clear();
  StartPy.clear();
  StartPz.clear(); 
  EndPointx.clear();
  EndPointy.clear(); 
  EndPointz.clear();  
  EndEnergy.clear();  
  EndPx.clear();  
  EndPy.clear();  
  EndPz.clear();  
  NDTrTrajPts.clear();
  DTrackId.clear();
  DPdgCode.clear();
  NTrTrajPts.clear();
  MidPosX.clear();
  MidPosY.clear();
  MidPosZ.clear();
  MidPx.clear();
  MidPy.clear();
  MidPz.clear();
  DMidPosX.clear();
  DMidPosY.clear();
  DMidPosZ.clear();
  DStartEnergy.clear();
  DStartP.clear();
  pixelIdVec.clear();
  pixelCountVec.clear();

  nPrimScint = 0;        
  nReadoutIncPhotons = 0;

  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  for (int i = 0; i<3; ++i){
    efield[i] = -99999;
  }

  no_primaries = -99999;
  geant_list_size=-999;
  primary_p=-999;
}

}

DEFINE_ART_MODULE(amselsim::AnaTree)
