namespace cheat{

  //--------------------------------------------------------------------
  template<typename Evt> //DO NOT USE THIS FUNCTION FROM WITHIN ART! The ParticleInventoryService is designed to impliment these methods as cleanly as possible within the art framework.
    void ParticleInventory::PrepEvent        (const Evt& evt ){
      if(!(this->CanRun(evt))){ 
        throw cet::exception("ParticleInventory") 
          << "Particle Inventory cannot function. "
          << "Is this file real data?";
      }
      fParticleList.clear();
      fMCTObj.fMCTruthList.clear();
      fMCTObj.fTrackIdToMCTruthIndex.clear();
      this->PrepParticleList(evt);
      this->PrepMCTruthList(evt);
      this->PrepTrackIdToMCTruthIndex(evt);
    }


  //--------------------------------------------------------------------
  template<typename Evt>
    void ParticleInventory::PrepParticleList(const Evt& evt ) const{

      if(this->ParticleListReady( )){ //The particle list already exists. Do nothing.
        return;
      }
      //The particle list needs to be built
      //We use auto so that we(the compiler) can determine which type we need for either art or gallery.
      const auto& pHandle = evt.template getValidHandle<std::vector<simb::MCParticle>>(fG4ModuleLabel); 
//      const auto& partVecIn = evt.template getValidHandle<std::vector<simb::MCParticle>>(fG4ModuleLabel); 
//      if(pHandle.failedToGet()){
//        /*mf::LogWarning("BackTracker") << "failed to get handle to simb::MCParticle from "
//          << fG4ModuleLabel
//          << ", return";*/ //Do this silently so we don't polute the logs. It is expected to fail for all gen and g4 jobs.
//        return;
//      }

      const auto& partVecIn = *pHandle;
      for(const auto& partIn : partVecIn){
        fParticleList.Add(new simb::MCParticle(partIn)); //Is this still doing a copy? If so, another method should be used.
      }
      if(fEveIdCalculator=="EmEveIdCalculator"){
        fParticleList.AdoptEveIdCalculator(new sim::EmEveIdCalculator);
      }else if(fEveIdCalculator=="EveIdCalculator"){
        fParticleList.AdoptEveIdCalculator(new sim::EveIdCalculator);
      }else{
        throw cet::exception("ParticleInventory3") 
          << "Particle Inventory cannot initialize the particle list.\n "
          << fEveIdCalculator <<" is not a known EveIdCalculator.\n";
      }
    }

  //--------------------------------------------------------------------
  template<typename Evt> //I may want to make this function private.
    void ParticleInventory::PrepMCTruthListAndTrackIdToMCTruthIndex(const Evt& evt ) const{
      if( this->TrackIdToMCTruthReady() && this->MCTruthListReady( ) ){ return;} 
      this->PrepParticleList( evt); //Make sure we have built the particle list for this event
      //const auto& mcpmctAssnsIn = *( evt.template getValidHandle<art::Assns<simb::MCParticle,simb::MCTruth>>(fG4ModuleLabel));
      const auto& mcpmctAssnsHandle =  evt.template getValidHandle<art::Assns<simb::MCParticle,simb::MCTruth>>(fG4ModuleLabel);
      const auto& mcpmctAssnsIn = *mcpmctAssnsHandle;
      for( const auto& mcpmctAssnIn : mcpmctAssnsIn){    //Assns are themselves a container. Loop over entries.
        const art::Ptr<simb::MCParticle>& part=mcpmctAssnIn.first;
        const art::Ptr<simb::MCTruth>&    mct =mcpmctAssnIn.second;
        unsigned short mctruth_idx = USHRT_MAX;
        for (size_t i = 0; i<fMCTObj.fMCTruthList.size(); ++i){
          if (fMCTObj.fMCTruthList[i] == mct){
            mctruth_idx = i;
            break;
          }
        }
        if (mctruth_idx == USHRT_MAX){
          fMCTObj.fMCTruthList.push_back(mct);
          fMCTObj.fTrackIdToMCTruthIndex.emplace(part->TrackId(), fMCTObj.fMCTruthList.size() - 1);
        }
        else{
          fMCTObj.fTrackIdToMCTruthIndex.emplace(part->TrackId(), mctruth_idx );
        }
      }
    }

  //--------------------------------------------------------------------
  template<typename Evt>
    void ParticleInventory::PrepMCTruthList             (const Evt& evt ) const{
      if(this->MCTruthListReady( ) ){ return;} //If the event is data or if the truth list is already built there is nothing for us to do.
      PrepMCTruthListAndTrackIdToMCTruthIndex( evt); //TrackIdToMCTruthIndex and MCTruthList are prepared at the same time. The access of information makes this the most convenient way to do so. It is only somewhat more expensive for the memory, but significantly less expensive for time.
    }

  //--------------------------------------------------------------------
  template<typename Evt>
    void ParticleInventory::PrepTrackIdToMCTruthIndex(const Evt& evt ) const{
      if(this->TrackIdToMCTruthReady()){ return;} //The list already exists. Do nothing.
      PrepMCTruthListAndTrackIdToMCTruthIndex( evt); //TrackIdToMCTruthIndex and MCTruthList are prepared at the same time. The access of information makes this the most convenient way to do so. It is only somewhat more expensive for the memory, but significantly less expensive for time.
    }

  template<typename Evt>
    bool ParticleInventory::CanRun(const Evt& evt) const{
      return !(evt.isRealData());
    }

}//end namespace
