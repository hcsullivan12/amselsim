////////////////////////////////////////////////////////////////////////
// DetectorPropertiesServiceAmSel.h
//
// Service interface for DetectorProperties functions
//
//  hsulliva@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETECTORPROPERTIESSERVICEAMSEL_H
#define DETECTORPROPERTIESSERVICEAMSEL_H

#include "amselsim/Services/DetectorPropertiesAmSel.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

///General LArSoft Utilities
namespace ldp{
  
  /**
   * @brief "AmSel" implementation of DetectorProperties service
   * 
   * This class wraps DetectorPropertiesAmSel provider into a art service.
   * It delivers the provider via the standard interface:
   *     
   *     detinfo::DetectorProperties const* detprop
   *       = art::ServiceHandle<detinfo::DetectorPropertiesAmSel>()
   *       ->provider();
   *     
   * or, using the standard interface in "CoreUtils/ServiceUtil.h":
   *     
   *     auto const* detprop
   *       = lar::providerFrom<detinfo::DetectorPropertiesAmSel>();
   *     
   * In addition to the functionality of the provider, this service allows
   * to read the configuration from the input file, inherited from a previous
   * run.
   * 
   * Configuration parameters
   * -------------------------
   * 
   * This service passes the whole configuration down to its service provider,
   * but it also reacts to:
   * - *InheritNumberTimeSamples* (boolean; default: false): if true, the
   *   configuration database in the ROOT input file is queried and if a
   *   configuration for this service is found, it's used instead of the
   *   one from the current FHiCL configuration
   * 
   */
  
  class DetectorPropertiesServiceAmSel : public detinfo::DetectorPropertiesService {

    public:
      
      // the following is currently not used for validation,
      // but only for documentation
      struct ServiceConfiguration_t {
        
        // service-specific configuration
        fhicl::Atom<bool> InheritNumberTimeSamples {
          fhicl::Name("InheritNumberTimeSamples"),
          fhicl::Comment(""),
          false /* default value */
        };
        
        // provider configuration
        ldp::DetectorPropertiesAmSel::Configuration_t ProviderConfiguration;
        
      }; // ServiceConfiguration_t
      
      
      // this enables art to print the configuration help:
      using Parameters = art::ServiceTable<ServiceConfiguration_t>;
      
      DetectorPropertiesServiceAmSel(fhicl::ParameterSet const& pset,
				art::ActivityRegistry& reg);

      virtual void   reconfigure(fhicl::ParameterSet const& pset) override;
      void   preProcessEvent(const art::Event& evt, art::ScheduleContext);
      void   postOpenFile(const std::string& filename);
      void   preBeginRun(const art::Run& run) { fGotElectronLifetimeFromDB=false; }
	
      virtual const provider_type* provider() const override { return fProp.get();}
      
    private:

      std::unique_ptr<ldp::DetectorPropertiesAmSel> fProp;
      fhicl::ParameterSet   fPS;       ///< Original parameter set.
      
      bool fInheritNumberTimeSamples; ///< Flag saying whether to inherit NumberTimeSamples
      
      bool isDetectorPropertiesServiceAmSel(const fhicl::ParameterSet& ps) const;

      bool fUseDatabaseForMC;
      bool fGotElectronLifetimeFromDB;
      
    }; // class DetectorPropertiesService
} //namespace detinfo
DECLARE_ART_SERVICE_INTERFACE_IMPL(ldp::DetectorPropertiesServiceAmSel, detinfo::DetectorPropertiesService, LEGACY)
#endif // DETECTORPROPERTIESSERVICEAMSEL_H
