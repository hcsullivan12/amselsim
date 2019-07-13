///////////////////////////////////////////////////////////////////////
///
/// \file   AmSelSignalShapingService.h
///
/// \brief  Service to provide microboone-specific signal shaping for
///         simulation (convolution) and reconstruction (deconvolution).
///
/// \author H. Greenlee 
///
/// This service inherits from SignalShaping and supplies
/// microboone-specific configuration.  It is intended that SimWire and
/// CalWire modules will access this service.
///
/// FCL parameters:
///
/// FieldBins       - Number of bins of field response.
/// Col3DCorrection - 3D path length correction for collection plane.
/// ColFieldRespAmp - Collection field response amplitude.
/// ShapeTimeConst  - Time constants for exponential shaping.
/// ColFilter       - Root parameterized collection plane filter function.
/// ColFilterParams - Collection filter function parameters.
///
////////////////////////////////////////////////////////////////////////

#ifndef AMSELSIM_SIGNALSHAPINGSERVICET_H
#define AMSELSIM_SIGNALSHAPINGSERVICET_H

#include <vector>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "lardata/Utilities/SignalShaping.h"
#include "TF1.h"
#include "TH1D.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace util {
  class AmSelSignalShapingService {
  public:

    // Constructor, destructor.

    AmSelSignalShapingService(const fhicl::ParameterSet& pset,
				   art::ActivityRegistry& reg);
    ~AmSelSignalShapingService();

    // Update configuration parameters.

    void reconfigure(const fhicl::ParameterSet& pset);
    const util::SignalShaping& SignalShaping() const;
    double ElectronsToCurrent(double const& el) const;

    // Do convolution calcution (for simulation).
    template <class T> void Convolute(std::vector<T>& func) const;

  private:

    void init() const{const_cast<AmSelSignalShapingService*>(this)->init();}
    void init();

    // Calculate response functions.
    // Copied from SimWireT1034.

    void SetElectResponse();

    bool fInit;               ///< Initialization flag.

    // Following attributes hold the convolution and deconvolution kernels
    util::SignalShaping fColSignalShaping;

    // Electronics response.
    std::vector<double> fElectResponse;
  };
}

//----------------------------------------------------------------------
// Do convolution.
template <class T> inline void util::AmSelSignalShapingService::Convolute(std::vector<T>& func) const
{ 
  for (auto& s : func) s = ElectronsToCurrent(s);
  SignalShaping().Convolute(func);
}


DECLARE_ART_SERVICE(util::AmSelSignalShapingService, LEGACY)
#endif
