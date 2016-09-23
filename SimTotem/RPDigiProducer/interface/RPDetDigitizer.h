#ifndef SimTotem_RPDigiProducer_RP_DET_DIGITIZER_H
#define SimTotem_RPDigiProducer_RP_DET_DIGITIZER_H

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include <vector>
#include <string>


//#include "SimGeneral/HepPDT/interface/HepPDTable.h"
#include "SimTracker/Common/interface/SiG4UniversalFluctuation.h"
#include "SimGeneral/NoiseGenerators/interface/GaussianTailNoiseGenerator.h"

#include "SimTotem/RPDigiProducer/interface/RPSimTypes.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "SimTotem/RPDigiProducer/interface/RPHitChargeConverter.h"
#include "SimTotem/RPDigiProducer/interface/RPVFATSimulator.h"
#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
#include "SimTotem/RPDigiProducer/interface/RPDisplacementGenerator.h"
#include "SimTotem/RPDigiProducer/interface/RPGaussianTailNoiseAdder.h"
#include "SimTotem/RPDigiProducer/interface/RPPileUpSignals.h"


namespace CLHEP{
        class HepRandomEngine;
}


class RPDetDigitizer
{
  public:
    RPDetDigitizer(const edm::ParameterSet &params, CLHEP::HepRandomEngine& eng, RPDetId det_id, const edm::EventSetup& iSetup);
    void run(const std::vector<PSimHit> &input, const std::vector<int> &input_links, 
        std::vector<RPStripDigi> &output_digi, std::vector<RPDetTrigger> &output_trig, 
        SimRP::DigiPrimaryMapType &output_digi_links, 
        SimRP::TriggerPrimaryMapType &output_trig_links);
    ~RPDetDigitizer();
      
  private:
    RPGaussianTailNoiseAdder *theRPGaussianTailNoiseAdder;
    RPPileUpSignals *theRPPileUpSignals;
    RPHitChargeConverter *theRPHitChargeConverter;
    RPVFATSimulator *theRPVFATSimulator;
    RPDisplacementGenerator *theRPDisplacementGenerator;

  private:
    const edm::ParameterSet &params_;

    int numStrips;
    double theNoiseInElectrons;   // Noise (RMS) in units of electrons.
    double theStripThresholdInE;  // Strip noise treshold in electorns.
    bool noNoise;                 //if the nos is included
    RPDetId det_id_;
    bool misalignment_simulation_on_;
    int verbosity_;
    bool  _links_persistence;
};

#endif  //SimTotem_RPDigiProducer_RP_DET_DIGITIZER_H
