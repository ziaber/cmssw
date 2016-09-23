#ifndef SimTotem_RPDigiProducer_RPDigiProducer_h
#define SimTotem_RPDigiProducer_RPDigiProducer_h

// -*- C++ -*-
//
// Package:    RPDigiProducer
// Class:      RPDigiProducer
// 
/**\class RPDigiProducer RPDigiProducer.cc SimTotem/RPDigiProducer/src/RPDigiProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Hubert NIEWIADOMSKI
//         Created:  Tue Mar 27 11:57:03 CEST 2007
// $Id: RPDigiProducer.h,v 1.6 2008/09/12 17:26:25 lgrzanka Exp $
//
//

#include "boost/shared_ptr.hpp"

// system include files
#include <memory>
#include <vector>
#include <map>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/TotemDigi/interface/TotemRPDigi.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
#include "SimTotem/RPDigiProducer/interface/RPSimTypes.h"

#include "SimTotem/RPDigiProducer/interface/RPDetDigitizer.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"

#include "DataFormats/Common/interface/DetSet.h"
#include "SimTotem/RPDigiProducer/interface/DeadChannelsManager.h"

//
// class decleration
//


namespace CLHEP {
  class HepRandomEngine;
}


class RPDigiProducer : public edm::EDProducer {
   public:
      explicit RPDigiProducer(const edm::ParameterSet&);
      ~RPDigiProducer();

   private:
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob();

      edm::DetSet<TotemRPDigi> convertRPStripDetSet(const edm::DetSet<RPStripDigi>&);

    // ----------member data ---------------------------
      std::vector<std::string> RP_hit_containers_;
      typedef std::map<unsigned int, std::vector<PSimHit> > simhit_map;
      typedef simhit_map::iterator simhit_map_iterator;
      simhit_map SimHitMap;
      
      edm::ParameterSet conf_;
      std::map<RPDetId, boost::shared_ptr<RPDetDigitizer> > theAlgoMap;
      std::vector<edm::DetSet<TotemRPDigi> > theDigiVector;
      std::vector<edm::DetSet<RPDetTrigger> > theTriggerVector;

      CLHEP::HepRandomEngine* rndEngine = nullptr;
      int verbosity_;

      /**
       * this variable answers the question whether given channel is dead or not
       */
      DeadChannelsManager deadChannelsManager;
      /**
       * this variable indicates whether we take into account dead channels or simulate as if all
       * channels work ok (by default we do not simulate dead channels)
       */
      bool simulateDeadChannels;

      edm::EDGetTokenT<CrossingFrame<PSimHit>> tokenCrossingFrameTotemRP;
};


#endif  //SimTotem_RPDigiProducer_RPDigiProducer_h
