/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
****************************************************************************/

#include "DataFormats/TotemRawDataFormats/interface/RawEvent.h"
#include "DataFormats/TotemRawDataFormats/interface/VFATFrame.h"
#include "DataFormats/TotemRawDataFormats/interface/FramePosition.h"
#include "DataFormats/TotemRawDataFormats/interface/VFATFrameCollection.h"
#include "DataFormats/TotemRawDataFormats/interface/SimpleVFATFrameCollection.h"
#include "DataFormats/TotemRawDataFormats/interface/OptoRxVFATFrameCollection.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include <vector>

namespace { namespace {
	edm::Wrapper<Totem::RawEvent> dummy0;

	Totem::VFATFrame dummy1;
	std::vector<Totem::VFATFrame> dummy2;

	Totem::FramePosition dummy4;

	// the class is abstract
	//Totem::VFATFrameCollection dummy5;
	//edm::Wrapper< Totem::VFATFrameCollection > dummy6;

	Totem::SimpleVFATFrameCollection dummy7;
	edm::Wrapper< Totem::SimpleVFATFrameCollection > dummy8;

	Totem::OptoRxVFATFrameCollection dummy9;
	edm::Wrapper< Totem::OptoRxVFATFrameCollection > dummy10;
	
	Totem::OptoRxMetaData dummy17;
	std::map<unsigned int, Totem::OptoRxMetaData> dummy18;

	Totem::TriggerData dummy19;

	std::map<unsigned int, long> dummy100;
}}
