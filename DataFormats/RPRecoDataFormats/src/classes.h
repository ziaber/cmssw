
#include "DataFormats/RPRecoDataFormats/interface/RPTrackCandidate.h"
#include "DataFormats/RPRecoDataFormats/interface/RPTrackCandidateCollection.h"
#include "DataFormats/RPRecoDataFormats/interface/RPMulTrackCandidateCollection.h"
#include "DataFormats/RPRecoDataFormats/interface/RPTrackCandidateDistinctCollectionsSet.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/RPRecoDataFormats/interface/RPFittedTrack.h"
#include "DataFormats/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "DataFormats/RPRecoDataFormats/interface/RPMulFittedTrackCollection.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDigCluster.h"
#include <map>
#include "TObject.h"
#include "DataFormats/RPRecoDataFormats/interface/RPReconstructedProton.h"
#include "DataFormats/RPRecoDataFormats/interface/RPReconstructedProtonCollection.h"
#include "DataFormats/RPRecoDataFormats/interface/RP2DHit.h"
#include "DataFormats/RPRecoDataFormats/interface/RP2DHitDebug.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "DataFormats/RPRecoDataFormats/interface/RPReconstructedProtonPair.h"
#include "DataFormats/RPRecoDataFormats/interface/RPReconstructedProtonPairCollection.h"
#include "DataFormats/RPRecoDataFormats/interface/RPRecoElasticEvent.h"

#include <vector>
#include "DataFormats/RPRecoDataFormats/interface/RPRecognizedPatterns.h"
#include "DataFormats/RPRecoDataFormats/interface/RPRecognizedPatternsCollection.h"
#include "DataFormats/RPRecoDataFormats/interface/CentralMassInfo.h"

#include "DataFormats/RPRecoDataFormats/interface/RPStationTrackFit.h"
#include "DataFormats/RPRecoDataFormats/interface/RPStationTrackFitCollection.h"



namespace {
  namespace {
    
    CentralMassInfo ctrinfo;
    edm::Wrapper<CentralMassInfo> wctrinfo;
    
    RPTrackCandidate tc;
    RPTrackCandidateCollection coll;
    RPMulTrackCandidateCollection mcoll;
    std::pair<RPId, RPTrackCandidate> p_rpid_rptrcand;
    std::pair<RPId, std::vector<RPTrackCandidate> > p_rpid_rptrcandvec;
    std::map<RPId, RPTrackCandidate> m_rpid_rptrcand;
    std::map<RPId, std::vector<RPTrackCandidate> > m_rpid_rptrcandvec;
    RPTrackCandidateDistinctCollectionsSet coll_set;
    std::vector<RPTrackCandidate> v_tr_cand;
    edm::Wrapper<RPTrackCandidateCollection> RPTrackCandidateCollectionWrapper;
    edm::Wrapper<RPMulTrackCandidateCollection> RPMulTrackCandidateCollectionWrapper;
    edm::Wrapper<RPTrackCandidateDistinctCollectionsSet> RPTrackCandidateDistinctCollectionsSetWrapper;
    
    std::vector<std::vector<RPFittedTrack> > vvrpfittrack;
    edm::Wrapper<std::vector<std::vector<RPFittedTrack> > > wrvvrpfittrack;
    std::map<RPId, std::vector<std::vector<RPFittedTrack> > > mapwrvvrpfittrack;
    edm::Wrapper<std::map<RPId, std::vector<std::vector<RPFittedTrack> > > > wrapmapwrvvrpfittrack;
    edm::Wrapper<RPMulFittedTrackSetsCollection> wraprpmulttracksetcol; 
    
    std::vector<std::vector<RPTrackCandidate> > vvrptrcand;
    std::map<RPId, std::vector<std::vector<RPTrackCandidate> > > mapvvrptrcand;
    RPMulTrackCandidateSetsCollection rpmultracsetcandcol;
    edm::Wrapper<std::vector<std::vector<RPTrackCandidate> > > wrapvvrptrcand;
    edm::Wrapper<std::map<RPId, std::vector<std::vector<RPTrackCandidate> > > > wrapmapvvrptrcand;
    edm::Wrapper<RPMulTrackCandidateSetsCollection> wraprpmultracsetcandcol;
    
    std::less<unsigned int> lui;
    std::allocator<std::pair<const unsigned int,RPTrackCandidate> > apcuirptc;
    std::binary_function<unsigned int,unsigned int,bool> bf;

    RPFittedTrack the_fitted_track;
    RPFittedTrackCollection the_track_cand_col;
    RPMulFittedTrackCollection the_mtrack_cand_col;
    std::pair<RPId, RPFittedTrack> p_rpid_rpfittedtr;
    std::pair<RPId, std::vector<RPFittedTrack> > p_rpid_rpfittedtrvec;
    std::vector<RPFittedTrack>  v_fitted_track;
    edm::Wrapper<RPFittedTrack> the_w_rpft;
    edm::Wrapper<RPFittedTrackCollection> the_w_rpftc;
    edm::Wrapper<RPMulFittedTrackCollection> the_w_mrpftc;
    
    RPDetHitPoint rpdhp;
    TVector3 tv3;
    std::vector<RPDetHitPoint> vrdhp;
    std::map<RPId, RPFittedTrack> mdeidfittrac;
    std::map<RPId, std::vector<RPFittedTrack> > vmdeidfittrac;
    edm::Wrapper<std::map<RPId, RPFittedTrack> > wmdeidfittrac;
    edm::Wrapper<std::map<RPId, std::vector<RPFittedTrack> > > vwmdeidfittrac;
    edm::Wrapper<RPDetHitPoint> wrpdhp;
    edm::Wrapper<TVector3> wtv3;
    edm::Wrapper<std::vector<RPDetHitPoint> > wvrdhp;
    
    TObject to;
    edm::Wrapper<TObject> wto;
    
    RPReconstructedProton rprecprot;
    edm::Wrapper<RPReconstructedProton> wraprprecprot;
    std::vector<RPReconstructedProton> prprecprotvec;
    
    RPReconstructedProtonCollection rprpcol;
    edm::Wrapper<RPReconstructedProtonCollection> wrrprpcol;
    
    RP2DHit rp2dhit;
    edm::Wrapper<RP2DHit> wrp2dhit;
    RP2DHitDebug rp2debugdhit;
    edm::Wrapper<RP2DHitDebug> rwp2debugdhit;
    
    std::map<RPId, RP2DHitDebug> vrp2debugdhit;
    edm::Wrapper<std::map<RPId, RP2DHitDebug> > wvrp2debugdhit;
    std::pair<unsigned int, RP2DHitDebug> vrp2debughitpair;
    
    RPReconstructedProtonPair p;
    edm::Wrapper<RPReconstructedProtonPair> wp;
    std::vector<RPReconstructedProtonPair> vp;
    
    RPReconstructedProtonPairCollection pc;
    edm::Wrapper<RPReconstructedProtonPairCollection> wpc;

	RPRecoElasticEvent rpree;
	RPRecoElasticEvent::road_type rpreert;
	std::vector<RPRecoElasticEvent::road_type> vrpreert;
	RPRecoElasticEvent::fit_type rpreerft;
	std::vector<RPRecoElasticEvent::fit_type> vrpreerft;
	edm::Wrapper<RPRecoElasticEvent> wrpree;
	
    RPRecognizedPatterns rprp;
    RPRecognizedPatterns::Line rprpl;
	std::vector<RPRecognizedPatterns::Line> vrprpl;
    RPRecognizedPatternsCollection rprpc;
	std::map<unsigned int, RPRecognizedPatterns> muirprp;
    edm::Wrapper<RPRecognizedPatterns> wrprp;
    edm::Wrapper<RPRecognizedPatternsCollection> wrprpc;
    
	std::allocator<std::pair<const unsigned int, RPFittedTrack> > aslpairinintrpfittrack;

	RPStationTrackFit rpstf;
	std::vector<RPStationTrackFit> vrpstf;
	std::map<unsigned int, std::vector<RPStationTrackFit> > m_ui_vrpstf;
	RPStationTrackFitCollection rpstfc;
	edm::Wrapper<RPStationTrackFitCollection> w_rpstfc;

  }
}
