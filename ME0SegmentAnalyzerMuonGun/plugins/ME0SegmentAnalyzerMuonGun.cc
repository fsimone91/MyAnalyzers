//Package:     ME0SegmentAnalyzerMuonGun
// Class:       ME0SegmentAnalyzerMuonGun
//
/**\class  ME0SegmentAnalyzerMuonGun  ME0SegmentAnalyzerMuonGun  ME0SegmentAnalyzerMuonGun/ ME0SegmentAnalyzerMuonGun/plugins/ ME0SegmentAnalyzerMuonGun.cc
 
 Description: [one line class summary]
 
 Implementation:
 [Notes on implementation]
 */
//
// Original Author:
//         Created:  Fri, 09 Oct 2015 10:34:53 GMT
// $Id$
//
//


// system include files
#include <memory>
#include <algorithm>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoMuon/MuonIdentification/interface/ME0MuonSelector.h"
#include <DataFormats/MuonReco/interface/ME0Muon.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/ME0MuonCollection.h>
//#include <DataFormats/MuonReco/interface/MuonME0Hits.h>

#include <DataFormats/GEMRecHit/interface/ME0RecHit.h>
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/GEMRecHit/interface/ME0RecHitCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/GeometrySurface/interface/LocalError.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//Geom
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include <Geometry/GEMGeometry/interface/ME0EtaPartition.h>
#include <DataFormats/MuonDetId/interface/ME0DetId.h>

#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"


#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/GEMDigi/interface/ME0DigiPreRecoCollection.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoMuon/MuonIdentification/plugins/ME0MuonSelector.cc"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include <memory>
#include <vector>
#include <cmath>
#include "TLorentzVector.h"

//
// class declaration
//
using namespace std;
using namespace edm;
class  ME0SegmentAnalyzerMuonGun : public edm::EDAnalyzer {
public:
    explicit  ME0SegmentAnalyzerMuonGun(const edm::ParameterSet&);
    ~ ME0SegmentAnalyzerMuonGun();
    void Initialize();
    
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 //   const edm::PSimHitContainer isTrackMatched(edm::PSimHitContainer::const_iterator, const edm::Event&,const edm::EventSetup&);
    
private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    edm::InputTag me0SegmentInputTag_;
    edm::InputTag me0DigiInputTag_;
    std::string wp_;
    double timeMin_;
    double timeMax_;
    double minEta_;
    double maxEta_;
    
  /*  InputTagToken_ = consumes<edm::PSimHitContainer>(cfg.getParameter<edm::InputTag>("simInputLabel"));
    InputTagToken_Digi = consumes<ME0DigiPreRecoCollection>(cfg.getParameter<edm::InputTag>("digiInputLabel"));
    InputTagToken_RecHit = consumes<ME0RecHitCollection>(cfg.getParameter<edm::InputTag>("recHitInputLabel"));*/
   
    edm::Service<TFileService> fs;
    
    TH1F *hNEvZmm;  TH1F * hDPhi; TH1F *  hme0machtMuonPt;  TH1F *  hme0machtMuonEta;  TH1F *  hme0machtMuonPhi;  TH1F *  hme0machtMuonCharge;
    TH1F *  hNME0Time ;  TH1F *  hNME0RecHits; TH1F *  hPtRes ;
  TH1F *  hSimPt; TH1F *  hSimEta;  TH1F *  hSimPhi; TH1F *  hSimPtDen; TH1F *  hSimEtaDen;  TH1F *  hSimPhiDen;
    TH1F *  hNzmm;  TH1F *  hNEv;  TH1F *hNGenMu;  TH1F * hNME0Mu;  TH1F *   hNMatchME0Mu;
    
    TH2F *  hPtVSDphi;TH2F *  hPtVSDEta; TH2F *  hPVSDphi ;  TH2F *   hPtVSDX;  TH2F *   hPtVSDY;  TH2F *   hPtVSDXLocal;  TH2F *   hPtVSDYLocal; TH2F *  hPtResVSDPhi;TH2F *   hSimPtVSDphi;
    
    
    TH1F *   hGenMuPt; TH1F *   hGenMuEta;
    TH1F *   hNumTight_Pt; TH1F *   hNumTight_Eta; TH1F *   hNumLoose_Pt; TH1F *   hNumLoose_Eta; TH1F *   hNumNoID_Eta; TH1F *   hNumNoID_Pt;
    TH1F * hAbsDPhi ; TH1F * hAbsDEta;  TH1F * hDEta;TH1F * hDX;
    TH1F *   hGenMass;  TH1F *   hRecoMass; TH1F *   hRecoMassIntime;	TH1F * 	hRecoMass_matchID;
    
    TH2F * hSimPtVSDeta;
    TH1F *   hSimDEta;
    TH1F *   hSimDPhi; TH1F *   hSimDX; TH1F *   hSimDY;
    
    TH1F * hSelectedSimTrack;
    TH1F * hGenMuonsME0;
    TH1F * hME0MuonsID; TH1F * hME0MuonsInMatchCone;
    
    TH1F * hNSimHitsME0;
    TH1F * hNRecHitsME0;
    TH1F * hRatioRecoToSimHits;
    TH1F * hNME0Segment;
    TH1F * hNBkgHitsME0;
    TH1F * hRecHitTime; TH1F *  hME0Segm_SimpT_Signal ;
    TH1F * hNDigiMatchedRH;	TH1F * hNPrompt;	TH1F * hN_noPrompt; TH1F * 	hNPromptMuHit;	TH1F * 	hNPromptNoMuHit;	TH1F * 	hNPromptHit_pdgId;	TH1F *  hME0SimTrackPdgID;
    TH1F *  hSimElePtME0;TH1F *   hSimMuPtME0;	TH1F *   hSimPionPtME0;	TH1F * 	hSimEleNHitsME0;	TH1F * 	hSimMuonNHitsME0;	TH1F * 	hSimPionNHitsME0; TH1F * hPdgIDCheck;
    TH2F *  hMuonDigiDPhiVsPT; TH2F * hNoEleDigiDPhiVsPT;
    
    TH1F * hDRME0SimTrack;	TH1F * hDRME0SimMuonEle;	TH1F * hN_noEleHit;	TH1F * hN_noPromptHit_pdgId;	TH1F *	hSegmentComposition;	TH1F * hNME0SegmentAfterDigiMatch;	TH1F * hMuEleinME0Segm;	TH1F * hMuEleOthersinME0Segm;TH1F * hMuOnlyinME0Segm;
    
    TH1F * hNEleBrem;	TH1F * hNEleDeltaRays;	TH1F * hNEle; TH1F * hMatchedSimEleTrack;
    
    TH1F * hRHDeltaPhiSameLayer;	TH1F * hRHDeltaEtaSameLayer;	TH1F * hRHDeltaTSameLayer;
    
    TH1F * hNoPromptRecHitTime;	TH1F * hMuonRecHitTime;	TH1F * hEleRecHitTime;	TH1F * hNMuonSameLayerTOF;	TH1F * hNEleSameLayerTOF;	TH1F * hNoPromptSameLayerTOF;
    
    TH1F * hDeltaPhiSimReco;	TH1F * hDeltaEtaSimReco;TH1F * hDeltaXSimReco; TH1F * hDeltaYSimReco;
    TH1F * hNoMuinME0Segm; TH2F * hMuonDigiDXVsPT;	TH2F * hMuonDigiLocalDPhiVsPT; 	TH2F * hMuonDigiLocalDXVsPT; TH1F * hDeltaXSimRecoLocal; TH1F * hDeltaYSimRecoLocal;
    
    TH1F * 	hBeamSpotX0; 	TH1F * 	hBeamSpotY0; 	TH2F * 	hBeamSpotX0Y0;	TH1F * 	hBeamSpotZ0; 	TH1F * 	hBeamSpotBkgmaZ;
    TH1F * 	hBeamSpotdxdz;	TH1F * 	hBeamSpotBeamWidthX;	TH1F * 	hBeamSpotBeamWidthY;	TH1F * 	hVertexMult ; 	TH1F * 	hverteX;	TH1F * 	hverteY;TH2F * 	hverteXY; TH1F * 	hverteZ ;
    
    TH2F * 	hRecDPhiVSimDphi;TH1F * hDiffRecDPhiVSimDphi; TH1F *  hRecDPhiOverSimDphi;
    
    TH1F *	hDeltaXSimRecoLocal_1;	TH1F *	hDeltaXSimRecoLocal_2;	TH1F *	hDeltaXSimRecoLocal_3;
    
    TH1F * hDeltaPhiSimReco_1; TH1F * hDeltaPhiSimReco_2; TH1F * hDeltaPhiSimReco_3;
    TH1F * hSimDPhiPos; TH1F * hSimDPhiNeg;
    TH1F * hDPhiPos; TH1F * hDPhiNeg;
    TH1F * hME0MuonsChargeNeg;
    TH1F * hME0MuonsChargePos;
    TH1F * hDqOverDphi;
    TH1F * hQSimQRecoOverDphi;
    TH1F * hGenMuNegEta;
    TH1F * hGenMuNegPhi;
    TH1F * hGenMuPosEta;
    TH1F * hGenMuPosPhi;
    
    TH1F * hSimDPhiPos_HighPt;
    TH1F * hSimDPhiPos_LowPt;
    
    TH1F * hSimDPhiNeg_HighPt;
    TH1F * hSimDPhiNeg_LowPt;
    TH1F * hDPhiPos_LowPt;
    TH1F * hDPhiPos_HighPt;
    
    TH1F * hDPhiNeg_LowPt;
    TH1F * hDPhiNeg_HighPt;
    
    TH2F * hPosMuonDigiDPhiVsPT;
    TH2F * hNegMuonDigiDPhiVsPT;
    TH1F * hDeltaPhiSimReco_pos;
    TH1F * hDeltaPhiSimReco_neg;
    
    TH1F * 	hSimLocalDX;	TH1F * 	hLocalDX;TH1F * 	hSimLocalDY;	TH1F * 	hLocalDY;
    TH1F * 	hSimLocalDXPos;	TH1F * 	hLocalDXPos;TH1F * 	hSimLocalDXNeg;TH1F * hLocalDXNeg;
    TH1F * 	hSimLocalDXPos_LowPt;	TH1F * 	hLocalDXPos_LowPt;TH1F * 	hSimLocalDXNeg_LowPt;TH1F * hLocalDXNeg_LowPt;
    TH1F * 	hSimLocalDXPos_HighPt;	TH1F * 	hLocalDXPos_HighPt;TH1F * 	hSimLocalDXNeg_HighPt;TH1F * hLocalDXNeg_HighPt;
    
    
    TH1F * 	hSimLocalDXPos_ME0Plus;TH1F * 	hLocalDXPos_ME0Plus;
    
    TH1F * 	hSimLocalDXNeg_ME0Plus;	TH1F * 	hLocalDXNeg_ME0Plus;
    
    TH1F * 	hSimLocalDXPos_LowPt_ME0Plus;	TH1F * 	hLocalDXPos_LowPt_ME0Plus;
    
    TH1F * 	hSimLocalDXNeg_LowPt_ME0Plus;	TH1F * 	hLocalDXNeg_LowPt_ME0Plus;
    
    TH1F * 	hSimLocalDXPos_HighPt_ME0Plus;	TH1F * 	hLocalDXPos_HighPt_ME0Plus;
    
    TH1F * 	hSimLocalDXNeg_HighPt_ME0Plus;	TH1F * 	hLocalDXNeg_HighPt_ME0Plus;
    
    TH1F * 	hSimDPhiPos_ME0Plus ;TH1F * 	hSimDPhiNeg_ME0Plus;
    TH1F * 	hSimDPhiPos_LowPt_ME0Plus; 	TH1F * 	hSimDPhiPos_HighPt_ME0Plus;
    TH1F * 	hSimDPhiNeg_LowPt_ME0Plus ;	TH1F * 	hSimDPhiNeg_HighPt_ME0Plus;
    
    TH1F * 	hDPhiPos_ME0Plus; 	TH1F * 	hDPhiNeg_ME0Plus ;
    TH1F * 	hDPhiPos_LowPt_ME0Plus; 	TH1F * 	hDPhiPos_HighPt_ME0Plus;
    TH1F * 	hDPhiNeg_LowPt_ME0Plus ; 	TH1F * 	hDPhiNeg_HighPt_ME0Plus;
    
    
    /////////////////////
    
    TH1F * 	hSimLocalDXPos_ME0Minus;TH1F * 	hLocalDXPos_ME0Minus;
    
    TH1F * 	hSimLocalDXNeg_ME0Minus;	TH1F * 	hLocalDXNeg_ME0Minus;
    
    TH1F * 	hSimLocalDXPos_LowPt_ME0Minus;	TH1F * 	hLocalDXPos_LowPt_ME0Minus;
    
    TH1F * 	hSimLocalDXNeg_LowPt_ME0Minus;	TH1F * 	hLocalDXNeg_LowPt_ME0Minus;
    
    TH1F * 	hSimLocalDXPos_HighPt_ME0Minus;	TH1F * 	hLocalDXPos_HighPt_ME0Minus;
    
    TH1F * 	hSimLocalDXNeg_HighPt_ME0Minus;	TH1F * 	hLocalDXNeg_HighPt_ME0Minus;
    
    TH1F * 	hSimDPhiPos_ME0Minus ;TH1F * 	hSimDPhiNeg_ME0Minus;
    TH1F * 	hSimDPhiPos_LowPt_ME0Minus; 	TH1F * 	hSimDPhiPos_HighPt_ME0Minus;
    TH1F * 	hSimDPhiNeg_LowPt_ME0Minus ;	TH1F * 	hSimDPhiNeg_HighPt_ME0Minus;
    
    TH1F * 	hDPhiPos_ME0Minus; 	TH1F * 	hDPhiNeg_ME0Minus ;
    TH1F * 	hDPhiPos_LowPt_ME0Minus; 	TH1F * 	hDPhiPos_HighPt_ME0Minus;
    TH1F * 	hDPhiNeg_LowPt_ME0Minus ; 	TH1F * 	hDPhiNeg_HighPt_ME0Minus;
    
    TH1F * 		hLocalDX_over_R;TH1F * 		hSimLocalDXoverR;
    TH2F * hverteZT;  TH1F *  hverteTime; TH1F *  hBX; TH1F *	hNPU; TH1F *	hTrueInt;
    TH1F * 	hDPhiNoMatch;  TH1F * 	hDPhiNoMatch_TimeWindow; TH1F * 	hDPhiNoMatch_TimeWindowTightID;
    
    TH1F * 	hTightME0SegmDigiHitTOF_noprompt; TH1F * 	hTightME0SegmDigiHitPdgID_noprompt;
    TH1F * 	hTightME0SegmDigiHitTOF_prompt; TH1F * 	hTightME0SegmDigiHitPdgID_prompt;    TH1F * hAllME0SegmentTime;  TH1F * hAllME0SegmentTimeErr;
    TH1F * hME0SegmentCollectionTime;
    
    TH1F * 	hDigiHitType;TH1F * 	hDigiHitPdgID;TH1F * 	hDigiHitToF_prompt;TH1F * 	hDigiHitToF_noprompt;
    
    TH1F *  hMuonInME0SegTOF;	TH1F *  hEleInME0SegTOF;	TH1F *  hPionInME0SegTOF;	TH1F *  hBoInME0SegTOF;	TH1F *  hProtonsInME0SegTOF;  TH1F * hBeamSpotSigmaZ;
    
    
    TH1F *  hME0SegmDigiHitME0Segm_prompt;	TH1F *  hME0SegmDigiHitPdgIDME0Segm_prompt;	TH1F *  hME0SegmDigiHitTOFME0Segm_noprompt;	TH1F *  hME0SegmDigiHitPdgIDME0Segm_noprompt; TH1F *  hDPhiNoMatchME0Seg;	TH1F *  hDPhiNoMatchME0Seg_TimeWindow;
    
    TH1F *  hMuonInME0SegTOF_L1;	TH1F *  hEleInME0SegTOF_L1;	TH1F *  hPionInME0SegTOF_L1;	TH1F *  hBoInME0SegTOF_L1;	TH1F *  hProtonsInME0SegTOF_L1;
    
    TH1F *  hMuonInME0SegTOF_inTime;	TH1F *  hEleInME0SegTOF_inTime;	TH1F *  hPionInME0SegTOF_inTime;	TH1F *  hBoInME0SegTOF_inTime;	TH1F *  hProtonsInME0SegTOF_inTime;
    
    TH1F * hDPhiMatchByHitsME0Seg_Signal; TH1F * hDPhiMatchByHitsME0Seg_Bkg;	TH1F * hTimeMatchByHitsME0Seg_Signal;	TH1F * hTimeMatchByHitsME0Seg_Bkg;
    TH2F * hDPhivsTOFMatchByHitsME0Seg_Bkg; TH2F * hDPhivsTOFMatchByHitsME0Seg_Signal;
    TH1F * 	hNME0SegmentSignal;TH1F * 	hNME0SegmentBkg;TH1F * 	hTimeDiffSigvsBkg;
    
    TH1F * 	hDPhiMatchByHitsME0Seg_Signal_timeCut;	TH1F * 	hDPhiMatchByHitsME0Seg_Bkg_timeCut	;	TH2F * 	hDPhivsTOFMatchByHitsME0Seg_Bkg_timeCut;	TH2F * 	hDPhivsTOFMatchByHitsME0Seg_Signal_timeCut;		TH1F * 	hNME0SegmentSignal_timeCut;	TH1F * 	hNME0SegmentBkg_timeCut;
    TH1F *  hSimMuonPt_DPhiCut30GeVAllCut;   TH1F *  hSimMuonPt_DPhiCut5GeVAllCut;   TH1F *  hEventPass5GeVCutAllCut;TH1F *  hEventPass30GeVCutAllCut;
    
    
    
    TH1F *	hSimMuonPt_DPhiCut5GeV ;	TH1F *	hSimMuonPt_DPhiCut30GeV;	TH1F *	hEventPass30GeVCut;	TH1F *	hEventPass5GeVCut;	TH1F *  hNEvME0;
    TH1F *	hSimVertZPosition; TH2F * hSimVertXYPosition; TH1F * hSimME0SegmentTime;TH2F * hSimVertZPositionVsSimSegmTime;
    
    TH1F *  hME0SegCompMuonOnly;TH1F * hME0SegCompMuonPlusEle;TH1F *hME0SegCompMuonPlusHad;TH1F * hME0SegCompMuonPlusElePlusHad;TH1F * hME0SegCompHadOnly;TH1F * hME0SegCompEleOnly;TH1F * hME0SegmNumber;
    
    TH1F * hME0SegCompMuonOnly_Sig; TH1F * hME0SegCompMuonPlusEle_Sig; TH1F *  hME0SegCompMuonPlusHad_Sig; TH1F * hME0SegCompMuonPlusElePlusHad_Sig; TH1F * hME0SegCompEleOnly_Sig;    TH1F * hME0SegCompHadOnly_Sig; TH1F * hME0SegmNumberSig;
    
    TH1F * hME0SegCompMuonOnly_Bkg; TH1F * hME0SegCompMuonPlusEle_Bkg; TH1F * hME0SegCompMuonPlusHad_Bkg;TH1F * hME0SegCompMuonPlusElePlusHad_Bkg; TH1F * hME0SegCompEleOnly_Bkg;TH1F * hME0SegCompHadOnly_Bkg; TH1F * hME0SegmNumberBkg;
    
    TH1F *  hPull_MyRecoSim_zVtx;   TH1F *  hPull_4DRecoSim_zVtx;
    TH1F * hME0Segm_nHit_bkg;    TH1F * hME0Segm_chi2_bkg;TH1F * hME0Segm_time_bkg;
    TH1F *  hME0Segm_nHit_sig;   TH1F *  hME0Segm_chi2_sig;   TH1F *  hME0Segm_time_sig;TH1F *  hME0Segm_nHit;   TH1F *  hME0Segm_chi2; TH1F * hME0SegmMyTOF;
    TH2F *  hME0Segm_DeltaPhiVsSimpT_Signal;
    TH2F *  hME0Segm_DeltaPhiVsSimpT_Signal_rebin;
    TH2F *  hME0Segm_DeltaPhiVsSimpT_Bkg;
    
    TH1F *   hME0SegmDigiHitTOF  ;
    TH1F *    hME0SegmDigiHitPdgID  ;
    TH1F *    hME0SegmDigiHitEta   ;
    
    TH1F *   hME0SegmDigiHitTOF_signal  ;
    TH1F *    hME0SegmDigiHitPdgID_signal  ;
    TH1F *    hME0SegmDigiHitEta_signal   ;
    
    TH1F *   hME0SegmDigiHitTOF_bkg;    TH1F *    hME0SegmDigiHitPdgID_bkg  ;
    TH1F *    hME0SegmDigiHitEta_bkg   ;TH1F *   hME0Segm_eta_bkg; TH1F *   hME0Segm_eta_sig;
    
    
  TH1F *  hME0Segm_phi_sig;   TH1F *  hME0Segm_X_sig;   TH1F *  hME0Segm_Y_sig;   TH1F *  hME0Segm_eta_sigCheck; TH1F *  hME0Segm_eta_sigCheckFail;
  TH1F *    hME0Segm_eta_bkgCheck;  TH1F *    hME0Segm_eta_bkgCheckFail;TH1F *    hME0Segm_eta;
  TH1F *    hDistME0Seg; TH1F *     hSegmentPerChamber; TH1F *     hSegmentPerChamber2;
  TH1F *  hME0Segm_DeltaPhi_pT0p5; TH1F *  hME0Segm_DeltaPhi_pT1;
  TH1F *  hME0Segm_DeltaPhi_pT3; TH1F *  hME0Segm_DeltaPhi_pT5; TH1F *  hME0Segm_DeltaPhi_pT10; TH1F *  hME0Segm_DeltaPhi_pT20;
  TH1F *  hME0Segm_DeltaPhi_pT50;
  TH1F *  hNumberOfMatchedHit; TH1F *    hME0SegmPullPhi; TH1F *    hME0SegmPullEta;
  TH1F *    hME0SegmPullPhi_fitFail; TH1F *    hME0SegmPullEta_fitFail;
  TH1F *    hME0SegmPullPhi_BestMatch;  TH1F *    hME0SegmPullEta_BestMatch;TH1F * hNumberOfDuplicatesPerSignalSimTrack;
    
   TH1F *  hNME0SegmentSignal_BestMatch;   TH1F *  hME0Segm_nHit_sig_BestMatch;
   TH1F *  hME0Segm_chi2_sig_BestMatch;   TH1F *  hME0Segm_time_sig_BestMatch;
    
   TH1F *  hnumAllME0SimHits;   TH1F *  hnumEleME0SimHits;   TH1F *  hnumMuonME0SimHits;   TH1F *  hnumPhotonME0SimHit;   TH1F *  hnumHadronME0SimHits;
   TH1F *  hEleSimHitEta;   TH1F *  hMuonSimHitEta;   TH1F *  hPhotonSimHitEta;   TH1F *  hHadronSimHitEta;TH1F *   hSimHitEta;
   TH1F *  hDuplicatesEta;   TH1F *   hDuplicatesPhi;
   TH1F *  hDeltaRoll;   TH2F *  hDeltaRollvsDeltaEta;


    // virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    


  //  edm::EDGetTokenT< std::vector<reco::ME0Muon>> OurMuons_;
    edm::EDGetTokenT<SimTrackContainer> simTracks_;
    edm::EDGetTokenT<ME0SegmentCollection> me0Segments_;
    edm::EDGetTokenT<ME0SegmentCollection> me0SegmentsRU_;
    edm::EDGetTokenT<reco::GenParticleCollection> genparticles_;
    edm::EDGetTokenT< std::vector<reco::Muon> > muons_;
    edm::EDGetTokenT< ME0DigiPreRecoCollection > me0_digis_;
    edm::EDGetTokenT< reco::BeamSpot > beamSpotHandle_;
    edm::EDGetTokenT< reco::VertexCollection > primaryVertices_;
    edm::EDGetTokenT< vector<PSimHit>   > ME0HitsCollection_;


    
  //  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
//    iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ME0SegmentAnalyzerMuonGun:: ME0SegmentAnalyzerMuonGun(const edm::ParameterSet& iConfig):
wp_(iConfig.getParameter<std::string>("wp")),
timeMin_(iConfig.getParameter<double>("timeMin")),
timeMax_(iConfig.getParameter<double>("timeMax")),
minEta_(iConfig.getParameter<double>("minEta")),
maxEta_(iConfig.getParameter<double>("maxEta")),
me0Segments_(consumes<ME0SegmentCollection>(iConfig.getParameter< edm::InputTag >("me0SegmentInputTag"))),
me0_digis_(consumes<ME0DigiPreRecoCollection>(iConfig.getParameter< edm::InputTag >("me0DigiInputTag")))
{
   // std::string labelMe0Muon_("me0SegmentMatching");
 //   OurMuons_ = consumes< std::vector<reco::ME0Muon> >(edm::InputTag("me0SegmentMatching"));
    simTracks_ = consumes< SimTrackContainer >(edm::InputTag("g4SimHits"));
    
    //me0Segments_ = consumes<ME0SegmentCollection  >(edm::InputTag("me0SegmentsPerfect","", "RECO"));
    //me0Segments_ = consumes<ME0SegmentCollection  >(edm::InputTag("me0SegmentsPerfectRU","", "RECO"));
    // me0Segments_ = consumes<ME0SegmentCollection  >(edm::InputTag("me0SegmentsSmeared","", "RECO"));
    //me0Segments_ = consumes<ME0SegmentCollection  >(edm::InputTag("me0SegmentsSmearedRU","", "RECO"));
    //vtx_Token(consumes<reco::VertexCollection>(parset.getParameter< edm::InputTag >("vtxTag"))),
    
    me0SegmentsRU_ = consumes<ME0SegmentCollection  >(edm::InputTag("me0SegmentsPerfectRU", "", "RECO"));
    //me0SegmentsRU_ = consumes<ME0SegmentCollection  >(edm::InputTag("me0SegmentsSmearedRU", "", "RECO"));
   
    
    genparticles_ = consumes< reco::GenParticleCollection >(edm::InputTag("genParticles")),
    muons_ = consumes< std::vector<reco::Muon> >(edm::InputTag("muons"));
    
    //me0_digis_ = consumes< ME0DigiPreRecoCollection >(edm::InputTag("simMuonME0Digis"));
    //me0_digis_ = consumes< ME0DigiPreRecoCollection >(edm::InputTag("simMuonME0ReDigis"));
    
    beamSpotHandle_ = consumes< reco::BeamSpot  >(edm::InputTag("offlineBeamSpot"));
    ME0HitsCollection_ = consumes< vector<PSimHit>  >(edm::InputTag("g4SimHits", "MuonME0Hits"));
    primaryVertices_ = consumes< reco::VertexCollection  >(edm::InputTag("primaryVertices"));


    /*ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
    iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
    edm::Handle<SimVertexContainer> simVertexCollection;
    iEvent.getByLabel("g4SimHits", simVertexCollection);
    const SimVertexContainer simVC = *(simVertexCollection.product());*/

}


ME0SegmentAnalyzerMuonGun::~ ME0SegmentAnalyzerMuonGun()
{
    
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
}

bool isSimMatched(edm::SimTrackContainer::const_iterator simTrack, edm::PSimHitContainer::const_iterator itHit)
{
    
    bool result = false;
    int trackId = simTrack->trackId();
    int trackId_sim = itHit->trackId();
    if(trackId == trackId_sim) result = true;
    //std::cout<<"ID: "<<trackId<<" "<<trackId_sim<<" "<<result<<std::endl;
    return result;
    
}

/*edm::PSimHitContainer isTrackMatched(SimTrackContainer::const_iterator simTrack, const Event & event, const EventSetup& eventSetup)
{
    edm::PSimHitContainer selectedME0Hits;
    
    edm::ESHandle<ME0Geometry> me0geom;
    eventSetup.get<MuonGeometryRecord>().get(me0geom);
    
    //edm::Handle<edm::PSimHitContainer> ME0Hits;
    // event.getByLabel(edm::InputTag("g4SimHits","MuonME0Hits"), ME0Hits);
    //   edm::EDGetTokenT< edm::PSimHitContainer  > ME0HitsCollection_;
       //ME0HitsCollection_= consumes< edm::PSimHitContainer >(edm::InputTag("ME0HitsCollection"));
    //  ME0HitsCollection_ = mayConsume<edm::PSimHitContainer >(edm::InputTag("ME0HitsCollection"));
    Handle<edm::PSimHitContainer> ME0HitsCollection;
    event.getByToken(consumes< edm::PSimHitContainer >edm::InputTag("ME0HitsCollection"),  ME0HitsCollection);
    
    
    ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
    eventSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
    
    for (edm::PSimHitContainer::const_iterator itHit =  ME0HitsCollection->begin(); itHit != ME0HitsCollection->end(); ++itHit){
        
        DetId id = DetId(itHit->detUnitId());
        if (!(id.subdetId() == MuonSubdetId::ME0)) continue;
        if(itHit->particleType() != (*simTrack).type()) continue;
        
        bool result = isSimMatched(simTrack, itHit);
        if(result) selectedME0Hits.push_back(*itHit);
        
    }
    //	std::cout<<"N simHit in ME0segm : "<<selectedME0Hits.size()<<std::endl;
    return selectedME0Hits;
    
}
*/


struct MyME0Digi
{
    Int_t detId, particleType;
    Short_t layer;
    Float_t g_eta, g_phi;
    Float_t tof;
    Float_t prompt;
};

void
ME0SegmentAnalyzerMuonGun::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //Initialize();
    using namespace edm;

  //  Handle<std::vector<reco::ME0Muon>> OurMuons;
  //  iEvent.getByToken( OurMuons_,  OurMuons);
    
    
    Handle<ME0SegmentCollection> me0Segments;
    iEvent.getByToken(me0Segments_,  me0Segments);
    
    Handle<ME0SegmentCollection> me0SegmentsRU;
    iEvent.getByToken(me0SegmentsRU_,  me0SegmentsRU);

   
    Handle<SimTrackContainer> simTracks;
    iEvent.getByToken(simTracks_,  simTracks);
    
    
    Handle<std::vector<reco::Muon>> muons;
    iEvent.getByToken(muons_,  muons);
    
    Handle<ME0DigiPreRecoCollection> me0_digis;
    iEvent.getByToken(me0_digis_,  me0_digis);
    
    Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(beamSpotHandle_,  beamSpotHandle);
    
    Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(primaryVertices_,  primaryVertices);
    
    Handle<edm::PSimHitContainer> ME0HitsCollection;
    iEvent.getByToken(ME0HitsCollection_,  ME0HitsCollection);
    
    Handle<reco::GenParticleCollection > genparticles;
    iEvent.getByToken(genparticles_,genparticles);
    
    

/*
    
    ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
    iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
    
    edm::Handle<SimVertexContainer> simVertexCollection;
    iEvent.getByLabel("g4SimHits", simVertexCollection);
    const SimVertexContainer simVC = *(simVertexCollection.product());
*/
    
    edm::ESHandle<ME0Geometry> me0geom;
    iSetup.get<MuonGeometryRecord>().get(me0geom);
    const ME0Geometry* me0Geom;
    me0Geom= &*me0geom;

    double c=29.97; //speed of light in cm/ns
    
    /*
     Run   = iEvent.id().run();
     Event = iEvent.id().event();
     Lumi  = iEvent.luminosityBlock();
     Bunch = iEvent.bunchCrossing();
     */
    
    std::cout<<"*********************************BeginEvent="<<iEvent.id().event()<<"******************************************** "<<std::endl;
    std::cout<</*" SizeME0SimHits="<< ME0HitsCollection->size()<<*/" ME0SegmentSize="<<me0Segments->size()<<" RUME0SegmentSize="<<me0SegmentsRU->size()<<std::endl;
    std::vector<int> indexgenmu;
    hNEv->Fill(1);
    SimTrackContainer::const_iterator simTrack;
    double numberOfSimTracks =0.;
    double nSegmPerChamber=0.;    
    //std::cout<<" num Simulated tracks: "<<simTracks->size()<<std::endl;
    std::vector<double> simHitPhi, simHitElePhi, simHitX, simHitLocalX;
    std::vector<double> simHitEta, simHitEleEta, simHitY, simHitLocalY;
    std::vector<double> simHitR, simHitTOF;
    std::vector<int> numME0SimHits;
     std::vector<int> numMuonME0SimHits, numEleME0SimHits, numPhotonME0SimHits, numHadronME0SimHits, numAllME0SimHits;
    SimTrackContainer ME0Tracks, ME0EleSimTracks,AllME0Tracks;
    edm::OwnVector<ME0Segment>  ME0SegBkg;
    edm::OwnVector<ME0Segment>  ME0SegSig;
    ME0SegBkg.clear();
    edm::PSimHitContainer selME0SimHits, allME0SimHits;
    double SimVtxZ=0.;double SimVtxX=0.;double SimVtxY=0.;
    
    for (edm::PSimHitContainer::const_iterator itHit =  ME0HitsCollection->begin(); itHit != ME0HitsCollection->end(); ++itHit){
        DetId id = DetId(itHit->detUnitId());
        if (!(id.subdetId() == MuonSubdetId::ME0)) continue;
        numAllME0SimHits.push_back(1);
        LocalPoint lp = itHit->entryPoint();
        GlobalPoint hitGP_sim( me0geom->idToDet(itHit->detUnitId())->surface().toGlobal(lp));
        hSimHitEta->Fill(hitGP_sim.eta());
        if(fabs(itHit->particleType())==11) {hEleSimHitEta->Fill(hitGP_sim.eta()); numEleME0SimHits.push_back(1);}
        if(fabs(itHit->particleType())==13) {hMuonSimHitEta->Fill(hitGP_sim.eta());numMuonME0SimHits.push_back(1);}
        if(fabs(itHit->particleType())==22) {hPhotonSimHitEta->Fill(hitGP_sim.eta()); numPhotonME0SimHits.push_back(1);}
        if(fabs(itHit->particleType())>200) {hHadronSimHitEta->Fill(hitGP_sim.eta()); numHadronME0SimHits.push_back(1);}
    }
    hnumAllME0SimHits->Fill(numAllME0SimHits.size());
    hnumEleME0SimHits->Fill(numEleME0SimHits.size());
    hnumMuonME0SimHits->Fill(numMuonME0SimHits.size());
    hnumPhotonME0SimHit->Fill(numPhotonME0SimHits.size());
    hnumHadronME0SimHits->Fill(numHadronME0SimHits.size());

    for (simTrack = simTracks->begin(); simTrack != simTracks->end(); ++simTrack){
        
        double simPt=(*simTrack).momentum().pt();
        double simEta = (*simTrack).momentum().eta();
        
        if (abs(simEta) > maxEta_ || abs(simEta) < minEta_) continue;

        //cout<<"Sim pT: "<<(*simTrack).momentum().pt()<<" carica="<<(*simTrack).charge()<<" pdgID="<<(*simTrack).type()<<endl;
	allME0SimHits.clear();

    for (edm::PSimHitContainer::const_iterator itHit =  ME0HitsCollection->begin(); itHit != ME0HitsCollection->end(); ++itHit){
	  DetId id = DetId(itHit->detUnitId());
	  if (!(id.subdetId() == MuonSubdetId::ME0)) continue;
	  if(itHit->particleType() != (*simTrack).type()) continue;
	  bool result = isSimMatched(simTrack, itHit);
	  if(result) allME0SimHits.push_back(*itHit);
        }

	if( allME0SimHits.size() !=0 ) {
	  AllME0Tracks.push_back(*simTrack);
	  hME0SimTrackPdgID->Fill( fabs((*simTrack).type()) );
	  if (fabs((*simTrack).type())==11)	{hSimElePtME0->Fill(simPt);	 ME0EleSimTracks.push_back(*simTrack);}
	  if (fabs((*simTrack).type())==13)	{hSimMuPtME0->Fill(simPt); }
	  if (fabs((*simTrack).type())==211)	{hSimPionPtME0->Fill(simPt); }
	}
        //if ((*simTrack).noVertex()) continue;
        if ((*simTrack).noGenpart()) continue;
        if (!(abs((*simTrack).type()) == 13)) continue;
        
        selME0SimHits.clear();
        
        for (edm::PSimHitContainer::const_iterator itHit =  ME0HitsCollection->begin(); itHit != ME0HitsCollection->end(); ++itHit){
            
            DetId id = DetId(itHit->detUnitId());
            if (!(id.subdetId() == MuonSubdetId::ME0)) continue;
            if(itHit->particleType() != (*simTrack).type()) continue;
            
            bool result = isSimMatched(simTrack, itHit);
            if(result) selME0SimHits.push_back(*itHit);
            
        }
        //	std::cout<<"N simHit in ME0segm : "<<selectedME0Hits.size()<<std::endl;
        
        
        //selME0SimHits=isTrackMatched(simTrack, iEvent , iSetup);
        //edm::PSimHitContainer selME0SimHits = isTrackMatched(simTrack, iEvent , iSetup);
        
        //int ME0SimHitsize = selME0SimHits.size();
        //std::cout<<"# me0 sim hits="<< selME0SimHits.size() <<" SimTrack vertex index "<<simTrack->vertIndex()<<std::endl;
        //uint TrkVtxID=simTrack->vertIndex();
        
        
        if( selME0SimHits.size() ==0 ) continue;
        hSimEta->Fill((*simTrack).momentum().eta());
        hSimPt->Fill((*simTrack).momentum().pt());
        hSimPhi->Fill((*simTrack).momentum().phi());
        ME0Tracks.push_back(*simTrack);
        numberOfSimTracks++;
        //std::cout<<"TrackID="<<simTrack->trackId()<<" track type="<<(*simTrack).type()<<"# me0 sim hits="<< selME0SimHits.size() <<" simPt="<<simPt<<std::endl;
        simHitPhi.clear(); simHitEta.clear();
        simHitX.clear();   simHitY.clear(); simHitLocalX.clear();  simHitLocalY.clear(); simHitR.clear();  simHitTOF.clear();
        int selhitcounter =0;
        
        for (edm::PSimHitContainer::const_iterator itHit =  selME0SimHits.begin(); itHit != selME0SimHits.end(); ++itHit){
            //ME0DetId idme0 = ME0DetId(itHit->detUnitId());
            //int layer_sim = idme0.layer();
            LocalPoint lp = itHit->entryPoint();
            GlobalPoint hitGP_sim( me0geom->idToDet(itHit->detUnitId())->surface().toGlobal(lp));
            simHitPhi.push_back(hitGP_sim.phi());
            simHitEta.push_back(hitGP_sim.eta());
            simHitX.push_back(hitGP_sim.x());
            simHitY.push_back(hitGP_sim.y());
            simHitLocalX.push_back(lp.x());
            simHitLocalY.push_back(lp.y());
            simHitR.push_back(hitGP_sim.perp());
            simHitTOF.push_back(itHit->tof());
            selhitcounter	++;
            //std::cout<<"track pt="<<(*simTrack).momentum().pt()<<" simHit eta="<<hitGP_sim.eta()<<" phi="<<hitGP_sim.phi()<<" simHit detID="<<idme0<<" layer sim="<<layer_sim<<" localX="<<lp.x()<<" trackID="<<itHit->trackId()<<" Q="<<(*simTrack).charge()<<" global point r="<<hitGP_sim.perp()<<" R="<< TMath::Sqrt( (hitGP_sim.x()*hitGP_sim.x()) + (hitGP_sim.y()*hitGP_sim.y()) )<<std::endl;
           // std::cout<<"Muon simHit eta="<<hitGP_sim.eta()<<" phi="<<hitGP_sim.phi()<<" layer sim="<<layer_sim<<" localX="<<lp.x()<<" localY="<<lp.y()<<" trackID="<<itHit->trackId()<<std::endl;
        }
        numME0SimHits.push_back(selhitcounter);
        int sizeSimPhi = simHitPhi.size()-1;
        double SimDeltaPhi = TMath::Abs( simHitPhi[0] - simHitPhi[sizeSimPhi] );
        double SimDeltaEta= TMath::Abs( simHitEta[0] - simHitEta[sizeSimPhi] );
        
        double sum=0.;
        for(uint k=0; k<simHitTOF.size(); k++){
            sum += simHitTOF.at(k);}
        
        double SimME0SegmentTime = sum/simHitTOF.size();
        //std::cout<<" sum="<<sum<<" mean="<<SimME0SegmentTime<<std::endl;
        hSimME0SegmentTime->Fill(SimME0SegmentTime);
        
        /*for(edm::SimVertexContainer::const_iterator v=simVertexCollection->begin(); v!=simVertexCollection->end(); ++v){
            //std::cout<<" -------sim vertex Time ---------"<<v->tof()<<std::endl;
            if(v->vertexId()==TrkVtxID )  {
                std::cout<<" simVtx pos along z= "<<v->position().z()<<std::endl;
                hSimVertZPosition->Fill(v->position().z());
                hSimVertXYPosition->Fill(v->position().x(), v->position().x());
                hSimVertZPositionVsSimSegmTime->Fill(v->position().z(), SimME0SegmentTime);
                SimVtxZ=v->position().z();
                SimVtxX=v->position().x();
                SimVtxY=v->position().y();
            }}*/
        
        cout<<"Selected Signal Track:"<<simTrack->trackId()<<" pdgID="<<(*simTrack).type()<<" SimDPhi="<<simHitPhi[0] - simHitPhi[sizeSimPhi]<<" pT="<<simPt<<" carica="<<(*simTrack).charge()<<" DX(l1,l6)"<< simHitLocalX[0] - simHitLocalX[sizeSimPhi]<<" DX/R"<< (simHitLocalX[0] - simHitLocalX[sizeSimPhi])/simHitR[0]<<endl;
        
        
        hSimPtVSDphi->Fill(SimDeltaPhi, simPt);
        hSimPtVSDeta->Fill(SimDeltaEta, simPt);
        
        
        hSimDEta->Fill( simHitEta[0] - simHitEta[sizeSimPhi]);
        hSimDPhi->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);
        hSimDX->Fill( simHitX[0] - simHitX[sizeSimPhi]);
        hSimDY->Fill( simHitY[0] - simHitX[sizeSimPhi]);
        
        hSimLocalDX->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);
        hSimLocalDY->Fill( simHitLocalY[0] - simHitLocalY[sizeSimPhi]);
        if (simHitR[0] !=0.) hSimLocalDXoverR->Fill( (simHitLocalX[0] - simHitLocalX[sizeSimPhi])/simHitR[0]);
        
        if ((*simTrack).charge() >0) {hSimDPhiPos->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]); hSimLocalDXPos->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
        if ((*simTrack).charge() <0) {hSimDPhiNeg->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]); hSimLocalDXNeg->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
        
        if (((*simTrack).charge() >0) && ((*simTrack).momentum().pt()) >20 ) {hSimDPhiPos_HighPt->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]); hSimLocalDXPos_HighPt->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
        if (((*simTrack).charge() >0) && ((*simTrack).momentum().pt()) <10 ) {hSimDPhiPos_LowPt->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);  hSimLocalDXPos_LowPt->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
        
        if (((*simTrack).charge() <0) && ((*simTrack).momentum().pt()) >20 ) {hSimDPhiNeg_HighPt->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);  hSimLocalDXNeg_HighPt->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
        if (((*simTrack).charge() <0) && ((*simTrack).momentum().pt()) <10 ) {hSimDPhiNeg_LowPt->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);   hSimLocalDXNeg_LowPt->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
        
        
        
        if((*simTrack).momentum().eta() > 0 ){
            
            
            if ((*simTrack).charge() >0) {hSimDPhiPos_ME0Plus->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]); hSimLocalDXPos_ME0Plus->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
            if ((*simTrack).charge() <0) {hSimDPhiNeg_ME0Plus->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]); hSimLocalDXNeg_ME0Plus->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
            
            if (((*simTrack).charge() >0) && ((*simTrack).momentum().pt()) >20 ) {hSimDPhiPos_HighPt_ME0Plus->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]); hSimLocalDXPos_HighPt_ME0Plus->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
            if (((*simTrack).charge() >0) && ((*simTrack).momentum().pt()) <10 ) {hSimDPhiPos_LowPt_ME0Plus->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);  hSimLocalDXPos_LowPt_ME0Plus->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
            
            if (((*simTrack).charge() <0) && ((*simTrack).momentum().pt()) >20 ) {hSimDPhiNeg_HighPt_ME0Plus->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);  hSimLocalDXNeg_HighPt_ME0Plus->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
            if (((*simTrack).charge() <0) && ((*simTrack).momentum().pt()) <10 ) {hSimDPhiNeg_LowPt_ME0Plus->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);   hSimLocalDXNeg_LowPt_ME0Plus->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
            
        }
        
        if((*simTrack).momentum().eta() < 0 ){
            
            if ((*simTrack).charge() >0) {hSimDPhiPos_ME0Minus->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]); hSimLocalDXPos_ME0Minus->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
            if ((*simTrack).charge() <0) {hSimDPhiNeg_ME0Minus->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]); hSimLocalDXNeg_ME0Minus->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
            
            if (((*simTrack).charge() >0) && ((*simTrack).momentum().pt()) >20 ) {hSimDPhiPos_HighPt_ME0Minus->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]); hSimLocalDXPos_HighPt_ME0Minus->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
            if (((*simTrack).charge() >0) && ((*simTrack).momentum().pt()) <10 ) {hSimDPhiPos_LowPt_ME0Minus->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);  hSimLocalDXPos_LowPt_ME0Minus->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
            
            if (((*simTrack).charge() <0) && ((*simTrack).momentum().pt()) >20 ) {hSimDPhiNeg_HighPt_ME0Minus->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);  hSimLocalDXNeg_HighPt_ME0Minus->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
            if (((*simTrack).charge() <0) && ((*simTrack).momentum().pt()) <10 ) {hSimDPhiNeg_LowPt_ME0Minus->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);   hSimLocalDXNeg_LowPt_ME0Minus->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
        }
        
        
    }
    std::cout<<" Num simTrack in ME0  "<<numberOfSimTracks<<" ME0Track.size()="<<ME0Tracks.size()<<" Num me0 segm="<<me0Segments->size()<<" tot num of simTrack in me0="<<AllME0Tracks.size()<<std::endl;
    //if (ME0Tracks.size()>0)
    hSelectedSimTrack->Fill(ME0Tracks.size());
   // if (me0Segments->size()>0)
    hNME0Segment->Fill(me0Segments->size());
   // if(ME0Tracks.size()>0)
    hNEvME0->Fill(1);
    
    ///////////////////////RECO VERTICES 4D Plots (by HGCal just for comparison)///////////////////////
 /*
    double reco4DVtx_x=0;double reco4DVtx_y=0;  double reco4DVtx_z=0;
    //double reco4DVtx_t=0;
    hVertexMult->Fill(primaryVertices->size());
    for (reco::VertexCollection::const_iterator it = primaryVertices->begin();  it != primaryVertices->end(); ++it) {
        hverteX->Fill( it->x());
        hverteY->Fill( it->y());
        hverteXY->Fill( it->x(), it->y());
        hverteZ->Fill( it->z());
       // hverteTime->Fill( it->t());
        
       // hverteZT->Fill( it->z(), it->t());
        //match the sim vertex with the reco vertex
        if( (fabs(SimVtxZ - it->z())< 0.1) && (fabs(SimVtxX - it->x())<0.1) && (fabs(SimVtxY-it->y()<0.1))){
            reco4DVtx_x=it->x();
            reco4DVtx_y=it->y();
            reco4DVtx_z=it->z();
            //reco4DVtx_t=it->t();
            std::cout<<"reco vtx z="<<reco4DVtx_z<<" x="<<reco4DVtx_x<<" y="<<reco4DVtx_y<<" sumpT="<<it->p4().pt()<<std::endl;
        }
    }
    */
    ///////////////////////////////////////////////////////Create ME0 simHit colletion with only muon hits////////////////////////////////////
    edm::PSimHitContainer muonSimME0Hits;
    for(uint r=0; r<ME0Tracks.size(); r++){
        for (edm::PSimHitContainer::const_iterator itHit =  ME0HitsCollection->begin(); itHit != ME0HitsCollection->end(); ++itHit){
            DetId id = DetId(itHit->detUnitId());
            if( (itHit->trackId() == ME0Tracks.at(r).trackId()) && (id.subdetId() == MuonSubdetId::ME0))	{
                muonSimME0Hits.push_back(*itHit);
            }}
    }
    
    for(ME0DigiPreRecoCollection::DigiRangeIterator cItr = me0_digis->begin(); cItr != me0_digis->end(); ++cItr){
        ME0DetId id = (*cItr).first;
        Short_t  me0_digiLayer =  id.layer();
        Short_t  me0_digiChamber =  id.chamber();
        Short_t  me0_digiRegion =  id.region();
        const GeomDet* gdet = me0Geom->idToDet(id);
        const BoundPlane & surface = gdet->surface();
        ME0DigiPreRecoCollection::const_iterator digiItr;
        for (digiItr = (*cItr ).second.first; digiItr != (*cItr ).second.second; ++digiItr){
            hDigiHitPdgID->Fill(digiItr->pdgid());
            hDigiHitToF_prompt->Fill(digiItr->tof());
            LocalPoint lp(digiItr->x(), digiItr->y(), 0);
            GlobalPoint gp = surface.toGlobal(lp);
            Float_t me0_digiX = gp.x();
            Float_t me0_digiY = gp.y();
            Float_t me0_digiEta = gp.eta();
            Float_t me0_digiPhi = gp.phi();
            Float_t me0_digiLocalX = lp.x();
            Float_t me0_digiLocalY = lp.y();
            //std::cout<<" DIGI global phi="<<me0_digiPhi<<" global eta="<<me0_digiEta<<" localX="<<me0_digiLocalX<< " localY="<<me0_digiLocalY<<" layer="<<me0_digiLayer<<" pdg="<<digiItr->pdgid()<<"  tof="<<digiItr->tof()<<" isPrompt?="<<digiItr->prompt()<<std::endl;
            //////////////////////////////Match digi in ME0Segment with simHit of muon simTtracks///////
            for (edm::PSimHitContainer::const_iterator itHit =muonSimME0Hits.begin(); itHit != muonSimME0Hits.end(); ++itHit){
                ME0DetId idme0 = ME0DetId(itHit->detUnitId());
                int layer_sim = idme0.layer();
                int chamber_sim = idme0.chamber();
                int region_sim = idme0.region();
                LocalPoint lp = itHit->entryPoint();
                GlobalPoint hitGP_sim( me0geom->idToDet( itHit->detUnitId())->surface().toGlobal(lp));
                if((itHit->particleType()==digiItr->pdgid()) && (layer_sim==me0_digiLayer) && (me0_digiChamber ==chamber_sim) && (region_sim==me0_digiRegion) ){
                    //std::cout<<"****signalMuon hits eta="<<gp.eta()<<"  phi="<<gp.phi()<<" localX="<<endl;
                    //  std::cout<<" SIMHIT global phi "<< hitGP_sim.phi()<<" global eta "<<hitGP_sim.eta()<<"  localX "<< lp.x() << " localY="<< lp.y() <<" layer "<<layer_sim<<std::endl;
                    //  std::cout<<" DIGI global phi "<<me0_digiPhi<<" global eta "<<me0_digiEta<<"  localX "<<me0_digiLocalX  << " localY="<<  me0_digiLocalY  <<" layer "<<me0_digiLayer<<" pdg="<<digiItr->pdgid()<<"  tof="<<digiItr->tof()<<std::endl;
                    hDeltaPhiSimReco->Fill(hitGP_sim.phi()-me0_digiPhi);
                    hDeltaEtaSimReco->Fill(hitGP_sim.eta()-me0_digiEta);
                    hDeltaXSimReco->Fill(hitGP_sim.x()-me0_digiX);
                    hDeltaYSimReco->Fill(hitGP_sim.y()-me0_digiY);
                    hDeltaXSimRecoLocal->Fill(lp.x()-me0_digiLocalX);
                    hDeltaYSimRecoLocal->Fill(lp.y()-me0_digiLocalY);
                    hDeltaRoll->Fill(idme0.roll()-id.roll());
                    hDeltaRollvsDeltaEta->Fill((idme0.roll()-id.roll()), (hitGP_sim.eta()-me0_digiEta));
                }
            } }
    }//loop over digi

    ////////////////////////////////////////////////////////Loop over ME0Segments////////////
    std::vector<int>   nMuonHitME0Segm;
    std::vector<int>   nEleHitME0Segm;
    std::vector<int>   nPionHitME0Segm;
    std::vector<int>   nProtonHitME0Segm;
    std::vector<int>   nBoHitME0Segm;
    std::vector<int>   nMatchedHit;
    std::vector<uint>   nTrackID;
    std::vector<double> SignalTime, BkgTime, ME0SegTime, recHitTOF   ;
    std::vector<int>  nSegmDPhiPt5GeV, nSegmDPhiPt5GeV_AllCut;
    std::vector<int>  nSegmDPhiPt30GeV, nSegmDPhiPt30GeV_AllCut;
    std::vector<int>   nAllHit, nAllHitSig, nAllHitBkg ;
    std::vector<int>  nMuonHitME0Segm_Sig, nEleHitME0Segm_Sig, nHadHitME0Segm_Sig;
    std::vector<int>  nMuonHitME0Segm_Bkg, nEleHitME0Segm_Bkg, nHadHitME0Segm_Bkg;
    int segIt=0;
    std::pair<int, int> SegIDMatchedHitPair;
    
    struct matchedSegm_t {
        int nHitMatched;
        int SegID;
        uint matchedTrackID;
        float SegEta;
        float SegPhi;
        int   SegNhits;
        float SegChi2;
        float SegTime;
    };
    matchedSegm_t matchedSegm;
    std::vector<matchedSegm_t>    matchedSegmDB;
    std::vector<std::pair<int, int>> SegIDMatchedHitPairVec; 
    if(selME0SimHits.size()>0){
        
        ME0SegmentCollection::const_iterator me0segIt;
        std::vector<int> nME0SegmentSignal, nME0SegmentBkg;
        std::vector<int> nME0SegmentSignal_timeCut, nME0SegmentBkg_timeCut;
        std::vector<int> nHitME0Segment;

  
        
        for(me0segIt=me0Segments->begin(); me0segIt<me0Segments->end(); me0segIt++ ){
            nHitME0Segment.clear(); nMatchedHit.clear(); nTrackID.clear();
            nMuonHitME0Segm.clear(); nEleHitME0Segm.clear(); nPionHitME0Segm.clear(); nProtonHitME0Segm.clear(); nBoHitME0Segm.clear(); nAllHit.clear();
            
            //cout<<"ME0Segment nHits="<<me0segIt->nRecHits()<<" tof="<<me0segIt->time()<<" chi2="<<me0segIt->chi2()<<endl;
            std::vector<double>   me0RHPhiME0Seg_noMatch;
            std::vector<double>   me0RHPhiME0Seg_noMatchTimeWindow;
            auto me0rhs =  me0segIt->specificRecHits();
            nHitME0Segment.push_back(me0rhs.size());
            hME0Segm_nHit->Fill(me0rhs.size());
            hME0Segm_chi2->Fill(me0segIt->chi2());
            //hME0Segm_time->Fill(me0segIt->time());
            hME0SegmentCollectionTime->Fill(me0segIt->time());
	   
            LocalPoint ME0SegmentLocP = me0segIt->localPosition();
            ME0DetId segME0idAll = me0segIt->me0DetId();
            //GlobalPoint ME0SegmentGPAll = (me0geom->chamber(segME0idAll))->toGlobal(ME0SegmentLocP);
            GlobalPoint ME0SegmentGPAll = (me0geom->chamber(segME0idAll))->toGlobal(ME0SegmentLocP);
            hME0Segm_eta->Fill(ME0SegmentGPAll.eta());
           
            for (auto rh = me0rhs.begin(); rh!= me0rhs.end(); rh++){
                auto rhLP = rh->localPosition();
                auto me0id = rh->me0Id();
                auto rhr = me0geom->chamber(me0id);
                auto rhGP = rhr->toGlobal(rhLP);
                double globalPhi = rhGP.phi();  double globalX = rhGP.x();			double globalY = rhGP.y();
                double localX = rhLP.x();			double localY = rhLP.y();
                const ME0EtaPartition* roll = me0Geom->etaPartition(me0id);
               // cout<<" RECHit global phi="<<globalPhi << "global eta="<<rhGP.eta()<<"  localX="<< rhLP.x() << " localY="<<  rhLP.y()  <<" layer="<<me0id.layer()<<" roll="<<me0id.roll()<<" time="<<rh->tof()<<endl;
                
                
                me0RHPhiME0Seg_noMatch.push_back(globalPhi);
                if( (me0segIt->time()<30.5) && (me0segIt->time()>5.5))	{
                    //	if((rh->tof()<30.5) && (rh->tof()>5.5)){
                    me0RHPhiME0Seg_noMatchTimeWindow.push_back(globalPhi);
                    
                }
                
                //for(ME0DigiPreRecoCollection::DigiRangeIterator cItr = me0_digis->begin(); cItr != me0_digis->end(); ++cItr){
                for(MuonDigiCollection<ME0DetId,ME0DigiPreReco>::DigiRangeIterator cItr = me0_digis->begin(); cItr != me0_digis->end(); ++cItr){
                 
                    ME0DetId id = (*cItr).first;
                    Short_t  me0_digiLayer =  id.layer();
                    Short_t  me0_digiChamber =  id.chamber();
                    const GeomDet* gdet = me0Geom->idToDet(id);
                    const BoundPlane & surface = gdet->surface();
                    //ME0DigiPreRecoCollection::const_iterator digiItr;
                    MuonDigiCollection<ME0DetId,ME0DigiPreReco>::const_iterator digiItr;
                    
                    for (digiItr = (*cItr ).second.first; digiItr != (*cItr ).second.second; ++digiItr){
                        
                        LocalPoint lp(digiItr->x(), digiItr->y(), 0);
                        GlobalPoint gp = surface.toGlobal(lp);
                        Float_t me0_digiX = gp.x();
                        Float_t me0_digiY = gp.y();
                        Float_t me0_digiEta = gp.eta();
                        Float_t me0_digiPhi = gp.phi();
                        Float_t me0_digiLocalX = lp.x();
                        Float_t me0_digiLocalY = lp.y();
                        //std::cout<<" DIGI global phi "<<me0_digiPhi<<" global eta "<<me0_digiEta<<"  localX "<< me0_digiLocalX << " localY="<<  me0_digiLocalY  <<" layer "<<me0_digiLayer<<" pdg="<<digiItr->pdgid()<<"  tof="<<digiItr->tof()<<std::endl;
                        if(localX == me0_digiLocalX && localY == me0_digiLocalY && me0_digiLayer==me0id.layer() && (me0_digiChamber ==me0id.chamber()) && (id.roll()==me0id.roll()) /*&& digiItr->prompt()*/  ){
                      //  if( me0_digiLayer==me0id.layer() && (me0_digiChamber ==me0id.chamber()) &&(id.roll()==me0id.roll())){
                        //  std::cout<<" DIGI global phi "<<me0_digiPhi<<" global eta "<<me0_digiEta<<"  localX "<< me0_digiLocalX << " localY="<<  me0_digiLocalY  <<" layer "<<me0_digiLayer<<" pdg="<<digiItr->pdgid()<<"  tof="<<digiItr->tof()<<std::endl;
                            //std::cout<<" digiHit: pdgID="<<digiItr->pdgid()<<" is prompt? "<<digiItr->prompt()<<"  tof="<<digiItr->tof()<<" global eta="<<gp.eta()<<" global phi="<<gp.phi()<<std::endl;
                            
                                recHitTOF.push_back(digiItr->tof());
                                hME0SegmDigiHitTOF->Fill(digiItr->tof());
                                hME0SegmDigiHitPdgID->Fill(digiItr->pdgid());
                                hME0SegmDigiHitEta->Fill(gp.eta());
                            
                            
                                nAllHit.push_back(1);
                                if(fabs(digiItr->pdgid())==13) {hMuonInME0SegTOF->Fill(digiItr->tof()) ;   nMuonHitME0Segm.push_back(1);}
                                if(fabs(digiItr->pdgid())==11) {hEleInME0SegTOF->Fill(digiItr->tof());    nEleHitME0Segm.push_back(1);}
                                if(fabs(digiItr->pdgid())==211) {hPionInME0SegTOF->Fill(digiItr->tof());  nPionHitME0Segm.push_back(1);}
                                if(fabs(digiItr->pdgid())==321) {hBoInME0SegTOF->Fill(digiItr->tof());    nBoHitME0Segm.push_back(1);}
                                if(fabs(digiItr->pdgid())==2212) {hProtonsInME0SegTOF->Fill(digiItr->tof()); nProtonHitME0Segm.push_back(1);}
                                
                                if(me0_digiLayer ==1){
                                    if(fabs(digiItr->pdgid())==13) hMuonInME0SegTOF_L1->Fill(digiItr->tof());
                                    if(fabs(digiItr->pdgid())==11) hEleInME0SegTOF_L1->Fill(digiItr->tof());
                                    if(fabs(digiItr->pdgid())==211) hPionInME0SegTOF_L1->Fill(digiItr->tof());
                                    if(fabs(digiItr->pdgid())==321) hBoInME0SegTOF_L1->Fill(digiItr->tof());
                                    if(fabs(digiItr->pdgid())==2212) hProtonsInME0SegTOF_L1->Fill(digiItr->tof());
                                }
                                if( (me0segIt->time()<30.5) && (me0segIt->time()>5.5))	{
                                    if(fabs(digiItr->pdgid())==13) hMuonInME0SegTOF_inTime->Fill(digiItr->tof());
                                    if(fabs(digiItr->pdgid())==11) hEleInME0SegTOF_inTime->Fill(digiItr->tof());
                                    if(fabs(digiItr->pdgid())==211) hPionInME0SegTOF_inTime->Fill(digiItr->tof());
                                    if(fabs(digiItr->pdgid())==321) hBoInME0SegTOF_inTime->Fill(digiItr->tof());
                                    if(fabs(digiItr->pdgid())==2212) hProtonsInME0SegTOF_inTime->Fill(digiItr->tof());
                                }
                                
                           

                                hME0SegmDigiHitTOFME0Segm_noprompt->Fill(digiItr->tof());
                                hME0SegmDigiHitPdgIDME0Segm_noprompt->Fill(digiItr->pdgid());
                            
                            
                            
				//////////////////////////////Match digi in ME0Segment with simHit of muon simTtracks///////
                            
                            for (edm::PSimHitContainer::const_iterator itHit =muonSimME0Hits.begin(); itHit != muonSimME0Hits.end(); ++itHit){
                                ME0DetId idme0 = ME0DetId(itHit->detUnitId());
                                int layer_sim = idme0.layer();
                                int chamber_sim = idme0.chamber();
                                LocalPoint lp = itHit->entryPoint();
                                GlobalPoint hitGP_sim( me0geom->idToDet( itHit->detUnitId())->surface().toGlobal(lp));
                                //cout<<"simhit pdgID="<<itHit->particleType()<<" eta="<<hitGP_sim.eta()<<" phi="<<hitGP_sim.phi()<<endl;
                                if( (itHit->particleType()==digiItr->pdgid()) && (layer_sim==me0_digiLayer) && /* (fabs(lp.x()-me0_digiLocalX)<=3*0.03) && (fabs(lp.y()-me0_digiLocalY)<=3*0.01) &&*/ (me0_digiChamber ==chamber_sim) &&  (id.roll()==idme0.roll())  ){
                                //if( (itHit->particleType()==digiItr->pdgid()) && (layer_sim==me0_digiLayer) && (fabs(lp.x()-me0_digiLocalX)<=0) && (fabs(lp.y()-me0_digiLocalY)<=0) && (me0_digiChamber ==chamber_sim)   ){
                                   // std::cout<<" SIMHIT global phi "<< hitGP_sim.phi()<<" global eta="<<hitGP_sim.eta()<<"  localX="<<lp.x() << " localY="<< lp.y() <<" layer="<<layer_sim<<" pdg="<<itHit->particleType()<<std::endl;
                                   // std::cout<<" DIGI global phi "<<me0_digiPhi<<" global eta="<<me0_digiEta<<"  localX="<<me0_digiLocalX  << " localY="<<  me0_digiLocalY  <<" layer="<<me0_digiLayer<<" pdg="<<digiItr->pdgid()<<" roll="<<id.roll()<<std::endl;
                                    nMatchedHit.push_back(1);
                                    nTrackID.push_back(itHit->trackId());
                                    /*hDeltaPhiSimReco->Fill(hitGP_sim.phi()-me0_digiPhi);
                                    hDeltaEtaSimReco->Fill(hitGP_sim.eta()-me0_digiEta);
                                    hDeltaXSimReco->Fill(hitGP_sim.x()-me0_digiX);
                                    hDeltaYSimReco->Fill(hitGP_sim.y()-me0_digiY);
                                    hDeltaXSimRecoLocal->Fill(lp.x()-me0_digiLocalX);
                                    hDeltaYSimRecoLocal->Fill(lp.y()-me0_digiLocalY);*/

                                }
                                //else{//std::cout<<"****BkgPart hits eta="<<gp.eta()<<"  phi="<<gp.phi()<<endl;}
                            }
                            
                        }//match digi-rechit
                        
                    }	}//loop over digi
                
            }//loop over rechit in me0segm

	    //	    SegIDMatchedHitPair = std::make_pair (segIt,nMatchedHit);
            double sumRHtof=0.;
            for(uint t=0; t<recHitTOF.size(); t++){
                sumRHtof += recHitTOF.at(t);}
                double ME0SegmMyTOF = sumRHtof/recHitTOF.size();
                hDPhiNoMatchME0Seg->Fill(me0RHPhiME0Seg_noMatch[0]- me0RHPhiME0Seg_noMatch[me0RHPhiME0Seg_noMatch.size()-1]);
                if(me0RHPhiME0Seg_noMatchTimeWindow.size())  hDPhiNoMatchME0Seg_TimeWindow->Fill(me0RHPhiME0Seg_noMatchTimeWindow[0] - me0RHPhiME0Seg_noMatchTimeWindow[me0RHPhiME0Seg_noMatchTimeWindow.size()-1]);
	    //std::cout<<" matched hit= "<<	nMatchedHit.size()<<" tof="<<ME0SegmMyTOF<<endl;
            hME0SegmMyTOF->Fill(ME0SegmMyTOF);
            hNumberOfMatchedHit->Fill(nMatchedHit.size());
            std::vector<double> me0RHPhiME0Seg_matchByHitsSignal,me0RHPhiME0Seg_matchByHitsSignal_timeCut ;
            std::vector<double> me0RHPhiME0Seg_matchByHitsBkg, me0RHPhiME0Seg_matchByHitsBkg_timeCut;
            ///////////////////////////////////SIGNAL////////////////////////////////
            //if(nMatchedHit.size()>=3 && nMatchedHit.size()<=6){
            if(nMatchedHit.size()>=3 && (nMatchedHit.size()/nHitME0Segment.size())>0.5){
                segIt++;
                SegIDMatchedHitPair = std::make_pair (segIt,nMatchedHit.size());
                SegIDMatchedHitPairVec.push_back(SegIDMatchedHitPair);
                matchedSegm.nHitMatched = nMatchedHit.size();
                matchedSegm.SegID = segIt;
                matchedSegm.matchedTrackID=nTrackID[0];
                matchedSegm.SegTime =me0segIt->time() ;
                matchedSegm.SegChi2 =me0segIt->chi2() ;
                matchedSegm.SegNhits =me0rhs.size() ;
                std::cout<<"Signal matched hit= "<<nMatchedHit.size()<<endl;
                // nMuonHitME0Segm_Sig.clear(); nEleHitME0Segm_Sig.clear(); nHadHitME0Segm_Sig.clear();
                SignalTime.push_back(me0segIt->time());
                auto me0rhs =  me0segIt->specificRecHits();
                hME0Segm_nHit_sig->Fill(me0rhs.size());
                hME0Segm_chi2_sig->Fill(me0segIt->chi2());
                hME0Segm_time_sig->Fill(me0segIt->time());
                ME0SegSig.push_back(*me0segIt);
               // LocalVector ME0SegmentLP_sig = me0segIt->localDirection();
               // ME0DetId segME0id_sig = me0segIt->me0DetId();
               // GlobalVector ME0SegmentGP_sig = (me0geom->chamber(segME0id_sig))->toGlobal(ME0SegmentLP_sig);
                
                LocalPoint ME0SegmentLP_sig = me0segIt->localPosition();
                ME0DetId segME0id_sig = me0segIt->me0DetId();
                GlobalPoint ME0SegmentGP_sig = (me0geom->chamber(segME0id_sig))->toGlobal(ME0SegmentLP_sig);
                matchedSegm.SegEta =ME0SegmentGP_sig.eta();
                matchedSegm.SegPhi =ME0SegmentGP_sig.phi();
                matchedSegmDB.push_back(matchedSegm);
                
                hME0Segm_eta_sig->Fill(ME0SegmentGP_sig.eta());
                hME0Segm_phi_sig->Fill(ME0SegmentGP_sig.phi());
                hME0Segm_X_sig->Fill(ME0SegmentGP_sig.x());
                hME0Segm_Y_sig->Fill(ME0SegmentGP_sig.y());
                if(me0segIt->chi2() !=0) hME0Segm_eta_sigCheck->Fill(ME0SegmentGP_sig.eta());
                if(me0segIt->chi2() ==0) hME0Segm_eta_sigCheckFail->Fill(ME0SegmentGP_sig.eta());
                for (auto rh = me0rhs.begin(); rh!= me0rhs.end(); rh++){
                    auto rhLP = rh->localPosition();
                    auto me0id = rh->me0Id();
                    auto rhr = me0geom->chamber(me0id);
                    auto rhGP = rhr->toGlobal(rhLP);
                    double globalPhi = rhGP.phi();   double globalX = rhGP.x() ;  double globalY = rhGP.y();
                    me0RHPhiME0Seg_matchByHitsSignal.push_back(globalPhi);
                    for(ME0DigiPreRecoCollection::DigiRangeIterator cItr = me0_digis->begin(); cItr != me0_digis->end(); ++cItr){
                        ME0DetId id = (*cItr).first;
                        Short_t  me0_digiLayer =  id.layer();
                        const GeomDet* gdet = me0Geom->idToDet(id);
                        const BoundPlane & surface = gdet->surface();
                        ME0DigiPreRecoCollection::const_iterator digiItr;
                        for (digiItr = (*cItr ).second.first; digiItr != (*cItr ).second.second; ++digiItr){
                            LocalPoint lp(digiItr->x(), digiItr->y(), 0);
                            GlobalPoint gp = surface.toGlobal(lp);
                            Float_t me0_digiX = gp.x();
                            Float_t me0_digiY = gp.y();
                            if(globalX == me0_digiX && globalY == me0_digiY && me0_digiLayer==me0id.layer()){
			      
			      hME0SegmDigiHitTOF_signal->Fill(digiItr->tof());
			      hME0SegmDigiHitPdgID_signal->Fill(digiItr->pdgid());
			      hME0SegmDigiHitEta_signal->Fill(gp.eta());
			      nAllHitSig.push_back(1);
			      if(fabs(digiItr->pdgid())==13) {nMuonHitME0Segm_Sig.push_back(1);}
			      if(fabs(digiItr->pdgid())==11) {nEleHitME0Segm_Sig.push_back(1);}
			      if(fabs(digiItr->pdgid())>200) {nHadHitME0Segm_Sig.push_back(1);}
			    }}}
                    
                }
                hDPhiMatchByHitsME0Seg_Signal->Fill(me0RHPhiME0Seg_matchByHitsSignal[0]- me0RHPhiME0Seg_matchByHitsSignal[me0RHPhiME0Seg_matchByHitsSignal.size()-1]);
                hDPhivsTOFMatchByHitsME0Seg_Signal->Fill(me0RHPhiME0Seg_matchByHitsSignal[0]- me0RHPhiME0Seg_matchByHitsSignal[me0RHPhiME0Seg_matchByHitsSignal.size()-1],me0segIt->time());
                nME0SegmentSignal.push_back(1);
                
                if( (me0segIt->time()>timeMin_)  && (me0segIt->time()<timeMax_) ){
                    hDPhiMatchByHitsME0Seg_Signal_timeCut->Fill(me0RHPhiME0Seg_matchByHitsSignal[0]- me0RHPhiME0Seg_matchByHitsSignal[me0RHPhiME0Seg_matchByHitsSignal.size()-1]);
                    hDPhivsTOFMatchByHitsME0Seg_Signal_timeCut->Fill(me0RHPhiME0Seg_matchByHitsSignal[0]- me0RHPhiME0Seg_matchByHitsSignal[me0RHPhiME0Seg_matchByHitsSignal.size()-1],me0segIt->time());
                    nME0SegmentSignal_timeCut.push_back(1);
                }
                double DeltaPhiSeg = me0RHPhiME0Seg_matchByHitsSignal[0]- me0RHPhiME0Seg_matchByHitsSignal[me0RHPhiME0Seg_matchByHitsSignal.size()-1];

                if( TMath::Abs(DeltaPhiSeg)<0.002  ){
                    for(uint r=0; r<ME0Tracks.size(); r++){
                        if(ME0Tracks.at(r).trackId() == nTrackID.at(0)) {
                            hSimMuonPt_DPhiCut30GeV->Fill(ME0Tracks.at(r).momentum().pt());
                            if(me0rhs.size()>3 && me0segIt->time()>17.5 && me0segIt->time()<19.1) hSimMuonPt_DPhiCut30GeVAllCut->Fill(ME0Tracks.at(r).momentum().pt());
                        }
                    }}// if 30 GeV
                
                if( TMath::Abs(DeltaPhiSeg)<0.01  ){
		  for(uint r=0; r<ME0Tracks.size(); r++){
		    if(ME0Tracks.at(r).trackId() == nTrackID.at(0)) {
		      hSimMuonPt_DPhiCut5GeV->Fill(ME0Tracks.at(r).momentum().pt());
		      if(me0rhs.size()>3 && me0segIt->time()>17.5 && me0segIt->time()<19.1) hSimMuonPt_DPhiCut5GeVAllCut->Fill(ME0Tracks.at(r).momentum().pt());
		    }
		  }}// if 5 GeV
		
		//////////////Denominator for a LocalReco Efficiency//////////////               
        for(uint r=0; r<ME0Tracks.size(); r++){
		  // cout<<"Muon ME0SimTrack Eta: "<<ME0Tracks.at(r).momentum().eta()<<" type="<<ME0Tracks.at(r).type()<<" ID="<<ME0Tracks.at(r).trackId()<<endl;
		  if(ME0Tracks.at(r).trackId() == nTrackID.at(0)) {
		    //cout<<"Matched ME0SimTrack pt="<<ME0Tracks.at(r).momentum().pt()<<" id="<<ME0Tracks.at(r).trackId()<<" deltaPhi="<<DeltaPhiSeg<<endl;
		
		    hME0Segm_DeltaPhiVsSimpT_Signal->Fill(DeltaPhiSeg, ME0Tracks.at(r).momentum().pt());
		    //hME0Segm_DeltaPhiVsSimpT_Signal_rebin->Fill(DeltaPhiSeg, ME0Tracks.at(r).momentum().pt());
		    if(ME0Tracks.at(r).momentum().pt()<=0.5) hME0Segm_DeltaPhi_pT0p5->Fill(DeltaPhiSeg);
			if(ME0Tracks.at(r).momentum().pt()>0.5 && ME0Tracks.at(r).momentum().pt()<=1)hME0Segm_DeltaPhi_pT1->Fill(DeltaPhiSeg);
			if(ME0Tracks.at(r).momentum().pt()>1 && ME0Tracks.at(r).momentum().pt()<=3)hME0Segm_DeltaPhi_pT3->Fill(DeltaPhiSeg);
			if(ME0Tracks.at(r).momentum().pt()>3 && ME0Tracks.at(r).momentum().pt()<=5)hME0Segm_DeltaPhi_pT5->Fill(DeltaPhiSeg);
			if(ME0Tracks.at(r).momentum().pt()>5 && ME0Tracks.at(r).momentum().pt()<=10)hME0Segm_DeltaPhi_pT10->Fill(DeltaPhiSeg);
			if(ME0Tracks.at(r).momentum().pt()>10 && ME0Tracks.at(r).momentum().pt()<=20)hME0Segm_DeltaPhi_pT20->Fill(DeltaPhiSeg);
			if(ME0Tracks.at(r).momentum().pt()>20 && ME0Tracks.at(r).momentum().pt()<=50)hME0Segm_DeltaPhi_pT50->Fill(DeltaPhiSeg);
			  hME0Segm_SimpT_Signal->Fill(ME0Tracks.at(r).momentum().pt());

			if(me0segIt->chi2() !=0){
			  hME0SegmPullPhi->Fill(ME0Tracks.at(r).momentum().phi()-ME0SegmentGP_sig.phi());
			  hME0SegmPullEta->Fill(ME0Tracks.at(r).momentum().eta()-ME0SegmentGP_sig.eta());
			}
			if(me0segIt->chi2()==0){
			  hME0SegmPullPhi_fitFail->Fill(ME0Tracks.at(r).momentum().phi()-ME0SegmentGP_sig.phi());
              hME0SegmPullEta_fitFail->Fill(ME0Tracks.at(r).momentum().eta()-ME0SegmentGP_sig.eta());
			}
                    }
                }
            
            
                
                ////////////////////////////////////VertexEstrapolation Attempt////////////////////////////////////////////////////////////////////////
                
                double Lenght_vtxToME0=c*me0segIt->time();
                double timeFlight= me0segIt->time();
                double ctau= c*timeFlight;
                LocalPoint ME0SegmentLP = me0segIt->localPosition();
                ME0DetId segME0id = me0segIt->me0DetId();
                GlobalPoint ME0SegmentGP = (me0geom->chamber(segME0id))->toGlobal(ME0SegmentLP);
                double myZrec = TMath::Sqrt(Lenght_vtxToME0*Lenght_vtxToME0 - (ME0SegmentGP.y()*ME0SegmentGP.y()) - (ME0SegmentGP.x()*ME0SegmentGP.x()) )  ;
                double myZrec2 = TMath::Sqrt(Lenght_vtxToME0*Lenght_vtxToME0 - (ME0SegmentGP.y() - SimVtxY)*(ME0SegmentGP.y() - SimVtxY) - (ME0SegmentGP.x()-SimVtxX)*(ME0SegmentGP.x()-SimVtxX))  ;
               // double myZrec3= TMath::Sqrt(ctau*ctau - (ME0SegmentGP.y()-reco4DVtx_y)*(ME0SegmentGP.y()-reco4DVtx_y) - (ME0SegmentGP.x()-reco4DVtx_x)*(ME0SegmentGP.x()-reco4DVtx_x) )  ;
                // double myZrec3  = ctau;
                
                //double zV=0.; double zV2=0.;  double zV3=0.;
                //if(ME0SegmentGP.z()>0) {zV=ME0SegmentGP.z()-myZrec; zV2=ME0SegmentGP.z()-myZrec2; ; zV3=ME0SegmentGP.z()/*-myZrec3*/;}
                //if(ME0SegmentGP.z()<0) {zV=ME0SegmentGP.z()+ myZrec;zV2=ME0SegmentGP.z()+ myZrec2 ; zV3=ME0SegmentGP.z()/*+ myZrec3*/ ;}
                
                //cout<<"  ME0Seg time="<<me0segIt->time()<<" Lenght_vtxToME0="<<Lenght_vtxToME0<<" ME0Seg localPosY="<<me0segIt->localPosition().y() <<" GlobalY="<<ME0SegmentGP.y()<<" GlobalX="<<ME0SegmentGP.x()<<" GlobalZ="<<ME0SegmentGP.z()<<endl;
		//                cout<<" myZrec="<<myZrec<<" Zv="<<zV<<" Zv2="<<zV2<<" ZV3="<<zV3<<" "<<std::endl;
                
                ME0SegTime.push_back(me0segIt->time());
                ////////////////////////////////////VertexEstrapolation////////////////////////////////////////////////////////////////////////
                
            }//if Signal
            //////////////////////////////////////
            
            ///compute zVertex and tVertex only for signal ME0Segments by crossing the 2 track info
            ///and solving the equation
            if(ME0SegTime.size()==2){
                double t1=ME0SegTime.at(0);
                double t2=ME0SegTime.at(1);
                double zME0=527;
                double tvertex=(t2+t1)/2 - (zME0/c);
                double zvertex=c*tvertex + zME0  -(c*t1);
                hPull_MyRecoSim_zVtx->Fill(zvertex-SimVtxZ);
               // hPull_4DRecoSim_zVtx->Fill(reco4DVtx_z-SimVtxZ);
            }
            
            ///////////////////////////////////BKG////////////////////////////////
    //        if(nMatchedHit.size()<=2){
            if(nMatchedHit.size()<=2 && (nMatchedHit.size()/nHitME0Segment.size())<0.5){
                //nMuonHitME0Segm_Bkg.clear(); nEleHitME0Segm_Bkg.clear(); nHadHitME0Segm_Bkg.clear();
                hTimeMatchByHitsME0Seg_Bkg->Fill( me0segIt->time());
                BkgTime.push_back(me0segIt->time());
                auto me0rhs =  me0segIt->specificRecHits();
                hME0Segm_nHit_bkg->Fill(me0rhs.size());
                hME0Segm_chi2_bkg->Fill(me0segIt->chi2());
                hME0Segm_time_bkg->Fill(me0segIt->time());
	            ME0SegBkg.push_back(*me0segIt);

                LocalPoint ME0SegmentLP_bkg = me0segIt->localPosition();
                ME0DetId segME0id_bkg = me0segIt->me0DetId();
                GlobalPoint ME0SegmentGP_bkg = (me0geom->chamber(segME0id_bkg))->toGlobal(ME0SegmentLP_bkg);
                hME0Segm_eta_bkg->Fill(ME0SegmentGP_bkg.eta());
                if(me0segIt->chi2() !=0) hME0Segm_eta_bkgCheck->Fill(ME0SegmentGP_bkg.eta());
                if(me0segIt->chi2() ==0) hME0Segm_eta_bkgCheckFail->Fill(ME0SegmentGP_bkg.eta());

                for (auto rh = me0rhs.begin(); rh!= me0rhs.end(); rh++){
                    auto rhLP = rh->localPosition();
                    auto me0id = rh->me0Id();
                    auto rhr = me0geom->chamber(me0id);
                    auto rhGP = rhr->toGlobal(rhLP);
                    double globalPhi = rhGP.phi();  double globalX = rhGP.x() ;  double globalY = rhGP.y();
                    me0RHPhiME0Seg_matchByHitsBkg.push_back(globalPhi);
                    for(ME0DigiPreRecoCollection::DigiRangeIterator cItr = me0_digis->begin(); cItr != me0_digis->end(); ++cItr){
                        ME0DetId id = (*cItr).first;
                        Short_t  me0_digiLayer =  id.layer();
                        const GeomDet* gdet = me0Geom->idToDet(id);
                        const BoundPlane & surface = gdet->surface();
                        ME0DigiPreRecoCollection::const_iterator digiItr;
                        for (digiItr = (*cItr ).second.first; digiItr != (*cItr ).second.second; ++digiItr){
                            LocalPoint lp(digiItr->x(), digiItr->y(), 0);
                            GlobalPoint gp = surface.toGlobal(lp);
                            Float_t me0_digiX = gp.x();
                            Float_t me0_digiY = gp.y();
                            if(globalX == me0_digiX && globalY == me0_digiY && me0_digiLayer==me0id.layer()){
                               
                                hME0SegmDigiHitTOF_bkg->Fill(digiItr->tof());
                                hME0SegmDigiHitPdgID_bkg->Fill(digiItr->pdgid());
                                hME0SegmDigiHitEta_bkg->Fill(gp.eta());
				if(fabs(digiItr->pdgid())==13) {nMuonHitME0Segm_Bkg.push_back(1);}
				if(fabs(digiItr->pdgid())==11) {nEleHitME0Segm_Bkg.push_back(1);}
				if(fabs(digiItr->pdgid())>200) {nHadHitME0Segm_Bkg.push_back(1);}
                                    
                                
                            }}}
                    
                }
                hDPhiMatchByHitsME0Seg_Bkg->Fill(me0RHPhiME0Seg_matchByHitsBkg[0]-me0RHPhiME0Seg_matchByHitsBkg[me0RHPhiME0Seg_matchByHitsBkg.size()-1]);
                hDPhivsTOFMatchByHitsME0Seg_Bkg->Fill(me0RHPhiME0Seg_matchByHitsBkg[0]-me0RHPhiME0Seg_matchByHitsBkg[me0RHPhiME0Seg_matchByHitsBkg.size()-1],me0segIt->time());
                nME0SegmentBkg.push_back(1);
                
                if( (me0segIt->time()>timeMin_)  && (me0segIt->time()<timeMax_) ){
                    hDPhiMatchByHitsME0Seg_Bkg_timeCut->Fill(me0RHPhiME0Seg_matchByHitsBkg[0]-me0RHPhiME0Seg_matchByHitsBkg[me0RHPhiME0Seg_matchByHitsBkg.size()-1]);
                    hDPhivsTOFMatchByHitsME0Seg_Bkg_timeCut->Fill(me0RHPhiME0Seg_matchByHitsBkg[0]-me0RHPhiME0Seg_matchByHitsBkg[me0RHPhiME0Seg_matchByHitsBkg.size()-1],me0segIt->time());
                    nME0SegmentBkg_timeCut.push_back(1);
                }
                
                double DeltaPhiBkg = me0RHPhiME0Seg_matchByHitsBkg[0]-me0RHPhiME0Seg_matchByHitsBkg[me0RHPhiME0Seg_matchByHitsBkg.size()-1];
                if ( TMath::Abs(DeltaPhiBkg)<0.002) { nSegmDPhiPt30GeV.push_back(1); }
                if ( TMath::Abs(DeltaPhiBkg)<0.01) { nSegmDPhiPt5GeV.push_back(1); }
                
                if ( TMath::Abs(DeltaPhiBkg)<0.002 && me0rhs.size()>3 && me0segIt->time()>17.5 && me0segIt->time()<19.1) nSegmDPhiPt30GeV_AllCut.push_back(1);
                if ( TMath::Abs(DeltaPhiBkg)<0.01 && me0rhs.size()>3 && me0segIt->time()>17.5 && me0segIt->time()<19.1) nSegmDPhiPt5GeV_AllCut.push_back(1);
                
        
                
            }//if background
            
	  
            ///Fill plot for all ME0Segment composition
            if(nAllHit.size()>0) hME0SegmNumber->Fill(1);
            
            int nHadAll= nPionHitME0Segm.size() + nBoHitME0Segm.size() +nProtonHitME0Segm.size();
            if((nMuonHitME0Segm.size()>0) && (nEleHitME0Segm.size()==0) && (nHadAll==0)   ) {hME0SegCompMuonOnly->Fill(1);}
            if((nMuonHitME0Segm.size()>0)  && (nEleHitME0Segm.size()>0) && (nHadAll==0)   ) {hME0SegCompMuonPlusEle->Fill(1);}
            if((nMuonHitME0Segm.size()>0)  && (nEleHitME0Segm.size()==0)&&  (nHadAll>0)  ) {hME0SegCompMuonPlusHad->Fill(1);}
            if((nMuonHitME0Segm.size()>0)  && (nEleHitME0Segm.size()>0) &&  (nHadAll>0)  ) {hME0SegCompMuonPlusElePlusHad->Fill(1);}
            if((nMuonHitME0Segm.size()==0)  && (nEleHitME0Segm.size()==0) && (nHadAll>0) ) {hME0SegCompHadOnly->Fill(1);}
            if((nMuonHitME0Segm.size()==0)  && (nEleHitME0Segm.size()>0) &&  (nHadAll==0) ) {hME0SegCompEleOnly->Fill(1);}
            
            ///Fill plot for all ME0Segment composition  (signal only)
            if(nAllHitSig.size()>0) hME0SegmNumberSig->Fill(1);
            if((nMuonHitME0Segm_Sig.size()>0) && (nEleHitME0Segm_Sig.size()==0) && (nHadHitME0Segm_Sig.size()==0)) hME0SegCompMuonOnly_Sig->Fill(1);
            if((nMuonHitME0Segm_Sig.size()>0) && (nEleHitME0Segm_Sig.size()>0) && (nHadHitME0Segm_Sig.size()==0)) hME0SegCompMuonPlusEle_Sig->Fill(1);
            if((nMuonHitME0Segm_Sig.size()>0) && (nEleHitME0Segm_Sig.size()==0) && (nHadHitME0Segm_Sig.size()>0)) hME0SegCompMuonPlusHad_Sig->Fill(1);
            if((nMuonHitME0Segm_Sig.size()>0) && (nEleHitME0Segm_Sig.size()>0) && (nHadHitME0Segm_Sig.size()>0)) hME0SegCompMuonPlusElePlusHad_Sig->Fill(1);
            if((nMuonHitME0Segm_Sig.size()==0) && (nEleHitME0Segm_Sig.size()>0) && (nHadHitME0Segm_Sig.size()==0)) hME0SegCompEleOnly_Sig->Fill(1);
            if((nMuonHitME0Segm_Sig.size()==0) && (nEleHitME0Segm_Sig.size()==0) && (nHadHitME0Segm_Sig.size()>0)) hME0SegCompHadOnly_Sig->Fill(1);
            
            ///Fill plot for all ME0Segment composition  (bkg only)
            if(nAllHitBkg.size()>0) hME0SegmNumberBkg->Fill(1);
            if((nMuonHitME0Segm_Bkg.size()>0) && (nEleHitME0Segm_Bkg.size()==0) && (nHadHitME0Segm_Bkg.size()==0)) hME0SegCompMuonOnly_Bkg->Fill(1);
            if((nMuonHitME0Segm_Bkg.size()>0) && (nEleHitME0Segm_Bkg.size()>0) && (nHadHitME0Segm_Bkg.size()==0)) hME0SegCompMuonPlusEle_Bkg->Fill(1);
            if((nMuonHitME0Segm_Bkg.size()>0) && (nEleHitME0Segm_Bkg.size()==0) && (nHadHitME0Segm_Bkg.size()>0)) hME0SegCompMuonPlusHad_Bkg->Fill(1);
            if((nMuonHitME0Segm_Bkg.size()>0) && (nEleHitME0Segm_Bkg.size()>0) && (nHadHitME0Segm_Bkg.size()>0)) hME0SegCompMuonPlusElePlusHad_Bkg->Fill(1);
            if((nMuonHitME0Segm_Bkg.size()==0) && (nEleHitME0Segm_Bkg.size()>0) && (nHadHitME0Segm_Bkg.size()==0)) hME0SegCompEleOnly_Bkg->Fill(1);
            if((nMuonHitME0Segm_Bkg.size()==0) && (nEleHitME0Segm_Bkg.size()==0) && (nHadHitME0Segm_Bkg.size()>0)) hME0SegCompHadOnly_Bkg->Fill(1);
            
            
            
        }//loop on me0segm

	    ////Loop on Bkg segments to compute the mean distance
	    cout<<"----------------- ME0BkgSegSize="<<ME0SegBkg.size()<<endl;
	    std::vector<int> BkgMult, BkgMult1,BkgMult2,BkgMult3,BkgMult4,BkgMult5,BkgMult6,BkgMult7,BkgMult8,BkgMult9; 
	    std::vector<int> BkgMult10, BkgMult11,BkgMult12,BkgMult13,BkgMult14,BkgMult15,BkgMult16,BkgMult17,BkgMult18;

	    for (uint b=0; b<ME0SegBkg.size();b++){
	      LocalPoint segmLP = ME0SegBkg[b].localPosition();
	      ME0DetId segmDetID = ME0SegBkg[b].me0DetId();
	      GlobalPoint segmGP = (me0geom->chamber(segmDetID))->toGlobal(segmLP);
	      int chB=segmDetID.chamber(); 
	      int reB=segmDetID.region();
	      if(reB>0){
	      if(chB==1)  BkgMult1.push_back(1);
	      if(chB==2)  BkgMult2.push_back(1);
	      if(chB==3)  BkgMult3.push_back(1);
	      if(chB==4)  BkgMult4.push_back(1);
	      if(chB==5)  BkgMult5.push_back(1);
	      if(chB==6)  BkgMult6.push_back(1);
	      if(chB==7)  BkgMult7.push_back(1);
	      if(chB==8)  BkgMult8.push_back(1);
	      if(chB==9)  BkgMult9.push_back(1);
	      if(chB==10) BkgMult10.push_back(1);
	      if(chB==11) BkgMult11.push_back(1);
	      if(chB==12) BkgMult12.push_back(1);
	      if(chB==13) BkgMult13.push_back(1);
	      if(chB==14) BkgMult14.push_back(1);
	      if(chB==15) BkgMult15.push_back(1);
	      if(chB==16) BkgMult16.push_back(1);
	      if(chB==17) BkgMult17.push_back(1);
	      if(chB==18) BkgMult18.push_back(1);
	      }

        for (uint a=b+1; a<ME0SegBkg.size();a++){
		if(a==b) continue;
		LocalPoint segmLPa = ME0SegBkg[a].localPosition();
		ME0DetId segmDetIDa = ME0SegBkg[a].me0DetId();
		GlobalPoint segmGPa = (me0geom->chamber(segmDetIDa))->toGlobal(segmLPa);
		int chA=segmDetIDa.chamber();
		//nSegmPerChamber.clear();

		if( (chB==chA) && (segmDetID.region()==segmDetIDa.region())){
		double dist = TMath::Sqrt((segmGPa.x()-segmGP.x())*(segmGPa.x()-segmGP.x()) + (segmGPa.y()-segmGP.y())*(segmGPa.y()-segmGP.y()) );
		hDistME0Seg->Fill(dist);
		BkgMult.push_back(1);
		nSegmPerChamber +=1;
		/*if((ME0SegBkg[a].tof<25 && ME0SegBkg[a].tof>0)  && (ME0SegBkg[b].tof<25 && ME0SegBkg[b].tof>0)){
		  cout<<"A:"<<a<<" chamber="<<chA<<" eta="<<segmGPa.eta()<<" phi="<<segmGPa.phi()<<endl;
		  cout<<"B:"<<b<<" chamber="<<chB<<" eta="<<segmGP.eta()<<" phi="<<segmGP.phi()<<endl;
		  }*/
		}
	      }
	    }
	    /*	    cout<<" nSegmPerChamber="<<nSegmPerChamber<<" bkg per ch="<<BkgMult.size()<<endl;
	    cout<<" bkg per ch 1="<<BkgMult1.size()<<endl;
	    cout<<" bkg per ch 2="<<BkgMult2.size()<<endl;
	    cout<<" bkg per ch 3="<<BkgMult3.size()<<endl;
	    cout<<" bkg per ch 4="<<BkgMult4.size()<<endl;
	    cout<<" bkg per ch 5="<<BkgMult5.size()<<endl;
	    cout<<" bkg per ch 6="<<BkgMult6.size()<<endl;
	    cout<<" bkg per ch 7="<<BkgMult7.size()<<endl;
	    cout<<" bkg per ch 8="<<BkgMult8.size()<<endl;
	    cout<<" bkg per ch 9="<<BkgMult9.size()<<endl;
	    cout<<" bkg per ch 10="<<BkgMult10.size()<<endl;
	    cout<<" bkg per ch 11="<<BkgMult11.size()<<endl;
	    cout<<" bkg per ch 12="<<BkgMult12.size()<<endl;
	    cout<<" bkg per ch 13="<<BkgMult13.size()<<endl;
	    cout<<" bkg per ch 14="<<BkgMult14.size()<<endl;
	    cout<<" bkg per ch 15="<<BkgMult15.size()<<endl;
	    cout<<" bkg per ch 16="<<BkgMult16.size()<<endl;
	    cout<<" bkg per ch 17="<<BkgMult17.size()<<endl;
	    cout<<" bkg per ch 18="<<BkgMult18.size()<<endl;*/
	    double nMeanBkgSegmPerCh=(BkgMult1.size() + BkgMult2.size() + BkgMult3.size() + BkgMult4.size() + BkgMult5.size()      + BkgMult6.size()+ BkgMult7.size()+ BkgMult8.size()+ BkgMult9.size()+ BkgMult10.size()+ BkgMult11.size() + BkgMult12.size() + BkgMult13.size() + BkgMult14.size() + BkgMult15.size()+ BkgMult16.size() + + BkgMult17.size() + + BkgMult18.size())/18;
	    hSegmentPerChamber->Fill(nMeanBkgSegmPerCh);
	    hSegmentPerChamber2->SetBinContent(1,BkgMult1.size());
	    hSegmentPerChamber2->SetBinContent(2,BkgMult2.size());
	    hSegmentPerChamber2->SetBinContent(3,BkgMult3.size());
	    hSegmentPerChamber2->SetBinContent(4,BkgMult4.size());
	    hSegmentPerChamber2->SetBinContent(5,BkgMult5.size());
	    hSegmentPerChamber2->SetBinContent(6,BkgMult6.size());
	    hSegmentPerChamber2->SetBinContent(7,BkgMult7.size());
	    hSegmentPerChamber2->SetBinContent(8,BkgMult8.size());
	    hSegmentPerChamber2->SetBinContent(9,BkgMult9.size());
	    hSegmentPerChamber2->SetBinContent(10,BkgMult10.size());
	    hSegmentPerChamber2->SetBinContent(11,BkgMult11.size());
	    hSegmentPerChamber2->SetBinContent(12,BkgMult12.size());
	    hSegmentPerChamber2->SetBinContent(13,BkgMult13.size());
	    hSegmentPerChamber2->SetBinContent(14,BkgMult14.size());
	    hSegmentPerChamber2->SetBinContent(15,BkgMult15.size());
	    hSegmentPerChamber2->SetBinContent(16,BkgMult16.size());
	    hSegmentPerChamber2->SetBinContent(17,BkgMult17.size());
	    hSegmentPerChamber2->SetBinContent(18,BkgMult18.size());
	    
	    ///////Check on matched signal segment::if there is more than one ME0Segment 
	    ///////associated to a given SimTrack. If yes, pick the one with the highest
	    ///////numbero of associated hits.
	    std::vector<uint> OneAssME0Seg;        std::vector<uint> DuplicatesCont;  std::vector<uint> DuplicatesID;
	     std::vector<uint> ME0SegIDSorted;  std::vector<uint> ME0SegIDSortedTrackID;
 

        
        for (uint s=0;s<matchedSegmDB.size();s++){
            cout<<"++++++++segID="<<matchedSegmDB[s].SegID<<" ,nMatchedHit="<<matchedSegmDB[s].nHitMatched<<" matchedSimTrackID="<<matchedSegmDB[s].matchedTrackID<<endl;
            for (uint t=s+1;t<matchedSegmDB.size();t++){
                if(matchedSegmDB[t].matchedTrackID==matchedSegmDB[s].matchedTrackID){
                    cout<<"------------->Found double match with id="<<matchedSegmDB[t].SegID<<endl;
                    if(matchedSegmDB[s].nHitMatched < matchedSegmDB[t].nHitMatched)  {
                        OneAssME0Seg.push_back(matchedSegmDB[s].SegID);
                        ME0SegIDSorted.push_back(matchedSegmDB[t].SegID);
                      //  ME0SegIDSortedTrackID.push_back(matchedSegmDB[t].matchedTrackID);
                    }
                    
                    if(matchedSegmDB[t].nHitMatched < matchedSegmDB[s].nHitMatched)  {
                        OneAssME0Seg.push_back(matchedSegmDB[t].SegID);
                        ME0SegIDSorted.push_back(matchedSegmDB[s].SegID);
                      //  ME0SegIDSortedTrackID.push_back(matchedSegmDB[s].matchedTrackID);
                    }
                    if(matchedSegmDB[t].nHitMatched == matchedSegmDB[s].nHitMatched){
                        OneAssME0Seg.push_back(matchedSegmDB[t].SegID);
                        ME0SegIDSorted.push_back(matchedSegmDB[t].SegID);
                       // ME0SegIDSortedTrackID.push_back(matchedSegmDB[t].matchedTrackID);
                    }
		  }
		}}
	      
        for(uint f=0;f<OneAssME0Seg.size();f++){cout<<"++++++++#Before: ToBe erased SegID="<<OneAssME0Seg[f]<<endl;}
        //for(uint f=0;f<ME0SegIDSorted.size();f++){cout<<"++++++++#Before: Sorted SegID="<<ME0SegIDSorted[f]<<endl;}
        cout<<"MatchedSegm="<<matchedSegmDB.size()<<" #Number of duplicates="<<OneAssME0Seg.size()+ DuplicatesCont.size()<<endl;
        //hNumberOfDuplicatesPerSignalSimTrack->Fill(OneAssME0Seg.size()+DuplicatesCont.size());
        cout<<"controllo i doppioni nella lista di segmenti da cancellare:"<<endl;
        for(uint d=0;d<OneAssME0Seg.size();d++){
            for(uint e=d+1;e<OneAssME0Seg.size();e++){
                if(OneAssME0Seg[d]==OneAssME0Seg[e])  {cout<<" doppione con segID="<<OneAssME0Seg[d]<<endl;  OneAssME0Seg.erase(OneAssME0Seg.begin()+d);}
            }
        }
        for(uint f=0;f<OneAssME0Seg.size();f++){cout<<"After cleaning doppione con segID="<<OneAssME0Seg[f]<<endl;}
        //for(uint f=0;f<DuplicatesCont.size();f++){cout<<"+++++++#After: Segment with same num of matched hitID="<<DuplicatesCont[f]<<endl;}

        
        for (uint h=0; h<matchedSegmDB.size();h++){
             cout<<"Before cleaning: segm in DB id="<<matchedSegmDB[h].SegID<<endl;
            for(uint f=0;f<OneAssME0Seg.size();f++){
                //if(b==OneAssME0Seg[f]) ME0SegSig.erase(ME0SegSig.begin()+b);
                
                int signedID = (int) OneAssME0Seg[f];
                if(matchedSegmDB[h].SegID==signedID) {
                   cout<<" duplicate id="<<signedID<<" segm in DB id="<<matchedSegmDB[h].SegID<<endl;
                    matchedSegmDB.erase(matchedSegmDB.begin()+h);
                }
		}
          /*  if(DuplicatesCont.size()==1 && matchedSegmDB[h].SegID== (int)DuplicatesCont[0]) matchedSegmDB.erase(matchedSegmDB.begin()+DuplicatesCont[0]);
            if(DuplicatesCont.size()>1 && matchedSegmDB[h].SegID== (int)DuplicatesCont[0]) matchedSegmDB.erase(matchedSegmDB.begin()+DuplicatesCont[0]);
            */
        }
        cout<<" ME0SegSig.size="<<ME0SegSig.size()<<" matchedSegmDB.size="<<matchedSegmDB.size()<<" ME0Tracks.size="<<ME0Tracks.size()<<endl;
        for (uint b=0; b<matchedSegmDB.size();b++){cout<<"ME0SegID after cleaning="<<matchedSegmDB[b].SegID<<" simTrackid="<<matchedSegmDB[b].matchedTrackID<<endl;}
        
        ///Fill histos for best matched me0segment-to-simTrack
        std::vector<int> NumberOfDuplicatesPerSignalSimTrack;
        std::vector<uint> ME0SegIDToSim;
        if(matchedSegmDB.size() !=0){
            for(uint r=0; r<ME0Tracks.size(); r++){
                NumberOfDuplicatesPerSignalSimTrack.clear();
                for (uint s=0;s<matchedSegmDB.size();s++){
                    if(ME0Tracks.at(r).trackId() == matchedSegmDB[s].matchedTrackID) {
                        NumberOfDuplicatesPerSignalSimTrack.push_back(1);
                        ME0SegIDToSim.push_back(matchedSegmDB[s].SegID);}
            }
            cout<<"SimTrackID="<<ME0Tracks.at(r).trackId()<<" #duplicates="<<NumberOfDuplicatesPerSignalSimTrack.size()<<endl;
                if(NumberOfDuplicatesPerSignalSimTrack.size()>0){
                    hSimPtDen->Fill(ME0Tracks.at(r).momentum().pt());
                    hSimEtaDen->Fill(ME0Tracks.at(r).momentum().eta());
                    hSimPhiDen->Fill(ME0Tracks.at(r).momentum().phi());
                    hME0SegmPullPhi_BestMatch->Fill(ME0Tracks.at(r).momentum().phi()-matchedSegmDB[0].SegPhi);
                    hME0SegmPullEta_BestMatch->Fill(ME0Tracks.at(r).momentum().eta()- matchedSegmDB[0].SegEta);
                    hME0Segm_nHit_sig_BestMatch->Fill(matchedSegmDB[0].SegNhits);
                    hME0Segm_chi2_sig_BestMatch->Fill(matchedSegmDB[0].SegChi2);
                    hME0Segm_time_sig_BestMatch->Fill(matchedSegmDB[0].SegTime);
                    
                }
            hNumberOfDuplicatesPerSignalSimTrack->Fill(NumberOfDuplicatesPerSignalSimTrack.size());
            
       }}
        ///Fill histos for the non-signal segment  for each simTrack
        for(uint d=0;d<NumberOfDuplicatesPerSignalSimTrack.size();d++){
            for (uint h=0; h<matchedSegmDB.size();h++){
                if(NumberOfDuplicatesPerSignalSimTrack.at(d)==matchedSegmDB[h].SegID){
                    hDuplicatesEta->Fill(matchedSegmDB[h].SegEta);
                    hDuplicatesPhi->Fill(matchedSegmDB[h].SegPhi);}
            }}
        
      /*  hNME0SegmentSignal_BestMatch->Fill(matchedSegmDB.size());
        for(uint i=0;i<matchedSegmDB.size();i++){
            if(matchedSegmDB[i].matchedTrackID){
                //cout<<"+++Check: matchedSimTrackID="<<matchedSegmDB[i].matchedTrackID<<" pt="<<ME0Tracks.at(matchedSegmDB[i].matchedTrackID).momentum().pt()<<endl;
                hSimPtDen->Fill(ME0Tracks.at(i).momentum().pt());
                hME0SegmPullPhi_BestMatch->Fill(ME0Tracks.at(i).momentum().phi()-matchedSegmDB[i].SegPhi);
                hME0SegmPullEta_BestMatch->Fill(ME0Tracks.at(i).momentum().eta()- matchedSegmDB[i].SegEta);
                
                hME0Segm_nHit_sig_BestMatch->Fill(matchedSegmDB[i].SegNhits);
                hME0Segm_chi2_sig_BestMatch->Fill(matchedSegmDB[i].SegChi2);
                hME0Segm_time_sig_BestMatch->Fill(matchedSegmDB[i].SegTime);
            }
	      }*/

        
        ////////////////////////////////////////Loop over ME0RU segments
   /*     ME0SegmentCollection::const_iterator me0segItRU;
        std::vector<int> nMatchedHitRU;
        std::vector<int> nTrackIDRU;
        int segItRU =0;
        
        matchedSegm_t matchedSegmRU;
        std::vector<matchedSegm_t>    matchedSegmDBRU;
        for(me0segItRU=me0SegmentsRU->begin(); me0segItRU<me0SegmentsRU->end(); me0segItRU++ ){
            nMatchedHitRU.clear(); nTrackIDRU.clear();
            auto me0rhs =  me0segItRU->specificRecHits();
            cout<<"RU: ME0Segment nHits="<<me0segItRU->nRecHits()<<" tof="<<me0segItRU->time()<<" chi2="<<me0segItRU->chi2()<<endl;
            for (auto rh = me0rhs.begin(); rh!= me0rhs.end(); rh++){
                auto rhLP = rh->localPosition();
                auto me0id = rh->me0Id();
                auto rhr = me0geom->chamber(me0id);
                auto rhGP = rhr->toGlobal(rhLP);
                double globalPhi = rhGP.phi();  double globalX = rhGP.x();			double globalY = rhGP.y();
                for(ME0DigiPreRecoCollection::DigiRangeIterator cItr = me0_digis->begin(); cItr != me0_digis->end(); ++cItr){
                    ME0DetId id = (*cItr).first;
                    Short_t  me0_digiLayer =  id.layer();
                    Short_t  me0_digiChamber =  id.chamber();
                    const GeomDet* gdet = me0Geom->idToDet(id);
                    const BoundPlane & surface = gdet->surface();
                    ME0DigiPreRecoCollection::const_iterator digiItr;
                    for (digiItr = (*cItr ).second.first; digiItr != (*cItr ).second.second; ++digiItr){
                        LocalPoint lp(digiItr->x(), digiItr->y(), 0);
                        GlobalPoint gp = surface.toGlobal(lp);
                        Float_t me0_digiX = gp.x();
                        Float_t me0_digiY = gp.y();
                        Float_t me0_digiEta = gp.eta();
                        Float_t me0_digiPhi = gp.phi();
                        Float_t me0_digiLocalX = lp.x();
                        Float_t me0_digiLocalY = lp.y();
                        if(globalX == me0_digiX && globalY == me0_digiY && me0_digiLayer==me0id.layer()){
                        std::cout<<" DIGI global phi "<<me0_digiPhi<<" global eta "<<me0_digiEta<<"  localX "<< me0_digiLocalX << " localY="<<  me0_digiLocalY  <<" layer "<<me0_digiLayer<<" pdg="<<digiItr->pdgid()<<"  tof="<<digiItr->tof()<<std::endl;
                            //////////////////////////////Match digi in ME0Segment with simHit of muon simTtracks///////
                            for (edm::PSimHitContainer::const_iterator itHit =muonSimME0Hits.begin(); itHit != muonSimME0Hits.end(); ++itHit)
                            {
                                ME0DetId idme0 = ME0DetId(itHit->detUnitId());
                                int layer_sim = idme0.layer();
                                int chamber_sim = idme0.chamber();
                                LocalPoint lp = itHit->entryPoint();
                                GlobalPoint hitGP_sim( me0geom->idToDet( itHit->detUnitId())->surface().toGlobal(lp));
                                //cout<<"simhit pdgID="<<itHit->particleType()<<" eta="<<hitGP_sim.eta()<<" phi="<<hitGP_sim.phi()<<endl;
                                if((itHit->particleType()==digiItr->pdgid()) && (layer_sim==me0_digiLayer) && (fabs(lp.x()-me0_digiLocalX)<=0) && (fabs(lp.y()-me0_digiLocalY)<=0.) && (me0_digiChamber ==chamber_sim)   ){
                                    //std::cout<<"****signalMuon hits eta="<<gp.eta()<<"  phi="<<gp.phi()<<" localX="<<endl;
                                   // std::cout<<"RU SIMHIT global phi "<< hitGP_sim.phi()<<" global eta "<<hitGP_sim.eta()<<"  localX "<< lp.x() << " localY="<< lp.y() <<" layer "<<layer_sim<<std::endl;
                                    //std::cout<<"RU DIGI global phi "<<me0_digiPhi<<" global eta "<<me0_digiEta<<"  localX "<<me0_digiLocalX  << " localY="<<  me0_digiLocalY  <<" layer "<<me0_digiLayer<<" pdg="<<digiItr->pdgid()<<"  tof="<<digiItr->tof()<<std::endl;
                                    nMatchedHitRU.push_back(1);
                                    nTrackIDRU.push_back(itHit->trackId());}
                            }
                            
                        }//match digi-rechit
                        
                        }	}//loop over digi
                
                }//loop over rechit in me0segm
            ///////////////////////////////////SIGNAL////////////////////////////////
            if(nMatchedHitRU.size()>=3 && nMatchedHitRU.size()<=6){
                segItRU++;
                matchedSegmRU.nHitMatched = nMatchedHitRU.size();
                matchedSegmRU.SegID = segItRU;
                matchedSegmRU.matchedTrackID=nTrackIDRU[0];
                std::cout<<"RU Signal matched hit= "<<nMatchedHitRU.size()<<endl;
                //matchedSegmRU.SegEta =ME0SegmentGP_sig.eta();
                //matchedSegmRU.SegPhi =ME0SegmentGP_sig.phi();
                matchedSegmDBRU.push_back(matchedSegmRU);

            }
            
        }
        ////////////////////////////////////////End loop over RUSegments//check duplicates in RUSegments
        std::vector<uint> OneAssME0SegRU;
        //for(uint r=0; r<ME0Tracks.size(); r++){
        //  cout<<"++++++++Matched ME0SimTrack pt="<<ME0Tracks.at(r).momentum().pt()<<" id="<<ME0Tracks.at(r).trackId()<<endl;
        //if(ME0Tracks.at(r).trackId() == matchedSegmDB[s].matchedTrackID) {
        
        for (uint s=0;s<matchedSegmDBRU.size();s++){
            cout<<"++++++++segID="<<matchedSegmDBRU[s].SegID<<" ,nMatchedHit="<<matchedSegmDBRU[s].nHitMatched<<" matchedSimTrackID="<<matchedSegmDBRU[s].matchedTrackID<<endl;
            for (uint t=s+1;t<matchedSegmDBRU.size();t++){
                if(matchedSegmDBRU[t].matchedTrackID==matchedSegmDBRU[s].matchedTrackID){
                    cout<<"------------->RU::Found double match with id="<<matchedSegmDBRU[t].SegID<<endl;
                    if(matchedSegmDBRU[s].nHitMatched < matchedSegmDBRU[t].nHitMatched)  OneAssME0SegRU.push_back(matchedSegmDBRU[s].SegID);
                    if(matchedSegmDBRU[t].nHitMatched < matchedSegmDBRU[s].nHitMatched)  OneAssME0SegRU.push_back(matchedSegmDBRU[t].SegID);
                    if(matchedSegmDBRU[t].nHitMatched == matchedSegmDBRU[s].nHitMatched)  OneAssME0SegRU.push_back(matchedSegmDBRU[s].SegID);
                }
            }}
        
        //for(uint f=0;f<OneAssME0Seg.size();f++){cout<<"++++++++Ass SegID="<<OneAssME0Seg[f]<<endl;}
        cout<<"RU #MatchedSegm="<<matchedSegmDBRU.size()<<" #Number of duplicates="<<OneAssME0SegRU.size()<<endl;
        for (uint h=0; h<matchedSegmDBRU.size();h++){
            cout<<"RU Before cleaning: segm in DB id="<<matchedSegmDBRU[h].SegID<<endl;
            for(uint f=0;f<OneAssME0SegRU.size();f++){
                //if(b==OneAssME0Seg[f]) ME0SegSig.erase(ME0SegSig.begin()+b);
                
                int signedIDRU = (int) OneAssME0SegRU[f];
                if(matchedSegmDBRU[h].SegID==signedIDRU) {
                    cout<<"RU duplicate id="<<signedIDRU<<" segm in DB id="<<matchedSegmDBRU[h].SegID<<endl;
                    matchedSegmDBRU.erase(matchedSegmDBRU.begin()+h);
                }
            }}
        cout<<"RU matchedSegmDB.size="<<matchedSegmDBRU.size()<<" ME0Tracks.size="<<ME0Tracks.size()<<endl;
        for (uint b=0; b<matchedSegmDBRU.size();b++)cout<<"RU ME0SegID after cleaning="<<matchedSegmDBRU[b].SegID<<endl;
        
 

        if(matchedSegmDBRU.size() < matchedSegmDB.size()) {
            cout<<"**************** FOUND EVENT with no matched RU segment! "<<" EVENT="<<iEvent.id().event()<<" run:"<<iEvent.id().run()<<endl;
            
        }
        */
        if( nSegmDPhiPt5GeV.size()>0) hEventPass5GeVCut->Fill(1);
        if( nSegmDPhiPt30GeV.size()>0) hEventPass30GeVCut->Fill(1);
        
        if( nSegmDPhiPt5GeV_AllCut.size()>0) hEventPass5GeVCutAllCut->Fill(1);
        if( nSegmDPhiPt30GeV_AllCut.size()>0) hEventPass30GeVCutAllCut->Fill(1);
        
        hNME0SegmentSignal->Fill(nME0SegmentSignal.size());
        hNME0SegmentBkg->Fill(nME0SegmentBkg.size());
        
        if(nME0SegmentSignal_timeCut.size()) hNME0SegmentSignal_timeCut->Fill(nME0SegmentSignal_timeCut.size());
        if(nME0SegmentBkg_timeCut.size())	hNME0SegmentBkg_timeCut->Fill(nME0SegmentBkg_timeCut.size());
        
        
        for(uint m=0; m<SignalTime.size(); m++ ){
            for(uint k=0; k<BkgTime.size(); k++ ){
                //	std::cout<<"(s="<<m<<" b="<<k<<") signal -bkg time="<<SignalTime.at(m)-BkgTime.at(k)<<endl;
                hTimeDiffSigvsBkg->Fill(SignalTime.at(m)-BkgTime.at(k));
            }
        }
        
        ///////////////////////////////////////////Loop over simTracks///////////////////////////
        
        //	SimTrackContainer ME0Tracks, ME0EleSimTracks;
        SimTrackContainer::const_iterator bkgSimTrackIt;
        SimTrackContainer::const_iterator simTrackIT;
        for (bkgSimTrackIt = simTracks->begin(); bkgSimTrackIt != simTracks->end(); ++bkgSimTrackIt){
            
            //double simPt=(*bkgSimTrackIt).momentum().pt();
            double simEta = (*bkgSimTrackIt).momentum().eta();
            if (abs(simEta) > maxEta_ || abs(simEta) < minEta_) continue;
            //edm::PSimHitContainer selAllME0SimHits = isTrackMatched(bkgSimTrackIt, iEvent , iSetup);
            
           // if( selAllME0SimHits.size() ==0 ) continue;
            //cout<<"All ME0 simTrackID:"<<bkgSimTrackIt->trackId()<<" pdgID="<<(*bkgSimTrackIt).type()<<" pT="<<simPt<<" carica="<<(*bkgSimTrackIt).charge()<<endl;
            
        }
        /*
         for(uint r=0; r<ME0Tracks.size(); r++){
         cout<<"Muon ME0SimTrack Eta: "<<ME0Tracks.at(r).momentum().eta()<<" type="<<ME0Tracks.at(r).type()<<" ID="<<ME0Tracks.at(r).trackId()<<endl;
         }*/
        
        ///////////////////////////////////////Loop over gen particles/////////////////////////////
        
        
        for(unsigned int i = 0; i < genparticles->size();i++) {
            if((abs(genparticles->at(i).pdgId()) == 13) && (genparticles->at(i).status() == 1) ) indexgenmu.push_back(i);
        }//filtro eventi
        
        //	std::cout<< "N gen muon " << indexgenmu.size() << " N me0 muon " << OurMuons->size() << std::endl;
        
        
        int counterZmm =0;
        if (indexgenmu.size()>0) {counterZmm ++; hNzmm->Fill(1);}
        
        std::vector<int>  indexGenMuInME0;
        std::vector<int>  indexGenMuElseWhere;
        std::vector<int>  indexRecoMuElseWhere;
        TLorentzVector genmu1, genmu2, genZ;
        TLorentzVector recomu1, recomu2, recoZ;
        TLorentzVector recomu1Intime, recomu2Intime, recoZIntime;
        
        
        for(uint i =0; i<indexgenmu.size(); i++){
            //std::cout<<i<<" particle= "<<genparticles->at(indexgenmu[i]).pdgId()<<" status="<<genparticles->at(indexgenmu[i]).status()<<" eta="<<genparticles->at(indexgenmu[i]).eta()<<" phi="<<genparticles->at(indexgenmu[i]).phi()<<" pt="<<genparticles->at(indexgenmu[i]).pt()<<std::endl;
            hGenMuPt->Fill( (genparticles->at(indexgenmu[i]).pt()));
            hGenMuEta->Fill( (genparticles->at(indexgenmu[i]).eta())  );
            
            if ( indexgenmu.size() ==2){
                genmu1.SetPtEtaPhiM(genparticles->at(indexgenmu[0]).pt(),genparticles->at(indexgenmu[0]).eta(), genparticles->at(indexgenmu[0]).phi(), 0.105   );
                genmu2.SetPtEtaPhiM(genparticles->at(indexgenmu[1]).pt(),genparticles->at(indexgenmu[1]).eta(), genparticles->at(indexgenmu[1]).phi(), 0.105   );
                genZ = genmu1 +genmu2;
                //std::cout<<" gen pt1="<<genmu1.Pt()<<" eta="<<genmu1.Eta()<<std::endl;
                //std::cout<<" gen pt2="<<genmu2.Pt()<<" eta="<<genmu2.Eta()<<std::endl;
                hGenMass->Fill(genZ.M());
            }
            
            if((abs(genparticles->at(indexgenmu[i]).eta())<maxEta_ ) && ( (genparticles->at(indexgenmu[i]).eta() > minEta_ )||(genparticles->at(indexgenmu[i]).eta()< (-minEta_ )) )  ){
                indexGenMuInME0.push_back(indexgenmu[i]);
                //std::cout<<i<<" particle= "<<genparticles->at(indexgenmu[i]).pdgId()<<" status="<<genparticles->at(indexgenmu[i]).status()<<" eta="<<genparticles->at(indexgenmu[i]).eta()<<" phi="<<genparticles->at(indexgenmu[i]).phi()<<" pt="<<genparticles->at(indexgenmu[i]).pt()<<std::endl;
            }else{
                indexGenMuElseWhere.push_back(indexgenmu[i]);
            }
        }
        
        hGenMuonsME0->Fill(indexGenMuInME0.size());
        hNGenMu->Fill(indexgenmu.size());
      //  hNME0Mu->Fill(OurMuons->size());
        
        //std::cout<< "N gen muon in me0 " << indexGenMuInME0.size() <<""<< std::endl;
        unsigned int k = 0;
        std::vector<double>   me0SegAbsDPhi;  std::vector<double>   me0SegDPhi;
        std::vector<double>   me0SegAbsDEta;  std::vector<double>   me0SegDEta;
        std::vector<int>  idxtmpreco;
        std::vector<uint>  idxtmpgen;
        //find all the me0 muon next to the genMuon in DR<0.1
        
        
        hNEvZmm->Fill(counterZmm);
        
        
        
        //////////////////////////////////////////////////////BeamSpot and VTX plot//////////////////////////////////////////////////////
       /* if ( beamSpotHandle.isValid() ){
            beamSpot = *beamSpotHandle;
            hBeamSpotX0->Fill(beamSpot.x0());
            hBeamSpotY0->Fill(beamSpot.y0());
            hBeamSpotX0Y0->Fill(beamSpot.x0(),beamSpot.y0());
            hBeamSpotZ0->Fill(beamSpot.z0());
            hBeamSpotSigmaZ->Fill( beamSpot.sigmaZ());
            hBeamSpotdxdz->Fill(beamSpot.dxdz());
            hBeamSpotBeamWidthX->Fill(beamSpot.BeamWidthX());
            hBeamSpotBeamWidthY->Fill(beamSpot.BeamWidthY());
            
        }*/
        //	if (primaryVertices.isValid()) {
        
        //}
        
       /* edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
        iEvent.getByLabel("addPileupInfo", PupInfo);
        
        std::vector<PileupSummaryInfo>::const_iterator PVI;
        for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
            hBX->Fill(PVI->getBunchCrossing());
            hNPU->Fill(PVI->getPU_NumInteractions());
            hTrueInt->Fill(PVI->getTrueNumInteractions());
        }*/
        
        //////////////////////////////////////////////////////BeamSpot and VTX plot//////////////////////////////////////////////////////
        //	}
        
        

        
    
    }
    
    

}




// ------------ method called once each job just before starting event loop  ------------
//void ME0SegmentAnalyzerMuonGun::beginRun(edm::Run const&, edm::EventSetup const&) {
void ME0SegmentAnalyzerMuonGun::beginJob()
{
    hNEvZmm = fs->make<TH1F>("hNEvZmm","hNEvZmm",10,0,10);
    hNEvME0 = fs->make<TH1F>("hNEvME0","hNEvME0",10,0,10);
    hNSimHitsME0= fs->make<TH1F>("hNSimHitsME0","hNSimHitsME0",20,0,20);
    hNRecHitsME0= fs->make<TH1F>("hNRecHitsME0","hNRecHitsME0",100,0,100);
    hNBkgHitsME0= fs->make<TH1F>("hNBkgHitsME0","hNBkgHitsME0",60,-10,50);
    
    hNDigiMatchedRH= fs->make<TH1F>("hNDigiMatchedRH","hNDigiMatchedRH",21,-0.5,20.5);
    hNPrompt= fs->make<TH1F>("hNPrompt","hNPrompt",21,-0.5,20.5);
    hN_noPrompt= fs->make<TH1F>("hN_noPrompt","hN_noPrompt",21,-0.5,20.5);
    
    hRatioRecoToSimHits= fs->make<TH1F>("hRatioRecoToSimHits","hRatioRecoToSimHits",1000,0,100);
    
    hNPromptMuHit= fs->make<TH1F>("hNPromptMuHit","hNPromptMuHit",21,-0.5,20.5);
    hNPromptNoMuHit= fs->make<TH1F>("hNPromptNoMuHit","hNPromptNoMuHit",21,-0.5,20.5);
    hNPromptHit_pdgId= fs->make<TH1F>("hNPromptHit_pdgId","hNPromptHit_pdgId",3001,0,3000.5);
    hN_noPromptHit_pdgId= fs->make<TH1F>("hN_noPromptHit_pdgId","hN_noPromptHit_pdgId",3001,0,3000.5);
    hN_noEleHit= fs->make<TH1F>("hN_noEleHit","hN_noEleHit",21,-0.5,20.5);
    
    hSimEta  = fs->make<TH1F>("hSimEta","hSimEta",400,-4,4);
    hSimPhi  = fs->make<TH1F>("hSimPhi","hSimPhi",400,-4,4);
    hSimPt  = fs->make<TH1F>("hSimPt","hSimPt",200, 0,200);
    hSimPtDen = fs->make<TH1F>("hSimPtDen","hSimPtDen",200, 0,200);
    hSimEtaDen  = fs->make<TH1F>("hSimEtaDen","hSimEtaDen",400,-4,4);
    hSimPhiDen  = fs->make<TH1F>("hSimPhiDen","hSimPhiDen",400,-4,4);
    hME0Segm_SimpT_Signal = fs->make<TH1F>(" hME0Segm_SimpT_Signal "," hME0Segm_SimpT_Signal ",200, 0,200);
    hPtVSDphi = fs->make<TH2F>("hPtVSDphi","hPtVSDphi",5000, 0, 0.05 , 200,0,200);

    hME0Segm_DeltaPhi_pT0p5=fs->make<TH1F>("hME0Segm_DeltaPhi_pT0p5","hME0Segm_DeltaPhi_pT0p5",5000, 0, 0.05);
    hME0Segm_DeltaPhi_pT1=fs->make<TH1F>("hME0Segm_DeltaPhi_pT1","hME0Segm_DeltaPhi_pT1",5000, 0, 0.05);
    hME0Segm_DeltaPhi_pT3=fs->make<TH1F>("hME0Segm_DeltaPhi_pT3","hME0Segm_DeltaPhi_pT3",5000, 0, 0.05);
    hME0Segm_DeltaPhi_pT5=fs->make<TH1F>("hME0Segm_DeltaPhi_pT5","hME0Segm_DeltaPhi_pT5",5000, 0, 0.05);
    hME0Segm_DeltaPhi_pT10=fs->make<TH1F>("hME0Segm_DeltaPhi_pT10","hME0Segm_DeltaPhi_pT10",5000, 0, 0.05);
    hME0Segm_DeltaPhi_pT20=fs->make<TH1F>("hME0Segm_DeltaPhi_pT20","hME0Segm_DeltaPhi_pT20",5000, 0, 0.05);
    hME0Segm_DeltaPhi_pT50=fs->make<TH1F>("hME0Segm_DeltaPhi_pT50","hME0Segm_DeltaPhi_pT50",5000, 0, 0.05);

    hSimPtVSDphi = fs->make<TH2F>("hSimPtVSDphi","hSimPtVSDphi",5000, 0, 0.05 , 200,0,200);
    hSimPtVSDeta = fs->make<TH2F>("hSimPtVSDeta","hSimPtVSDeta",10000, 0, 0.1 , 200,0,200);
    hSimDEta = fs->make<TH1F>("hSimDEta","hSimDEta",1000,-0.5,0.5);
    hSimDY = fs->make<TH1F>("hSimDY","hSimDY",1000,-20,20);
    hSimDX = fs->make<TH1F>("hSimDX","hSimDX",1000,-20,20);
    hSimDPhi = fs->make<TH1F>("hSimDPhi","hSimDPhi",1000,-0.05,0.05);
    hSimDPhiPos = fs->make<TH1F>("hSimDPhiPos","hSimDPhiPos",1000,-0.05,0.05);
    hSimDPhiNeg = fs->make<TH1F>("hSimDPhiNeg","hSimDPhiNeg",1000,-0.05,0.05);
    
    hSimDPhiPos_LowPt = fs->make<TH1F>("hSimDPhiPos_LowPt","hSimDPhiPos_LowPt",1000,-0.05,0.05);
    hSimDPhiPos_HighPt = fs->make<TH1F>("hSimDPhiPos_HighPt","hSimDPhiPos_HighPt",1000,-0.05,0.05);
    
    hSimDPhiNeg_LowPt = fs->make<TH1F>("hSimDPhiNeg_LowPt","hSimDPhiNeg_LowPt",1000,-0.05,0.05);
    hSimDPhiNeg_HighPt = fs->make<TH1F>("hSimDPhiNeg_HighPt","hSimDPhiNeg_HighPt",1000,-0.05,0.05);
    
    hDPhi = fs->make<TH1F>("hDPhi","hDPhi",1000,-0.05,0.05);
    hDPhiPos = fs->make<TH1F>("hDPhiPos","hDPhiPos",1000,-0.05,0.05);
    hDPhiNeg = fs->make<TH1F>("hDPhiNeg","hDPhiNeg",1000,-0.05,0.05);
    
    hDPhiNoMatch= fs->make<TH1F>("hDPhiNoMatch","hDPhiNoMatch",1000,-0.05,0.05);
    hDPhiNoMatch_TimeWindow= fs->make<TH1F>("hDPhiNoMatch_TimeWindow","hDPhiNoMatch_TimeWindow",1000,-0.05,0.05);
    hDPhiNoMatch_TimeWindowTightID= fs->make<TH1F>("hDPhiNoMatch_TimeWindowTightID","hDPhiNoMatch_TimeWindowTightID",1000,-0.05,0.05);
    
    hDPhiPos_LowPt = fs->make<TH1F>("hDPhiPos_LowPt","hDPhiPos_LowPt",1000,-0.05,0.05);
    hDPhiPos_HighPt = fs->make<TH1F>("hDPhiPos_HighPt","hDPhiPos_HighPt",1000,-0.05,0.05);
    
    hDPhiNeg_LowPt = fs->make<TH1F>("hDPhiNeg_LowPt","hDPhiPos_NegPt",1000,-0.05,0.05);
    hDPhiNeg_HighPt = fs->make<TH1F>("hDPhiNeg_HighPt","hDPhiNeg_HighPt",1000,-0.05,0.05);
    
    hDPhiMatchByHitsME0Seg_Signal= fs->make<TH1F>("hDPhiMatchByHitsME0Seg_Signal","hDPhiMatchByHitsME0Seg_Signal",1000,-0.05,0.05);
    hDPhiMatchByHitsME0Seg_Bkg= fs->make<TH1F>("DPhiMatchByHitsME0Seg_Bkg","DPhiMatchByHitsME0Seg_Bkg",1000,-0.05,0.05);
    
    hDPhiMatchByHitsME0Seg_Signal_timeCut = fs->make<TH1F>("hDPhiMatchByHitsME0Seg_Signal_timeCut","hDPhiMatchByHitsME0Seg_Signal_timeCut",1000,-0.05,0.05);
    hDPhiMatchByHitsME0Seg_Bkg_timeCut = fs->make<TH1F>("DPhiMatchByHitsME0Seg_Bkg_timeCut","DPhiMatchByHitsME0Seg_Bkg_timeCut",1000,-0.05,0.05);
    
    
    
    hTimeMatchByHitsME0Seg_Signal=fs->make<TH1F>("hTimeMatchByHitsME0Seg_Signal","hTimeMatchByHitsME0Seg_Signal",4000,0,100);
    hTimeMatchByHitsME0Seg_Bkg=fs->make<TH1F>("hTimeMatchByHitsME0Seg_Bkg","hTimeMatchByHitsME0Seg_Bkg",4000,0,100);
    
    hDEta = fs->make<TH1F>("hDEta","hDEta",1000,-0.5,0.5);
    hDX = fs->make<TH1F>("hDX","hDX",1000,-20,20);
    hAbsDPhi = fs->make<TH1F>("hAbsDPhi","hAbsDPhi",1000,0.,0.5);
    hAbsDEta = fs->make<TH1F>("hAbsDEta","hAbsDEta",1000,0.,0.5);
    
    
    hMuonDigiDPhiVsPT= fs->make<TH2F>("hMuonDigiDPhiVsPT","hMuonDigiDPhiVsPT", 5000, 0, 0.05 , 200,0,200);
    hPosMuonDigiDPhiVsPT= fs->make<TH2F>("hPosMuonDigiDPhiVsPT","hPosMuonDigiDPhiVsPT", 5000, 0, 0.05 , 200,0,200);
    hNegMuonDigiDPhiVsPT= fs->make<TH2F>("hNegMuonDigiDPhiVsPT","hNegMuonDigiDPhiVsPT", 5000, 0, 0.05 , 200,0,200);
    
    hMuonDigiDXVsPT= fs->make<TH2F>("hMuonDigiDXVsPT","hMuonDigiDXVsPT", 5000, 0, 20. , 200,0,200);
    hMuonDigiLocalDPhiVsPT= fs->make<TH2F>("hMuonDigiLocalDPhiVsPT","hMuonDigiLocalDPhiVsPT", 5000, 0, 0.1 , 200,0,200);
    hMuonDigiLocalDXVsPT= fs->make<TH2F>("hMuonDigiLocalDXVsPT","hMuonDigiLocalDXVsPT", 5000, 0, 20. , 200,0,200);
    hNoEleDigiDPhiVsPT= fs->make<TH2F>("hNoEleDigiDPhiVsPT","hNoEleDigiDPhiVsPT", 5000, 0, 0.05 , 200,0,200);
    hPtVSDEta = fs->make<TH2F>("hPtVSDeta","hPtVSDeta",5000, 0, 0.1 , 200,0,200);
    
    
    hme0machtMuonPt = fs->make<TH1F>("hme0machtMuonPt","hme0machtMuonPt",200,0,200);
    hme0machtMuonEta= fs->make<TH1F>("hme0machtMuonEta","hme0machtMuonEta",200,-4,4);
    hme0machtMuonPhi= fs->make<TH1F>("hme0machtMuonPhi","hme0machtMuonPhi",200,-4,4);
    hme0machtMuonCharge= fs->make<TH1F>("hme0machtMuonCharge","hme0machtMuonCharge",10,-5,5);
    hNME0Time = fs->make<TH1F>("hNME0Time","hNME0Time",2000,-200,200);
    hNME0RecHits = fs->make<TH1F>("hNME0RecHits","hNME0RecHits",100,0,100);
    
    hNzmm  = fs->make<TH1F>("hNzmm","hNzmm",10,0,10);
    hNEv  = fs->make<TH1F>("hNEv","hNEv",10,0,10);
    hNGenMu = fs->make<TH1F>("hNGenMu","hNGenMu",10,0,10);
    hNME0Mu = fs->make<TH1F>("hNME0Mu","hNME0Mu",10,0,10);
    hNMatchME0Mu = fs->make<TH1F>("hNMatchME0Mu","hNMatchME0Mu",10,0,10);
    
    hPVSDphi  = fs->make<TH2F>("hPVSDPhi","hPVSDPhi",1000, 0, 0.01 , 200,0,200);
    
    
    hPtVSDX = fs->make<TH2F>("hPtVSDX","hPtVSDX",1000, 0, 10 , 200,0,200);
    hPtVSDY = fs->make<TH2F>("hPtVSDY","hPtVSDY",1000, 0, 10 , 200,0,200);
    
    
    hPtVSDXLocal = fs->make<TH2F>("hPtVSDXLocal","hPtVSDXLocal",1000, 0, 10 , 200,0,200);
    hPtVSDYLocal =  fs->make<TH2F>("hPtVSDYLocal","hPtVSDYLocal",1000, 0, 10 , 200,0,200);
    
    
    
    hGenMuPt = fs->make<TH1F>("hGenMuPt","hGenMuPt",200,0,200);
    hGenMuEta = fs->make<TH1F>("hGenMuEta","hGenMuEta",200,0,4);
    
    hNumTight_Pt = fs->make<TH1F>("hNumTight_Pt","hNumTight_Pt",200,0,200);
    hNumTight_Eta = fs->make<TH1F>("hNumTight_Eta","hNumTight_Eta",200,0,4);
    hNumLoose_Pt = fs->make<TH1F>("hNumLoose_Pt","hNumLoose_Pt",200,0,200);
    hNumLoose_Eta = fs->make<TH1F>("hNumLoose_Eta","hNumLoose_Eta",200,0,4);
    hNumNoID_Eta = fs->make<TH1F>("hNumNoID_Eta","hNumNoID_Eta",200,0,4);
    hNumNoID_Pt = fs->make<TH1F>("hNumNoID_Pt","hNumNoID_Pt",200,0,200);
    
    hGenMass =  fs->make<TH1F>("hGenMass","hGenMass",200,0,200);
    hRecoMass =  fs->make<TH1F>("hRecoMass","hRecoMass",200,0,200);
    hRecoMassIntime =  fs->make<TH1F>("hRecoMassIntime","hRecoMassIntime",200,0,200);
    hRecoMass_matchID = fs->make<TH1F>("hRecoMass_matchID","hRecoMass_matchID",200,0,200);
    
    hSelectedSimTrack =  fs->make<TH1F>("hSelectedSimTrack","hSelectedSimTrack",20,-0.5,20.5);
    hNME0Segment =  fs->make<TH1F>("hNME0Segment","hNME0Segment ",500,0,500);
    hGenMuonsME0= fs->make<TH1F>("hGenMuonsME0","hGenMuonsME0",20,0,20);
    hME0MuonsInMatchCone = fs->make<TH1F>("hME0MuonsInMatchCone","hME0MuonsInMatchCone",20,0,20);
    hME0MuonsID = fs->make<TH1F>("hME0MuonsID","hME0MuonsID",20,0,20);
    hME0SimTrackPdgID = fs->make<TH1F>("hME0SimTrackPdgID","hME0SimTrackPdgID",3000,0,3000);
    
    hSimElePtME0  = fs->make<TH1F>("hSimElePtME0","hSimElePtME0",100,0,10);
    hSimMuPtME0   = fs->make<TH1F>("hSimMuPtME0","hSimMuPtME0",100,0,100);
    hSimPionPtME0 = fs->make<TH1F>("hSimPionPtME0","hSimPionPtME0",100,0,10);
    
    hSimEleNHitsME0 = fs->make<TH1F>("hSimEleNHitsME0","hSimEleNHitsME0",10,0,10);
    hSimMuonNHitsME0 = fs->make<TH1F>("hSimMuonNHitsME0","hSimMuonNHitsME0",10,0,10);
    hSimPionNHitsME0 = fs->make<TH1F>("hSimPionNHitsME0","hSimPionNHitsME0",10,0,10);
    hPdgIDCheck= fs->make<TH1F>("hPdgIDCheck","hPdgIDCheck",220,0,220);
    hDRME0SimTrack  = fs->make<TH1F>("hDRME0SimTrack","hDRME0SimTrack",200,0,10);
    hDRME0SimMuonEle= fs->make<TH1F>("hDRME0SimMuonEle","hDRME0SimMuonEle",200,0,10);
    
    hNME0SegmentAfterDigiMatch= fs->make<TH1F>("hNME0SegmentAfterDigiMatch","hNME0SegmentAfterDigiMatch",100,0,500);
    hSegmentComposition= fs->make<TH1F>("hSegmentComposition","hSegmentComposition",3,0,3);
    hSegmentComposition->GetXaxis()->SetBinLabel(1,"13 || 11");
    hSegmentComposition->GetXaxis()->SetBinLabel(2,"13 || 11 || >200");
    hSegmentComposition->GetXaxis()->SetBinLabel(3,"13 only");
    
    hMuEleinME0Segm = fs->make<TH1F>("hMuEleinME0Segm","hMuEleinME0Segm",3,0,3);
    hMuEleOthersinME0Segm = fs->make<TH1F>("hMuEleOthersinME0Seg","hMuEleOthersinME0Seg",3,0,3);
    hMuOnlyinME0Segm = fs->make<TH1F>("hMuOnlyinME0Segm","hMuOnlyinME0Segm",3,0,3);
    hNoMuinME0Segm= fs->make<TH1F>("hNoMuinME0Segm","hNoMuinME0Segm",3,0,3);
    
    hNEleBrem= fs->make<TH1F>("hNEleBrem","hNEleBrem",3,0,3);
    hNEleDeltaRays= fs->make<TH1F>("hNEleDeltaRays","hNEleDeltaRays",3,0,3);
    hNEle= fs->make<TH1F>("hNEle","hNEle",3,0,3);
    hMatchedSimEleTrack= fs->make<TH1F>("hMatchedSimEleTrack","hMatchedSimEleTrack",100,0,5);
    
    hRHDeltaPhiSameLayer= fs->make<TH1F>("hRHDeltaPhiSameLayer","hRHDeltaPhiSameLayer",1000,-0.1,0.1);
    hRHDeltaEtaSameLayer= fs->make<TH1F>("hRHDeltaEtaSameLayer","hRHDeltaEtaSameLayer",1000,-0.5,0.5);
    hRHDeltaTSameLayer= fs->make<TH1F>("hRHDeltaTSameLayer","hRHDeltaTSameLayer",1000,-5,5);
    
    hNMuonSameLayerTOF = fs->make<TH1F>("hNMuonSameLayerTOF","hNMuonSameLayerTOF",1000,0,100);
    hNEleSameLayerTOF = fs->make<TH1F>("hNEleSameLayerTOF","hNEleSameLayerTOF",1000,0,100);
    hNoPromptSameLayerTOF = fs->make<TH1F>("hNoPromptLayerTOF","hNoPromptLayerTOF",1000,0,100);
    
    hNoPromptRecHitTime= fs->make<TH1F>("hNoPromptRecHitTime","hNoPromptRecHitTime",1000,0,100);
    hMuonRecHitTime= fs->make<TH1F>("hMuonRecHitTime","hMuonRecHitTime",1000,0,100);
    hEleRecHitTime= fs->make<TH1F>("hEleRecHitTime","hEleRecHitTime",1000,0,100);
    
    
    hDeltaPhiSimReco =  fs->make<TH1F>("hDeltaPhiSimReco","hDeltaPhiSimReco",1000,-0.1,0.1);
    
    hDeltaPhiSimReco_1 =  fs->make<TH1F>("hDeltaPhiSimReco_1","hDeltaPhiSimReco_1",1000,-0.1,0.1);
    hDeltaPhiSimReco_2 =  fs->make<TH1F>("hDeltaPhiSimReco_2","hDeltaPhiSimReco_2",1000,-0.1,0.1);
    hDeltaPhiSimReco_3 =  fs->make<TH1F>("hDeltaPhiSimReco_3","hDeltaPhiSimReco_3",1000,-0.1,0.1);
    
    hDeltaPhiSimReco_pos =  fs->make<TH1F>("hDeltaPhiSimReco_pos","hDeltaPhiSimReco_pos",1000,-0.1,0.1);
    hDeltaPhiSimReco_neg =  fs->make<TH1F>("hDeltaPhiSimReco_neg","hDeltaPhiSimReco_neg",1000,-0.1,0.1);
    
    hDeltaEtaSimReco =  fs->make<TH1F>("hDeltaEtaSimReco","hDeltaEtaSimReco",1000,-1,1);
    
    hDeltaXSimReco =  fs->make<TH1F>("hDeltaXSimReco","hDeltaXSimReco",1000,-1,1);
    hDeltaYSimReco =  fs->make<TH1F>("hDeltaYSimReco","hDeltaYSimReco",1000,-50,50);
    
    hDeltaXSimRecoLocal = fs->make<TH1F>("hDeltaXSimRecoLocal","hDeltaXSimRecoLocal",1000,-1,1);
    hDeltaYSimRecoLocal =  fs->make<TH1F>("hDeltaYSimRecoLocal","hDeltaYSimRecoLocal",1000,-50,50);
    
    hDeltaXSimRecoLocal_1 = fs->make<TH1F>("hDeltaXSimRecoLocal_1","hDeltaXSimRecoLocal_1",1000,-1,1);
    hDeltaXSimRecoLocal_2 = fs->make<TH1F>("hDeltaXSimRecoLocal_2","hDeltaXSimRecoLocal_2",1000,-1,1);
    hDeltaXSimRecoLocal_3 = fs->make<TH1F>("hDeltaXSimRecoLocal_3","hDeltaXSimRecoLocal_3",1000,-1,1);
    
    hRecDPhiVSimDphi=fs->make<TH2F>("hRecDPhiVSimDphi","hRecDPhiVSimDphi",1000, -0.01, 0.01 , 1000, -0.01, 0.01);
    hDiffRecDPhiVSimDphi =  fs->make<TH1F>("hDiffRecDPhiVSimDphi","hDiffRecDPhiVSimDphi",500,-0.05,0.05);
    hRecDPhiOverSimDphi =  fs->make<TH1F>("hRecDPhiOverSimDphi","hRecDPhiOverSimDphi",1000,-4,4);
    hDqOverDphi=  fs->make<TH1F>("hDqOverDphi","hDqOverDphi",1000,-0.1,0.1);
    hQSimQRecoOverDphi =  fs->make<TH1F>("hQSimQRecoOverDphi","hQSimQRecoOverDphi",1000,-0.1,0.1);
    
    hBeamSpotX0 =  fs->make<TH1F>("hBeamSpotX0","hBeamSpotX0",1000,-10,10);
    hBeamSpotY0 =  fs->make<TH1F>("hBeamSpotY0","hBeamSpotY0",1000,-10,10);
    hBeamSpotX0Y0 =  fs->make<TH2F>("hBeamSpotX0Y0","hBeamSpotX0Y0",1000,-10,10,1000,-10,10);
    hBeamSpotZ0 = fs->make<TH1F>("hBeamSpotZ0","hBeamSpotZ0",1000,-30,30);
    hBeamSpotSigmaZ = fs->make<TH1F>("hBeamSpotSigmaZ","hBeamSpotSigmaZ",1000,-30,30);
    hBeamSpotdxdz  = fs->make<TH1F>("hBeamSpotdxdz","hBeamSpotdxdz",1000,-30,30);
    hBeamSpotBeamWidthX = fs->make<TH1F>("hBeamSpotBeamWidthX","hBeamSpotBeamWidthX",1000,-30,30);
    hBeamSpotBeamWidthY = fs->make<TH1F>("hBeamSpotBeamWidthY","hBeamSpotBeamWidthY",1000,-30,30);
    hVertexMult  = fs->make<TH1F>("hVertexMult","hVertexMult",100,0,300);
    hverteX  = fs->make<TH1F>("hverteX","hverteX",1000,-5,5);
    hverteY  = fs->make<TH1F>("hverteY","hverteY",1000,-5,5);
    hverteXY  = fs->make<TH2F>("hverteXY","hverteXY",1000,-0.5,0.5,1000, -0.5,0.5);
    hverteZT  = fs->make<TH2F>("hverteZT","hverteZT",1000,-30,30,1000, -1,1 );
    hverteZ = fs->make<TH1F>("hverteZ","hverteZ",1000,-30,30);
    hverteTime = fs->make<TH1F>("hverteTime","hverteTime",500,-1,1);
    hBX = fs->make<TH1F>("hBX","hBX",30,-15,15);
    hNPU = fs->make<TH1F>("hNPU","hNPU",400,0,400);
    hTrueInt= fs->make<TH1F>("hTrueInt","hTrueInt",400,0,400);
    
    hSimVertZPosition = fs->make<TH1F>("hSimVertZPosition","hSimVertZPosition",100,-30,30);
    hSimVertXYPosition = fs->make<TH2F>("hSimVertXYPosition","hSimVertXYPosition",100,-5,5,100, -5,5 );
    
    hME0MuonsChargeNeg = fs->make<TH1F>("hME0MuonsChargeNeg","hME0MuonsChargeNeg",5,-2.5,2.5);
    hME0MuonsChargePos = fs->make<TH1F>("hME0MuonsChargePos","hME0MuonsChargePos",5,-2.5,2.5);
    
    hGenMuNegEta= fs->make<TH1F>("hGenMuNegEta","hGenMuNegEta",100,-3,3);
    hGenMuNegPhi= fs->make<TH1F>("hGenMuNegPhi","hGenMuNegPhi",100,-4,4);
    hGenMuPosEta= fs->make<TH1F>("hGenMuPosEta","hGenMuPosEta",100,-3,3);
    hGenMuPosPhi= fs->make<TH1F>("hGenMuPosPhi","hGenMuPosPhi",100,-4,4);
    
    hSimLocalDX= fs->make<TH1F>("hSimLocalDX","hSimLocalDX",1000,-50,50);
    hLocalDX= fs->make<TH1F>("hLocalDX","hLocalDX",1000,-50,50);
    
    hSimLocalDXPos= fs->make<TH1F>("hSimLocalDXPos","hSimLocalDXPos",1000,-20,20);
    hLocalDXPos= fs->make<TH1F>("hLocalDXPos","hLocalDXPos",1000,-20,20);
    
    hSimLocalDXNeg= fs->make<TH1F>("hSimLocalDXNeg","hSimLocalDXNeg",1000,-20,20);
    hLocalDXNeg= fs->make<TH1F>("hLocalDXNeg","hLocalDXNeg",1000,-20,20);
    
    hSimLocalDXPos_LowPt= fs->make<TH1F>("hSimLocalDXPos_LowPt","hSimLocalDXPos_LowPt",1000,-20,20);
    hLocalDXPos_LowPt= fs->make<TH1F>("hLocalDXPos_LowPt","hLocalDXPos_LowPt",1000,-20,20);
    
    hSimLocalDXNeg_LowPt= fs->make<TH1F>("hSimLocalDXNeg_LowPt","hSimLocalDXNeg_LowPt",1000,-20,20);
    hLocalDXNeg_LowPt= fs->make<TH1F>("hLocalDXNeg_LowPt","hLocalDXNeg_LowPt",1000,-20,20);
    
    hSimLocalDXPos_HighPt= fs->make<TH1F>("hSimLocalDXPos_HighPt","hSimLocalDXPos_HighPt",1000,-20,20);
    hLocalDXPos_HighPt= fs->make<TH1F>("hLocalDXPos_HighPt","hLocalDXPos_HighPt",1000,-20,20);
    
    hSimLocalDXNeg_HighPt= fs->make<TH1F>("hSimLocalDXNeg_HighPt","hSimLocalDXNeg_HighPt",1000,-20,20);
    hLocalDXNeg_HighPt= fs->make<TH1F>("hLocalDXNeg_HighPt","hLocalDXNeg_HighPt",1000,-20,20);
    
    //////////////////////
    //	hSimLocalDX= fs->make<TH1F>("hSimLocalDX","hSimLocalDX",1000,-50,50);
    //	hLocalDX= fs->make<TH1F>("hLocalDX","hLocalDX",1000,-50,50);
    
    hSimLocalDXPos_ME0Plus= fs->make<TH1F>("hSimLocalDXPos_ME0Plus","hSimLocalDXPos_ME0Plus",1000,-20,20);
    hLocalDXPos_ME0Plus= fs->make<TH1F>("hLocalDXPos_ME0Plus","hLocalDXPos_ME0Plus",1000,-20,20);
    
    hSimLocalDXNeg_ME0Plus= fs->make<TH1F>("hSimLocalDXNeg_ME0Plus","hSimLocalDXNeg_ME0Plus",1000,-20,20);
    hLocalDXNeg_ME0Plus= fs->make<TH1F>("hLocalDXNeg_ME0Plus","hLocalDXNeg_ME0Plus",1000,-20,20);
    
    hSimLocalDXPos_LowPt_ME0Plus= fs->make<TH1F>("hSimLocalDXPos_LowPt_ME0Plus","hSimLocalDXPos_LowPt_ME0Plus",1000,-20,20);
    hLocalDXPos_LowPt_ME0Plus= fs->make<TH1F>("hLocalDXPos_LowPt_ME0Plus","hLocalDXPos_LowPt_ME0Plus",1000,-20,20);
    
    hSimLocalDXNeg_LowPt_ME0Plus= fs->make<TH1F>("hSimLocalDXNeg_LowPt_ME0Plus","hSimLocalDXNeg_LowPt_ME0Plus",1000,-20,20);
    hLocalDXNeg_LowPt_ME0Plus= fs->make<TH1F>("hLocalDXNeg_LowPt_ME0Plus","hLocalDXNeg_LowPt_ME0Plus",1000,-20,20);
    
    hSimLocalDXPos_HighPt_ME0Plus= fs->make<TH1F>("hSimLocalDXPos_HighPt_ME0Plus","hSimLocalDXPos_HighPt_ME0Plus",1000,-20,20);
    hLocalDXPos_HighPt_ME0Plus= fs->make<TH1F>("hLocalDXPos_HighPt_ME0Plus","hLocalDXPos_HighPt_ME0Plus",1000,-20,20);
    
    hSimLocalDXNeg_HighPt_ME0Plus= fs->make<TH1F>("hSimLocalDXNeg_HighPt_ME0Plus","hSimLocalDXNeg_HighPt_ME0Plus",1000,-20,20);
    hLocalDXNeg_HighPt_ME0Plus= fs->make<TH1F>("hLocalDXNeg_HighPt_ME0Plus","hLocalDXNeg_HighPt_ME0Plus",1000,-20,20);
    
    hSimDPhiPos_ME0Plus = fs->make<TH1F>("hSimDPhiPos_ME0Plus","hSimDPhiPos_ME0Plus",1000,-0.05,0.05);
    hSimDPhiNeg_ME0Plus = fs->make<TH1F>("hSimDPhiNeg_ME0Plus","hSimDPhiNeg_ME0Plus",1000,-0.05,0.05);
    
    hSimDPhiPos_LowPt_ME0Plus = fs->make<TH1F>("hSimDPhiPos_LowPt_ME0Plus","hSimDPhiPos_LowPt_ME0Plus",1000,-0.05,0.05);
    hSimDPhiPos_HighPt_ME0Plus = fs->make<TH1F>("hSimDPhiPos_HighPt_ME0Plus","hSimDPhiPos_HighPt_ME0Plus",1000,-0.05,0.05);
    
    hSimDPhiNeg_LowPt_ME0Plus = fs->make<TH1F>("hSimDPhiNeg_LowP_ME0Plust","hSimDPhiNeg_LowPt_ME0Plus",1000,-0.05,0.05);
    hSimDPhiNeg_HighPt_ME0Plus = fs->make<TH1F>("hSimDPhiNeg_HighPt_ME0Plus","hSimDPhiNeg_HighPt_ME0Plus",1000,-0.05,0.05);
    
    hDPhiPos_ME0Plus = fs->make<TH1F>("hDPhiPos_ME0Plus","hDPhiPos_ME0Plus",1000,-0.05,0.05);
    hDPhiNeg_ME0Plus = fs->make<TH1F>("hDPhiNeg_ME0Plus","hDPhiNeg_ME0Plus",1000,-0.05,0.05);
    
    hDPhiPos_LowPt_ME0Plus = fs->make<TH1F>("hDPhiPos_LowPt_ME0Plus","hDPhiPos_LowPt_ME0Plus",1000,-0.05,0.05);
    hDPhiPos_HighPt_ME0Plus = fs->make<TH1F>("hDPhiPos_HighPt_ME0Plus","hDPhiPos_HighPt_ME0Plus",1000,-0.05,0.05);
    
    hDPhiNeg_LowPt_ME0Plus = fs->make<TH1F>("hDPhiNeg_LowPt_ME0Plus","hDPhiPos_NegPt_ME0Plus",1000,-0.05,0.05);
    hDPhiNeg_HighPt_ME0Plus = fs->make<TH1F>("hDPhiNeg_HighPt_ME0Plus","hDPhiNeg_HighPt_ME0Plus",1000,-0.05,0.05);
    
    
    
    /////////////////////
    
    hSimLocalDXPos_ME0Minus= fs->make<TH1F>("hSimLocalDXPos_ME0Minus","hSimLocalDXPos_ME0Minus",1000,-20,20);
    hLocalDXPos_ME0Minus= fs->make<TH1F>("hLocalDXPos_ME0Minus","hLocalDXPos_ME0Minus",1000,-20,20);
    
    hSimLocalDXNeg_ME0Minus= fs->make<TH1F>("hSimLocalDXNeg_ME0Minus","hSimLocalDXNeg_ME0Minus",1000,-20,20);
    hLocalDXNeg_ME0Minus= fs->make<TH1F>("hLocalDXNeg_ME0Minus","hLocalDXNeg_ME0Minus",1000,-20,20);
    
    hSimLocalDXPos_LowPt_ME0Minus= fs->make<TH1F>("hSimLocalDXPos_LowPt_ME0Minus","hSimLocalDXPos_LowPt_ME0Minus",1000,-20,20);
    hLocalDXPos_LowPt_ME0Minus= fs->make<TH1F>("hLocalDXPos_LowPt_ME0Minus","hLocalDXPos_LowPt_ME0Minus",1000,-20,20);
    
    hSimLocalDXNeg_LowPt_ME0Minus= fs->make<TH1F>("hSimLocalDXNeg_LowPt_ME0Minus","hSimLocalDXNeg_LowPt_ME0Minus",1000,-20,20);
    hLocalDXNeg_LowPt_ME0Minus= fs->make<TH1F>("hLocalDXNeg_LowPt_ME0Minus","hLocalDXNeg_LowPt_ME0Minus",1000,-20,20);
    
    hSimLocalDXPos_HighPt_ME0Minus= fs->make<TH1F>("hSimLocalDXPos_HighPt_ME0Minus","hSimLocalDXPos_HighPt_ME0Minus",1000,-20,20);
    hLocalDXPos_HighPt_ME0Minus= fs->make<TH1F>("hLocalDXPos_HighPt_ME0Minus","hLocalDXPos_HighPt_ME0Minus",1000,-20,20);
    
    hSimLocalDXNeg_HighPt_ME0Minus= fs->make<TH1F>("hSimLocalDXNeg_HighPt_ME0Minus","hSimLocalDXNeg_HighPt_ME0Minus",1000,-20,20);
    hLocalDXNeg_HighPt_ME0Minus= fs->make<TH1F>("hLocalDXNeg_HighPt_ME0Minus","hLocalDXNeg_HighPt_ME0Minus",1000,-20,20);
    
    hSimDPhiPos_ME0Minus = fs->make<TH1F>("hSimDPhiPos_ME0Minus","hSimDPhiPos_ME0Minus",1000,-0.05,0.05);  
    hSimDPhiNeg_ME0Minus = fs->make<TH1F>("hSimDPhiNeg_ME0Minus","hSimDPhiNeg_ME0Minus",1000,-0.05,0.05);  
    
    hSimDPhiPos_LowPt_ME0Minus = fs->make<TH1F>("hSimDPhiPos_LowPt_ME0Minus","hSimDPhiPos_LowPt_ME0Minus",1000,-0.05,0.05);  
    hSimDPhiPos_HighPt_ME0Minus = fs->make<TH1F>("hSimDPhiPos_HighPt_ME0Minus","hSimDPhiPos_HighPt_ME0Minus",1000,-0.05,0.05);  
    
    hSimDPhiNeg_LowPt_ME0Minus = fs->make<TH1F>("hSimDPhiNeg_LowP_ME0Minust","hSimDPhiNeg_LowPt_ME0Minus",1000,-0.05,0.05);  
    hSimDPhiNeg_HighPt_ME0Minus = fs->make<TH1F>("hSimDPhiNeg_HighPt_ME0Minus","hSimDPhiNeg_HighPt_ME0Minus",1000,-0.05,0.05);  
    
    hDPhiPos_ME0Minus = fs->make<TH1F>("hDPhiPos_ME0Minus","hDPhiPos_ME0Minus",1000,-0.05,0.05); 
    hDPhiNeg_ME0Minus = fs->make<TH1F>("hDPhiNeg_ME0Minus","hDPhiNeg_ME0Minus",1000,-0.05,0.05); 
    
    hDPhiPos_LowPt_ME0Minus = fs->make<TH1F>("hDPhiPos_LowPt_ME0Minus","hDPhiPos_LowPt_ME0Minus",1000,-0.05,0.05); 
    hDPhiPos_HighPt_ME0Minus = fs->make<TH1F>("hDPhiPos_HighPt_ME0Minus","hDPhiPos_HighPt_ME0Minus",1000,-0.05,0.05); 
    
    hDPhiNeg_LowPt_ME0Minus = fs->make<TH1F>("hDPhiNeg_LowPt_ME0Minus","hDPhiPos_NegPt_ME0Minus",1000,-0.05,0.05); 
    hDPhiNeg_HighPt_ME0Minus = fs->make<TH1F>("hDPhiNeg_HighPt_ME0Minus","hDPhiNeg_HighPt_ME0Minus",1000,-0.05,0.05); 
    
    hLocalDX_over_R  = fs->make<TH1F>("hLocalDX_over_R","hLocalDX_over_R",1000,-0.05,0.05); 
    hSimLocalDXoverR  = fs->make<TH1F>("hSimLocalDXoverR","hSimLocalDXoverR",1000,-0.05,0.05); 
    
    ///////////////////
    
    hSimLocalDY= fs->make<TH1F>("hSimLocalDY","hSimLocalDY",1000,-50,50);
    hLocalDY= fs->make<TH1F>("hLocalDY","hLocalDY",1000,-50,50);
    
    hTightME0SegmDigiHitTOF_noprompt = fs->make<TH1F>("hTightME0SegmDigiHitTOF_noprompt","hTightME0SegmDigiHitTOF_noprompt",2000,-200,200) ;
    hTightME0SegmDigiHitPdgID_noprompt = fs->make<TH1F>("hTightME0SegmDigiHitPdgID_noprompt","hTightME0SegmDigiHitPdgID_noprompt",2214,-0.5,2213.5)  ;
    
    hTightME0SegmDigiHitTOF_prompt = fs->make<TH1F>("hTightME0SegmDigiHitTOF_prompt","hTightME0SegmDigiHitTOF_prompt",2000,-200,200) ;
    hTightME0SegmDigiHitPdgID_prompt = fs->make<TH1F>("hTightME0SegmDigiHitPdgID_prompt","hTightME0SegmDigiHitPdgID_prompt",2214,-0.5,2213.5)  ;
    
    hAllME0SegmentTime = fs->make<TH1F>("hAllME0SegmentTime","hAllME0SegmentTime",200,-100,100)  ;
    hAllME0SegmentTimeErr= fs->make<TH1F>("hAllME0SegmentTimeErr","hAllME0SegmentTimeErr",200,-100,100)  ;
    hME0SegmentCollectionTime= fs->make<TH1F>("hME0SegmentCollectionTime","hME0SegmentCollectionTime",200,-100,100)  ;
    
    hDigiHitType = fs->make<TH1F>("hDigiHitType","hDigiHitType",5,-0.5,4)  ;
    hDigiHitPdgID = fs->make<TH1F>("hDigiHitPdgID","hDigiHitPdgID",2214,-0.5,2213.5)  ;
    hDigiHitToF_prompt= fs->make<TH1F>("hDigiHitToF_prompt","hDigiHitToF_prompt",200,-100,100)  ;
    hDigiHitToF_noprompt= fs->make<TH1F>("hDigiHitToF_noprompt","hDigiHitToF_noprompt",200,-100,100)  ;
    
    hMuonInME0SegTOF = fs->make<TH1F>("hMuonInME0SegTOF","hMuonInME0SegTOF",20000,-100,100)  ;	
    hEleInME0SegTOF= fs->make<TH1F>("hEleInME0SegTOF","hEleInME0SegTOF",20000,-100,100)  ;	
    hPionInME0SegTOF=	fs->make<TH1F>("hPionInME0SegTOF","hPionInME0SegTOF",20000,-100,100)  ;	
    hBoInME0SegTOF=	fs->make<TH1F>("hBoInME0SegTOF","hBoInME0SegTOF",20000,-100,100) ;	
    hProtonsInME0SegTOF=	fs->make<TH1F>("hProtonsInME0SegTO","hProtonsInME0SegTOF",20000,-100,100) ;
    
    hMuonInME0SegTOF_L1 = fs->make<TH1F>("hMuonInME0SegTOF_L1","hMuonInME0SegTOF_L1",20000,-100,100)  ;	
    hEleInME0SegTOF_L1= fs->make<TH1F>("hEleInME0SegTOF_L1","hEleInME0SegTOF_L1",20000,-100,100)  ;	
    hPionInME0SegTOF_L1=	fs->make<TH1F>("hPionInME0SegTOF_L1","hPionInME0SegTOF_L1",20000,-100,100)  ;	
    hBoInME0SegTOF_L1=	fs->make<TH1F>("hBoInME0SegTOF_L1","hBoInME0SegTOF",20000,-100,100) ;	
    hProtonsInME0SegTOF_L1=	fs->make<TH1F>("hProtonsInME0SegTOF_L1","hProtonsInME0SegTOF_L1",20000,-100,100) ;
    
    hMuonInME0SegTOF_inTime		= fs->make<TH1F>("hMuonInME0SegTOF_inTime","hMuonInME0SegTOF_inTime",20000,-100,100)  ;	
    hEleInME0SegTOF_inTime		= fs->make<TH1F>("hEleInME0SegTOF_inTime","hEleInME0SegTOF_inTime",20000,-100,100)  ;	
    hPionInME0SegTOF_inTime     = fs->make<TH1F>("hPionInME0SegTOF_inTime","hPionInME0SegTOF_inTime",20000,-100,100)  ;	
    hBoInME0SegTOF_inTime		= fs->make<TH1F>("hBoInME0SegTOF_inTime","hBoInME0SegTOF_inTime",20000,-100,100) ;	
    hProtonsInME0SegTOF_inTime	= fs->make<TH1F>("hProtonsInME0SegTOF_inTime","hProtonsInME0SegTOF_inTime",20000,-100,100) ;
    
    
    hME0SegmDigiHitME0Segm_prompt = fs->make<TH1F>("hME0SegmDigiHitTOFME0Segm_prompt","hME0SegmDigiHitTOFME0Segm_prompt",2000,-200,200) ;
    hME0SegmDigiHitPdgIDME0Segm_prompt = fs->make<TH1F>("hME0SegmDigiHitPdgIDME0Segm_prompt","hME0SegmDigiHitPdgIDME0Segm_prompt",2214,-0.5,2213.5)  ;
    hME0SegmDigiHitTOFME0Segm_noprompt = fs->make<TH1F>("hME0SegmDigiHitTOFME0Segm_noprompt", "hME0SegmDigiHitTOFME0Segm_noprompt",2000,-200,200) ;
    hME0SegmDigiHitPdgIDME0Segm_noprompt=  fs->make<TH1F>("hME0SegmDigiHitPdgIDME0Segm_noprompt", "hME0SegmDigiHitPdgIDME0Segm_noprompt", 2214,-0.5,2213.5)  ;
    

    
    hDPhiNoMatchME0Seg = fs->make<TH1F>("hDPhiNoMatchME0Seg","hDPhiNoMatchME0Seg",1000,-0.05,0.05);
    hDPhiNoMatchME0Seg_TimeWindow = fs->make<TH1F>("hDPhiNoMatchME0Seg_TimeWindow","hDPhiNoMatchME0Seg_TimeWindow",1000,-0.05,0.05);
    
    hDPhivsTOFMatchByHitsME0Seg_Bkg  =  fs->make<TH2F>("hDPhivsTOFMatchByHitsME0Seg_Bkg","hDPhivsTOFMatchByHitsME0Seg_Bkg", 1000, -0.05, 0.05, 2000,0,100)  ; 
    hDPhivsTOFMatchByHitsME0Seg_Signal =  fs->make<TH2F>("hDPhivsTOFMatchByHitsME0Seg_Signal","hDPhivsTOFMatchByHitsME0Seg_Signal", 1000, -0.05, 0.05, 2000,0,100)  ; 
    
    hDPhivsTOFMatchByHitsME0Seg_Bkg_timeCut  =  fs->make<TH2F>("hDPhivsTOFMatchByHitsME0Seg_Bkg_timeCut","hDPhivsTOFMatchByHitsME0Seg_Bkg_timeCut", 1000, -0.05, 0.05, 2000,0,100)  ; 
    hDPhivsTOFMatchByHitsME0Seg_Signal_timeCut =  fs->make<TH2F>("hDPhivsTOFMatchByHitsME0Seg_Signal_timeCut","hDPhivsTOFMatchByHitsME0Seg_Signal_timeCut", 1000, -0.05, 0.05, 2000,0,100)  ; 
    
    hNME0SegmentSignal =  fs->make<TH1F>("hNME0SegmentSignal","hNME0SegmentSignal", 50, -0.5, 50.5)  ;
    hNME0SegmentBkg   =  fs->make<TH1F>("hNME0SegmentBkg","hNME0SegmentBkg", 500, 0, 500)  ;
    
    hNME0SegmentSignal_timeCut  =  fs->make<TH1F>("hNME0SegmentSignal_timeCut ","hNME0SegmentSignal_timeCut ", 50, 0, 50)  ; 
    hNME0SegmentBkg_timeCut    =  fs->make<TH1F>("hNME0SegmentBkg_timeCut ","hNME0SegmentBkg_timeCut ", 100, 0, 500)  ;
    
    hTimeDiffSigvsBkg  =  fs->make<TH1F>("hTimeDiffSigvsBkg","hTimeDiffSigvsBkg", 2000, -10, 10)  ; 
    
    hSimMuonPt_DPhiCut5GeV = fs->make<TH1F>("hSimMuonPt_DPhiCut5GeV","hSimMuonPt_DPhiCut5GeV", 200, 0, 200)  ; 
    hSimMuonPt_DPhiCut30GeV = fs->make<TH1F>("hSimMuonPt_DPhiCut30GeV","hSimMuonPt_DPhiCut30GeV", 200, 0, 200)  ; 
    hEventPass30GeVCut = fs->make<TH1F>("hEventPass30GeVCut","hEventPass30GeVCut", 10, 0, 10)  ; 
    hEventPass5GeVCut = fs->make<TH1F>("hEventPass5GeVCut","hEventPass5GeVCut", 10, 0, 10)  ;
    
    hSimMuonPt_DPhiCut30GeVAllCut = fs->make<TH1F>("hSimMuonPt_DPhiCut30GeVAllCut","hSimMuonPt_DPhiCut30GeVAllCut", 200, 0, 200)  ;
    hSimMuonPt_DPhiCut5GeVAllCut  = fs->make<TH1F>("hSimMuonPt_DPhiCut5GeVAllCut","hSimMuonPt_DPhiCut5GeVAllCut", 200, 0, 200)  ;
    hEventPass5GeVCutAllCut = fs->make<TH1F>("hEventPass5GeVCutAllCut","hEventPass5GeVCutAllCut", 10, 0, 10)  ;
    hEventPass30GeVCutAllCut = fs->make<TH1F>("hEventPas30GeVCutAllCut","hEventPass30GeVCutAllCut", 10, 0, 10)  ;
    
    
    
    hSimME0SegmentTime = fs->make<TH1F>("hSimME0SegmentTime","hSimME0SegmentTime", 4000, 0, 100)  ;
    hSimVertZPositionVsSimSegmTime =  fs->make<TH2F>("hSimVertZPositionVsSimSegmTime","hSimVertZPositionVsSimSegmTime", 1000, -30, 30, 1000,-30,30)  ;
    
    hME0SegmNumber  = fs->make<TH1F>("hME0SegmNumber","hME0SegmNumber", 3, 0, 3)  ;
    hME0SegCompMuonOnly = fs->make<TH1F>("hME0SegCompMuonOnly","hME0SegCompMuonOnly", 3, 0, 3)  ;
    hME0SegCompMuonPlusEle = fs->make<TH1F>("hME0SegCompMuonPlusEle","hME0SegCompMuonPlusEle", 3, 0, 3) ;
    hME0SegCompMuonPlusHad = fs->make<TH1F>("hME0SegCompMuonPlusHad","hME0SegCompMuonPlusHad", 3, 0, 3) ;
    hME0SegCompMuonPlusElePlusHad = fs->make<TH1F>("hME0SegCompMuonPlusElePlusHad","hME0SegCompMuonPlusElePlusHad", 3, 0, 3) ;
    hME0SegCompHadOnly = fs->make<TH1F>("hME0SegCompHadOnly","hME0SegCompHadOnly", 3, 0, 3) ;
    hME0SegCompEleOnly = fs->make<TH1F>("hME0SegCompEleOnly","hME0SegCompEleOnly", 3, 0, 3) ;
    
    hME0SegmNumberSig  = fs->make<TH1F>("hME0SegmNumberSig","hME0SegmNumberSig", 3, 0, 3)  ;
    hME0SegCompMuonOnly_Sig = fs->make<TH1F>("hME0SegCompMuonOnly_Sig","hME0SegCompMuonOnly_Sig", 3, 0, 3)  ;
    hME0SegCompMuonPlusEle_Sig = fs->make<TH1F>("hME0SegCompMuonPlusEle_Sig","hME0SegCompMuonPlusEle_Sig", 3, 0, 3)  ;
    hME0SegCompMuonPlusHad_Sig = fs->make<TH1F>("hME0SegCompMuonPlusHad_Sig","hME0SegCompMuonPlusHad_Sig", 3, 0, 3)  ;
    hME0SegCompMuonPlusElePlusHad_Sig = fs->make<TH1F>("hME0SegCompMuonPlusElePlusHad_Sig","hME0SegCompMuonPlusElePlusHad_Sig", 3, 0, 3)  ;
    hME0SegCompEleOnly_Sig = fs->make<TH1F>("hME0SegCompEleOnly_Sig","hME0SegCompEleOnly_Sig", 3, 0, 3)  ;
    hME0SegCompHadOnly_Sig = fs->make<TH1F>("hME0SegCompHadOnly_Sig","hME0SegCompHadOnly_Sig", 3, 0, 3)  ;
    
    hME0SegmNumberBkg  = fs->make<TH1F>("hME0SegmNumberBkg","hME0SegmNumberBkg", 3, 0, 3)  ;
    hME0SegCompMuonOnly_Bkg = fs->make<TH1F>("hME0SegCompMuonOnly_Bkg","hME0SegCompMuonOnly_Bkg", 3, 0, 3)  ;
    hME0SegCompMuonPlusEle_Bkg = fs->make<TH1F>("hME0SegCompMuonPlusEle_Bkg","hME0SegCompMuonPlusEle_Bkg", 3, 0, 3)  ;
    hME0SegCompMuonPlusHad_Bkg = fs->make<TH1F>("hME0SegCompMuonPlusHad_Bkg","hME0SegCompMuonPlusHad_Bkg", 3, 0, 3)  ;
    hME0SegCompMuonPlusElePlusHad_Bkg = fs->make<TH1F>("hME0SegCompMuonPlusElePlusHad_Bkg","hME0SegCompMuonPlusElePlusHad_Bkg", 3, 0, 3)  ;
    hME0SegCompEleOnly_Bkg = fs->make<TH1F>("hME0SegCompEleOnly_Bkg","hME0SegCompEleOnly_Bkg", 3, 0, 3)  ;
    hME0SegCompHadOnly_Bkg = fs->make<TH1F>("hME0SegCompHadOnly_Bkg","hME0SegCompHadOnly_Bkg", 3, 0, 3)  ;
    
    hPull_MyRecoSim_zVtx = fs->make<TH1F>("hPull_MyRecoSim_zVtx","hPull_MyRecoSim_zVtx", 2000, -20, 20)  ;
    hPull_4DRecoSim_zVtx = fs->make<TH1F>("hPull_4DRecoSim_zVtx"," hPull_4DRecoSim_zVtx", 2000, -20, 20)  ;
    
    
    hME0Segm_nHit_bkg = fs->make<TH1F>("hME0Segm_nHit_bkg","hME0Segm_nHit_bkg", 20, -0.5, 19.5)  ;
    hME0Segm_chi2_bkg = fs->make<TH1F>("hME0Segm_chi2_bkg","hME0Segm_chi2_bkg", 100, 0, 100)  ;
    hME0Segm_time_bkg = fs->make<TH1F>("hME0Segm_time_bkg","hME0Segm_time_bkg", 2000, -200, 200)  ;
    hME0Segm_eta_bkg  =  fs->make<TH1F>("hME0Segm_eta_bkg","hME0Segm_eta_bkg", 600, -3, 3)  ;    

    hME0Segm_nHit_sig = fs->make<TH1F>("hME0Segm_nHit_sig","hME0Segm_nHit_sig", 20, -0.5, 19.5)  ;
    hME0Segm_chi2_sig = fs->make<TH1F>("hME0Segm_chi2_sig","hME0Segm_chi2_sig", 100, 0, 100)  ;
    hME0Segm_time_sig = fs->make<TH1F>("hME0Segm_time_sig","hME0Segm_time_sig", 2000, -200, 200)  ;
    hME0Segm_eta_sig  =  fs->make<TH1F>("hME0Segm_eta_sig","hME0Segm_eta_sig", 600, -3, 3)  ;

    hME0Segm_nHit = fs->make<TH1F>("hME0Segm_nHit","hME0Segm_nHit", 20, -0.5, 19.5)  ;
    hME0Segm_chi2 = fs->make<TH1F>("hME0Segm_chi2","hME0Segm_chi2", 100, 0, 100)  ;

    hME0SegmMyTOF = fs->make<TH1F>("hME0SegmMyTOF","hME0SegmMyTOF", 400, -200, 200)  ;
    hME0Segm_DeltaPhiVsSimpT_Signal = fs->make<TH2F>("hME0Segm_DeltaPhiVsSimpT_Signal"," hME0Segm_DeltaPhiVsSimpT_Signal",5000, 0, 0.05 , 200,0,200);
    //float binDiv[] = { 0.5, 1., 3., 5., 10., 20., 50., 200. };
    //int   binNum   = sizeof(binDiv)/sizeof(float) - 1; // or just = 9
    //    hME0Segm_DeltaPhiVsSimpT_Signal_rebin = fs->make<TH2F>("hME0Segm_DeltaPhiVsSimpT_Signal_rebin"," hME0Segm_DeltaPhiVsSimpT_Signal_rebin", binNum , binDiv, 200,0,200);
    
    hME0SegmDigiHitTOF_signal = fs->make<TH1F>("hME0SegmDigiHitTOF_signal","hME0SegmDigiHitTOF_signal", 4000, -200, 200)  ;
    hME0SegmDigiHitPdgID_signal = fs->make<TH1F>("hME0SegmDigiHitPdgID_signal","hME0SegmDigiHitPdgID_signal", 2221, 0, 2221)  ;
    hME0SegmDigiHitEta_signal= fs->make<TH1F>("hME0SegmDigiHitEta_signal","hME0SegmDigiHitEta_signal", 600, -3, 3)  ;
    
    hME0SegmDigiHitTOF_bkg = fs->make<TH1F>("hME0SegmDigiHitTOF_bkg","hME0SegmDigiHitTOF_bkg", 4000, -200, 200)  ;
    hME0SegmDigiHitPdgID_bkg = fs->make<TH1F>("hME0SegmDigiHitPdgID_bkg","hME0SegmDigiHitPdgID_bkg", 2221, 0, 2221)  ;
    hME0SegmDigiHitEta_bkg= fs->make<TH1F>("hME0SegmDigiHitEta_bkg","hME0SegmDigiHitEta_bkg", 600, -3, 3)  ;
    
    hME0SegmDigiHitTOF = fs->make<TH1F>("hME0SegmDigiHitTOF","hME0SegmDigiHitTOF", 4000, -200, 200)  ;
    hME0SegmDigiHitPdgID = fs->make<TH1F>("hME0SegmDigiHitPdgID","hME0SegmDigiHitPdgID", 2221, 0, 2221)  ;
    hME0SegmDigiHitEta= fs->make<TH1F>("hME0SegmDigiHitEta","hME0SegmDigiHitEta", 600, -3, 3)  ;
    
    hME0Segm_phi_sig = fs->make<TH1F>("hME0Segm_phi_sig","hME0Segm_phi_sig", 800, -4, 4)  ;
    hME0Segm_X_sig   = fs->make<TH1F>("hME0Segm_X_sig","hME0Segm_X_sig", 800, -400, 400)  ;
    hME0Segm_Y_sig   = fs->make<TH1F>("hME0Segm_Y_sig","hME0Segm_Y_sig", 800, -400, 400)  ;
    hME0Segm_eta_sigCheck =  fs->make<TH1F>("hME0Segm_eta_sigCheck","hME0Segm_eta_sigCheck", 600, -3, 3)  ;
    hME0Segm_eta_sigCheckFail=  fs->make<TH1F>("hME0Segm_eta_sigCheckFail","hME0Segm_eta_sigCheckFail", 600, -3, 3)  ;
    
     hME0Segm_eta_bkgCheck =  fs->make<TH1F>("hME0Segm_eta_bkgCheck","hME0Segm_eta_bkgCheck", 600, -3, 3)  ;
     hME0Segm_eta_bkgCheckFail=  fs->make<TH1F>("hME0Segm_eta_bkgCheckFail","hME0Segm_eta_bkgCheckFail", 600, -3, 3)  ;
    
    
    hME0Segm_eta=fs->make<TH1F>("hME0Segm_eta","hME0Segm_eta", 600, -3, 3)  ;
    hDistME0Seg =fs->make<TH1F>("hDistME0Seg","hDistME0Seg", 500, 0, 100)  ;
    hSegmentPerChamber=fs->make<TH1F>("hSegmentPerChamber","hSegmentPerChamber", 31, -0.5, 30.5)  ;
    hSegmentPerChamber2=fs->make<TH1F>("hSegmentPerChamber2","hSegmentPerChamber2", 21, -0.5, 20.5)  ;
    hNumberOfMatchedHit =fs->make<TH1F>("hNumberOfMatchedHit","hNumberOfMatchedHit", 21, -0.5, 20.5)  ;
    hME0SegmPullPhi =fs->make<TH1F>("hME0SegmPullPhi","hME0SegmPullPhi", 1000, -1., 1.)  ;
    hME0SegmPullEta =fs->make<TH1F>("hME0SegmPullEta","hME0SegmPullEta", 1000, -1., 1.)  ;
    hME0SegmPullPhi_fitFail =fs->make<TH1F>("hME0SegmPullPhi_fitFail","hME0SegmPullPhi_fitFail", 1000, -1., 1.)  ;
    hME0SegmPullEta_fitFail =fs->make<TH1F>("hME0SegmPullEta_fitFail","hME0SegmPullEta_fitFail", 1000, -1., 1.)  ;
    
    hME0SegmPullPhi_BestMatch =fs->make<TH1F>("hME0SegmPullPhi_BestMatch","hME0SegmPullPhi_BestMatch", 1000, -1., 1.)  ;
    hME0SegmPullEta_BestMatch =fs->make<TH1F>("hME0SegmPullEta_BestMatch","hME0SegmPullEta_BestMatch", 1000, -1., 1.)  ;
    hNumberOfDuplicatesPerSignalSimTrack =fs->make<TH1F>("hNumberOfDuplicatesPerSignalSimTrack","hNumberOfDuplicatesPerSignalSimTrack", 20, -0.5, 19.5)  ;

    hNME0SegmentSignal_BestMatch =fs->make<TH1F>("hNME0SegmentSignal_BestMatch","hNME0SegmentSignal_BestMatch", 20, -0.5, 19.5)  ;
    hME0Segm_nHit_sig_BestMatch =fs->make<TH1F>("hME0Segm_nHit_sig_BestMatch","hME0Segm_nHit_sig_BestMatch", 20, -0.5, 19.5)  ;
    hME0Segm_chi2_sig_BestMatch  =fs->make<TH1F>("hME0Segm_chi2_sig_BestMatch","hME0Segm_chi2_sig_BestMatch", 100, -0.5, 99.5)  ;
    hME0Segm_time_sig_BestMatch  =fs->make<TH1F>("hME0Segm_time_sig_BestMatch","hME0Segm_time_sig_BestMatch", 4000, -200, 200)  ;

    
    hnumAllME0SimHits   =fs->make<TH1F>("hnumAllME0SimHits","hnumAllME0SimHits", 200, -0.5, 199.5)  ;
    hnumEleME0SimHits   =fs->make<TH1F>("hnumEleME0SimHits","hnumEleME0SimHits", 200, -0.5, 199.5)  ;
    hnumMuonME0SimHits  =fs->make<TH1F>("hnumMuonME0SimHits","hnumMuonME0SimHits", 200, -0.5, 199.5)  ;
    hnumPhotonME0SimHit =fs->make<TH1F>("hnumPhotonME0SimHit","hnumPhotonME0SimHit", 200, -0.5, 199.5)  ;
    hnumHadronME0SimHits =fs->make<TH1F>("hnumHadronME0SimHits","hnumHadronME0SimHits", 200, -0.5, 199.5)  ;
    hEleSimHitEta     =fs->make<TH1F>("hEleSimHitEta","hEleSimHitEta", 600, -3, 3)  ;
    hMuonSimHitEta    =fs->make<TH1F>("hMuonSimHitEta","hMuonSimHitEta", 600, -3, 3)  ;
    hPhotonSimHitEta  =fs->make<TH1F>("hPhotonSimHitEta","hPhotonSimHitEta", 600, -3, 3)  ;
    hHadronSimHitEta  =fs->make<TH1F>("hHadronSimHitEta","hHadronSimHitEta", 600, -3, 3)  ;
    hSimHitEta        =fs->make<TH1F>("hSimHitEta","hSimHitEta", 600, -3, 3)  ;

    hDuplicatesEta  =fs->make<TH1F>("hDuplicatesEta"," hDuplicatesEta", 600, -3, 3)  ;
    hDuplicatesPhi  =fs->make<TH1F>("hDuplicatesPhi"," hDuplicatesPhi", 800, -4, 4)  ;
    hDeltaRoll  =fs->make<TH1F>("hDeltaRoll","hDeltaRoll", 11 , -5.5, 5.5)  ;
    hDeltaRollvsDeltaEta =fs->make<TH2F>("hDeltaRollvsDeltaEta","hDeltaRollvsDeltaEta", 11 , -5.5, 5.5, 100,-0.5,0.5)  ;

}


// ------------ method called once each job just after ending the event loop  ------------
void ME0SegmentAnalyzerMuonGun::endJob() {
    /*
     rootfile->cd();
     mytree->Write();
     rootfile->Close();*/
    
}

// ------------ method called when starting to processes a run  ------------

//void  ME0SegmentAnalyzerMuonGun::beginRun(edm::Run const&, edm::EventSetup const&)
//{



// ------------ method called when ending the processing of a run  ------------

    /*
 void 
 ME0SegmentAnalyzerMuonGun::endRun(edm::Run const&, edm::EventSetup const&)
 {
 }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
 void 
 ME0SegmentAnalyzerMuonGun::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
 void 
 ME0SegmentAnalyzerMuonGun::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ME0SegmentAnalyzerMuonGun::fillDescriptions(edm::ConfigurationDescriptions& descriptions){
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE( ME0SegmentAnalyzerMuonGun);

