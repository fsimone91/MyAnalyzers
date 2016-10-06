//Package:     ME0SegmentAnalyzerSim
// Class:       ME0SegmentAnalyzerSim
// 
/**\class  ME0SegmentAnalyzerSim  ME0SegmentAnalyzerSim  ME0SegmentAnalyzerSim/ ME0SegmentAnalyzerSim/plugins/ ME0SegmentAnalyzerSim.cc
 
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

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoMuon/MuonIdentification/plugins/ME0MuonSelector.cc"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
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
class  ME0SegmentAnalyzerSim: public edm::EDAnalyzer {
public:
	explicit  ME0SegmentAnalyzerSim(const edm::ParameterSet&);
	~ ME0SegmentAnalyzerSim();
	void Initialize(); 
	
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
	
private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;
	std::string wp_;
	double timeMin_;
	double timeMax_;
	double minEta_;
	double maxEta_;
	
	edm::Service<TFileService> fs;
	
	TH1F *hNEvZmm;  TH1F * hDPhi; TH1F *  hme0machtMuonPt;  TH1F *  hme0machtMuonEta;  TH1F *  hme0machtMuonPhi;  TH1F *  hme0machtMuonCharge;
	TH1F *  hNME0Time ;  TH1F *  hNME0RecHits; TH1F *  hPtRes ;  
	TH1F *  hSimPt; TH1F *  hSimEta;
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
	TH1F * hRecHitTime;
	
	TH1F * hNDigiMatchedRH;	TH1F * hNPrompt;	TH1F * hN_noPrompt; TH1F * 	hNPromptMuHit;	TH1F * 	hNPromptNoMuHit;	TH1F * 	hNPromptHit_pdgId;	TH1F *  hME0SimTrackPdgID;	
	TH1F *  hSimElePtME0;TH1F *   hSimMuPtME0;	TH1F *   hSimPionPtME0;	TH1F * 	hSimEleNHitsME0;	TH1F * 	hSimMuonNHitsME0;	TH1F * 	hSimPionNHitsME0; TH1F * hPdgIDCheck; 
	TH2F *  hMuonDigiDPhiVsPT; TH2F * hNoEleDigiDPhiVsPT;
	
	TH1F * hDRME0SimTrack;	TH1F * hDRME0SimMuonEle;	TH1F * hN_noEleHit;	TH1F * hN_noPromptHit_pdgId;	TH1F *	hSegmentComposition;	TH1F * hNME0SegmentAfterDigiMatch;	TH1F * hMuEleinME0Segm;	TH1F * hMuEleOthersinME0Segm;TH1F * hMuOnlyinME0Segm;
	
	TH1F * hNEleBrem;	TH1F * hNEleDeltaRays;	TH1F * hNEle; TH1F * hMatchedSimEleTrack;
	
	TH1F * hRHDeltaPhiSameLayer;	TH1F * hRHDeltaEtaSameLayer;	TH1F * hRHDeltaTSameLayer;
	
	TH1F * hNoPromptRecHitTime;	TH1F * hMuonRecHitTime;	TH1F * hEleRecHitTime;	TH1F * hNMuonSameLayerTOF;	TH1F * hNEleSameLayerTOF;	TH1F * hNoPromptSameLayerTOF;
	
	TH1F * hDeltaPhiSimReco;	TH1F * hDeltaEtaSimReco;TH1F * hDeltaXSimReco; TH1F * hDeltaYSimReco; 
	TH1F * hNoMuinME0Segm; TH2F * hMuonDigiDXVsPT;	TH2F * hMuonDigiLocalDPhiVsPT; 	TH2F * hMuonDigiLocalDXVsPT; TH1F * hDeltaXSimRecoLocal; TH1F * hDeltaYSimRecoLocal;
	
	TH1F * 	hBeamSpotX0; 	TH1F * 	hBeamSpotY0; 	TH2F * 	hBeamSpotX0Y0;	TH1F * 	hBeamSpotZ0; 	TH1F * 	hBeamSpotSigmaZ;
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
	
	TH1F * 	hSimLocalDX;	TH1F * 	hLocalDX;TH1F * 	hSimLocalDY;	TH1F * 	hLocalDY;TH1F * 	hSimLocalDXPos;	TH1F * 	hLocalDXPos;TH1F * 	hSimLocalDXNeg;TH1F * hLocalDXNeg;
    TH1F *    hSimMuonME0_Pt;
    TH1F *    hSimMuonME0_Eta; TH1F *    hSimMuEtaME0;
  TH1F *    hGenMuME0Pt;TH1F *    hGenMuME0Energy; TH1F * hGenMuME0Pt_genSimMatch;
	
	
	// virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	
	// ----------member data ---------------------------
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
 ME0SegmentAnalyzerSim:: ME0SegmentAnalyzerSim(const edm::ParameterSet& iConfig):
wp_(iConfig.getParameter<std::string>("wp")),
timeMin_(iConfig.getParameter<double>("timeMin")),
timeMax_(iConfig.getParameter<double>("timeMax")),
minEta_(iConfig.getParameter<double>("minEta")),
maxEta_(iConfig.getParameter<double>("maxEta"))

{
	
}


 ME0SegmentAnalyzerSim::~ ME0SegmentAnalyzerSim()
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


edm::PSimHitContainer isTrackMatched(SimTrackContainer::const_iterator simTrack, const Event & event, const EventSetup& eventSetup)
{
	edm::PSimHitContainer selectedME0Hits;
	
	edm::ESHandle<ME0Geometry> me0geom;
	eventSetup.get<MuonGeometryRecord>().get(me0geom);
	
	edm::Handle<edm::PSimHitContainer> ME0Hits;
	event.getByLabel(edm::InputTag("g4SimHits","MuonME0Hits"), ME0Hits);
	
	ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
	eventSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
	
	for (edm::PSimHitContainer::const_iterator itHit =  ME0Hits->begin(); itHit != ME0Hits->end(); ++itHit){
		
		DetId id = DetId(itHit->detUnitId());
		if (!(id.subdetId() == MuonSubdetId::ME0)) continue;
		if(itHit->particleType() != (*simTrack).type()) continue;
		
		bool result = isSimMatched(simTrack, itHit);
		if(result) selectedME0Hits.push_back(*itHit);
		
	}
	
	//	std::cout<<"N simHit in ME0segm : "<<selectedME0Hits.size()<<std::endl;
	return selectedME0Hits;
	
}

struct MyME0Digi
{
	Int_t detId, particleType;
	Short_t layer;
	Float_t g_eta, g_phi;
	Float_t tof;
	Float_t prompt;
};

void
 ME0SegmentAnalyzerSim::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	//Initialize();
	using namespace edm;
	
	edm::Handle<SimTrackContainer> simTracks;
	iEvent.getByLabel("g4SimHits",simTracks);
	
	
	
	edm::Handle <reco::GenParticleCollection> genparticles;
	iEvent.getByLabel("genParticles",genparticles);
	
	
	
	edm::Handle<ME0DigiPreRecoCollection> me0_digis;
	iEvent.getByLabel("simMuonME0Digis",  me0_digis); 
	
	edm::Handle<edm::PSimHitContainer>  ME0HitsCollection;
	iEvent.getByLabel(edm::InputTag("g4SimHits","MuonME0Hits"), ME0HitsCollection);
	
	ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
	iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
	
	
	edm::ESHandle<ME0Geometry> me0geom;
	iSetup.get<MuonGeometryRecord>().get(me0geom);
//	const ME0Geometry* me0Geom;
//	me0Geom= &*me0geom;
	
	
	
	/*
	 Run   = iEvent.id().run();
	 Event = iEvent.id().event();
	 Lumi  = iEvent.luminosityBlock();
	 Bunch = iEvent.bunchCrossing();
	 */
	
	std::cout<<"********************************************BeginEvent="<<iEvent.id().event()<<"******************************************** "<<std::endl;
	std::vector<int> indexgenmu;
	hNEv->Fill(1);  
		
	SimTrackContainer::const_iterator simTrack;
	double numberOfSimTracks =0.;
	//std::cout<<" num Simulated tracks: "<<simTracks->size()<<std::endl;
	std::vector<double> simHitPhi, simHitElePhi, simHitX, simHitLocalX, simTrackEta;
	std::vector<double> simHitEta, simHitEleEta, simHitY, simHitLocalY, simTrackPt;
	std::vector<int> numME0SimHits;
	SimTrackContainer ME0Tracks, ME0EleSimTracks;
	for (simTrack = simTracks->begin(); simTrack != simTracks->end(); ++simTrack){
		
		
		//if ((*simTrack).noVertex()) continue;
		if ((*simTrack).noGenpart()) continue;
		if (!(abs((*simTrack).type()) == 13)) continue;
        
		hSimEta->Fill((*simTrack).momentum().eta());
		hSimPt->Fill((*simTrack).momentum().pt());

		cout<<"Sim pT: "<<(*simTrack).momentum().pt()<<" carica="<<(*simTrack).charge()<<" energia="<<(*simTrack).momentum().energy()<<endl;
		double simPt=(*simTrack).momentum().pt();
		
		//cout<<"Sim Eta: "<<(*simTrack).momentum().eta()<<endl;
		double simEta = (*simTrack).momentum().eta();
		
		//cout<<"Sim Phi: "<<(*simTrack).momentum().phi()<<endl;
		//double simPhi = (*simTrack).momentum().phi();
		
		if (abs(simEta) > maxEta_ || abs(simEta) < minEta_) continue;
		
		
		if (fabs((*simTrack).type())==11)	{hSimElePtME0->Fill(simPt);	 ME0EleSimTracks.push_back(*simTrack);}
		if (fabs((*simTrack).type())==13)	{
                hSimMuPtME0->Fill(simPt); hSimMuEtaME0->Fill(simEta);  }
		if (fabs((*simTrack).type())==211)	{hSimPionPtME0->Fill(simPt); } 
		hME0SimTrackPdgID->Fill( fabs((*simTrack).type()) );
		
		edm::PSimHitContainer selME0SimHits = isTrackMatched(simTrack, iEvent , iSetup);
		
		//int ME0SimHitsize = selME0SimHits.size();
		//std::cout<<"# me0 sim hits="<< selME0SimHits.size() <<std::endl;
		
		if( selME0SimHits.size() ==0 ) continue;
		simTrackEta.push_back(simEta); simTrackPt.push_back(simPt);
		
					
		ME0Tracks.push_back(*simTrack);
		numberOfSimTracks++;
		//std::cout<<"TrackID="<<simTrack->trackId()<<" track type="<<(*simTrack).type()<<"# me0 sim hits="<< selME0SimHits.size() <<" simPt="<<simPt<<std::endl;
		simHitPhi.clear(); simHitEta.clear();
		simHitX.clear();   simHitY.clear(); simHitLocalX.clear();  simHitLocalY.clear(); 
		int selhitcounter =0;
	
		for (edm::PSimHitContainer::const_iterator itHit =  selME0SimHits.begin(); itHit != selME0SimHits.end(); ++itHit){
			ME0DetId idme0 = ME0DetId(itHit->detUnitId());	
			int layer_sim = idme0.layer();		
			LocalPoint lp = itHit->entryPoint();
			GlobalPoint hitGP_sim( me0geom->idToDet(itHit->detUnitId())->surface().toGlobal(lp));		
			simHitPhi.push_back(hitGP_sim.phi());
			simHitEta.push_back(hitGP_sim.eta());
			simHitX.push_back(hitGP_sim.x());
			simHitY.push_back(hitGP_sim.y());
			simHitLocalX.push_back(lp.x());
			simHitLocalY.push_back(lp.y());
			selhitcounter	++;
			std::cout<<"track pt="<<(*simTrack).momentum().pt()<<" simHit eta="<<hitGP_sim.eta()<<" phi="<<hitGP_sim.phi()<<" simHit detID="<<idme0<<" layer sim="<<layer_sim<<" localX="<<lp.x()<<" trackID="<<itHit->trackId()<<" Q="<<(*simTrack).charge()<<std::endl;
		}
		numME0SimHits.push_back(selhitcounter);
		int sizeSimPhi = simHitPhi.size()-1;
		double SimDeltaPhi = TMath::Abs( simHitPhi[0] - simHitPhi[sizeSimPhi] );
		double SimDeltaEta= TMath::Abs( simHitEta[0] - simHitEta[sizeSimPhi] );

		 
		cout<<"Track:"<<simTrack->trackId()<<" pdgID="<<(*simTrack).type()<<" SimDPhi="<<simHitPhi[0] - simHitPhi[sizeSimPhi]<<" pT="<<simPt<<" carica="<<(*simTrack).charge()<<" DX(l1,l6)"<< simHitLocalX[0] - simHitLocalX[sizeSimPhi]<<endl;
		
		
			hSimPtVSDphi->Fill(SimDeltaPhi, simPt);
			hSimPtVSDeta->Fill(SimDeltaEta, simPt);

		
			hSimDEta->Fill( simHitEta[0] - simHitEta[sizeSimPhi]);
			hSimDPhi->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);
			hSimDX->Fill( simHitX[0] - simHitX[sizeSimPhi]);
			hSimDY->Fill( simHitY[0] - simHitX[sizeSimPhi]);
		
			hSimLocalDX->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);
			hSimLocalDY->Fill( simHitLocalY[0] - simHitLocalY[sizeSimPhi]);
		
			if ((*simTrack).charge() >0) {hSimDPhiPos->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]); hSimLocalDXPos->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
			if ((*simTrack).charge() <0) {hSimDPhiNeg->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]); hSimLocalDXNeg->Fill( simHitLocalX[0] - simHitLocalX[sizeSimPhi]);}
		
			if (((*simTrack).charge() >0) && ((*simTrack).momentum().pt()) >20 ) hSimDPhiPos_HighPt->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);
			if (((*simTrack).charge() >0) && ((*simTrack).momentum().pt()) <10 ) hSimDPhiPos_LowPt->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);
		
			if (((*simTrack).charge() <0) && ((*simTrack).momentum().pt()) >20 ) hSimDPhiNeg_HighPt->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);
			if (((*simTrack).charge() <0) && ((*simTrack).momentum().pt()) <10 ) hSimDPhiNeg_LowPt->Fill( simHitPhi[0] - simHitPhi[sizeSimPhi]);
	
		
	}
//	std::cout<<" Num simTrack in ME0  "<<numberOfSimTracks<<" ME0Track.size()="<<ME0Tracks.size()<<" Num me0 segm="<<me0Segments->size()<<std::endl;
	if (ME0Tracks.size()>0) hSelectedSimTrack->Fill(ME0Tracks.size());
	for(uint i =0; i< ME0Tracks.size(); i++){
        hSimMuonME0_Pt->Fill(ME0Tracks.at(i).momentum().pt());
        hSimMuonME0_Eta->Fill(ME0Tracks.at(i).momentum().eta());
    }
 //   if (ME0Tracks.size() ==3) hSelectedSimTrack->Fill(ME0Tracks.size());

	
	//////////////////////////////////////////////////////////////////////////////////////////Loop over gen particles//////////////////////////////////////////////////////////////////////////////////////////
	
	
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
	  std::cout<<i<<" particle= "<<genparticles->at(indexgenmu[i]).pdgId()<<" status="<<genparticles->at(indexgenmu[i]).status()<<" eta="<<genparticles->at(indexgenmu[i]).eta()<<" phi="<<genparticles->at(indexgenmu[i]).phi()<<" pt="<<genparticles->at(indexgenmu[i]).pt()<<std::endl;
	  hGenMuPt->Fill( (genparticles->at(indexgenmu[i]).pt()));
	  hGenMuEta->Fill( (genparticles->at(indexgenmu[i]).eta())  );
		
		
		if((abs(genparticles->at(indexgenmu[i]).eta())<maxEta_ ) && ( (genparticles->at(indexgenmu[i]).eta() > minEta_ )||(genparticles->at(indexgenmu[i]).eta()< (-minEta_ )) )  ){
			indexGenMuInME0.push_back(indexgenmu[i]);
			std::cout<<i<<" --------particle in ME0= "<<genparticles->at(indexgenmu[i]).pdgId()<<" status="<<genparticles->at(indexgenmu[i]).status()<<" eta="<<genparticles->at(indexgenmu[i]).eta()<<" phi="<<genparticles->at(indexgenmu[i]).phi()<<" pt="<<genparticles->at(indexgenmu[i]).pt()<<std::endl;
			hGenMuME0Pt->Fill(genparticles->at(indexgenmu[i]).pt());
			hGenMuME0Energy->Fill(genparticles->at(indexgenmu[i]).energy());


			for(uint i =0; i< ME0Tracks.size(); i++){
			  double DEtaGenSim = ME0Tracks.at(i).momentum().eta() -genparticles->at(indexgenmu[i]).eta() ;
			  double DPhiGenSim = ME0Tracks.at(i).momentum().phi() -genparticles->at(indexgenmu[i]).phi() ;
			  double DRGenSim = TMath::Sqrt(DEtaGenSim*DEtaGenSim + DPhiGenSim*DPhiGenSim);
			  if(DRGenSim<0.001){
			    hGenMuME0Pt_genSimMatch->Fill(genparticles->at(indexgenmu[i]).pt());
			    //			    std::cout<<"  "<<
			  }

			}

			}else{
			indexGenMuElseWhere.push_back(indexgenmu[i]);
			}
	}
	
	hGenMuonsME0->Fill(indexGenMuInME0.size());
	hNGenMu->Fill(indexgenmu.size());
	}



// ------------ method called once each job just before starting event loop  ------------
void 
 ME0SegmentAnalyzerSim::beginJob()
{

  hGenMuME0Pt= fs->make<TH1F>("hGenMuME0Pt","hGenMuME0Pt",100,0,5);
  hGenMuME0Energy= fs->make<TH1F>("hGenMuME0Energy","hGenMuME0Energy",500,0,50);
  hGenMuME0Pt_genSimMatch = fs->make<TH1F>("hGenMuME0Pt_genSimMatch","hGenMuME0Pt_genSimMatch",100,0,5);
  hNEvZmm = fs->make<TH1F>("hNEvZmm","hNEvZmm",10,0,10); 
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
	
	hSimEta  = fs->make<TH1F>("hSimEta","hSimEta",100,-4,4);  
	hSimPt  = fs->make<TH1F>("hSimPt","hSimPt",200, 0,200);  
	hPtVSDphi = fs->make<TH2F>("hPtVSDphi","hPtVSDphi",5000, 0, 0.05 , 200,0,200); 
	hSimPtVSDphi = fs->make<TH2F>("hSimPtVSDphi","hSimPtVSDphi",5000, 0, 0.05 , 200,0,200); 
	hSimPtVSDeta = fs->make<TH2F>("hSimPtVSDeta","hSimPtVSDeta",10000, 0, 0.1 , 200,0,200); 
	hSimDEta = fs->make<TH1F>("hSimDEta","hSimDEta",1000,-0.5,0.5);  
	hSimDY = fs->make<TH1F>("hSimDY","hSimDY",1000,-20,20);  
	hSimDX = fs->make<TH1F>("hSimDX","hSimDX",1000,-20,20); 
	hSimDPhi = fs->make<TH1F>("hSimDPhi","hSimDPhi",1000,-0.1,0.1);
	hSimDPhiPos = fs->make<TH1F>("hSimDPhiPos","hSimDPhiPos",1000,-0.05,0.05);  
	hSimDPhiNeg = fs->make<TH1F>("hSimDPhiNeg","hSimDPhiNeg",1000,-0.05,0.05);  
	
	hSimDPhiPos_LowPt = fs->make<TH1F>("hSimDPhiPos_LowPt","hSimDPhiPos_LowPt",1000,-0.05,0.05);  
	hSimDPhiPos_HighPt = fs->make<TH1F>("hSimDPhiPos_HighPt","hSimDPhiPos_HighPt",1000,-0.05,0.05);  

	hSimDPhiNeg_LowPt = fs->make<TH1F>("hSimDPhiNeg_LowPt","hSimDPhiNeg_LowPt",1000,-0.05,0.05);  
	hSimDPhiNeg_HighPt = fs->make<TH1F>("hSimDPhiNeg_HighPt","hSimDPhiNeg_HighPt",1000,-0.05,0.05);  
	
	hDPhi = fs->make<TH1F>("hDPhi","hDPhi",1000,-0.05,0.05); 
	hDPhiPos = fs->make<TH1F>("hDPhiPos","hDPhiPos",1000,-0.05,0.05); 
	hDPhiNeg = fs->make<TH1F>("hDPhiNeg","hDPhiNeg",1000,-0.05,0.05); 
	
	hDPhiPos_LowPt = fs->make<TH1F>("hDPhiPos_LowPt","hDPhiPos_LowPt",1000,-0.05,0.05); 
	hDPhiPos_HighPt = fs->make<TH1F>("hDPhiPos_HighPt","hDPhiPos_HighPt",1000,-0.05,0.05); 
	
	hDPhiNeg_LowPt = fs->make<TH1F>("hDPhiNeg_LowPt","hDPhiPos_NegPt",1000,-0.05,0.05); 
	hDPhiNeg_HighPt = fs->make<TH1F>("hDPhiNeg_HighPt","hDPhiNeg_HighPt",1000,-0.05,0.05); 
	
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
	
	
	
	hNME0Time = fs->make<TH1F>("hNME0Time","hNME0Time",300,0,300);  
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
	
	
	hSelectedSimTrack =  fs->make<TH1F>("hSelectedSimTrack","hSelectedSimTrack",20,0,20);
	hNME0Segment =  fs->make<TH1F>("hNME0Segment","hNME0Segment ",20,0,20);
	hGenMuonsME0= fs->make<TH1F>("hGenMuonsME0","hGenMuonsME0",20,0,20);
	hME0MuonsInMatchCone = fs->make<TH1F>("hME0MuonsInMatchCone","hME0MuonsInMatchCone",20,0,20);
	hME0MuonsID = fs->make<TH1F>("hME0MuonsID","hME0MuonsID",20,0,20);
	hME0SimTrackPdgID = fs->make<TH1F>("hME0SimTrackPdgID","hME0SimTrackPdgID",3000,0,3000);
	
	hSimElePtME0  = fs->make<TH1F>("hSimElePtME0","hSimElePtME0",100,0,10);			
	hSimMuPtME0   = fs->make<TH1F>("hSimMuPtME0","hSimMuPtME0",100,0,5);
    hSimMuEtaME0   = fs->make<TH1F>("hSimMuEtaME0","hSimMuEtaME0",100,-5,5);
	hSimPionPtME0 = fs->make<TH1F>("hSimPionPtME0","hSimPionPtME0",100,0,10);
	
	hSimEleNHitsME0 = fs->make<TH1F>("hSimEleNHitsME0","hSimEleNHitsME0",10,0,10);	
	hSimMuonNHitsME0 = fs->make<TH1F>("hSimMuonNHitsME0","hSimMuonNHitsME0",10,0,10);
	hSimPionNHitsME0 = fs->make<TH1F>("hSimPionNHitsME0","hSimPionNHitsME0",10,0,10);
	hPdgIDCheck= fs->make<TH1F>("hPdgIDCheck","hPdgIDCheck",220,0,220);
	hDRME0SimTrack  = fs->make<TH1F>("hDRME0SimTrack","hDRME0SimTrack",200,0,10);
	hDRME0SimMuonEle= fs->make<TH1F>("hDRME0SimMuonEle","hDRME0SimMuonEle",200,0,10);
	
	hNME0SegmentAfterDigiMatch= fs->make<TH1F>("hNME0SegmentAfterDigiMatch","hNME0SegmentAfterDigiMatch",3,0,3);
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
	hverteX  = fs->make<TH1F>("hverteX","hverteX",100,-10,10);
	hverteY  = fs->make<TH1F>("hverteY","hverteY",100,-10,10);
	hverteXY  = fs->make<TH2F>("hverteXY","hverteXY",100,-10,10,100, -10,10 );
	hverteZ = fs->make<TH1F>("hverteZ","hverteZ",100,-30,30);
								   
	
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


	
	hSimLocalDY= fs->make<TH1F>("hSimLocalDY","hSimLocalDY",1000,-50,50);
	hLocalDY= fs->make<TH1F>("hLocalDY","hLocalDY",1000,-50,50);
    
    hSimMuonME0_Pt  = fs->make<TH1F>("hSimMuonME0_Pt","hSimMuonME0_Pt ",100,0,5);
    hSimMuonME0_Eta  = fs->make<TH1F>("hSimMuonME0_Eta "," hSimMuonME0_Eta",100,-5,5);


}



// ------------ method called once each job just after ending the event loop  ------------
void 
 ME0SegmentAnalyzerSim::endJob() 
{
	/*
	 rootfile->cd();
	 mytree->Write();
	 rootfile->Close();*/
	
}

// ------------ method called when starting to processes a run  ------------

//void  ME0SegmentAnalyzerSim::beginRun(edm::Run const&, edm::EventSetup const&)
//{



// ------------ method called when ending the processing of a run  ------------
/*
 void 
  ME0SegmentAnalyzerSim::endRun(edm::Run const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method called when starting to processes a luminosity block  ------------
/*
 void 
  ME0SegmentAnalyzerSim::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
 void 
  ME0SegmentAnalyzerSim::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
 ME0SegmentAnalyzerSim::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE( ME0SegmentAnalyzerSim);
