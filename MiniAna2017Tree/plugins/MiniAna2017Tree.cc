// -*- C++ -*-
//
// Package:    MiniAna2017/MiniAna2017Tree
// Class:      MiniAna2017Tree
// 
/**\class MiniAna2017Tree MiniAna2017Tree.cc MiniAna2017/MiniAna2017Tree/plugins/MiniAna2017Tree.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  venditti
//         Created:  Tue, 18 Dec 2018 09:30:06 GMT
//
//


// system include files
#include <memory>
#include <algorithm> 
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TFile.h"
#include "TH1.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"


#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicTree.h"


#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"


#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

class MiniAna2017Tree : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MiniAna2017Tree(const edm::ParameterSet&);
      ~MiniAna2017Tree();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  float dR(float eta1, float eta2, float phi1, float phi2);


  private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;  
  virtual void endJob() override;
  edm::EDGetTokenT<edm::View<pat::Muon> > muons_; 
  edm::EDGetTokenT<edm::View<reco::Vertex> > vertex_; 
  edm::EDGetTokenT<edm::View<reco::Track> > trackToken_;
  edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > Cand3Mu_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticles_;
  bool isMc;
  //  edm::EDGetTokenT<edm::TriggerResults> trigResultsToken;
  //  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigObjCollToken;
  const TransientTrackBuilder* theTransientTrackBuilder_; 


  edm::Service<TFileService> fs;

  TH1F *hEvents;
  
  /*
  TH1F *hEvents_3Mu, *hEvents_MuFilter, *hEvents_DiMuonDz, * hEvent_TauCharge,  * hEvent_3MuonVtx, *  hEvent_3MuonMass, *hEvent_Valid3MuonVtx;

  TH1F *hGenMuonPt;  TH1F *hGenMuonEta;
  TH1F *hMuonPt;TH1F *hMuonP; TH1F *hMuonEta, *hGlobalMuonEta, *hTrackerMuonEta, *hLooseMuonEta, *hSoftMuonEta;
  TH1F *  hGlobalMuonPt,  *  hSoftMuonPt,  *  hLooseMuonPt, *  hTrackerMuonPt;
  TH1F *hMuonNumberOfValidHits;
  TH1F * hMuonDB;   TH1F * hMuonTime;   TH1F * hMuonTimeErr;
  TH1F  *hGenTauPt, *hDiMuonDz, *hDiMuonDR, *hGoodMuSize, *hMuonSize_GoodDiMu;
  TH1F  *hThreeMuonInvMass, *hThreeMuonCharge, * hVtxSize, * hTau_vFit_chi2 ,  *h3MuonMassVtxFit  ;
  TH2F * hMuonTimeVsP;
  */
  //tree 
  TTree*      tree_;
  std::vector<float>  MuonPt, MuonEta, MuonPhi, MuonChi2P, MuonChi2LocalPosition, MuonGlbTrackProbability, MuonTrkRelChi2, MuonTrkKink;
  std::vector<double> MuonEnergy,  MuonCharge;

//Vtx position                                                                                                                                                                               
  std::vector<double>  Muon_vx,  Muon_vy,  Muon_vz;

  //MuonID                                                                                                                                                                                    
  std::vector<double>  Muon_isGlobal,  Muon_isTracker,  Muon_isSoft,  Muon_isLoose,  Muon_isPF,  Muon_isRPCMuon,  Muon_isStandAloneMuon,  Muon_isTrackerMuon,  Muon_isCaloMuon,  Muon_isQualityValid,  Muon_isTimeValid,  Muon_isIsolationValid,  Muon_numberOfMatchedStations,  Muon_numberOfMatches;

  //MuonTime                                                                                                                                                                                   
  std::vector<double>  Muon_timeAtIpInOut,Muon_timeAtIpInOutErr;

  //Muon inner + outer track                                                                                                                                                                   
  std::vector<double>  Muon_GLnormChi2, Muon_GLhitPattern_numberOfValidMuonHits,  Muon_trackerLayersWithMeasurement,  Muon_Numberofvalidpixelhits,  Muon_outerTrack_p,  Muon_outerTrack_eta,
    Muon_outerTrack_phi,  Muon_outerTrack_normalizedChi2,  Muon_outerTrack_muonStationsWithValidHits,  Muon_innerTrack_p,  Muon_innerTrack_eta,  Muon_innerTrack_phi,  Muon_innerTrack_normalizedChi2,  Muon_QInnerOuter;

  std::vector<double>   Muon_combinedQuality_updatedSta,  Muon_combinedQuality_trkKink,  Muon_combinedQuality_glbKink,  Muon_combinedQuality_trkRelChi2,  Muon_combinedQuality_staRelChi2,  Muon_combinedQuality_chi2LocalPosition,  Muon_combinedQuality_chi2LocalMomentum,  Muon_combinedQuality_localDistance,  Muon_combinedQuality_globalDeltaEtaPhi,  Muon_combinedQuality_tightMatch,  Muon_combinedQuality_glbTrackProbability,  Muon_calEnergy_em,  Muon_calEnergy_emS9,  Muon_calEnergy_emS25,  Muon_calEnergy_had,  Muon_calEnergy_hadS9,  Muon_segmentCompatibility,  Muon_caloCompatibility,  Muon_ptErrOverPt;

  std::vector<double>  Mu1_Pt,  Mu1_Eta,  Mu1_Phi,  Mu2_Pt,  Mu2_Eta,  Mu2_Phi,  Mu3_Pt,  Mu3_Eta,  Mu3_Phi,  Mu1_SimPt,  Mu1_SimEta,  Mu1_SimPhi,  Mu2_SimPt,  Mu2_SimEta,  Mu2_SimPhi,
    Mu3_SimPt,  Mu3_SimEta,  Mu3_SimPhi;

  int TripletCollectionSize, PVCollection_Size, MuonCollectionSize;
  std::vector<double>  TripletVtx_x,  TripletVtx_y,  TripletVtx_z,  TripletVtx_Chi2,  TripletVtx_NDOF,  Triplet_Mass,  Triplet_Pt,  Triplet_Eta,  Triplet_Phi, Triplet_Charge;


  std::vector<double>  RefittedPV_x;
  std::vector<double>  RefittedPV_y;
  std::vector<double>  RefittedPV_z;
  std::vector<double>  RefittedPV_NTracks;
  //RefittedPV_Chi2.push_back(PVertex.);                                                                                                                                                   

  std::vector<double>  FlightDistPVSV;
  std::vector<double>  FlightDistPVSV_Err;
  std::vector<double>  FlightDistPVSV_Significance;
  std::vector<double>  FlightDistPVSV_chi2;

  double PV_x,  PV_y,  PV_z,  PV_NTracks;

  //SyncTree
  /*  TTree*      SyncTree_;
  std::vector<float>  allmuons_pt, leadmuon_pt, leadmuon_phi, leadmuon_eta;
  std::vector<float>  alltracks_pt, leadtrack_pt,  leadtrack_eta,  leadtrack_phi;
  uint nprimevtxs, nmuons, evt, run, lumi;
  */
};



MiniAna2017Tree::MiniAna2017Tree(const edm::ParameterSet& iConfig){
  isMc = iConfig.getUntrackedParameter<bool>("isMcLabel");
  muons_ = consumes<edm::View<pat::Muon> >  (iConfig.getParameter<edm::InputTag>("muonLabel"));
  vertex_ = consumes<edm::View<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("VertexLabel"));
  trackToken_ = consumes<edm::View<reco::Track> > (edm::InputTag("generalTracks"));
  genParticles_ = consumes<edm::View<reco::GenParticle>  > (iConfig.getParameter<edm::InputTag>("genParticleLabel"));
  Cand3Mu_ = consumes<edm::View<reco::CompositeCandidate> > (iConfig.getParameter<edm::InputTag>("Cand3MuLabel"));
  //tauToken_(consumes(iConfig.getParameter("taus"))),
  //metToken_(consumes(iConfig.getParameter("mets")))
  //tree_(0);
  //MuonPt(0);
  }
     MiniAna2017Tree::~MiniAna2017Tree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


float MiniAna2017Tree::dR(float eta1, float eta2, float phi1, float phi2){
  float dphi=(phi1-phi2);
  float deta=(eta1-eta2);
  float deltaR= TMath::Sqrt(dphi*dphi + deta*deta);
  return deltaR;
}

bool isGoodTrack(const reco::Track &track) {
  if(track.pt()>1){
    if(std::fabs(track.eta())<2.4){
      if(track.hitPattern().trackerLayersWithMeasurement()>5){
	if(track.hitPattern().pixelLayersWithMeasurement()>1) return true;
      }
    }
  }
  return false;
}

typedef std::map<const reco::Track*, reco::TransientTrack> TransientTrackMap;
// auxiliary function to exclude tracks associated to tau lepton decay "leg"
// from primary event vertex refit                                                                                                                                          
bool tracksMatchByDeltaR(const reco::Track* trk1, const reco::Track* trk2)
{
  if ( reco::deltaR(*trk1, *trk2) < 1.e-2 && trk1->charge() == trk2->charge() ) return true;
  else return false;
}

void removeTracks(TransientTrackMap& pvTracks_toRefit, const std::vector<reco::Track*> svTracks)
{
  for ( std::vector<reco::Track*>::const_iterator svTrack = svTracks.begin(); svTrack != svTracks.end(); ++svTrack ){
    //--- remove track from list of tracks included in primary event vertex refit  
    //    if track matches by reference or in eta-phi                                                                                                                      
    //    any of the tracks associated to tau lepton decay "leg"
    for ( TransientTrackMap::iterator pvTrack = pvTracks_toRefit.begin(); pvTrack != pvTracks_toRefit.end(); ++pvTrack ) {
      if ( tracksMatchByDeltaR(pvTrack->first, *svTrack) ) {
        pvTracks_toRefit.erase(pvTrack);
        break;
      }
    }
}
}




void
MiniAna2017Tree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using std::vector;



  edm::Handle< edm::View<reco::Vertex> >vertices;
  iEvent.getByToken(vertex_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vertices->front();

  edm::Handle< edm::View<pat::Muon> > muons;
  iEvent.getByToken(muons_, muons);


  edm::Handle<edm::View<reco::CompositeCandidate> > Cand3Mu;
  iEvent.getByToken(Cand3Mu_, Cand3Mu);
  
  
  edm::Handle< edm::View<reco::GenParticle> > genParticles;
  iEvent.getByToken(genParticles_, genParticles);
  

  edm::Handle<edm::View<reco::Track> > trackCollection;
  iEvent.getByToken(trackToken_, trackCollection);

  edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);
  theTransientTrackBuilder_ = theTransientTrackBuilder.product();


  hEvents->Fill(1);

  
  if(isMc){
    uint j=0; 
    uint ngenP=genParticles->size();
    std::vector<int> genPidx;
    //uint tauRaw=-999; //int DsRaw=-999;

    cout<<"****************GenLevel Info Begin********************"<<endl;
    for(edm::View<reco::GenParticle>::const_iterator gp=genParticles->begin(); gp!=genParticles->end(), j<ngenP; ++gp , ++j){
      //if(fabs(gp->pdgId())==15) tauRaw = j;
    
      if(fabs(gp->pdgId())==13) {
	for (uint i=0; i<gp->numberOfMothers();i++){
	  if(fabs(gp->mother(i)->pdgId())==15) {
	    /*hGenMuonPt->Fill(gp->p()); 
	    hGenMuonEta->Fill(gp->eta());
	    hGenTauPt->Fill(gp->mother(i)->pt());*/
	    std::cout<<j<<"--genMu pt="<<gp->pt()<<" eta="<<gp->eta()<<" phi="<<gp->phi()<<" pdgID="<<gp->pdgId()<<" tau pt="<<gp->mother(i)->pt()<<" mu vtx_x="<<gp->vx()<<" mu vtx_y="<<gp->vy()<<" mu vtx_z="<<gp->vz()<<endl;
	    //cout<<tauRaw<<"--mother pdgID="<<gp->mother(i)->pdgId()<<" mother vtx_x="<<gp->mother(i)->vx()<<" vy="<<gp->mother(i)->vy()<<" vz="<<gp->mother(i)->vz()<<endl;
	    //cout<<"TauMother pdgId="<<gp->mother(i)->mother(0)->pdgId()<<" vx="<<gp->mother(i)->mother(0)->vx()<<" vy="<<gp->mother(i)->mother(0)->vy()<<" vz="<<gp->mother(i)->mother(0)->vz()<<endl;
	    genPidx.push_back(j);
	  }
	}
      }
    }
  
    if(genPidx.size()==3){
      float genDR12=dR((*genParticles)[genPidx[0]].eta(), (*genParticles)[genPidx[1]].eta(), (*genParticles)[genPidx[0]].phi(), (*genParticles)[genPidx[1]].phi()); 
      float genDR23=dR((*genParticles)[genPidx[0]].eta(), (*genParticles)[genPidx[2]].eta(), (*genParticles)[genPidx[0]].phi(), (*genParticles)[genPidx[2]].phi());
      float genDR13=dR((*genParticles)[genPidx[1]].eta(), (*genParticles)[genPidx[2]].eta(), (*genParticles)[genPidx[1]].phi(), (*genParticles)[genPidx[2]].phi());
      //cout<<"--genDR12="<<genDR12<<" genDR23="<<genDR23<<" genDR13="<<genDR13<<endl;
    }
    cout<<"****************GenLevel Info End ********************"<<endl;
  }
 

  //Primary Vtx
  std::vector<reco::TransientTrack> pvTracks_original;
  TransientTrackMap pvTrackMap_refit;
  //  const reco::Vertex* eventVertex;

  PVCollection_Size = vertices->size();

  for ( reco::Vertex::trackRef_iterator pvTrack = (*vertices)[0].tracks_begin(); pvTrack != (*vertices)[0].tracks_end(); ++pvTrack ) {
    reco::TransientTrack pvTrack_transient =theTransientTrackBuilder_->build(pvTrack->get());
    pvTracks_original.push_back(pvTrack_transient);
    pvTrackMap_refit.insert(std::make_pair(pvTrack->get(), pvTrack_transient));
  }
  //  cout<<" Number of tracks associated to the PV="<<pvTracks_original.size()<<endl;

  PV_x = (*vertices)[0].x();
  PV_y = (*vertices)[0].y();
  PV_z = (*vertices)[0].z();
  PV_NTracks = pvTracks_original.size();
  // PV_Chi2.push_back((*vertices)[0].Chi2());

  //Triplets  Loop
  cout<<"Number Of Triplets="<<Cand3Mu->size()<<endl;
  TripletCollectionSize = Cand3Mu->size() ;

  for(edm::View<reco::CompositeCandidate>::const_iterator TauIt=Cand3Mu->begin(); TauIt!=Cand3Mu->end(); ++TauIt){

    //    if (!(TauIt->vertexChi2() < 20)) continue ;
    //    TripletsCounter.push_back(1);

    //Daughter Kinematics at reco+gen level
    const Candidate * c1 = TauIt->daughter(0)->masterClone().get();
    const pat::Muon *mu1 = dynamic_cast<const pat::Muon *>(c1);

    const Candidate * c2 = TauIt->daughter(1)->masterClone().get();
    const pat::Muon *mu2 = dynamic_cast<const pat::Muon *>(c2);

    const Candidate * c3 = TauIt->daughter(2)->masterClone().get();
    const pat::Muon *mu3 = dynamic_cast<const pat::Muon *>(c3);

    Mu1_Pt.push_back(mu1->pt());  
    Mu1_Eta.push_back(mu1->eta());  
    Mu1_Phi.push_back(mu1->phi());  

    Mu2_Pt.push_back(mu2->pt());  
    Mu2_Eta.push_back(mu2->eta());  
    Mu2_Phi.push_back(mu2->phi());  
    
    Mu3_Pt.push_back(mu3->pt());  
    Mu3_Eta.push_back(mu3->eta());  
    Mu3_Phi.push_back(mu3->phi());  

    //    cout<<"Reco mu1 pt="<<mu1->pt()<<" mu2 pt="<<mu2->pt()<<" mu3 pt="<<mu3->pt()<<endl;

    reco::Vertex TripletVtx = reco::Vertex(TauIt->vertex(), TauIt->vertexCovariance(), TauIt->vertexChi2(), TauIt->vertexNdof(), TauIt->numberOfDaughters() );
    //TransientVertex TransientTripletVtx = reco::Vertex(TauIt->vertex(), TauIt->vertexCovariance(), TauIt->vertexChi2(), TauIt->vertexNdof(), TauIt->numberOfDaughters() );

    if(isMc){    
      bool isMatch1=false; bool isMatch2=false; bool isMatch3=false;
    if( (mu1->simType() == reco::MatchedMuonFromHeavyFlavour) && (fabs(mu1->simMotherPdgId()) == 15) ){
      isMatch1=true;
      Mu1_SimPt.push_back(mu1->simPt());
      Mu1_SimEta.push_back(mu1->simEta());
      Mu1_SimPhi.push_back(mu1->simPhi());
    }
    if( (mu2->simType() == reco::MatchedMuonFromHeavyFlavour) && (fabs(mu2->simMotherPdgId()) == 15) ){
      isMatch2=true;
      Mu2_SimPt.push_back(mu2->simPt());
      Mu2_SimEta.push_back(mu2->simEta());
      Mu2_SimPhi.push_back(mu2->simPhi());
    }
    if( (mu3->simType() == reco::MatchedMuonFromHeavyFlavour) && (fabs(mu3->simMotherPdgId()) == 15) ){
      isMatch3=true;
      Mu3_SimPt.push_back(mu3->simPt());
      Mu3_SimEta.push_back(mu3->simEta());
      Mu3_SimPhi.push_back(mu3->simPhi());
    }
    if( isMatch1 && isMatch2 && isMatch3){
      cout<<"Triplet Mass:"<<TauIt->mass()<<" pt="<<TauIt->pt()<<" vtx.x="<<TauIt->vx()<<" vtx x="<<TripletVtx.x()<<" chi2="<<TauIt->vertexChi2()<<" ndof="<<TauIt->vertexNdof()<<endl;
    }
    
    //GenVtx vars to be added
    }
    
    //Triplets Vars


    TripletVtx_x.push_back(TauIt->vx());
    TripletVtx_y.push_back(TauIt->vy());
    TripletVtx_z.push_back(TauIt->vz());

    TripletVtx_Chi2.push_back(TauIt->vertexChi2());
    TripletVtx_NDOF.push_back(TauIt->vertexNdof());

    Triplet_Mass.push_back(TauIt->mass());
    Triplet_Pt.push_back(TauIt->pt());
    Triplet_Eta.push_back(TauIt->eta());
    Triplet_Phi.push_back(TauIt->phi());
    Triplet_Charge.push_back(TauIt->charge());
    //Matrix covariance to be added!!!!

    //Refitted Vars
    //vector < TransientTrack > ttrks = TripletVtx.refittedTracks();

    ////
    TrackRef trk1, trk2, trk3;
    if (mu1->isGlobalMuon()) { trk1 = mu1->get<TrackRef,reco::CombinedMuonTag>();}
    else { trk1 = mu1->get<TrackRef>();}
    if (mu2->isGlobalMuon()) { trk2 = mu2->get<TrackRef,reco::CombinedMuonTag>();}
    else{ trk2 = mu2->get<TrackRef>();}
    if (mu3->isGlobalMuon()) { trk3 = mu3->get<TrackRef,reco::CombinedMuonTag>();}
    else{  trk3 = mu3->get<TrackRef>();}
    //cout<<" trk1 id="<<trk1.id()<<" tr2:"<<trk2.id()<<" trk3="<<trk3.id()<<endl;                                                                                        
    const reco::TransientTrack transientTrack1=theTransientTrackBuilder_->build( trk1 );
    const reco::TransientTrack transientTrack2=theTransientTrackBuilder_->build( trk2 );
    const reco::TransientTrack transientTrack3=theTransientTrackBuilder_->build( trk3 );
    reco::Track Track1 =transientTrack1.track();
    reco::Track Track2 =transientTrack2.track();
    reco::Track Track3 =transientTrack3.track();
    reco::Track* TrackRef1=&Track1;
    reco::Track* TrackRef2=&Track2;
    reco::Track* TrackRef3=&Track3;
    vector<reco::Track*> SVTrackRef;
    SVTrackRef.push_back(TrackRef1);
    SVTrackRef.push_back(TrackRef2);
    SVTrackRef.push_back(TrackRef3);
    removeTracks(pvTrackMap_refit,  SVTrackRef);

    std::vector<reco::TransientTrack> pvTracks_refit;
    for ( TransientTrackMap::iterator pvTrack = pvTrackMap_refit.begin();  pvTrack != pvTrackMap_refit.end(); ++pvTrack ) {
      pvTracks_refit.push_back(pvTrack->second);}
    cout<<" PV Tracks after refit="<<pvTracks_refit.size()<<endl;
    /*for(uint i=0; i<pvTracks_refit.size(); i++){                                                                                                                    
      TrackRef tr = TrackRef(pvTracks_refit, i);                                                                                                                      
      //reco::Track pvTr=pvTracks_refit.at(i).track();                                                                                                                
      //TrackRef pvTrRef = pvTr.get<TrackRef>();                                                                                                                      
      cout<<i<<"PV track ID="<<tr.id()<<endl;                                                                                                                         
      }*/

    if(pvTracks_refit.size() >0){
    KalmanVertexFitter PV_fitter (true);
    TransientVertex PVertex = PV_fitter.vertex(pvTracks_refit);
    //CachingVertex<5> fittedVertex = vertexFitter.vertex(tracksToVertex);                                                                                            
    GlobalPoint PVertexPos  (PVertex.position());
    GlobalPoint SVertexPos  (TripletVtx.x(), TripletVtx.y(), TripletVtx.z());
    double FlightDist = TMath::Sqrt( pow(( PVertexPos.x() -SVertexPos.x()),2)+ pow(( PVertexPos.y() -SVertexPos.y()),2) + pow(( PVertexPos.z() -SVertexPos.z()),2));

    VertexDistance3D vertTool;
    VertexState PVstate(PVertex.position(),PVertex.positionError());
    //VertexState SVstate(SVertexPos,TripletVtx.position());
    double distance = vertTool.distance(PVstate, TripletVtx).value();
    double dist_err = vertTool.distance(PVstate, TripletVtx).error();
    double dist_sign =vertTool.distance(PVstate, TripletVtx).significance();
    double chi2 = vertTool.compatibility(PVstate, TripletVtx);

    ////

    RefittedPV_x.push_back(PVertexPos.x());
    RefittedPV_y.push_back(PVertexPos.y());
    RefittedPV_z.push_back(PVertexPos.z());
    RefittedPV_NTracks.push_back(pvTracks_refit.size());
    //RefittedPV_Chi2.push_back(PVertex.);

    FlightDistPVSV.push_back(distance);
    FlightDistPVSV_Err.push_back(dist_err);
    FlightDistPVSV_Significance.push_back(dist_sign);
    FlightDistPVSV_chi2.push_back(chi2);
    }
  }

   cout<<"***Number of Muons="<<muons->size()<<endl; uint k=0;


   std::vector<int> MuFilter;
   vector<pat::Muon>    MyMu, MyMu2, SyncMu;
   double AllMuPt =0;
   MuonCollectionSize = muons->size();

   for(edm::View<pat::Muon>::const_iterator mu=muons->begin(); mu!=muons->end(), k<muons->size(); ++mu, ++k){
    
     /*hMuonPt->Fill(mu->pt());
    hMuonP->Fill(mu->p());
    hMuonEta->Fill(mu->eta());
    hMuonNumberOfValidHits->Fill(mu->numberOfValidHits());
    hMuonDB->Fill(mu->dB());
    hMuonTime->Fill(mu->time().timeAtIpInOut);
    hMuonTimeErr->Fill(mu->time().timeAtIpInOutErr);
    hMuonTimeVsP->Fill(mu->p(), mu->time().timeAtIpInOut);    
    
    if (mu->isSoftMuon(PV)) {  hSoftMuonPt->Fill(mu->pt());     hSoftMuonEta->Fill(mu->eta());  }
    if (mu->isTrackerMuon()){  hTrackerMuonPt->Fill(mu->pt());  hTrackerMuonEta->Fill(mu->eta());}
    if (mu->isGlobalMuon()) {  hGlobalMuonPt->Fill(mu->pt());   hGlobalMuonEta->Fill(mu->eta());}
    if (mu->isLooseMuon())  {  hLooseMuonPt->Fill(mu->pt());    hLooseMuonEta->Fill(mu->eta());}

    if(!(mu->isSoftMuon(PV) || mu->isTrackerMuon() || mu->isGlobalMuon() || mu->isLooseMuon())) continue;
    cout<<"Reco loose mu pt="<<mu->pt()<<endl;
    
    if((mu->pt() > 2) &&  (fabs(mu->eta()) < 2.4) && (mu->isPFMuon()) && (mu->isGlobalMuon()) ){
      SyncMu.push_back(*mu);
      AllMuPt +=mu->pt();
      }  

    
    if (!(mu->track().isNonnull()))continue;
    if (!(mu->track()->quality(reco::TrackBase::highPurity))) continue;
    if ((!(mu->track()->hitPattern().numberOfValidPixelHits()>0))  ) continue;
    if (!(mu->pt() > 0.5) ) continue;
    */
    

    //mu.muonBestTrack()->dz(PV.position()), mu.isTightMuon(PV));
    //std::cout<<"RecoMu pt="<<mu->pt()<<" eta="<<mu->eta()<<" phi="<<mu->phi()<<" segm compatibility="<<mu->segmentCompatibility()<<" trkRelChi2="<<mu->combinedQuality().trkRelChi2<<"  vz="<<mu->vz()<<" mu charge="<<mu->charge()<<" mu ECAL energy="<<mu->calEnergy().em<<std::endl;
    //" SymType="<<mu->simExtType()<<" SimPt="<<mu->simPt()<<
    //cout<<" pixel hits="<<mu->track()->hitPattern().numberOfValidPixelHits()<<" pt="<<mu->track()->pt()<<" mu ECAL energy="<<mu->calEnergy().em<<" mu.vx="<<mu->vx()<<endl;
    

    MuFilter.push_back(1);  
    MyMu.push_back(*mu);

    //Basic Kinematics
    MuonPt.push_back(mu->pt());
    MuonEta.push_back(mu->eta());
    MuonPhi.push_back(mu->phi());
    MuonEnergy.push_back(mu->energy());
    MuonCharge.push_back(mu->charge());

    //Vtx position
    Muon_vx.push_back(mu->vx());
    Muon_vy.push_back(mu->vy());
    Muon_vz.push_back(mu->vz());
    
    //MuonID
    Muon_isGlobal.push_back(mu->isGlobalMuon());
    Muon_isTracker.push_back(mu->isTrackerMuon());
    Muon_isSoft.push_back(mu->isSoftMuon(PV));
    Muon_isLoose.push_back(mu->isLooseMuon());
    Muon_isPF.push_back(mu->isPFMuon());
    Muon_isRPCMuon.push_back(mu->isRPCMuon());
    Muon_isStandAloneMuon.push_back(mu->isStandAloneMuon());
    Muon_isTrackerMuon.push_back(mu->isTrackerMuon());
    Muon_isCaloMuon.push_back(mu->isCaloMuon());
    Muon_isQualityValid.push_back(mu->isQualityValid());
    Muon_isTimeValid.push_back(mu->isTimeValid());
    Muon_isIsolationValid.push_back(mu->isIsolationValid());
    Muon_numberOfMatchedStations.push_back(mu->numberOfMatchedStations());
    Muon_numberOfMatches.push_back(mu->numberOfMatches(reco::Muon::SegmentArbitration));

    MuonChi2P.push_back(mu->combinedQuality().chi2LocalMomentum);
    MuonChi2LocalPosition.push_back(mu->combinedQuality().chi2LocalPosition);
    MuonGlbTrackProbability.push_back(mu->combinedQuality().glbTrackProbability);
    MuonTrkRelChi2.push_back(mu->combinedQuality().trkRelChi2);
    MuonTrkKink.push_back(mu->combinedQuality().trkKink);

    Muon_timeAtIpInOut.push_back(mu->time().timeAtIpInOut);
    Muon_timeAtIpInOutErr.push_back(mu->time().timeAtIpInOutErr);

    if (mu->isGlobalMuon()) {
      Muon_GLnormChi2.push_back(mu->globalTrack()->normalizedChi2());
      Muon_GLhitPattern_numberOfValidMuonHits.push_back(mu->globalTrack()->hitPattern().numberOfValidMuonHits());
    }else
      {
      Muon_GLnormChi2.push_back(-999);
      Muon_GLhitPattern_numberOfValidMuonHits.push_back(-999);
      }

    if (mu->innerTrack().isNonnull()){ 
    Muon_trackerLayersWithMeasurement.push_back(mu->innerTrack()->hitPattern().trackerLayersWithMeasurement());
    Muon_Numberofvalidpixelhits.push_back(mu->innerTrack()->hitPattern().numberOfValidPixelHits());
    Muon_innerTrack_p.push_back(mu->innerTrack()->p());
    Muon_innerTrack_eta.push_back(mu->innerTrack()->eta());
    Muon_innerTrack_phi.push_back(mu->innerTrack()->phi());
    Muon_innerTrack_normalizedChi2.push_back(mu->innerTrack()->normalizedChi2());
    }else
      {
      Muon_trackerLayersWithMeasurement.push_back(-999);
      Muon_Numberofvalidpixelhits.push_back(-999);
      Muon_innerTrack_p.push_back(-999);
      Muon_innerTrack_eta.push_back(-999);
      Muon_innerTrack_phi.push_back(-999);
      Muon_innerTrack_normalizedChi2.push_back(-999);
    }
    if (mu->outerTrack().isNonnull()){
    Muon_outerTrack_p.push_back(mu->outerTrack()->p());
    Muon_outerTrack_eta.push_back(mu->outerTrack()->eta());
    Muon_outerTrack_phi.push_back(mu->outerTrack()->phi());
    Muon_outerTrack_normalizedChi2.push_back(mu->outerTrack()->normalizedChi2());
    Muon_outerTrack_muonStationsWithValidHits.push_back(mu->outerTrack()->hitPattern().muonStationsWithValidHits());
    }else{
      Muon_outerTrack_p.push_back(-999);
      Muon_outerTrack_eta.push_back(-999);
      Muon_outerTrack_phi.push_back(-999);
      Muon_outerTrack_normalizedChi2.push_back(-999);
      Muon_outerTrack_muonStationsWithValidHits.push_back(-999);
    }
    if (mu->innerTrack().isNonnull() && mu->outerTrack().isNonnull()){
    Muon_QInnerOuter.push_back(mu->outerTrack()->charge()*mu->innerTrack()->charge());
    }else{
      Muon_QInnerOuter.push_back(-999);
    }


    Muon_combinedQuality_updatedSta.push_back(mu->combinedQuality().updatedSta);
    Muon_combinedQuality_trkKink.push_back(mu->combinedQuality().trkKink);
    Muon_combinedQuality_glbKink.push_back(mu->combinedQuality().glbKink);
    Muon_combinedQuality_trkRelChi2.push_back(mu->combinedQuality().trkRelChi2);
    Muon_combinedQuality_staRelChi2.push_back(mu->combinedQuality().staRelChi2);
    Muon_combinedQuality_chi2LocalPosition.push_back(mu->combinedQuality().chi2LocalPosition);
    Muon_combinedQuality_chi2LocalMomentum.push_back(mu->combinedQuality().chi2LocalMomentum);
    Muon_combinedQuality_localDistance.push_back(mu->combinedQuality().localDistance);
    Muon_combinedQuality_globalDeltaEtaPhi.push_back(mu->combinedQuality().globalDeltaEtaPhi);
    Muon_combinedQuality_tightMatch.push_back(mu->combinedQuality().tightMatch);
    Muon_combinedQuality_glbTrackProbability.push_back(mu->combinedQuality().glbTrackProbability);
        
    Muon_calEnergy_em.push_back(mu->calEnergy().em);
    Muon_calEnergy_emS9.push_back(mu->calEnergy().emS9);
    Muon_calEnergy_emS25.push_back(mu->calEnergy().emS25);
    Muon_calEnergy_had.push_back(mu->calEnergy().had);
    Muon_calEnergy_hadS9.push_back(mu->calEnergy().hadS9);
        
    Muon_segmentCompatibility.push_back(muon::segmentCompatibility(*mu));
    Muon_caloCompatibility.push_back(muon::caloCompatibility(*mu));
        
    Muon_ptErrOverPt.push_back( (mu->muonBestTrack()->ptError()/mu->muonBestTrack()->pt()) );

   }

   ////Synch Tree//////
   /*   double maxPt =0; double maxPhi=0, maxEta=0; vector<pat::Muon> SyncSortedMu ;
   double maxTrPt =0; double maxTrPhi=0, maxTrEta=0; vector<reco::Track> SyncSortedTr ;
   for(uint i=0; i<SyncMu.size();i++){
     if(SyncMu.at(i).pt() > maxPt){
       maxPt  = SyncMu.at(i).pt();
       maxPhi = SyncMu.at(i).phi();
       maxEta = SyncMu.at(i).eta();
       SyncSortedMu.push_back(SyncMu.at(i));
     }
   }

   allmuons_pt.push_back(AllMuPt);

   leadmuon_pt.push_back(maxPt);
   leadmuon_eta.push_back(maxEta);
   leadmuon_phi.push_back(maxPhi);
   nmuons = SyncMu.size();
   nprimevtxs =vertices->size();
   
   edm::View<reco::Track>::const_iterator trIt  = trackCollection->begin();
   edm::View<reco::Track>::const_iterator trEnd = trackCollection->end();

   double AllTrPt=0;
   for (; trIt != trEnd; ++trIt)
     {
      
       const reco::Track track = (*trIt);
       if(  (track.pt()>1) && (fabs(track.eta())<2.4) && (track.hitPattern().trackerLayersWithMeasurement()>5) && (track.hitPattern().pixelLayersWithMeasurement()>1)  ){
	 AllTrPt +=trIt->pt();
	 if(track.pt() > maxTrPt){
	   maxTrPt = track.pt();
	   maxTrEta = track.eta();
	   maxTrPhi= track.phi();
	 }
       }
     }

   alltracks_pt.push_back(AllTrPt);
   leadtrack_pt.push_back(maxTrPt);
   leadtrack_eta.push_back(maxTrEta);
   leadtrack_phi.push_back(maxTrPhi);

   evt   = iEvent.id().event();
   run = iEvent.id().run();
   lumi = iEvent.luminosityBlock();
   */
   ///////SyncTree

   //hGoodMuSize->Fill(MyMu.size());    
   //if (MuFilter.size()>0)   hEvents_MuFilter->Fill(1);
   //cout<<"***Number of Muons after preliminary requ="<<MyMu.size()<<endl;

   ///DiMuon Selections
   /*
   for (uint mm=0; mm<MyMu.size();mm++){
     for (uint kk=mm+1; kk<MyMu.size();kk++){
       float DiMuonDz = fabs(MyMu.at(mm).vz()-MyMu.at(kk).vz());
       float DiMuonDR = dR(MyMu.at(mm).eta(), MyMu.at(kk).eta(), MyMu.at(mm).phi(), MyMu.at(kk).phi());
       //cout<<mm<<","<<kk<<"  Dz="<<DiMuonDz<<" dr="<<DiMuonDR<<endl;
       hDiMuonDz->Fill(DiMuonDz);
       hDiMuonDR->Fill(DiMuonDR);
       //if (DiMuonDz<1) MyMu2.push_back(MyMu.at(mm));   //Questo taglio fa morire il 20% del segnale--> Guardare la distrubuzione di Dz nelle genParticle o nelle simTrack
     }
     MyMu2.push_back(MyMu.at(mm)); //skippiamo il taglio sul DZ
   }

   hMuonSize_GoodDiMu->Fill(MyMu2.size());
  cout<<"***Number of Muons after Dz="<<MyMu2.size()<<endl;
  if ( MyMu2.size()>0) hEvents_DiMuonDz->Fill(1);
   */
  //Selezione su TriMuonMass + vertex

  //std::vector<reco::TransientTrack> Ttracks1, Ttracks2, Ttracks3;
  //vector<TransientVertex> MyVertexC;
  // vector<RefCountedKinematicVertex>  VertexContainer;
  //  std::vector<reco::TransientTrack> Ttracks;
  //vector<TransientVertex> KVContainer, KV_ValidContainer;

  //vector<RefCountedKinematicTree> MyKinVtxTree;
  //Creating a KinematicParticleFactory                                                                                                                            
  /*  KinematicParticleFactoryFromTransientTrack pFactory;
  vector<RefCountedKinematicParticle> muonParticles;
  //initial chi2 and ndf before kinematic fits.                                                                                                                    
  ParticleMass muon_mass = 0.10565837; //pdg mass
  float muon_sigma = muon_mass*1.e-6;
  float chi = 0.;
  float ndf = 0.;  
  */

  //TLorentzVector Mu1(0.,0.,0.,0.), Mu2(0.,0.,0.,0.), Mu3(0.,0.,0.,0.),  ThreeMuons; 
   //   vector<vector<int>> My3Mu;
   /* struct My3Mu{
     int Mu1Idx;
     int Mu2Idx;
     int Mu3Idx;
     } ;*/
   //   My3Mu My3MuC;
   //vector<My3Mu> My3MuV;

  /*
  vector<int> VtxNDofSel, MassSel, ValidVtx,TauCharge;
   if(MyMu2.size()>2){
     hEvents_3Mu->Fill(1);
     for (uint a=0; a<MyMu2.size();a++){
       const reco::TransientTrack transientTrack1=theTransientTrackBuilder_->build( MyMu2.at(a).track() ); 
       Ttracks1.push_back(transientTrack1);
       Mu1.SetPtEtaPhiE(MyMu2.at(a).pt(), MyMu2.at(a).eta(), MyMu2.at(a).phi(), MyMu2.at(a).energy());
       //transientTrack.setBeamSpot(vertexBeamSpot);
       //FreeTrajectoryState mu1State = transientTrack1.impactPointTSCP().theState();
       //muonParticles.push_back(pFactory.particle(transientTrack1,muon_mass,chi,ndf,muon_sigma));


       for(uint b=a+1; b<MyMu2.size();b++){
	 const reco::TransientTrack transientTrack2 = theTransientTrackBuilder_->build( MyMu2.at(b).track() ); 
	 Ttracks2.push_back(transientTrack2);
	 Mu2.SetPtEtaPhiE(MyMu2.at(b).pt(), MyMu2.at(b).eta(), MyMu2.at(b).phi(), MyMu2.at(b).energy());
	 //	 FreeTrajectoryState mu2State = transientTrack2.impactPointTSCP().theState();
	 //if( !transientTrack2.impactPointTSCP().isValid()) continue;
	 //muonParticles.push_back(pFactory.particle(transientTrack2,muon_mass,chi,ndf,muon_sigma));
	 
	 for(uint c=b+1; c<MyMu2.size();c++){
	   const reco::TransientTrack transientTrack3 = theTransientTrackBuilder_->build( MyMu2.at(c).track() ); 
	   Ttracks3.push_back(transientTrack3);
	   Mu3.SetPtEtaPhiE(MyMu2.at(c).pt(), MyMu2.at(c).eta(), MyMu2.at(c).phi(), MyMu2.at(c).energy());
	   // FreeTrajectoryState mu3State = transientTrack3.impactPointTSCP().theState();
	   //if( !transientTrack3.impactPointTSCP().isValid()) continue;
	   ThreeMuons = Mu1+Mu2+Mu3;


	   //Creating a KinematicParticleFactory                                                                                                                            
	   KinematicParticleFactoryFromTransientTrack pFactory;
	   vector<RefCountedKinematicParticle> muonParticles;
	   //initial chi2 and ndf before kinematic fits.                                                                                                                    
	   ParticleMass muon_mass = 0.10565837; //pdg mass
	   //float muon_sigma = muon_mass*1.e-6;
	   float muon_sigma = muon_mass*1.e-3;
	   float chi = 0.;
	   float ndf = 0.;  
	   
	   muonParticles.push_back( pFactory.particle(transientTrack1,muon_mass,chi,ndf,muon_sigma));
	   muonParticles.push_back( pFactory.particle(transientTrack2,muon_mass,chi,ndf,muon_sigma));
	   muonParticles.push_back( pFactory.particle(transientTrack3,muon_mass,chi,ndf,muon_sigma));
	   
	   float ThreeMuonCharge = MyMu2.at(a).charge()+MyMu2.at(b).charge()+MyMu2.at(c).charge();
	   hThreeMuonCharge->Fill(ThreeMuonCharge);

	   if (!(fabs(ThreeMuonCharge)==1)) continue;
	   TauCharge.push_back(1);
	   hThreeMuonInvMass->Fill(ThreeMuons.M());

	   //  cout<<a<<","<<b<<","<<c<<" 3Muon mass="<<ThreeMuons.M()<<" charge="<<ThreeMuonCharge<<endl;
	   std::vector<reco::TransientTrack> Ttracks;
	   Ttracks.push_back(transientTrack1);
	   Ttracks.push_back(transientTrack2);
	   Ttracks.push_back(transientTrack3);
	   //cout<<" Tracks used to refit the SV with Kalman="<<Ttracks.size()<<" with Kin="<<muonParticles.size()<<endl;
	   KalmanVertexFitter KVfitter (true);
	   TransientVertex KVertex = KVfitter.vertex(Ttracks);
	   KVContainer.push_back(KVertex);
	   if (KVertex.isValid()) KV_ValidContainer.push_back(KVertex);

	   KinematicParticleVertexFitter fitter;
	   RefCountedKinematicTree TauVertexFitTree;
	   //TauVertexFitTree = fitter.fit(muonParticles);
	   try {
	     TauVertexFitTree = fitter.fit(muonParticles);
	   }catch (...) {
	     std::cout<<" Exception caught ... continuing 2 "<<std::endl;
	     continue;
	   }
	   if (!(TauVertexFitTree->isValid())) continue;
	   ValidVtx.push_back(1);
	   TauVertexFitTree->movePointerToTheTop();

	   RefCountedKinematicParticle Tau_vFit_noMC = TauVertexFitTree->currentParticle();
	   RefCountedKinematicVertex Tau_vFit_vertex_noMC = TauVertexFitTree->currentDecayVertex();
	   //cout<<"3muon vtx: Kin Fit mass="<<Tau_vFit_noMC->currentState().mass()<<" mass"<<ThreeMuons.M()<<endl;


	   hTau_vFit_chi2->Fill(Tau_vFit_vertex_noMC->chiSquared());
	   if(Tau_vFit_vertex_noMC->chiSquared()>50.) continue;
	   VtxNDofSel.push_back(1);
	   h3MuonMassVtxFit->Fill(Tau_vFit_noMC->currentState().mass());
	   if(Tau_vFit_noMC->currentState().mass()>1. && Tau_vFit_noMC->currentState().mass()<4.) {
	     MassSel.push_back(1); 
	     VertexContainer.push_back(Tau_vFit_vertex_noMC);
	     MyKinVtxTree.push_back(TauVertexFitTree);
	   }
	   
	   
	 }

	 //Ttracks.clear();   
	 //muonParticles.clear();

       }
     }

   }//Eventi con almeno 3muoni
   if ( TauCharge.size()>0)  hEvent_TauCharge->Fill(1);
   if ( ValidVtx.size()>0)   hEvent_Valid3MuonVtx->Fill(1);
   if ( VtxNDofSel.size()>0) hEvent_3MuonVtx->Fill(1);
   if ( MassSel.size()>0)    hEvent_3MuonMass->Fill(1);

   //cout<<" #of 3Mu cand="<<My3Mu.size()<<endl;
   //cout<<" Size of Transient Tracks="<<Ttracks.size()<<endl;
   hVtxSize->Fill(VertexContainer.size());
   cout<<" Number of recontructed vertices: KinFit="<<VertexContainer.size()<<" KalmanFit="<<KVContainer.size()<<" valid="<<KV_ValidContainer.size()<<endl;
  */
   ///////////////////////////////////////////metto le tracce del PV in una mappa/////////////////////////////////////////////////////////////                     
   
   //std::vector<reco::TransientTrack> pvTracks_original;
  //std::vector<reco::Track*> tauTrk;
   /*
   TransientTrackMap pvTrackMap_refit;
   const reco::Vertex* eventVertex;

   if ((!vertices->empty()) && (*vertices)[0].isValid()){
     cout<<" PV is Valid: "<< (*vertices)[0].isValid()<<endl;
     for ( reco::Vertex::trackRef_iterator pvTrack = (*vertices)[0].tracks_begin(); pvTrack != (*vertices)[0].tracks_end(); ++pvTrack ) {
   //for ( reco::Vertex::trackRef_iterator pvTrack = (*vertices).tracks_begin(); pvTrack != (*vertices)[0].tracks_end(); ++pvTrack ) {
       reco::TransientTrack pvTrack_transient = theTransientTrackBuilder_->build(pvTrack->get());
     pvTracks_original.push_back(pvTrack_transient);
     pvTrackMap_refit.insert(std::make_pair(pvTrack->get(), pvTrack_transient)); //a questo punto nella mappa entrambi gli elementi sono le tracce associate al PV  
   }
   }
   */

   //accedere alle tracce del vertice del tau
   
   //if(VertexContainer.size()>0){
  //if(! (MyKinVtxTree.at(0)->isEmpty()) ){
  //   std::vector<RefCountedKinematicParticle> daughters = MyKinVtxTree.at(0)->daughterParticles();
       //cout<<"number of tracks in the PV="<< pvTracks_original.size()<<" number of daughters="<<daughters.size()<<endl;
  //   cout<< " number of daughters="<<daughters.size()<<endl; 
       /*
       for (std::vector<RefCountedKinematicParticle>::const_iterator i = daughters.begin();	i != daughters.end(); ++i) {
	 const TransientTrackKinematicParticle * ttkp = dynamic_cast<const TransientTrackKinematicParticle * >(&(**i));
	 const reco::TransientTrack * svTrack_transient = dynamic_cast<const reco::TransientTrack*>(ttkp->initialTransientTrack()->basicTransientTrack());
	 reco::Track  tmp_svtrack = (svTrack_transient->track());
	 reco::Track* tmp_ptr= &tmp_svtrack;
	 tauTrk.push_back(tmp_ptr); 
       }*/
       //removeTracks(pvTrackMap_refit, tauTrk);
  // }
  //   }
      

    /*
   for (reco::Vertex::trackRef_iterator svit = VertexContainer[0]->tracks_begin(); svit!=VertexContainer[0].tracks_end(); ++svit){
     std::vector<RefCountedKinematicParticle> daughters = tree->daughterParticles();
     //cout<<" Tau Vtx position "<<svit.position()<<endl;
     reco::TransientTrack svTrack_transient = theTransientTrackBuilder_->build(svit->get());
     tauTrk.push_back(svTrack_transient);
     //svTrackMap_refit.insert(std::make_pair(svit->get(), svtrack_transient);
     }*/

   


     ////////////////definisco un vettore con le tracce del PV da cui ho tolto le tracce del SV///////////////////////////////////////////                            
     
   /*
   std::vector<reco::TransientTrack> pvTracks_refit;
     for ( TransientTrackMap::iterator pvTrack = pvTrackMap_refit.begin();  pvTrack != pvTrackMap_refit.end(); ++pvTrack ) {
       pvTracks_refit.push_back(pvTrack->second);}
     cout<<"number of tracks in the PV="<<pvTracks_refit.size()<<endl;
   */

  /*
 vector<TransientVertex> MyVertexC;
  cout<<" tracks to be fitted="<<Ttracks.size()<<endl;
  if( Ttracks.size()>2 ){
    AdaptiveVertexFitter  theFitter;
    //TransientVertex myVertex = theFitter.vertex(mytracks, vertexBeamSpot);  // if you want the beam constraint
    TransientVertex myVertex = theFitter.vertex(Ttracks);  // if you don't want the beam constraint
    // now you have a new vertex, can e.g. be compared with the original
    MyVertexC.push_back(myVertex);
    std::cout << "TEST   PV_z=" << PV.position().z() << " ADV_z=" << myVertex.position().z() << std::endl;
    std::cout << "TEST   PV_x=" << PV.position().x() << " ADV_x=" << myVertex.position().x() << std::endl;
    std::cout << "TEST   PV_y=" << PV.position().y() << " ADV_y=" << myVertex.position().y() << std::endl;
    float dxy=TMath::Sqrt(myVertex.position().x()*myVertex.position().x()+PV.position().y()*PV.position().y());
    cout<<" dxy="<<dxy<<" chi2="<<myVertex.totalChiSquared()<<" ndof="<<myVertex.degreesOfFreedom()<<" valid="<<myVertex.isValid()<<endl;
  }else{
    std::cout << "not enough tracks left" <<std::endl;

  }

  cout<<" Number of 3Muons vertices ="<<MyVertexC.size()<<endl;

  
    edm::Handle mets;
   iEvent.getByToken(metToken_, mets);
   const pat::MET &met = mets->front();
   printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
	  met.pt(), met.phi(), met.sumEt(),
	  met.genMET()->pt(),
	  met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));

   printf("\n");
   */



   //  SyncTree_->Fill();
  tree_->Fill();

  /*
  allmuons_pt.clear();
  alltracks_pt.clear();
  leadmuon_pt.clear();
  leadmuon_phi.clear();
  leadmuon_eta.clear(); 
  allmuons_pt.clear(); 
  leadtrack_pt.clear();
  leadtrack_eta.clear();
  leadtrack_phi.clear();
  evt= -999;
  run= -999;
  lumi= -999;
  nmuons = -999;  
  */
  MuonCollectionSize =0;
  MuonPt.clear();
  MuonEta.clear();
  MuonPhi.clear();
  MuonChi2P.clear();
  MuonChi2LocalPosition.clear();
  MuonGlbTrackProbability.clear();
  MuonTrkRelChi2.clear();
  MuonTrkKink.clear();

  MuonEnergy.clear();
  MuonCharge.clear();

  //Vtx position                                                                                                                                                                             
  Muon_vx.clear();
  Muon_vy.clear();
  Muon_vz.clear();

  //MuonID                                                                                                                                                                                   
  Muon_isGlobal.clear();
  Muon_isTracker.clear();
  Muon_isSoft.clear();
  Muon_isLoose.clear();
  Muon_isPF.clear();
  Muon_isRPCMuon.clear();
  Muon_isStandAloneMuon.clear();
  Muon_isTrackerMuon.clear();
  Muon_isCaloMuon.clear();
  Muon_isQualityValid.clear();
  Muon_isTimeValid.clear();
  Muon_isIsolationValid.clear();
  Muon_numberOfMatchedStations.clear();
  Muon_numberOfMatches.clear();

  //MuonTime
  Muon_timeAtIpInOut.clear();
  Muon_timeAtIpInOutErr.clear();

  //Muon inner + outer track
  Muon_GLnormChi2.clear();
  Muon_GLhitPattern_numberOfValidMuonHits.clear();

  Muon_trackerLayersWithMeasurement.clear();
  Muon_Numberofvalidpixelhits.clear();

  Muon_outerTrack_p.clear();
  Muon_outerTrack_eta.clear();
  Muon_outerTrack_phi.clear();
  Muon_outerTrack_normalizedChi2.clear();
  Muon_outerTrack_muonStationsWithValidHits.clear();
  Muon_innerTrack_p.clear();
  Muon_innerTrack_eta.clear();
  Muon_innerTrack_phi.clear();
  Muon_innerTrack_normalizedChi2.clear();
  Muon_QInnerOuter.clear();


  Muon_combinedQuality_updatedSta.clear();
  Muon_combinedQuality_trkKink.clear();
  Muon_combinedQuality_glbKink.clear();
  Muon_combinedQuality_trkRelChi2.clear();
  Muon_combinedQuality_staRelChi2.clear();
  Muon_combinedQuality_chi2LocalPosition.clear();
  Muon_combinedQuality_chi2LocalMomentum.clear();
  Muon_combinedQuality_localDistance.clear();
  Muon_combinedQuality_globalDeltaEtaPhi.clear();
  Muon_combinedQuality_tightMatch.clear();
  Muon_combinedQuality_glbTrackProbability.clear();

  Muon_calEnergy_em.clear();
  Muon_calEnergy_emS9.clear();
  Muon_calEnergy_emS25.clear();
  Muon_calEnergy_had.clear();
  Muon_calEnergy_hadS9.clear();

  Muon_segmentCompatibility.clear();
  Muon_caloCompatibility.clear();

  Muon_ptErrOverPt.clear();

  PV_x=-99;
  PV_y=-99;
  PV_z=-99;
  PV_NTracks=-99; 

  Mu1_Pt.clear();
  Mu1_Eta.clear();
  Mu1_Phi.clear();

  Mu2_Pt.clear();
  Mu2_Eta.clear();
  Mu2_Phi.clear();

  Mu3_Pt.clear();
  Mu3_Eta.clear();
  Mu3_Phi.clear();


  Mu1_SimPt.clear();
  Mu1_SimEta.clear();
  Mu1_SimPhi.clear();
  
  Mu2_SimPt.clear();
  Mu2_SimEta.clear();
  Mu2_SimPhi.clear();
  
  Mu3_SimPt.clear();
  Mu3_SimEta.clear();
  Mu3_SimPhi.clear();
  TripletCollectionSize = -99;

  TripletVtx_x.clear();
  TripletVtx_y.clear();
  TripletVtx_z.clear();

  TripletVtx_Chi2.clear();
  TripletVtx_NDOF.clear();

  Triplet_Mass.clear();
  Triplet_Pt.clear();
  Triplet_Eta.clear();
  Triplet_Phi.clear();
  Triplet_Charge.clear();

  RefittedPV_x.clear();
  RefittedPV_y.clear();
  RefittedPV_z.clear();
  RefittedPV_NTracks.clear();
  //RefittedPV_Chi2.push_back(PVertex.);                                                                                                                                                   

  FlightDistPVSV.clear();
  FlightDistPVSV_Err.clear();
  FlightDistPVSV_Significance.clear();
  FlightDistPVSV_chi2.clear();
  PVCollection_Size =0;
}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAna2017Tree::beginJob()
{
  hEvents = fs->make<TH1F>("hEvents","hEvents",10,0,10);
  /*
  hEvents_MuFilter = fs->make<TH1F>("hEvents_MuFilter","hEvents_MuFilter",10,0,10);
  hEvents_DiMuonDz = fs->make<TH1F>("hEvents_DiMuonDz","hEvents_DiMuonDz",10,0,10);
  hEvents_3Mu = fs->make<TH1F>("hEvents_3Mu","hEvents_3Mu",10,0,10);
  hEvent_TauCharge = fs->make<TH1F>("hEvent_TauCharge", "hEvent_TauCharge",10,0,10);
  hEvent_Valid3MuonVtx= fs->make<TH1F>("hEvent_Valid3MuonVtx", "hEvent_Valid3MuonVtx", 10, 0, 10);
  hEvent_3MuonVtx = fs->make<TH1F>("hEvents_3MuonVtx","hEvents_3MuonVtx",10,0,10);
  hEvent_3MuonMass = fs->make<TH1F>("hEvents_3MuonMass","hEvents_3MuonMass",10,0,10);
  

  hGenMuonPt = fs->make<TH1F>("hGenMuonPt","hGenMuonPt",200,0,100);
  hGenMuonEta = fs->make<TH1F>("hGenMuonEta","hGenMuonEta",60,-3,3);
  hGenTauPt = fs->make<TH1F>("hGenTauPt","hGenTauPt",200,0,100);
  hMuonPt = fs->make<TH1F>("hMuonPt","hMuonPt",200,0,100);

  hGlobalMuonPt = fs->make<TH1F>("hGlobalMuonPt","hMuonPt",200,0,100);
  hSoftMuonPt = fs->make<TH1F>("hSoftMuonPt","hMuonPt",200,0,100);
  hTrackerMuonPt = fs->make<TH1F>("hTrackerMuonPt","Tracker Muon Pt",200,0,100);
  hLooseMuonPt = fs->make<TH1F>("hLooseMuonPt","hMuonPt",200,0,100);

  hMuonP = fs->make<TH1F>("hMuonP","hMuonPt",200,0,100);
  hMuonEta = fs->make<TH1F>("hMuonEta","hMuonEta",60,-3,3);

  hTrackerMuonEta = fs->make<TH1F>("hTrackerMuonEta","Muon Eta",60,-3,3);
  hSoftMuonEta = fs->make<TH1F>("hSoftMuonEta","Muon Eta",60,-3,3);
  hGlobalMuonEta = fs->make<TH1F>("hGlobalMuonEta","Muon Eta",60,-3,3);
  hLooseMuonEta = fs->make<TH1F>("hLooseMuonEta","Muon Eta",60,-3,3);

  hGoodMuSize = fs->make<TH1F>("hGoodMuSize","Number of selected muons per event", 20,0,20);
  hMuonSize_GoodDiMu = fs->make<TH1F>("hMuonSize_GoodDiMu","Number of selected muons per event", 20,0,20);
  hMuonNumberOfValidHits = fs->make<TH1F>("hMuonNumberOfValidHits","hMuonNumberOfValidHits",50,0,50);
  hMuonDB= fs->make<TH1F>("hMuonDB","Muon Impact Parameter",100,-5,5);
  hMuonTime = fs->make<TH1F>("hMuonTime","Time of arrival at the IP",100,-0.005,0.005);
  hMuonTimeErr = fs->make<TH1F>("hMuonTimeErr","Error on Time of arrival at the IP",100,-0.005,0.005);
  hMuonTimeVsP = fs->make<TH2F>("hMuonTimeVsP","muon time vs p",100,-0.005,0.005, 500, 0, 50);
  hDiMuonDz = fs->make<TH1F>("hDiMuonDz","Distance in z between 2 muon pairs",200,0,10);
  hDiMuonDR = fs->make<TH1F>("hDiMuonDR","Distance in Eta,Phi between 2 muon pairs",200,0,10);
  hThreeMuonInvMass = fs->make<TH1F>("hThreeMuonInvMass","Invariant mass of the 3Muons Candidate",400,0,10);
  hThreeMuonCharge = fs->make<TH1F>("hThreeMuonCharge", "Charge of the 3Muons Candidate", 6, -3.25, 3.25);
  hVtxSize = fs->make<TH1F>("hVtxSize", "Number of ReFitted 3Muon Vtx",200, 0, 200);


  hTau_vFit_chi2 =  fs->make<TH1F>("hTau_vFit_chi2", "hTau_vFit_chi2", 100, 0, 100);
  h3MuonMassVtxFit = fs->make<TH1F>("h3MuonMassVtxFit", "h3MuonMassVtxFit", 100, 0, 50);
  */

  tree_ = fs->make<TTree>("ntuple","LFVTau ntuple");
 
  tree_->Branch("MuonCollectionSize",&MuonCollectionSize);
  tree_->Branch("MuonPt",&MuonPt);
  tree_->Branch("MuonEnergy", &MuonEnergy);
  tree_->Branch("MuonCharge", &MuonCharge);
  tree_->Branch("MuonEta",&MuonEta);
  tree_->Branch("MuonPhi",&MuonPhi);

  tree_->Branch("MuonChi2P", &MuonChi2P);
  tree_->Branch("MuonChi2LocalPosition", &MuonChi2LocalPosition);
  tree_->Branch("MuonGlbTrackProbability", &MuonGlbTrackProbability);
  tree_->Branch("MuonTrkRelChi2", &MuonTrkRelChi2);
  tree_->Branch("MuonTrkKink", &MuonTrkKink);



  //Vtx position                                                                                                                                                                               
  tree_->Branch("Muon_vx", &Muon_vx);
  tree_->Branch("Muon_vy", &Muon_vy);
  tree_->Branch("Muon_vz", &Muon_vz);

  //MuonID                                                                                                                                                                                     
  tree_->Branch("Muon_isGlobal", &Muon_isGlobal);
  tree_->Branch("Muon_isTracker", &Muon_isTracker);
  tree_->Branch("Muon_isSoft", &Muon_isSoft);
  tree_->Branch("Muon_isLoose", &Muon_isLoose);
  tree_->Branch("Muon_isPF", &Muon_isPF);
  tree_->Branch("Muon_isRPCMuon", &Muon_isRPCMuon);
  tree_->Branch("Muon_isStandAloneMuon", &Muon_isStandAloneMuon);
  tree_->Branch("Muon_isTrackerMuon", &Muon_isTrackerMuon);
  tree_->Branch("Muon_isCaloMuon", &Muon_isCaloMuon);
  tree_->Branch("Muon_isQualityValid", &Muon_isQualityValid);
  tree_->Branch("Muon_isTimeValid", &Muon_isTimeValid);
  tree_->Branch("Muon_isIsolationValid", &Muon_isIsolationValid);
  tree_->Branch("Muon_numberOfMatchedStations", &Muon_numberOfMatchedStations);
  tree_->Branch("Muon_numberOfMatches", &Muon_numberOfMatches);


  tree_->Branch("Muon_timeAtIpInOut",&Muon_timeAtIpInOut);

  //Muon inner + outer track                                                                                                                                                                   
  tree_->Branch("Muon_GLnormChi2", &Muon_GLnormChi2);
  tree_->Branch("Muon_GLhitPattern_numberOfValidMuonHits", &Muon_GLhitPattern_numberOfValidMuonHits);

  tree_->Branch("Muon_trackerLayersWithMeasurement", &Muon_trackerLayersWithMeasurement);
  tree_->Branch("Muon_Numberofvalidpixelhits", &Muon_Numberofvalidpixelhits);

  tree_->Branch("Muon_outerTrack_p", &Muon_outerTrack_p);
  tree_->Branch("Muon_outerTrack_eta", &Muon_outerTrack_eta);
  tree_->Branch("Muon_outerTrack_phi", &Muon_outerTrack_phi);
  tree_->Branch("Muon_outerTrack_normalizedChi2", &Muon_outerTrack_normalizedChi2);
  tree_->Branch("Muon_outerTrack_muonStationsWithValidHits", &Muon_outerTrack_muonStationsWithValidHits);
  tree_->Branch("Muon_innerTrack_p", &Muon_innerTrack_p);
  tree_->Branch("Muon_innerTrack_eta", &Muon_innerTrack_eta);
  
  tree_->Branch("Muon_innerTrack_phi", &Muon_innerTrack_phi);
  tree_->Branch("Muon_innerTrack_normalizedChi2", &Muon_innerTrack_normalizedChi2);
  tree_->Branch("Muon_QInnerOuter", &Muon_QInnerOuter);


  tree_->Branch("Muon_combinedQuality_updatedSta", &Muon_combinedQuality_updatedSta);
  tree_->Branch("Muon_combinedQuality_trkKink", &Muon_combinedQuality_trkKink);
  tree_->Branch("Muon_combinedQuality_glbKink", &Muon_combinedQuality_glbKink);
  tree_->Branch("Muon_combinedQuality_trkRelChi2", &Muon_combinedQuality_trkRelChi2);
  tree_->Branch("Muon_combinedQuality_staRelChi2", &Muon_combinedQuality_staRelChi2);
  tree_->Branch("Muon_combinedQuality_chi2LocalPosition", &Muon_combinedQuality_chi2LocalPosition);
  tree_->Branch("Muon_combinedQuality_chi2LocalMomentum", &Muon_combinedQuality_chi2LocalMomentum);
  tree_->Branch("Muon_combinedQuality_localDistance", &Muon_combinedQuality_localDistance);
  tree_->Branch("Muon_combinedQuality_globalDeltaEtaPhi", &Muon_combinedQuality_globalDeltaEtaPhi);
  tree_->Branch("Muon_combinedQuality_tightMatch", &Muon_combinedQuality_tightMatch); 
  tree_->Branch("Muon_combinedQuality_glbTrackProbability", &Muon_combinedQuality_glbTrackProbability);
  
  tree_->Branch("Muon_calEnergy_em", &Muon_calEnergy_em);
  tree_->Branch("Muon_calEnergy_emS9", &Muon_calEnergy_emS9);
  tree_->Branch("Muon_calEnergy_emS25", &Muon_calEnergy_emS25);
  tree_->Branch("Muon_calEnergy_had", &Muon_calEnergy_had);
  tree_->Branch("Muon_calEnergy_hadS9", &Muon_calEnergy_hadS9);
  
  tree_->Branch("Muon_segmentCompatibility", &Muon_segmentCompatibility);
  tree_->Branch("Muon_caloCompatibility", &Muon_caloCompatibility);
  
  tree_->Branch("Muon_ptErrOverPt", &Muon_caloCompatibility);

  tree_->Branch("PVCollection_Size", &PVCollection_Size);
  tree_->Branch("PV_x", &PV_x);
  tree_->Branch("PV_y", &PV_y);
  tree_->Branch("PV_z", &PV_z);
  tree_->Branch("PV_NTracks", &PV_NTracks);

  tree_->Branch("TripletCollectionSize", &TripletCollectionSize);  
  tree_->Branch("Mu1_Pt",&Mu1_Pt);
  tree_->Branch("Mu1_Eta", &Mu1_Eta);
  tree_->Branch("Mu1_Phi", &Mu1_Phi);
  tree_->Branch("Mu2_Pt", &Mu2_Pt);
  tree_->Branch("Mu2_Eta", &Mu2_Eta);
  tree_->Branch("Mu2_Phi", &Mu2_Phi);
  tree_->Branch("Mu3_Pt", &Mu3_Pt);
  tree_->Branch("Mu3_Eta",&Mu3_Eta);
  tree_->Branch("Mu3_Phi", &Mu3_Phi);

  
  tree_->Branch("Mu1_SimPt", &Mu1_SimPt);
  tree_->Branch("Mu1_SimEta", &Mu1_SimEta);
  tree_->Branch("Mu1_SimPhi", &Mu1_SimPhi);
  tree_->Branch("Mu2_SimPt", &Mu2_SimPt);
  tree_->Branch("Mu2_SimEta", &Mu2_SimEta);
  tree_->Branch("Mu2_SimPhi", &Mu2_SimPhi);
  tree_->Branch("Mu3_SimPt", &Mu3_SimPt);
  tree_->Branch("Mu3_SimEta", &Mu3_SimEta);
  tree_->Branch("Mu3_SimPhi", &Mu3_SimPhi);

  tree_->Branch("TripletVtx_x", &TripletVtx_x);
  tree_->Branch("TripletVtx_y", &TripletVtx_y);
  tree_->Branch("TripletVtx_z", &TripletVtx_z);

  tree_->Branch("TripletVtx_Chi2", &TripletVtx_Chi2);
  tree_->Branch("TripletVtx_NDOF", &TripletVtx_NDOF);

  tree_->Branch("Triplet_Mass", &Triplet_Mass);
  tree_->Branch("Triplet_Pt", &Triplet_Pt);
  tree_->Branch("Triplet_Eta", &Triplet_Eta);
  tree_->Branch("Triplet_Phi", &Triplet_Phi);
  tree_->Branch("Triplet_Charge", &Triplet_Charge);
  
  tree_->Branch("RefittedPV_x", &RefittedPV_x);
  tree_->Branch("RefittedPV_y", &RefittedPV_y);
  tree_->Branch("RefittedPV_z", &RefittedPV_z);
  tree_->Branch("RefittedPV_NTracks", &RefittedPV_NTracks);
  //RefittedPV_Chi2.push_back(PVertex.);                                                                                                                                                     

  tree_->Branch("FlightDistPVSV", &FlightDistPVSV);
  tree_->Branch("FlightDistPVSV_Err", &FlightDistPVSV_Err);
  tree_->Branch("FlightDistPVSV_Significance", &FlightDistPVSV_Significance);
  tree_->Branch("FlightDistPVSV_chi2", &FlightDistPVSV_chi2);

  /*  SyncTree_ = fs->make<TTree>("t","Sync ntuple");
  SyncTree_ ->Branch("allmuons_pt",&allmuons_pt);
  SyncTree_->Branch("leadmuon_pt",&leadmuon_pt);
  SyncTree_->Branch("leadmuon_phi",&leadmuon_phi);
  SyncTree_->Branch("leadmuon_eta",&leadmuon_eta);
  SyncTree_->Branch("nmuons",&nmuons);
  SyncTree_->Branch("nprimevtxs",&nprimevtxs); 

  SyncTree_->Branch("leadtrack_pt", &leadtrack_pt);
  SyncTree_->Branch("leadtrack_eta", &leadtrack_eta);
  SyncTree_->Branch("leadtrack_phi", &leadtrack_phi);
  SyncTree_->Branch("alltracks_pt", &alltracks_pt);
  SyncTree_->Branch("evt", &evt);
  SyncTree_->Branch("run", &run);
  SyncTree_->Branch("lumi", &lumi);
  */



}
 

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAna2017Tree::endJob() 
{
  //  tree_->GetDirectory()->cd();
  tree_->Write();

  //  SyncTree_->GetDirectory()->cd();
  //  SyncTree_->Write();

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------


void MiniAna2017Tree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
  }
//define this as a plug-in
DEFINE_FWK_MODULE(MiniAna2017Tree);
