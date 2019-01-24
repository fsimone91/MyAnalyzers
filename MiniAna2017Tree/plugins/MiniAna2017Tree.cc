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

#include "/lustrehome/venditti/TestMiniAOD2017/CMSSW_9_4_5/src/MiniAna2017/JPsiKsPAT/src/miniAODmuons.h"

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
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > genParticles_;
  bool isMc;
  //  edm::EDGetTokenT<edm::TriggerResults> trigResultsToken;
  //  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigObjCollToken;
  const TransientTrackBuilder* theTransientTrackBuilder_; 


  edm::Service<TFileService> fs;
  TH1F *hEvents, *hEvents_3Mu, *hEvents_MuFilter, *hEvents_DiMuonDz, * hEvent_TauCharge,  * hEvent_3MuonVtx, *  hEvent_3MuonMass, *hEvent_Valid3MuonVtx;

  TH1F *hGenMuonPt;  TH1F *hGenMuonEta;
  TH1F *hMuonPt;TH1F *hMuonP; TH1F *hMuonEta, *hGlobalMuonEta, *hTrackerMuonEta, *hLooseMuonEta, *hSoftMuonEta;
  TH1F *  hGlobalMuonPt,  *  hSoftMuonPt,  *  hLooseMuonPt, *  hTrackerMuonPt;
  TH1F *hMuonNumberOfValidHits;
  TH1F * hMuonDB;   TH1F * hMuonTime;   TH1F * hMuonTimeErr;
  TH1F  *hGenTauPt, *hDiMuonDz, *hDiMuonDR, *hGoodMuSize, *hMuonSize_GoodDiMu;
  TH1F  *hThreeMuonInvMass, *hThreeMuonCharge, * hVtxSize, * hTau_vFit_chi2 ,  *h3MuonMassVtxFit  ;
  TH2F * hMuonTimeVsP;
 

  TTree*      tree_;
  std::vector<float>  MuonPt, MuonEta, MuonPhi, MuonChi2P, MuonChi2LocalPosition, MuonGlbTrackProbability, MuonTrkRelChi2, MuonTrkKink;

};



MiniAna2017Tree::MiniAna2017Tree(const edm::ParameterSet& iConfig){
  isMc = iConfig.getUntrackedParameter<bool>("isMcLabel");
  muons_ = consumes<edm::View<pat::Muon> >  (iConfig.getParameter<edm::InputTag>("muonLabel"));
  vertex_ = consumes<edm::View<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("VertexLabel"));
  genParticles_ = consumes<edm::View<pat::PackedGenParticle>  > (iConfig.getParameter<edm::InputTag>("genParticleLabel"));
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

  //isMc = true;  

  
  edm::Handle< edm::View<pat::PackedGenParticle> > genParticles;
  iEvent.getByToken(genParticles_, genParticles);
 
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
    for(edm::View<pat::PackedGenParticle>::const_iterator gp=genParticles->begin(); gp!=genParticles->end(), j<ngenP; ++gp , ++j){
      //if(fabs(gp->pdgId())==15) tauRaw = j;
    
      if(fabs(gp->pdgId())==13) {
	for (uint i=0; i<gp->numberOfMothers();i++){
	  if(fabs(gp->mother(i)->pdgId())==15) {
	    hGenMuonPt->Fill(gp->p()); 
	    hGenMuonEta->Fill(gp->eta());
	    hGenTauPt->Fill(gp->mother(i)->pt());
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
      cout<<"--genDR12="<<genDR12<<" genDR23="<<genDR23<<" genDR13="<<genDR13<<endl;
    }
    cout<<"****************GenLevel Info End ********************"<<endl;
  }
 


   cout<<"***Number of Muons="<<muons->size()<<endl; uint k=0;


   std::vector<int> MuFilter;
   vector<pat::Muon>    MyMu, MyMu2;
   for(edm::View<pat::Muon>::const_iterator mu=muons->begin(); mu!=muons->end(), k<muons->size(); ++mu, ++k){
    hMuonPt->Fill(mu->pt());
    hMuonP->Fill(mu->p());
    hMuonEta->Fill(mu->eta());
    hMuonNumberOfValidHits->Fill(mu->numberOfValidHits());

    if (!(mu->track().isNonnull()))continue;
    if(!(mu->track()->quality(reco::TrackBase::highPurity))) continue;
    if ( (!(mu->track()->hitPattern().numberOfValidPixelHits()>0))  ) continue;
    if ( !(mu->pt() > 0.5) ) continue;
    
    //mu.muonBestTrack()->dz(PV.position()), mu.isTightMuon(PV));
    std::cout<<"RecoMu pt="<<mu->pt()<<" eta="<<mu->eta()<<" phi="<<mu->phi()<<" segm compatibility="<<mu->segmentCompatibility()<<" trkRelChi2="<<mu->combinedQuality().trkRelChi2<<"  vz="<<mu->vz()<<" mu charge="<<mu->charge()<<" mu ECAL energy="<<mu->calEnergy().em<<std::endl;
    //" SymType="<<mu->simExtType()<<" SimPt="<<mu->simPt()<<
    //cout<<" pixel hits="<<mu->track()->hitPattern().numberOfValidPixelHits()<<" pt="<<mu->track()->pt()<<" mu ECAL energy="<<mu->calEnergy().em<<" mu.vx="<<mu->vx()<<endl;
    
    hMuonDB->Fill(mu->dB());
    hMuonTime->Fill(mu->time().timeAtIpInOut);
    hMuonTimeErr->Fill(mu->time().timeAtIpInOutErr);
    
    if (mu->isSoftMuon(PV)) {  hSoftMuonPt->Fill(mu->pt());     hSoftMuonEta->Fill(mu->eta());  }
    if (mu->isTrackerMuon()){  hTrackerMuonPt->Fill(mu->pt());  hTrackerMuonEta->Fill(mu->eta());}
    if (mu->isGlobalMuon()) {  hGlobalMuonPt->Fill(mu->pt());   hGlobalMuonEta->Fill(mu->eta());}
    if (mu->isLooseMuon())  {  hLooseMuonPt->Fill(mu->pt());    hLooseMuonEta->Fill(mu->eta());}

    if(!(mu->isSoftMuon(PV) || mu->isTrackerMuon() || mu->isGlobalMuon() || mu->isLooseMuon())) continue;

    MuFilter.push_back(1);  
    MyMu.push_back(*mu);


    MuonPt.push_back(mu->pt());
    MuonEta.push_back(mu->eta());
    MuonPhi.push_back(mu->phi());
    MuonChi2P.push_back(mu->combinedQuality().chi2LocalMomentum);
    MuonChi2LocalPosition.push_back(mu->combinedQuality().chi2LocalPosition);
    MuonGlbTrackProbability.push_back(mu->combinedQuality().glbTrackProbability);
    MuonTrkRelChi2.push_back(mu->combinedQuality().trkRelChi2);
    MuonTrkKink.push_back(mu->combinedQuality().trkKink);


    hMuonTimeVsP->Fill(mu->p(), mu->time().timeAtIpInOut);
      
    /*
    const reco::TransientTrack transientTrack = theTransientTrackBuilder_->build( *(mu->track()) );
    //transientTrack.setBeamSpot(vertexBeamSpot);
    Ttracks.push_back(transientTrack);
    */
   }

   hGoodMuSize->Fill(MyMu.size());    

   if (MuFilter.size()>0)   hEvents_MuFilter->Fill(1);
   cout<<"***Number of Muons after preliminary requ="<<MyMu.size()<<endl;

   ///DiMuon Selections
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

  //Selezione su TriMuonMass + vertex

  std::vector<reco::TransientTrack> Ttracks1, Ttracks2, Ttracks3;
  //vector<TransientVertex> MyVertexC;
  vector<RefCountedKinematicVertex>  VertexContainer;
  //  std::vector<reco::TransientTrack> Ttracks;
  vector<TransientVertex> KVContainer, KV_ValidContainer;

  vector<RefCountedKinematicTree> MyKinVtxTree;
  //Creating a KinematicParticleFactory                                                                                                                            
  /*  KinematicParticleFactoryFromTransientTrack pFactory;
  vector<RefCountedKinematicParticle> muonParticles;
  //initial chi2 and ndf before kinematic fits.                                                                                                                    
  ParticleMass muon_mass = 0.10565837; //pdg mass
  float muon_sigma = muon_mass*1.e-6;
  float chi = 0.;
  float ndf = 0.;  
  */

  TLorentzVector Mu1(0.,0.,0.,0.), Mu2(0.,0.,0.,0.), Mu3(0.,0.,0.,0.),  ThreeMuons; 
   //   vector<vector<int>> My3Mu;
   /* struct My3Mu{
     int Mu1Idx;
     int Mu2Idx;
     int Mu3Idx;
     } ;*/
   //   My3Mu My3MuC;
   //vector<My3Mu> My3MuV;


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
	   cout<<" Tracks used to refit the SV with Kalman="<<Ttracks.size()<<" with Kin="<<muonParticles.size()<<endl;
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
	   cout<<"3muon vtx: Kin Fit mass="<<Tau_vFit_noMC->currentState().mass()<<" mass"<<ThreeMuons.M()<<endl;


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
   ///////////////////////////////////////////metto le tracce del PV in una mappa/////////////////////////////////////////////////////////////                     
   
   std::vector<reco::TransientTrack> pvTracks_original;
   std::vector<reco::Track*> tauTrk;
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
   
   if(VertexContainer.size()>0){
     if(! (MyKinVtxTree.at(0)->isEmpty()) ){
       std::vector<RefCountedKinematicParticle> daughters = MyKinVtxTree.at(0)->daughterParticles();
       //cout<<"number of tracks in the PV="<< pvTracks_original.size()<<" number of daughters="<<daughters.size()<<endl;
       cout<< " number of daughters="<<daughters.size()<<endl; 
       /*
       for (std::vector<RefCountedKinematicParticle>::const_iterator i = daughters.begin();	i != daughters.end(); ++i) {
	 const TransientTrackKinematicParticle * ttkp = dynamic_cast<const TransientTrackKinematicParticle * >(&(**i));
	 const reco::TransientTrack * svTrack_transient = dynamic_cast<const reco::TransientTrack*>(ttkp->initialTransientTrack()->basicTransientTrack());
	 reco::Track  tmp_svtrack = (svTrack_transient->track());
	 reco::Track* tmp_ptr= &tmp_svtrack;
	 tauTrk.push_back(tmp_ptr); 
       }*/
       //removeTracks(pvTrackMap_refit, tauTrk);
     }
   }
      

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




  tree_->Fill();

  MuonPt.clear();
  MuonEta.clear();
  MuonPhi.clear();
  MuonChi2P.clear();
  MuonChi2LocalPosition.clear();
  MuonGlbTrackProbability.clear();
  MuonTrkRelChi2.clear();
  MuonTrkKink.clear();

}

// ------------ method called once each job just before starting event loop  ------------
void 
MiniAna2017Tree::beginJob()
{
  hEvents = fs->make<TH1F>("hEvents","hEvents",10,0,10);
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

  tree_ = fs->make<TTree>("ntuple","LFVTau ntuple");
  tree_->Branch("MuonPt",&MuonPt);
  tree_->Branch("MuonPt",&MuonPt);
  tree_->Branch("MuonEta",&MuonEta);
  tree_->Branch("MuonPhi",&MuonPhi);

  tree_->Branch("MuonChi2P", &MuonChi2P);
  tree_->Branch("MuonChi2LocalPosition", &MuonChi2LocalPosition);
  tree_->Branch("MuonGlbTrackProbability", &MuonGlbTrackProbability);
  tree_->Branch("MuonTrkRelChi2", &MuonTrkRelChi2);
  tree_->Branch("MuonTrkKink", &MuonTrkKink);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAna2017Tree::endJob() 
{
  tree_->GetDirectory()->cd();
  tree_->Write();


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
