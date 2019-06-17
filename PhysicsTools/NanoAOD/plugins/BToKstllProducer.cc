#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include <TLorentzVector.h>
#include <TVector.h>
#include <TMatrix.h>

////
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
////

#include <map>
#include <vector>


//
// class declaration
//

using namespace std;

class BToKstllProducer : public edm::EDProducer {
  
public:
  
  explicit BToKstllProducer(const edm::ParameterSet &iConfig);
  
  ~BToKstllProducer() override {};
  
    
private:

  void Init();
  void identifyTriggeringMuons(edm::Event&, const edm::EventSetup&, const reco::Vertex &);
  void selectPatElectrons(edm::Event&, const edm::EventSetup&);
  void selectPatMuons(edm::Event&, const edm::EventSetup&);
  void selectPackedCandidate(edm::Event&, const edm::EventSetup&, pat::PackedCandidateCollection&);
  void selectLostTracks(edm::Event&, const edm::EventSetup&, pat::PackedCandidateCollection&);
  void selectLowPtGsfTracks(edm::Event&, const edm::EventSetup&, 
			    edm::Handle<std::vector<reco::GsfTrack>>& lowPtGsfTracksHandle,                    
			    std::vector<edm::Handle<edm::ValueMap<float>>>& mvaSeedsHandle,                    
			    edm::Handle<edm::Association<std::vector<pat::PackedCandidate>>>& gsfLinkLTHandle, 
			    edm::Handle<edm::Association<std::vector<pat::PackedCandidate>>>& gsfLinkPCHandle, 
			    pat::CompositeCandidateCollection&);

  bool getCandLepX_PFele(unsigned int loop_index,
			 const reco::Candidate*& candLepX,
			 reco::TransientTrack& leptonX_TT,
			 float candLepX_Dxy,
			 float candLepX_Dz,
			 float candLepX_DxyS,
			 float candLepX_DzS,
			 bool isSubleading, unsigned int muonTrg_index, unsigned int iRun);

  bool getCandLepX_LowPtEle(unsigned int loop_index,
			    reco::TransientTrack& leptonX_TT,
			    float candLepX_Dxy,
			    float candLepX_Dz,
			    float candLepX_DxyS,
			    float candLepX_DzS,
			    bool isSubleading, unsigned int muonTrg_index, unsigned int iRun,
			    edm::Handle<std::vector<reco::GsfTrack>>& lowPtGsfTracksHandle,                     
			    std::vector<edm::Handle<edm::ValueMap<float>>>& mvaSeedsHandle,                     
			    edm::Handle<edm::Association<std::vector<pat::PackedCandidate>>>& gsfLinkLTHandle,  
			    edm::Handle<edm::Association<std::vector<pat::PackedCandidate>>>& gsfLinkPCHandle); 
    
  bool getCandLepX_PFmuon(unsigned int loop_index,
			  const reco::Candidate*& candLepX,
			  reco::TransientTrack& leptonX_TT,
			  float candLepX_Dxy,
			  float candLepX_Dz,
			  float candLepX_DxyS,
			  float candLepX_DzS,
			  bool isSubleading, unsigned int muonTrg_index, unsigned int iRun,
			  const reco::Vertex& PV);

  bool getCandX_PackedCand(unsigned int loop_index,
			   const reco::Candidate*& candLepX,
			   reco::TransientTrack& candX_TT,
			   float candX_Dxy,
			   float candX_Dz,
			   float candX_DxyS,
			   float candX_DzS,
			   bool isLeading, bool isSubleading, 
			   bool isKloop, unsigned int muonTrg_index, unsigned int iRun,
			   const math::XYZPoint &leplepRefitVertex=math::XYZPoint(0,0,0), 
			   float kaon_dxyFromRefitllVtx=0);

  bool getCandX_LostTracks(unsigned int loop_index,
			   const reco::Candidate*& candLepX,
			   reco::TransientTrack& candX_TT,
			   float candX_Dxy,
			   float candX_Dz,
			   float candX_DxyS,
			   float candX_DzS,
			   bool isLeading, bool isSubleading, 
			   bool isKloop, unsigned int muonTrg_index, unsigned int iRun,
			   const math::XYZPoint &leplepRefitVertex=math::XYZPoint(0,0,0), 
			   float kaon_dxyFromRefitllVtx=0);

  virtual void produce(edm::Event&, const edm::EventSetup&);
    
  bool LepLepVertexRefitting(const reco::TransientTrack &lep1TT,
			     const reco::TransientTrack &lep2TT,
                             RefCountedKinematicVertex &refitVertex,
                             math::XYZVector &refitLepLep, unsigned int iRun);
    
  bool KstVertexRefitting(const reco::TransientTrack &kaonTT,
			  const reco::TransientTrack &pionTT,
			  RefCountedKinematicVertex &refitVertex,
			  RefCountedKinematicParticle &refitKst);			 

  bool BToKLepLepVertexRefitting(const reco::TransientTrack &lep1TT,
				 const reco::TransientTrack &lep2TT,
				 const reco::TransientTrack &kaonTT,
                                 RefCountedKinematicVertex &refitVertex,
				 unsigned int iRun);

  bool BToKstLepLepVertexRefitting(const reco::TransientTrack &lep1TT,
				   const reco::TransientTrack &lep2TT,
				   const RefCountedKinematicParticle &refitKst,
				   RefCountedKinematicVertex &refitVertex,
				   math::XYZVector &refitKPi, unsigned int iRun);
    
  pair<double,double> computeLS(RefCountedKinematicVertex refitVertex,
				reco::BeamSpot beamSpot);

  pair<double,double> computeLS(RefCountedKinematicVertex refitVertex,
				reco::Vertex beamSpot);
  
  double computeCosAlpha(math::XYZVector &refitBToMLepLep,
			 RefCountedKinematicVertex vertexFitTree,
			 reco::BeamSpot beamSpot);


  pair<double,double> computeDCA(const reco::TransientTrack &hadronTT,
				 reco::BeamSpot beamSpot);
    

  // ----------member data ---------------------------
    
  edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
  edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
  edm::EDGetTokenT<std::vector<pat::Electron>> electronSrc_;  

  edm::EDGetTokenT<std::vector<reco::GsfTrack>> lowPtGsfTracksSrc_;
  std::vector<edm::InputTag> mvaSeedTags_;
  edm::EDGetTokenT<edm::Association<std::vector<pat::PackedCandidate>>> lowPtGsfLinksLT_;
  edm::EDGetTokenT<edm::Association<std::vector<pat::PackedCandidate>>> lowPtGsfLinksPC_;

  std::vector<edm::EDGetTokenT<edm::ValueMap<float>>>  mvaSeeds_;
  edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> PFCandSrc_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> lostSubLeadLepTrackSrc_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> lostChHadrTrackSrc_;

  /////////
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  
  float trgExclusiondR_;
  float trgAcceptdz_;

  int nSelectedTriplets_;

  //vectors to handle several configuration runs
  std::vector<unsigned> isLepEle_;
  std::vector<unsigned> isLowPtEle_;
  std::vector<unsigned> isLowPtAndPfEle_;
  std::vector<unsigned> isChKst_;
  std::vector<unsigned> useLostLeadMuonTracks_;
  std::vector<unsigned> useLostSubLeadLepTracks_;
  std::vector<unsigned> useLostChHadrTracks_;

  float ptMinLeadEle_;
  float etaMaxLeadEle_;
  float ptMinSubLeadEle_;
  float etaMaxSubLeadEle_;
  float bdtUnbiasedLeadLep_;
  float bdtUnbiasedSubLeadLep_;

  float ptMinLeadMuon_;
  float etaMaxLeadMuon_;
  float ptMinSubLeadMuon_;
  float etaMaxSubLeadMuon_;

  float ptMinKaon_;
  float etaMaxKaon_;
  float DCASigMinKaon_;

  float ptMinPion_;
  float etaMaxPion_;
  float DCASigMinPion_;

  bool diLepCharge_;
  bool KstCharge_;
  float JPsiMassConstraint_;
  float KstMassConstraint_;
  bool save2TrkRefit_;

  float vtxCL_min_;
  float Bmass_min_;
  float Bmass_max_;
  float Bmass_Kst_min_;
  float Bmass_Kst_max_;

  float diLepton_dz_max_;
  float lepKaon_dz_max_;
  float lepPion_dz_max_;
  float kaonPion_dz_max_;
  float kaonRefitllVertex_dxy_max_;
  float kll_dxyPV_min_;
  float IPPV_llRefitVtx_min_;

  float ElectronMass_ = 0.5109989e-3;
  float ElectronMassErr_ = 3.1*1e-12;
  float MuonMass_ = 0.10565837;
  float MuonMassErr_ = 3.5*1e-9;
  float KaonMass_ = 0.493677;
  float KaonMassErr_ = 1.6e-5;
  float PionMass_ = 0.139570;
  float PionMassErr_ = 3.5e-7;
  
  std::vector<float> lep1Mass_;
  std::vector<float> lep2Mass_;
  std::vector<float> lep1MassErr_;
  std::vector<float> lep2MassErr_;


  //float JPsiMass_ = 3.096916;  //Configurable parameter
  float JPsiMassErr_ = 0.011;
  
  //float KstMass_ = 0.89176;  //Configurable parameter  
  float KstMassErr_ = 0.25e-3;


  int nEvent;
  bool debug;


  edm::Handle<std::vector<pat::Muon>> muonForTrgHandle;
  edm::Handle<std::vector<pat::Electron>> electronHandle;
  edm::Handle<std::vector<pat::Muon>> muonHandle;
  edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle;
  edm::Handle<edm::View<pat::PackedCandidate>> lostSubLeadLepTrackHandle;
  edm::Handle<edm::View<pat::PackedCandidate>> lostChHadrTrackHandle;

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;

  TLorentzVector lepton1, lepton2, kaonTL, pionTL;
  float chargeL1, chargeL2, chargeK, chargeP;
  float vzL1, vzL2, vzK, vzP;
  std::vector<float> vzTrgMuo;
  float gsfTrk1_seedBDT_unbiased, gsfTrk1_seedBDT_ptbiased;
  float gsfTrk2_seedBDT_unbiased, gsfTrk2_seedBDT_ptbiased;

  //for the analysis
  std::vector<int> triggeringRecoMuons;
  std::map<int, std::vector<int>> selectedElectrons;
  std::map<int, std::vector<int>> selectedMuons;
  std::map<int, std::vector<int>> selectedPackedCandidates;
  std::map<int, std::vector<int>> selectedLostTracks;
  std::map<int, std::vector<int>> selectedLowPtGsfTracks;
  std::map<int, std::map<std::pair<edm::ProductID,int>, unsigned int>> checkLeptonsDuplicate;
  std::map<int, std::map<std::pair<edm::ProductID,int>, unsigned int>> checkLeptonsDuplicateTriplet;
};



BToKstllProducer::BToKstllProducer(const edm::ParameterSet &iConfig):
  beamSpotSrc_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
  vertexSrc_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ),
  electronSrc_( consumes<std::vector<pat::Electron>> ( iConfig.getParameter<edm::InputTag>( "electronCollection" ) ) ),
  lowPtGsfTracksSrc_( consumes<std::vector<reco::GsfTrack>>(iConfig.getParameter<edm::InputTag>("lowPtGsfTrackCollection") ) ),
  mvaSeedTags_( iConfig.getParameter<std::vector<edm::InputTag>>("mvaSeeds") ),
  lowPtGsfLinksLT_( consumes<edm::Association<std::vector<pat::PackedCandidate>>>(iConfig.getParameter<edm::InputTag>("lowPtGsfLinksLT")) ),
  lowPtGsfLinksPC_( consumes<edm::Association<std::vector<pat::PackedCandidate>>>(iConfig.getParameter<edm::InputTag>("lowPtGsfLinksPC")) ),
  muonSrc_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
  PFCandSrc_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
  lostSubLeadLepTrackSrc_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "lostSubLeadLepTrackCollection" ) ) ),
  lostChHadrTrackSrc_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "lostChHadrTrackCollection" ) ) ),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
   ////////
  trgExclusiondR_((float)iConfig.getParameter<double>("trgExclusiondR")),
  trgAcceptdz_((float)iConfig.getParameter<double>("trgAcceptdz")),
  nSelectedTriplets_( iConfig.getParameter<int>( "nSelectedTriplets" ) ),

  isLepEle_( iConfig.getParameter<std::vector<unsigned>>( "isLeptonElectron" ) ),
  isLowPtEle_( iConfig.getParameter<std::vector<unsigned>>( "isLowPtEle" ) ),
  isLowPtAndPfEle_( iConfig.getParameter<std::vector<unsigned>>( "isLowPtAndPfEle" ) ),
  isChKst_( iConfig.getParameter<std::vector<unsigned>>( "isChannelKst" ) ),
  useLostLeadMuonTracks_( iConfig.getParameter<std::vector<unsigned>>( "useLostLeadMuonTracks" ) ),
  useLostSubLeadLepTracks_( iConfig.getParameter<std::vector<unsigned>>( "useLostSubLeadLepTracks" ) ),
  useLostChHadrTracks_( iConfig.getParameter<std::vector<unsigned>>( "useLostChHadrTracks" ) ),

  ptMinLeadEle_( (float)iConfig.getParameter<double>( "LeadEleMinPt" ) ),
  etaMaxLeadEle_( (float)iConfig.getParameter<double>( "LeadEleMaxEta" ) ),
  ptMinSubLeadEle_( (float)iConfig.getParameter<double>( "SubLeadEleMinPt" ) ),
  etaMaxSubLeadEle_( (float)iConfig.getParameter<double>( "SubLeadEleMaxEta" ) ),
  bdtUnbiasedLeadLep_( (float)iConfig.getParameter<double>( "LeadBDTUnbiased") ),
  bdtUnbiasedSubLeadLep_( (float)iConfig.getParameter<double>( "SubLeadBDTUnbiased") ),

  ptMinLeadMuon_( (float)iConfig.getParameter<double>( "LeadMuonMinPt" ) ),
  etaMaxLeadMuon_( (float)iConfig.getParameter<double>( "LeadMuonMaxEta") ),
  ptMinSubLeadMuon_( (float)iConfig.getParameter<double>( "SubLeadMuonMinPt" ) ),
  etaMaxSubLeadMuon_( (float)iConfig.getParameter<double>( "SubLeadMuonMaxEta" ) ),

  ptMinKaon_( (float)iConfig.getParameter<double>( "KaonMinPt" ) ),
  etaMaxKaon_( (float)iConfig.getParameter<double>( "KaonMaxEta" ) ),
  DCASigMinKaon_( (float)iConfig.getParameter<double>( "KaonMinDCASig" ) ),
  ptMinPion_( (float)iConfig.getParameter<double>( "PionMinPt" ) ),
  etaMaxPion_( (float)iConfig.getParameter<double>( "PionMaxEta" ) ),
  DCASigMinPion_( (float)iConfig.getParameter<double>( "PionMinDCASig" ) ),
  diLepCharge_( iConfig.getParameter<bool>( "DiLeptonChargeCheck" ) ),
  KstCharge_( iConfig.getParameter<bool>( "KstarChargeCheck" ) ),

  JPsiMassConstraint_( (float)iConfig.getParameter<double>( "JPsiMassConstraint" ) ),
  KstMassConstraint_( (float)iConfig.getParameter<double>( "KstMassConstraint" ) ),

  save2TrkRefit_( iConfig.getParameter<bool>( "save2TrackRefit" ) ),
//  save4TrkRefit_( iConfig.getParameter<bool>( "save4TrackRefit" ) ),

  vtxCL_min_( (float)iConfig.getParameter<double>( "vtxCL_min" ) ),
  Bmass_min_( (float)iConfig.getParameter<double>( "Bmass_min" ) ),
  Bmass_max_( (float)iConfig.getParameter<double>( "Bmass_max" ) ),
  Bmass_Kst_min_( (float)iConfig.getParameter<double>( "Bmass_Kst_min" ) ),
  Bmass_Kst_max_( (float)iConfig.getParameter<double>( "Bmass_Kst_max" ) ),
  diLepton_dz_max_( (float)iConfig.getParameter<double>( "diLepton_dz_max" ) ),
  lepKaon_dz_max_( (float)iConfig.getParameter<double>( "lepKaon_dz_max" ) ),
  lepPion_dz_max_( (float)iConfig.getParameter<double>( "lepPion_dz_max" ) ),
  kaonPion_dz_max_( (float)iConfig.getParameter<double>( "kaonPion_dz_max" ) ),
  kaonRefitllVertex_dxy_max_( (float)iConfig.getParameter<double>( "kaonRefitllVertex_dxy_max" ) ),
  kll_dxyPV_min_( (float)iConfig.getParameter<double>( "kll_dxyPV_min" ) ),
  IPPV_llRefitVtx_min_( (float)iConfig.getParameter<double>( "IPPV_llRefitVtx_min" ) )
{
  nEvent = 0;
  debug = false;

  if(debug) std::cout << " costruttore " << std::endl;

  if(debug){
    std::cout << " in costruttore isLepEle_.size() = " << isLepEle_.size() << std::endl;
  }

  for(unsigned int iRun = 0; iRun<isLepEle_.size(); ++iRun){

    if(isLowPtAndPfEle_[iRun] && isLowPtEle_[iRun]){
      throw cms::Exception("configError") << " ERROR isLowPtAndPfEle_ and isLowPtEle_ not compatible at the moment in config \n";
    }

    if(useLostLeadMuonTracks_[iRun] && (isLepEle_[iRun] || isLowPtEle_[iRun] || isLowPtAndPfEle_[iRun])){
      throw cms::Exception("configError") << " ERROR useLostLeadMuonTracks_ and isLepEle_/isLowPtEle_/isLowPtAndPfEle_ not compatible at the moment in config \n";
    }

    if(useLostSubLeadLepTracks_[iRun] && (isLowPtEle_[iRun] || isLowPtAndPfEle_[iRun])){
      //std::cout << " ERROR useLostSubLeadLepTracks_ and isLowPtEle_ not compatible at the moment " << std::endl;
      throw cms::Exception("configError") << " ERROR useLostSubLeadLepTracks_ and isLowPtEle_/isLowPtAndPfEle_ not compatible at the moment in config \n";
    }

    if(!isLepEle_[iRun] && (isLowPtEle_[iRun] || isLowPtAndPfEle_[iRun])){
      throw cms::Exception("configError") << " ERROR isLepEle_ should be true with  isLowPtEle_/isLowPtAndPfEle_  in config \n";
    }

    lep1Mass_.push_back( (isLepEle_[iRun] == 1) ? ElectronMass_ : MuonMass_);
    lep2Mass_.push_back(lep1Mass_[iRun]);
    lep1MassErr_.push_back( (isLepEle_[iRun] == 1) ? ElectronMassErr_ : MuonMassErr_);
    lep2MassErr_.push_back(lep1MassErr_[iRun]);
  }


  for ( const auto& tag : mvaSeedTags_ ){
    mvaSeeds_.push_back( consumes<edm::ValueMap<float>>(tag) );
  }

  produces<pat::CompositeCandidateCollection>();
  produces<pat::PackedCandidateCollection>("PFCand");
  produces<pat::PackedCandidateCollection>("LostTracks");
  produces<pat::CompositeCandidateCollection>("LowPtGsfTrack");

}

void BToKstllProducer::Init(){

  if(debug) std::cout << " BToKstllProducer::Init " << std::endl;


  muonForTrgHandle.clear();
  electronHandle.clear();
  muonHandle.clear();
  pfCandHandle.clear();
  lostSubLeadLepTrackHandle.clear();
  lostChHadrTrackHandle.clear();

  vzTrgMuo.clear();
  triggeringRecoMuons.clear();
  selectedElectrons.clear();
  selectedMuons.clear();
  selectedPackedCandidates.clear();
  selectedLostTracks.clear();
  selectedLowPtGsfTracks.clear();
  return;
}


void BToKstllProducer::identifyTriggeringMuons(edm::Event& iEvent, const edm::EventSetup& iSetup, 
					       const reco::Vertex & PV){

  if(debug) std::cout << " BToKstllProducer::identifyTriggeringMuons " << std::endl;

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

  std::vector<std::vector<float>> triggeringMuons;
  //taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Trigger
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  if(debug) std::cout << "\n TRIGGER OBJECTS " << std::endl;

  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
    obj.unpackFilterLabels(iEvent, *triggerBits);
    obj.unpackPathNames(names);

    bool isTriggerMuon = false;
    for (unsigned h = 0; h < obj.filterIds().size(); ++h)
      if(obj.filterIds()[h] == 83){ 
	isTriggerMuon = true; 
	if(debug) std::cout << "\t   Type IDs:   " << 83;  //83 = muon
	break;
      } 

    if(!isTriggerMuon) continue; 
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
      std::string filterName = obj.filterLabels()[h];
      if(filterName.find("hltL3") != std::string::npos  && filterName.find("Park") != std::string::npos){
	isTriggerMuon = true;
	if(debug) std::cout << "\t   Filters:   " << filterName; 
	break;
      }
      else{ isTriggerMuon = false; }
    }

    if(!isTriggerMuon) continue;
    std::vector<float> localMuon;
    localMuon.push_back(obj.pt());
    localMuon.push_back(obj.eta());
    localMuon.push_back(obj.phi());
    triggeringMuons.push_back(localMuon);
    if(debug){ std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      // Print trigger object collection and type
      std::cout << "\t   Collection: " << obj.collection() << std::endl;
    }
  }//trigger objects

  if(debug){
    std::cout << "\n total n of triggering muons = " << triggeringMuons.size() << std::endl;
    for(auto ij : triggeringMuons){
      for(auto il : ij) std::cout << " >>> components (pt, eta, phi) = " << il << std::endl;
    }
  }

  triggeringRecoMuons.resize(triggeringMuons.size(), -1);
  vzTrgMuo.resize(triggeringMuons.size(), -1);
  std::vector<float> dRMuonMatching(triggeringMuons.size(), -1);

  //  edm::Handle<std::vector<pat::Muon>> muonForTrgHandle;
  iEvent.getByToken(muonSrc_, muonForTrgHandle);

  for(unsigned int iTrg=0; iTrg<muonForTrgHandle->size(); ++iTrg){
    const pat::Muon & muon1 = (*muonForTrgHandle)[iTrg];

    if(!(muon1.isLooseMuon() && muon1.isSoftMuon(PV))) continue;
    TLorentzVector recoMuon(muon1.pt(), muon1.eta(), muon1.phi(), MuonMass_);

    for(unsigned int ij=0; ij<triggeringMuons.size(); ++ij){
      //std::cout << " >>> comp = " << (triggeringMuons.at(ij))[0] << std::endl;
      TLorentzVector trgMuon((triggeringMuons.at(ij))[0], (triggeringMuons.at(ij))[1], (triggeringMuons.at(ij))[2], MuonMass_);
      float dR = recoMuon.DeltaR(trgMuon);
      if((dR < dRMuonMatching[ij] || dRMuonMatching[ij] == -1) && dR < 0.01){
	dRMuonMatching[ij] = dR;
	triggeringRecoMuons[ij] = iTrg;
	vzTrgMuo[ij] = muon1.vz();
       	if(debug) std::cout << " dR = " << dR 
			    << " reco = " << recoMuon.Pt() << " " << recoMuon.Eta() << " " << recoMuon.Phi() << " " 
			    << " HLT = " << trgMuon.Pt() << " " << trgMuon.Eta() << " " << trgMuon.Phi()
			    << std::endl;
      }
    }
  }

  return;
}

void BToKstllProducer::selectPatElectrons(edm::Event& iEvent, const edm::EventSetup& iSetup){

  if(debug) std::cout << " BToKstllProducer::selectPatElectrons " << std::endl;

  //  edm::Handle<std::vector<pat::Electron>> electronHandle;
  iEvent.getByToken(electronSrc_, electronHandle);

  //loop over triggering muons
  for (auto iTM : triggeringRecoMuons){
    if(iTM == -1) continue;
    const pat::Muon &trgMuon = (*muonForTrgHandle)[iTM];

    if(debug)  std::cout << " electronHandle->size() = " << electronHandle->size() << std::endl;
    for(unsigned int iEle=0; iEle<electronHandle->size(); ++iEle){
      const pat::Electron & ele = (*electronHandle)[iEle];
      float dR = deltaR(trgMuon.eta(), trgMuon.phi(), ele.eta(), ele.phi());
      
      if(debug) std::cout << " dR wrt trg muon = " << dR
                          << " pt = " << ele.pt() << " eta = " << ele.eta() << std::endl;

      if(dR < trgExclusiondR_ || std::fabs(trgMuon.vz() - ele.vz()) > trgAcceptdz_) continue;
      if(ele.pt()<std::min(ptMinLeadEle_, ptMinSubLeadEle_) || abs(ele.eta()) > std::max(etaMaxLeadEle_, etaMaxSubLeadEle_)) continue;

      selectedElectrons[iTM].push_back(iEle);

      //just save ref a PFCand
      for(unsigned int ic = 0; ic < ele.numberOfSourceCandidatePtrs(); ++ic){
	if(!ele.sourceCandidatePtr(ic).isNonnull() || !ele.sourceCandidatePtr(ic).isAvailable()) continue;
	reco::CandidatePtr dummyCandLep = ele.sourceCandidatePtr(ic);
	checkLeptonsDuplicate[iTM][std::pair<edm::ProductID,int>(dummyCandLep.id(), dummyCandLep.key())] = iEle;
      }
    }
    if(debug) std::cout << " post selectPatElectrons size = " << selectedElectrons[iTM].size() << std::endl;
  } //loop over trigger muons
  return;
}


void BToKstllProducer::selectPatMuons(edm::Event& iEvent, const edm::EventSetup& iSetup){

  if(debug) std::cout << " BToKstllProducer::selectPatMuons " << std::endl;

  //  edm::Handle<std::vector<pat::Muon>> muonHandle;
  iEvent.getByToken(muonSrc_, muonHandle);

  //loop over triggering muons
  for (auto iTM : triggeringRecoMuons){
    if(iTM == -1) continue;
    const pat::Muon &trgMuon = (*muonForTrgHandle)[iTM];
    
    if(debug)  std::cout << " muonHandle->size() = " << muonHandle->size() << std::endl;
    for(unsigned int iMuo=0; iMuo<muonHandle->size(); ++iMuo){
      if(debug)  std::cout << " iMuo = " << iMuo << " iTM = " << iTM << std::endl;
      if((unsigned int)iTM == iMuo) continue;
      const pat::Muon & muon = (*muonHandle)[iMuo];
      //selections for 1st lepton
      //if(!(muon1.isLooseMuon() && muon1.isSoftMuon(PV))) continue;
      float dR = deltaR(trgMuon.eta(), trgMuon.phi(), muon.eta(), muon.phi());
     
      if(debug) std::cout << " dR wrt trg muon = " << dR
			  << " pt = " << muon.pt() << " eta = " << muon.eta() << std::endl;

      if(dR < trgExclusiondR_ || std::fabs(trgMuon.vz() - muon.vz()) > trgAcceptdz_) continue;
      if(muon.pt()<std::min(ptMinLeadMuon_, ptMinSubLeadMuon_) || abs(muon.eta()) > std::max(etaMaxLeadMuon_, etaMaxSubLeadMuon_)) continue;
      
      selectedMuons[iTM].push_back(iMuo);
      
      //just save ref a PFCand
      for(unsigned int ic = 0; ic < muon.numberOfSourceCandidatePtrs(); ++ic){
	if(!muon.sourceCandidatePtr(ic).isNonnull() || !muon.sourceCandidatePtr(ic).isAvailable()) continue;
	reco::CandidatePtr dummyCandLep = muon.sourceCandidatePtr(ic);
	checkLeptonsDuplicate[iTM][std::pair<edm::ProductID,int>(dummyCandLep.id(), dummyCandLep.key())] = iMuo;
      }
    }
    if(debug) std::cout << " post selectPatMuons size = " << selectedMuons[iTM].size() << std::endl;
  }// loop over trigger muons

  return;
}

void BToKstllProducer::selectPackedCandidate(edm::Event& iEvent, const edm::EventSetup& iSetup, 
					     pat::PackedCandidateCollection& resultPF){

  if(debug) std::cout << " BToKstllProducer::selectPackedCandidate " << std::endl;

  //  edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle;
  iEvent.getByToken(PFCandSrc_, pfCandHandle);

  //loop over triggering muons                                                        
  for (auto iTM : triggeringRecoMuons){
    if(iTM == -1) continue;
    const pat::Muon &trgMuon = (*muonForTrgHandle)[iTM];

    for(unsigned int iPcand=0; iPcand<pfCandHandle->size(); ++iPcand){

      edm::ProductID pair1 = pfCandHandle.id();
      int pair2 = iPcand;
      std::pair<edm::ProductID, int> dummyPair(pair1, pair2);
      std::map<std::pair<edm::ProductID,int>, unsigned int>::iterator checkLD_it = checkLeptonsDuplicate[iTM].find(dummyPair);
      if(checkLD_it != checkLeptonsDuplicate[iTM].end()) {
	if(debug) std::cout << " found duplicate track" << std::endl;
	continue;
      }
      //      if(debug) std::cout << " NON duplicate track " << std::endl;
      
      const pat::PackedCandidate & pfC = (*pfCandHandle)[iPcand];
      float dR = deltaR(trgMuon.eta(), trgMuon.phi(), pfC.eta(), pfC.phi());
      
      if(dR < trgExclusiondR_ || std::fabs(trgMuon.vz() - pfC.vz()) > trgAcceptdz_) continue;

      //ID criteria
      //      if(debug) std::cout << " >>> pfC.hasTrackDetails() = " << pfC.hasTrackDetails() << " pfC.pdgId() = " << pfC.pdgId()  << std::endl;
      if(!pfC.hasTrackDetails()) continue;
      //exclude neutral should be safe do not ask too much ID                                                                       
      if(abs(pfC.pdgId()) == 0) continue;

      float minPt = ptMinKaon_;
      float maxEta = etaMaxKaon_;
      for(auto ij : isChKst_){
	if(ij){
	  minPt = std::min(ptMinKaon_, ptMinPion_);
	  maxEta = std::max(etaMaxKaon_, etaMaxPion_);
	  break;
	}
      }
      if(pfC.pt()<std::min(minPt, ptMinSubLeadMuon_) || abs(pfC.eta())>std::max(maxEta, etaMaxSubLeadMuon_)) continue;

      // in the analysis consider pt > 3 if subleading lepton all otherwise
      //if(pfC.pt() > 3.) continue;
      
      selectedPackedCandidates[iTM].push_back(iPcand);

      //prepare skimmed PFCandidate collection
      resultPF.push_back(pfC);
      //      checkLeptonsDuplicateCollection[iTM][dummyPair] = iMuo;
    }
    if(debug) std::cout << " post selectPackedCandidate size = " << selectedPackedCandidates[iTM].size() << std::endl;
  }//triggering muon

  return;
}

void BToKstllProducer::selectLostTracks(edm::Event& iEvent, const edm::EventSetup& iSetup, 
					pat::PackedCandidateCollection& resultLT){

  if(debug) std::cout << " BToKstllProducer::selectLostTracks " << std::endl;

  //  edm::Handle<edm::View<pat::PackedCandidate>> lostChHadrTrackHandle;
  iEvent.getByToken(lostChHadrTrackSrc_, lostChHadrTrackHandle);

  //loop over triggering muons
  for (auto iTM : triggeringRecoMuons){
    if(iTM == -1) continue;
    const pat::Muon &trgMuon = (*muonForTrgHandle)[iTM];
    
    for(unsigned int iLtrk=0; iLtrk<lostChHadrTrackHandle->size(); ++iLtrk){

      edm::ProductID pair1 = lostChHadrTrackHandle.id();
      int pair2 = iLtrk;
      std::pair<edm::ProductID, int> dummyPair(pair1, pair2);
      std::map<std::pair<edm::ProductID,int>, unsigned int>::iterator checkLD_itk = checkLeptonsDuplicate[iTM].find(dummyPair);
      if(checkLD_itk != checkLeptonsDuplicate[iTM].end()){
	if(debug) std::cout << " found duplicate lost track " << std::endl;
	continue;
      }
      //      if(debug) std::cout << " NON duplicate lost track " << std::endl;

      const pat::PackedCandidate & pfC = (*lostChHadrTrackHandle)[iLtrk];
      float dR = deltaR(trgMuon.eta(), trgMuon.phi(), pfC.eta(), pfC.phi());

      if(dR < trgExclusiondR_ || std::fabs(trgMuon.vz() - pfC.vz()) > trgAcceptdz_) continue;

      if(!pfC.hasTrackDetails()) continue;
      if(abs(pfC.pdgId()) != 211) continue; // && abs(kaon.pdgId()) != 11 && abs(kaon.pdgId()) != 13) continue; //Charged hadrons

      float minPt = ptMinKaon_;
      float maxEta = etaMaxKaon_;
      for(auto ij : isChKst_){
        if(ij){
          minPt = std::min(ptMinKaon_, ptMinPion_);
          maxEta = std::max(etaMaxKaon_, etaMaxPion_);
          break;
        }
      }
      if(pfC.pt()<std::min(minPt, ptMinSubLeadMuon_) || abs(pfC.eta())>std::max(maxEta, etaMaxSubLeadMuon_)) continue;
	 
      selectedLostTracks[iTM].push_back(iLtrk);

      //prepare skimmed PFCandidate collection          
      resultLT.push_back(pfC);
    }
    if(debug) std::cout << " post selectLostTracks size = " << selectedLostTracks[iTM].size() << std::endl;
  }//triggering muon
  return;
}


void BToKstllProducer::selectLowPtGsfTracks(edm::Event& iEvent, const edm::EventSetup& iSetup, 
					    edm::Handle<std::vector<reco::GsfTrack>>& lowPtGsfTracksHandle,                    
					    std::vector<edm::Handle<edm::ValueMap<float>>>& mvaSeedsHandle,                    
					    edm::Handle<edm::Association<std::vector<pat::PackedCandidate>>>& gsfLinkLTHandle, 
					    edm::Handle<edm::Association<std::vector<pat::PackedCandidate>>>& gsfLinkPCHandle, 
					    pat::CompositeCandidateCollection& resultGsf){

  if(debug) std::cout << " BToKstllProducer::selectLowPtGsfTracks " << std::endl;

  // edm::Handle<std::vector<reco::GsfTrack>> lowPtGsfTracksHandle;
  // std::vector<edm::Handle<edm::ValueMap<float>>> mvaSeedsHandle;
  // edm::Handle<edm::Association<std::vector<pat::PackedCandidate>>> gsfLinkLTHandle;
  // edm::Handle<edm::Association<std::vector<pat::PackedCandidate>>> gsfLinkPCHandle;

  // iEvent.getByToken(lowPtGsfTracksSrc_, lowPtGsfTracksHandle);
  // for (const auto& token : mvaSeeds_){
  //   edm::Handle<edm::ValueMap<float> > h;
  //   iEvent.getByToken(token, h);
  //   mvaSeedsHandle.push_back(h);
  // }
  // iEvent.getByToken(lowPtGsfLinksLT_, gsfLinkLTHandle);
  // iEvent.getByToken(lowPtGsfLinksPC_, gsfLinkPCHandle);


  //loop over triggering muons                                                      
  for (auto iTM : triggeringRecoMuons){
    if(iTM == -1) continue;
    const pat::Muon &trgMuon = (*muonForTrgHandle)[iTM];


    if(debug)  std::cout << " lowPtGsfTracksHandle->size() = " << lowPtGsfTracksHandle->size() << std::endl;    
    for(unsigned int iLgsf=0; iLgsf<lowPtGsfTracksHandle->size(); ++iLgsf){
      
      if(mvaSeedsHandle.size() != 2 || !mvaSeedsHandle[0].isValid() || mvaSeedsHandle[0]->empty() ||
	 !mvaSeedsHandle[1].isValid() || mvaSeedsHandle[1]->empty() || !lowPtGsfTracksHandle.isValid()) continue;
      
      const reco::GsfTrackRef gsfTrk(lowPtGsfTracksHandle, iLgsf);
      float gsfTrk_seedBDT_unbiased = float((*mvaSeedsHandle[0])[gsfTrk]);
      float gsfTrk_seedBDT_ptbiased = float((*mvaSeedsHandle[1])[gsfTrk]);
      
      if(gsfTrk_seedBDT_unbiased < std::min(bdtUnbiasedLeadLep_, bdtUnbiasedSubLeadLep_)) continue;
      
      float dR = deltaR(trgMuon.eta(), trgMuon.phi(), gsfTrk->etaMode(), gsfTrk->phiMode());

      //      if(debug) std::cout << " dR wrt trg muon = " << dR << " dz = " <<  std::fabs(trgMuon.vz() - gsfTrk->vz())
      //			  << " pt = " << gsfTrk->ptMode() << " eta = " << gsfTrk->etaMode() << std::endl;

      if(dR < trgExclusiondR_ || std::fabs(trgMuon.vz() - gsfTrk->vz()) > trgAcceptdz_) continue;
      
      pat::PackedCandidateRef ltRef = (*gsfLinkLTHandle)[gsfTrk];
      pat::PackedCandidateRef pcRef = (*gsfLinkPCHandle)[gsfTrk];
      if(debug){
	std::cout << " ltRef.isNonnull() = " << ltRef.isNonnull() << " ltRef.isAvailable() = " << ltRef.isAvailable()
		  << " pcRef.isNonnull() = " << pcRef.isNonnull() << " pcRef.isAvailable() = " << pcRef.isAvailable()
		  << std::endl;
      }

      for(auto ij : isLowPtAndPfEle_){
	if(ij){
	edm::ProductID pair1lt = ltRef.id();
	int pair2lt = ltRef.key();
	std::pair<edm::ProductID, int> dummyPairlt(pair1lt, pair2lt);
	std::map<std::pair<edm::ProductID,int>, unsigned int>::iterator checkLT_itk = checkLeptonsDuplicate[iTM].find(dummyPairlt);
	edm::ProductID pair1pc = pcRef.id();
	int pair2pc = pcRef.key();
	std::pair<edm::ProductID, int> dummyPairpc(pair1pc, pair2pc);
	std::map<std::pair<edm::ProductID,int>, unsigned int>::iterator checkPC_itk = checkLeptonsDuplicate[iTM].find(dummyPairpc);
	
	if(checkLT_itk != checkLeptonsDuplicate[iTM].end() || checkPC_itk != checkLeptonsDuplicate[iTM].end()){
	  if(debug) std::cout << " found duplicate low pt gsf track " << std::endl;
	  continue;
	}
	if(debug) std::cout << " NON duplicate low pt gsf track " << std::endl;

	if(gsfTrk->ptMode()<std::min(ptMinLeadEle_, ptMinSubLeadEle_) || abs(gsfTrk->etaMode()) > std::max(etaMaxLeadEle_, etaMaxSubLeadEle_)) continue;
	break;
	}
      }

      selectedLowPtGsfTracks[iTM].push_back(iLgsf);

      /*
      //just save ref a PFCand                                                                                                              
      if(ltRef.isNonnull() && ltRef.isAvailable()){
	checkLeptonsDuplicateCollection[std::pair<edm::ProductID,int>(ltRef.id(), ltRef.key())] = iLgsf;
      }
      if(pcRef.isNonnull() && pcRef.isAvailable()){
	checkLeptonsDuplicateCollection[std::pair<edm::ProductID,int>(pcRef.id(), pcRef.key())] = iLgsf;
	}
      */

      pat::CompositeCandidate resultGsfpcNew;
      resultGsfpcNew.addUserFloat("pt", gsfTrk->ptMode());
      resultGsfpcNew.addUserFloat("eta", gsfTrk->etaMode());
      resultGsfpcNew.addUserFloat("phi", gsfTrk->phiMode());
      resultGsfpcNew.addUserFloat("charge", gsfTrk->chargeMode());

      resultGsfpcNew.addUserFloat("seedBDT_unbiased", gsfTrk_seedBDT_unbiased);
      resultGsfpcNew.addUserFloat("seedBDT_ptbiased", gsfTrk_seedBDT_ptbiased);
      resultGsfpcNew.addUserFloat("dz", gsfTrk->dz());
      resultGsfpcNew.addUserFloat("dxy", gsfTrk->dxy());
      resultGsfpcNew.addUserFloat("dzS", gsfTrk->dz()/gsfTrk->dzError());
      resultGsfpcNew.addUserFloat("dxyS", gsfTrk->dxy()/gsfTrk->dxyError());
      resultGsfpcNew.addUserFloat("vz", gsfTrk->vz());
      resultGsf.push_back(resultGsfpcNew);
    }
    if(debug) std::cout << " post selectLowPtGsfTracks size = " << selectedLowPtGsfTracks[iTM].size() << std::endl;
  }//triggering muon
  return;
}


bool BToKstllProducer::getCandLepX_PFele(unsigned int loop_index,
					 const reco::Candidate*& candLepX,
					 reco::TransientTrack& leptonX_TT,
					 float candLepX_Dxy,
					 float candLepX_Dz,
					 float candLepX_DxyS,
					 float candLepX_DzS,
					 bool isSubleading, unsigned int muonTrg_index, unsigned int iRun){

  if(debug) std::cout << " BToKstllProducer::getCandLepX_PFele " << std::endl;

  const pat::Electron & eleX = (*electronHandle)[loop_index];
  candLepX = &eleX;

  if(!isSubleading){
    lepton1.SetPtEtaPhiM(candLepX->pt(), candLepX->eta(), candLepX->phi(), (isLepEle_[iRun] == 1) ? ElectronMass_ : MuonMass_);
    chargeL1 = candLepX->charge();
    vzL1 = candLepX->vz();
  }
  else{
    lepton2.SetPtEtaPhiM(candLepX->pt(), candLepX->eta(), candLepX->phi(), (isLepEle_[iRun] == 1) ? ElectronMass_ : MuonMass_);
    chargeL2 = candLepX->charge();
    vzL2 = candLepX->vz();
  }

  leptonX_TT = theTTBuilder->build(eleX.gsfTrack());
  candLepX_Dxy = eleX.dB(pat::Electron::PV2D);
  candLepX_Dz = eleX.dB(pat::Electron::PVDZ);
  candLepX_DxyS = candLepX_Dxy / eleX.edB(pat::Electron::PV2D);
  candLepX_DzS = candLepX_Dz / eleX.edB(pat::Electron::PVDZ);

  if(!isSubleading && (lepton1.Pt()<ptMinLeadEle_ || abs(lepton1.Eta()) > etaMaxLeadEle_)) return false;
  if(isSubleading){
    if(lepton2.Pt()<ptMinSubLeadEle_ || abs(lepton2.Eta()) > etaMaxSubLeadEle_) return false;

    //Lepton 1 is always saved as the leading one
    //should be already ok by looping on subleading from leading+1
    if(lepton1.Pt() < lepton2.Pt()) return false; 
    if(debug && chargeL1 != chargeL2) std::cout << " lepton1Charge = " << chargeL1 << " lepton2Charge = " << chargeL2 << std::endl;
    if(diLepCharge_ && chargeL1*chargeL2 > 0) return false;

  }
  

  //just save ref a PFCand
  for(unsigned int ic = 0; ic < eleX.numberOfSourceCandidatePtrs(); ++ic){
    if(!eleX.sourceCandidatePtr(ic).isNonnull() || !eleX.sourceCandidatePtr(ic).isAvailable()) continue;
    reco::CandidatePtr dummyCandLep = eleX.sourceCandidatePtr(ic);
    std::pair<edm::ProductID,int> keyVal(dummyCandLep.id(), dummyCandLep.key());
    checkLeptonsDuplicateTriplet[muonTrg_index][keyVal] = loop_index;
  }


  if(debug) std::cout << " passed lepton subleading = " << isSubleading << " - pfEle - " << std::endl;
  return true;
}


bool BToKstllProducer::getCandLepX_LowPtEle(unsigned int loop_index,
					    reco::TransientTrack& leptonX_TT,
					    float candLepX_Dxy,
					    float candLepX_Dz,
					    float candLepX_DxyS,
					    float candLepX_DzS,
					    bool isSubleading, unsigned int muonTrg_index, unsigned int iRun,
					    edm::Handle<std::vector<reco::GsfTrack>>& lowPtGsfTracksHandle,                    
					    std::vector<edm::Handle<edm::ValueMap<float>>>& mvaSeedsHandle,                    
					    edm::Handle<edm::Association<std::vector<pat::PackedCandidate>>>& gsfLinkLTHandle, 
					    edm::Handle<edm::Association<std::vector<pat::PackedCandidate>>>& gsfLinkPCHandle){

  if(debug) std::cout << " BToKstllProducer::getCandLepX_LowPtEle " << std::endl;

  if(mvaSeedsHandle.size() != 2 || !mvaSeedsHandle[0].isValid() || mvaSeedsHandle[0]->empty() ||
     !mvaSeedsHandle[1].isValid() || mvaSeedsHandle[1]->empty() || !lowPtGsfTracksHandle.isValid()) return false;

  const reco::GsfTrackRef gsfTrk(lowPtGsfTracksHandle, loop_index);
  if(!isSubleading){
    gsfTrk1_seedBDT_unbiased = float((*mvaSeedsHandle[0])[gsfTrk]);
    gsfTrk1_seedBDT_ptbiased = float((*mvaSeedsHandle[1])[gsfTrk]);

    lepton1.SetPtEtaPhiM(gsfTrk->ptMode(), gsfTrk->etaMode(), gsfTrk->phiMode(), ElectronMass_);
    chargeL1 = gsfTrk->chargeMode();
    vzL1 = gsfTrk->vz();
  }
  else{
    gsfTrk2_seedBDT_unbiased = float((*mvaSeedsHandle[0])[gsfTrk]);
    gsfTrk2_seedBDT_ptbiased = float((*mvaSeedsHandle[1])[gsfTrk]);

    lepton2.SetPtEtaPhiM(gsfTrk->ptMode(), gsfTrk->etaMode(), gsfTrk->phiMode(), ElectronMass_);
    chargeL2 = gsfTrk->chargeMode();
    vzL2 = gsfTrk->vz();
  }

  leptonX_TT = theTTBuilder->build(gsfTrk);
  candLepX_Dxy = gsfTrk->dxy();            
  candLepX_Dz = gsfTrk->dz(); 
  candLepX_DxyS = candLepX_Dxy / gsfTrk->dxyError();
  candLepX_DzS = candLepX_Dz / gsfTrk->dzError();


  if(!isSubleading){
    if(lepton1.Pt()<ptMinLeadEle_ || abs(lepton1.Eta()) > etaMaxLeadEle_) return false;

    // std::cout << "pre !isSubleading gsfTrk1_seedBDT_unbiased = " << gsfTrk1_seedBDT_unbiased 
    // 	      << " bdtUnbiasedLeadLep_ = " << bdtUnbiasedLeadLep_ << std::endl;

    if(float(gsfTrk1_seedBDT_unbiased) < float(bdtUnbiasedLeadLep_)) return false;

    // std::cout << "post !isSubleading gsfTrk1_seedBDT_unbiased = " << gsfTrk1_seedBDT_unbiased 
    // 	      << " bdtUnbiasedLeadLep_ = " << bdtUnbiasedLeadLep_ << std::endl;

    //pfEle + lowptGsf high pT objects already covered by pfEle
    if(isLowPtAndPfEle_[iRun] && lepton1.Pt() > 5.) return false;
  }
  if(isSubleading){
    if(lepton2.Pt()<ptMinSubLeadEle_ || abs(lepton2.Eta()) > etaMaxSubLeadEle_) return false;
    if(gsfTrk2_seedBDT_unbiased < bdtUnbiasedSubLeadLep_) return false;

    if(lepton1.Pt() < lepton2.Pt()) return false;
    if(diLepCharge_ && chargeL1*chargeL2 > 0) return false;

    if(isLowPtAndPfEle_[iRun] && deltaR(lepton1.Eta(), lepton1.Phi(), lepton2.Eta(), lepton2.Phi()) < 0.01) return false;

    //pfEle + lowptGsf high pT objects already covered by pfEle
    if(isLowPtAndPfEle_[iRun] && lepton2.Pt() > 5.) return false;
  }

  //add reference to low pt gsf used to clean tracks in the kaon/pion search
  pat::PackedCandidateRef ltRef = (*gsfLinkLTHandle)[gsfTrk];
  pat::PackedCandidateRef pcRef = (*gsfLinkPCHandle)[gsfTrk];

  //just save ref a PFCand                                                                        
  if(ltRef.isNonnull() && ltRef.isAvailable()){
    checkLeptonsDuplicateTriplet[muonTrg_index][std::pair<edm::ProductID,int>(ltRef.id(), ltRef.key())] = loop_index;
  }
  if(pcRef.isNonnull() && pcRef.isAvailable()){
    checkLeptonsDuplicateTriplet[muonTrg_index][std::pair<edm::ProductID,int>(pcRef.id(), pcRef.key())] = loop_index;
  }

  if(debug) std::cout << " passed lepton subleading = " << isSubleading << " - low pt gsf - " << std::endl;
  return true;

}


bool BToKstllProducer::getCandLepX_PFmuon(unsigned int loop_index,
					  const reco::Candidate*& candLepX,
					  reco::TransientTrack& leptonX_TT,
					  float candLepX_Dxy,
					  float candLepX_Dz,
					  float candLepX_DxyS,
					  float candLepX_DzS,
					  bool isSubleading, unsigned int muonTrg_index, unsigned int iRun, const reco::Vertex& PV){

  if(debug) std::cout << " BToKstllProducer::getCandLepX_PFmuon " << std::endl;

  const pat::Muon & muonX = (*muonHandle)[loop_index];
  candLepX = &muonX;

  if(!isSubleading){
    lepton1.SetPtEtaPhiM(candLepX->pt(), candLepX->eta(), candLepX->phi(), (isLepEle_[iRun] == 1) ? ElectronMass_ : MuonMass_);
    chargeL1 = candLepX->charge();
    vzL1 = candLepX->vz();
  }
  else{
    lepton2.SetPtEtaPhiM(candLepX->pt(), candLepX->eta(), candLepX->phi(), (isLepEle_[iRun] == 1) ? ElectronMass_ : MuonMass_);
    chargeL2 = candLepX->charge();
    vzL2 = candLepX->vz();
  }
  if(debug) std::cout << " post1 " << std::endl;

  leptonX_TT = theTTBuilder->build(muonX.bestTrack());
  candLepX_Dxy = muonX.dB(pat::Muon::PV2D);
  candLepX_Dz = muonX.dB(pat::Muon::PVDZ);
  candLepX_DxyS = candLepX_Dxy/muonX.edB(pat::Muon::PV2D);
  candLepX_DzS = candLepX_Dz/muonX.edB(pat::Muon::PVDZ);

  if(debug) std::cout << " post2 " << std::endl;

  if(!isSubleading){
    if(lepton1.Pt()<ptMinLeadMuon_ || abs(lepton1.Eta()) > etaMaxLeadMuon_) return false;
    if( !muonX.isLooseMuon() || !muonX.isSoftMuon(PV)) return false;
  }
  if(isSubleading){
    if(lepton2.Pt()<ptMinSubLeadMuon_ || abs(lepton2.Eta()) > etaMaxSubLeadMuon_) return false;

    //Lepton 1 is always saved as the leading one
    //should be already ok by looping on subleading from leading+1
    if(lepton1.Pt() < lepton2.Pt()) return false; 
    if(debug && chargeL1 != chargeL2) std::cout << " lepton1Charge = " << chargeL1 << " lepton2Charge = " << chargeL2 << std::endl;
    if(diLepCharge_ && chargeL1*chargeL2 > 0) return false;
  }
  
  //just save ref a PFCand
  for(unsigned int ic = 0; ic < muonX.numberOfSourceCandidatePtrs(); ++ic){
    if(!muonX.sourceCandidatePtr(ic).isNonnull() || !muonX.sourceCandidatePtr(ic).isAvailable()) continue;
    reco::CandidatePtr dummyCandLep = muonX.sourceCandidatePtr(ic);
    std::pair<edm::ProductID,int> keyVal(dummyCandLep.id(), dummyCandLep.key());
    checkLeptonsDuplicateTriplet[muonTrg_index][keyVal] = loop_index;
    //removableL2.push_back(keyVal);                                                                       
  }
  if(debug) std::cout << " passed lepton subleading = " << isSubleading << " - pfMuon - " << std::endl;
  return true;
}



bool BToKstllProducer::getCandX_PackedCand(unsigned int loop_index,
					   const reco::Candidate*& candX,
					   reco::TransientTrack& candX_TT,
					   float candX_Dxy,
					   float candX_Dz,
					   float candX_DxyS,
					   float candX_DzS,
					   bool isLeading, bool isSubleading, 
					   bool isKloop, unsigned int muonTrg_index, unsigned int iRun,
					   const math::XYZPoint &leplepRefitVertex, float kaon_dxyFromRefitllVtx){
  
  if(debug) std::cout << " BToKstllProducer::getCandLepX_PackedCand " << std::endl;

  const pat::PackedCandidate &cand = (*pfCandHandle)[loop_index];
  candX = &cand;

  edm::ProductID pair1 = pfCandHandle.id();
  int pair2 = loop_index;
  std::pair<edm::ProductID, int> dummyPair(pair1, pair2);
  std::map<std::pair<edm::ProductID,int>, unsigned int>::iterator checkLD_it = checkLeptonsDuplicateTriplet[muonTrg_index].find(dummyPair);
  if(checkLD_it != checkLeptonsDuplicateTriplet[muonTrg_index].end()) {
    if(debug) std::cout << " found duplicate L2" << std::endl;
    return false;
  }
  if(debug) std::cout << " NON duplicate L2 " << std::endl;

  if(isLeading){
    lepton1.SetPtEtaPhiM(candX->pt(), candX->eta(), candX->phi(), (isLepEle_[iRun] == 1) ? ElectronMass_ : MuonMass_);
    chargeL1 = candX->charge();
    vzL1 = candX->vz();
  }
  else if(isSubleading){
    lepton2.SetPtEtaPhiM(candX->pt(), candX->eta(), candX->phi(), (isLepEle_[iRun] == 1) ? ElectronMass_ : MuonMass_);
    chargeL2 = candX->charge();
    vzL2 = candX->vz();
  }
  else if(isKloop){
    kaonTL.SetPtEtaPhiM(candX->pt(), candX->eta(), candX->phi(), KaonMass_);
    chargeK = candX->charge();
    vzK = candX->vz();
  }
  else{
    pionTL.SetPtEtaPhiM(candX->pt(), candX->eta(), candX->phi(), PionMass_);
    chargeP = candX->charge();
    vzP = candX->vz();
  }

  candX_TT = theTTBuilder->build(cand.bestTrack());
  candX_Dxy = cand.dxy();
  candX_Dz = cand.dz();
  candX_DxyS = candX_Dxy/cand.dxyError();
  candX_DzS = candX_Dz/cand.dzError();

  if(isLeading){
    if(lepton1.Pt() < ptMinLeadMuon_ || abs(lepton1.Eta()) > etaMaxLeadMuon_) return false;
    if(lepton1.Pt() > 5.) return false;  
  }
  else if(isSubleading){
    if(lepton2.Pt()<ptMinSubLeadMuon_ || abs(lepton2.Eta()) > etaMaxSubLeadMuon_) return false;

    //Lepton 1 is always saved as the leading one
    //should be already ok by looping on subleading from leading+1
    if(lepton1.Pt() < lepton2.Pt()) return false; 
    if(debug && chargeL1 != chargeL2) std::cout << " lepton1Charge = " << chargeL1 << " lepton2Charge = " << chargeL2 << std::endl;
    if(diLepCharge_ && chargeL1*chargeL2 > 0) return false;

    //assume in mixed config you want to cover only low pt with pf candidates (tracks)
    if(lepton2.Pt() > 5.) return false;

    //lepton1 and lepton2 belong to different collections need to check they are different
    if(deltaR(lepton1.Eta(), lepton1.Phi(), lepton2.Eta(), lepton2.Phi()) < 0.01) return false;

    if(diLepton_dz_max_ > -1. && std::abs(vzL1 - vzL2) > diLepton_dz_max_) return false;
  }
  else if(isKloop){
    if(abs(candX->pdgId()) != 211) return false;
    if(kaonTL.Pt() < ptMinKaon_ || abs(kaonTL.Eta()) > etaMaxKaon_) return false;

    if(deltaR(lepton1.Eta(), lepton1.Phi(), kaonTL.Eta(), kaonTL.Phi()) < 0.01 ||
       deltaR(lepton2.Eta(), lepton2.Phi(), kaonTL.Eta(), kaonTL.Phi()) < 0.01 ) return false;

    if(save2TrkRefit_){
      kaon_dxyFromRefitllVtx = cand.dxy(leplepRefitVertex);
      if(kaonRefitllVertex_dxy_max_ != -1 &&
	 std::abs(kaon_dxyFromRefitllVtx) > kaonRefitllVertex_dxy_max_) return false;
    }
  }
  else{
    if(abs(candX->pdgId())!=211) return false;
    if(pionTL.Pt() < ptMinPion_ || abs(pionTL.Eta()) > etaMaxPion_) return false;

    if(KstCharge_ && chargeK * chargeP>0) return false;

    if(deltaR(lepton1.Eta(), lepton1.Phi(), pionTL.Eta(), pionTL.Phi()) < 0.01 ||
       deltaR(lepton2.Eta(), lepton2.Phi(), pionTL.Eta(), pionTL.Phi()) < 0.01 ||
       deltaR(kaonTL.Eta(), kaonTL.Phi(), pionTL.Eta(), pionTL.Phi()) < 0.01) return false;
  }

  if(isLeading || isSubleading)
    checkLeptonsDuplicateTriplet[muonTrg_index][dummyPair] = loop_index;

  if(debug) std::cout << " passed pfcandidate subleading = " << isSubleading  << " isKloop = " << isKloop << " - pfCand - " << std::endl;
  return true;
}



bool BToKstllProducer::getCandX_LostTracks(unsigned int loop_index,
					   const reco::Candidate*& candX,
					   reco::TransientTrack& candX_TT,
					   float candX_Dxy,
					   float candX_Dz,
					   float candX_DxyS,
					   float candX_DzS,
					   bool isLeading, bool isSubleading, 
					   bool isKloop, unsigned int muonTrg_index, unsigned int iRun,
					   const math::XYZPoint &leplepRefitVertex, float kaon_dxyFromRefitllVtx){
  

  if(debug) std::cout << " BToKstllProducer::getCandLepX_LostTracks " << std::endl;

  const pat::PackedCandidate &cand = (*lostChHadrTrackHandle)[loop_index];
  candX = &cand;


  edm::ProductID pair1 = lostChHadrTrackHandle.id();
  int pair2 = loop_index;
  std::pair<edm::ProductID, int> dummyPair(pair1, pair2);
  std::map<std::pair<edm::ProductID,int>, unsigned int>::iterator checkLD_it = checkLeptonsDuplicateTriplet[muonTrg_index].find(dummyPair);
  if(checkLD_it != checkLeptonsDuplicateTriplet[muonTrg_index].end()) {
    if(debug) std::cout << " found duplicate L2" << std::endl;
    return false;
  }
  if(debug) std::cout << " NON duplicate L2 " << std::endl;

  
  if(isLeading){
    lepton1.SetPtEtaPhiM(candX->pt(), candX->eta(), candX->phi(), (isLepEle_[iRun] == 1) ? ElectronMass_ : MuonMass_);
    chargeL1 = candX->charge();
    vzL1 = candX->vz();
  }
  else if(isSubleading){
    lepton2.SetPtEtaPhiM(candX->pt(), candX->eta(), candX->phi(), (isLepEle_[iRun] == 1) ? ElectronMass_ : MuonMass_);
    chargeL2 = candX->charge();
    vzL2 = candX->vz();
  }
  else if(isKloop){
    kaonTL.SetPtEtaPhiM(candX->pt(), candX->eta(), candX->phi(), KaonMass_);
    chargeK = candX->charge();
    vzK = candX->vz();
  }
  else{
    pionTL.SetPtEtaPhiM(candX->pt(), candX->eta(), candX->phi(), PionMass_);
    chargeP = candX->charge();
    vzP = candX->vz();
  }

  candX_TT = theTTBuilder->build(cand.bestTrack());
  candX_Dxy = cand.dxy();
  candX_Dz = cand.dz();
  candX_DxyS = candX_Dxy/cand.dxyError();
  candX_DzS = candX_Dz/cand.dzError();

  if(isLeading){
    if(lepton1.Pt() < ptMinLeadMuon_ || abs(lepton1.Eta()) > etaMaxLeadMuon_) return false;
    if(lepton1.Pt() > 5.) return false;
  }  
  else if(isSubleading){
    if(lepton2.Pt()<ptMinSubLeadMuon_ || abs(lepton2.Eta()) > etaMaxSubLeadMuon_) return false;
    
    //Lepton 1 is always saved as the leading one
    //should be already ok by looping on subleading from leading+1
    if(lepton1.Pt() < lepton2.Pt()) return false; 
    if(debug && chargeL1 != chargeL2) std::cout << " lepton1Charge = " << chargeL1 << " lepton2Charge = " << chargeL2 << std::endl;
    if(diLepCharge_ && chargeL1*chargeL2 > 0) return false;

    //assume in mixed config you want to cover only low pt with pf candidates (tracks)
    if(lepton2.Pt() > 5.) return false;

    //lepton1 and lepton2 belong to different collections need to check they are different      
    if(deltaR(lepton1.Eta(), lepton1.Phi(), lepton2.Eta(), lepton2.Phi()) < 0.01) return false;
    if(diLepton_dz_max_ > -1. && std::abs(vzL1 - vzL2) > diLepton_dz_max_) return false;
  }
  else if(isKloop){
    if(abs(candX->pdgId()) != 211) return false;
    if(kaonTL.Pt() < ptMinKaon_ || abs(kaonTL.Eta()) > etaMaxKaon_) return false;

    if(deltaR(lepton1.Eta(), lepton1.Phi(), kaonTL.Eta(), kaonTL.Phi()) < 0.01 ||
       deltaR(lepton2.Eta(), lepton2.Phi(), kaonTL.Eta(), kaonTL.Phi()) < 0.01 ) return false;

    if(save2TrkRefit_){
      kaon_dxyFromRefitllVtx = cand.dxy(leplepRefitVertex);
      if(kaonRefitllVertex_dxy_max_ != -1 &&
	 std::abs(kaon_dxyFromRefitllVtx) > kaonRefitllVertex_dxy_max_) return false;
    }
  }
  else{
    if(abs(candX->pdgId())!=211) return false;
    if(pionTL.Pt() < ptMinPion_ || abs(pionTL.Eta()) > etaMaxPion_) return false;

    if(KstCharge_ && chargeK * chargeP>0) return false;

    if(deltaR(lepton1.Eta(), lepton1.Phi(), pionTL.Eta(), pionTL.Phi()) < 0.01 ||
       deltaR(lepton2.Eta(), lepton2.Phi(), pionTL.Eta(), pionTL.Phi()) < 0.01 ||
       deltaR(kaonTL.Eta(), kaonTL.Phi(), pionTL.Eta(), pionTL.Phi()) < 0.01) return false;

  }

  if(isLeading || isSubleading)
    checkLeptonsDuplicateTriplet[muonTrg_index][dummyPair] = loop_index;

  if(debug) std::cout << " passed losttrack subleading = " << isSubleading  << " isKloop = " << isKloop << " - lostTrack - " << std::endl;
  return true;
}



void BToKstllProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
  ++nEvent;
  if(debug) std::cout << "\n new event = " << nEvent << std::endl;

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vertexSrc_, vertexHandle);
  const reco::Vertex & PV = vertexHandle->front();

  Init();

  identifyTriggeringMuons(iEvent, iSetup, PV);

  // Output collection
  std::unique_ptr<pat::CompositeCandidateCollection> tripletResults (new pat::CompositeCandidateCollection);
  std::unique_ptr<pat::PackedCandidateCollection> resultPF( new pat::PackedCandidateCollection );
  std::unique_ptr<pat::PackedCandidateCollection> resultLT( new pat::PackedCandidateCollection );
  std::unique_ptr<pat::CompositeCandidateCollection> resultGsf( new pat::CompositeCandidateCollection );


  if(triggeringRecoMuons.size() == 0){
    if(debug) std::cout << " NO reco muon natched to triggering muon " << std::endl;
    iEvent.put(std::move(tripletResults));
    iEvent.put(std::move(resultPF), "PFCand");
    iEvent.put(std::move(resultLT), "LostTracks");
    iEvent.put(std::move(resultGsf), "LowPtGsfTrack");
    return;
  }

  selectPatElectrons(iEvent, iSetup);
  selectPatMuons(iEvent, iSetup);

  edm::Handle<std::vector<reco::GsfTrack>> lowPtGsfTracksHandle;                    
  std::vector<edm::Handle<edm::ValueMap<float>>> mvaSeedsHandle;                    
  edm::Handle<edm::Association<std::vector<pat::PackedCandidate>>> gsfLinkLTHandle; 
  edm::Handle<edm::Association<std::vector<pat::PackedCandidate>>> gsfLinkPCHandle; 

  //for the moment because standard MC does not have lowptgsftracks
  //to be removed then
  for(auto ij: isLepEle_){
    if(ij){
      iEvent.getByToken(lowPtGsfTracksSrc_, lowPtGsfTracksHandle);
      for (const auto& token : mvaSeeds_){
	edm::Handle<edm::ValueMap<float> > h;
	iEvent.getByToken(token, h);
	mvaSeedsHandle.push_back(h);
      }
      iEvent.getByToken(lowPtGsfLinksLT_, gsfLinkLTHandle);
      iEvent.getByToken(lowPtGsfLinksPC_, gsfLinkPCHandle);

      selectLowPtGsfTracks(iEvent, iSetup, lowPtGsfTracksHandle, mvaSeedsHandle, gsfLinkLTHandle, gsfLinkPCHandle, *resultGsf);
      break;
    }
  }
  selectPackedCandidate(iEvent, iSetup, *resultPF);
  for(auto ij: useLostChHadrTracks_){
    if(ij){
      selectLostTracks(iEvent, iSetup, *resultLT);
      break;
    }
  }


  edm::ESHandle<MagneticField> bFieldHandle;
  
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
    
  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("BToKstllProducer") << "No beam spot available from EventSetup" ;
  }
  
  reco::BeamSpot beamSpot = *beamSpotHandle;

  // Output collection
  std::vector<std::unique_ptr<pat::CompositeCandidateCollection>> result;
  for(unsigned int iRun=0; iRun<isLepEle_.size(); ++iRun)
    result.push_back(std::unique_ptr<pat::CompositeCandidateCollection> (new pat::CompositeCandidateCollection));
  std::vector<float> resultTag[2];
  std::pair<float, int> worstTag_val_idx[2];

  //RAFIX loop twice: 1 with ele config & 2 with muon config

  for(unsigned int iRun=0; iRun<isLepEle_.size(); ++iRun){

  //init ad ogni run sul config
  checkLeptonsDuplicate.clear();
  checkLeptonsDuplicateTriplet.clear();

  //loop over triggering muons
  for (auto iTM : triggeringRecoMuons){
    if(iTM == -1) continue;
    /*
    const pat::Muon &trgMuon = (*muonForTrgHandle)[iTM];
    TLorentzVector trgMuonTL(trgMuon.pt(), trgMuon.eta(), trgMuon.phi(), MuonMass_);
    float trgMuonVZ = trgMuon.vz();
    */

    unsigned int leptonNumber = (isLepEle_[iRun] && isLowPtEle_[iRun]) ? (selectedLowPtGsfTracks[iTM].size()) : (isLepEle_[iRun] ? selectedElectrons[iTM].size() : selectedMuons[iTM].size());
    if(isLepEle_[iRun] && isLowPtAndPfEle_[iRun]) leptonNumber = selectedElectrons[iTM].size()+selectedLowPtGsfTracks[iTM].size();
    unsigned int pfCandNumber = selectedPackedCandidates[iTM].size();    
    unsigned int LeadMuonTrackNumber = useLostLeadMuonTracks_[iRun] ? (pfCandNumber + selectedLostTracks[iTM].size()) : 0;
    unsigned int subLeadLeptonTrackNumber = useLostSubLeadLepTracks_[iRun] ? (leptonNumber + pfCandNumber + selectedLostTracks[iTM].size()) : leptonNumber;
    unsigned int lostChHadrTrackNumber = useLostChHadrTracks_[iRun] ? selectedLostTracks[iTM].size() : 0;
    
    if(debug){
      std::cout << " leptonNumber = " << leptonNumber << " pfCandNumber = " << pfCandNumber
		<< " subLeadLeptonTrackNumber = " << subLeadLeptonTrackNumber << " lostChHadrTrackNumber = " << lostChHadrTrackNumber << std::endl; 
      std::cout << " [iRun] = " << iRun << " isLepEle_ = " << isLepEle_[iRun] << " isLowPtEle_ = " << isLowPtEle_[iRun]
		<< " isLowPtAndPfEle_ = " << isLowPtAndPfEle_[iRun] << " isChKst_ = " << isChKst_[iRun] << std::endl;
    }
    
    // loop on leading lepton from pat::Lepton 
    // subleading from pat::Lepton or PFcandidate+LT
    // and then kaon (+pion) for K or Kstar
    // to build triplets
    for (unsigned int i = 0; i < leptonNumber + LeadMuonTrackNumber; ++i) {
      //checkLeptonsDuplicate.clear();

      float candLep1Dxy = -99.;
      float candLep1Dz = -99.;      
      float candLep1DxyS = -99.;
      float candLep1DzS = -99.;
      const reco::Candidate* candLep1;
      reco::TransientTrack lepton1TT; 

      bool isLep1LowPt = bool((isLowPtEle_[iRun]) || (isLepEle_[iRun] && isLowPtAndPfEle_[iRun] && i >= selectedElectrons[iTM].size()));
      if(debug) std::cout << " isLep1LowPt = " << isLep1LowPt << std::endl;
      bool isLep1PFL = bool(i < leptonNumber && isLep1LowPt==false);
      bool isLep1PFC = bool( i < (pfCandNumber+leptonNumber) && isLep1LowPt == false && isLep1PFL == false );
      bool isLep1LT = bool(i >= (pfCandNumber+leptonNumber));

      bool objectPassed = false;

      if(isLepEle_[iRun] || isLowPtAndPfEle_[iRun]){


	if(isLep1LowPt)
	  objectPassed = getCandLepX_LowPtEle(selectedLowPtGsfTracks[iTM].at(i-selectedElectrons[iTM].size()), lepton1TT, 
					      candLep1Dxy, candLep1Dz,candLep1DxyS, candLep1DzS,
					      false, iTM, iRun,
					      lowPtGsfTracksHandle, mvaSeedsHandle, gsfLinkLTHandle, gsfLinkPCHandle);
	
	else
	  objectPassed = getCandLepX_PFele(selectedElectrons[iTM].at(i), candLep1, lepton1TT, 
					   candLep1Dxy, candLep1Dz,candLep1DxyS, candLep1DzS,
					   false, iTM, iRun);
	
      }
      else if(isLepEle_[iRun] && !isLowPtAndPfEle_[iRun] && isLowPtEle_[iRun])
	objectPassed = getCandLepX_LowPtEle(selectedLowPtGsfTracks[iTM].at(i), lepton1TT, 
					    candLep1Dxy, candLep1Dz,candLep1DxyS, candLep1DzS,
					    false, iTM, iRun,
					    lowPtGsfTracksHandle, mvaSeedsHandle, gsfLinkLTHandle, gsfLinkPCHandle);

      else if(!isLepEle_[iRun]){

        if(i < leptonNumber)
            objectPassed = getCandLepX_PFmuon(selectedMuons[iTM].at(i), candLep1, lepton1TT,
                                              candLep1Dxy, candLep1Dz, candLep1DxyS, candLep1DzS,false, iTM, iRun, PV);
        else if(isLep1PFC)
            objectPassed = getCandX_PackedCand(selectedPackedCandidates[iTM].at(i-leptonNumber), candLep1, lepton1TT,
                                               candLep1Dxy, candLep1Dz,candLep1DxyS, candLep1DzS, true, false, false, iTM, iRun);
        else if(isLep1LT)
            objectPassed = getCandX_LostTracks(selectedLostTracks[iTM].at(i-leptonNumber-pfCandNumber), candLep1, lepton1TT, 
					       candLep1Dxy, candLep1Dz, candLep1DxyS, candLep1DzS, true, false, false, iTM, iRun);

      }

      if(!objectPassed) continue;
      if(debug) std::cout << " passed lepton 1 " << std::endl;
      //      if(debug && !isLowPtAndPfEle_ && !isLowPtEle_) std::cout << " candLep1->charge() = " << candLep1ref->charge() << std::endl;

      for (unsigned int j = i+1; j < subLeadLeptonTrackNumber; ++j) {
	bool isLep2LowPt = bool((isLowPtEle_[iRun] && !isLowPtAndPfEle_[iRun]) || (isLepEle_[iRun] && isLowPtAndPfEle_[iRun] && j >= selectedElectrons[iTM].size()));
	bool isLep2PFL = bool(j < leptonNumber && isLep2LowPt == false);
	bool isLep2PFC = bool(j < (pfCandNumber+leptonNumber) && isLep2LowPt == false && isLep2PFL == false );
	bool isLep2LT = bool(j >= (pfCandNumber+leptonNumber));

	if(debug) std::cout << " i = " << i << " j = " << j << " leptonNumber = " << leptonNumber <<std::endl;

	if(debug) std::cout << " isLep2LowPt = " << isLep2LowPt << " isLep2PFL = " << isLep2PFL
			    << " isLep2PFC = " << isLep2PFC << " isLep2LT = " << isLep2LT << std::endl;

 
	float candLep2Dxy = -99.;
	float candLep2Dz = -99.;
	float candLep2DxyS = -99.;
	float candLep2DzS = -99.;
	const reco::Candidate* candLep2;
	reco::TransientTrack lepton2TT; 

	if(isLep2PFL || isLep2LowPt){
	  if(isLepEle_[iRun] || isLowPtAndPfEle_[iRun]){

	    if(isLep2LowPt)
	      objectPassed = getCandLepX_LowPtEle(selectedLowPtGsfTracks[iTM].at(j-selectedElectrons[iTM].size()), lepton2TT,
						  candLep2Dxy, candLep2Dz,candLep2DxyS, candLep2DzS,
						  true, iTM, iRun,
						  lowPtGsfTracksHandle, mvaSeedsHandle, gsfLinkLTHandle, gsfLinkPCHandle);
	    
	    else
	      objectPassed = getCandLepX_PFele(selectedElectrons[iTM].at(j), candLep2, lepton2TT,
					       candLep2Dxy, candLep2Dz,candLep2DxyS, candLep2DzS,
					       true, iTM, iRun);

	  }
	  else if(isLepEle_[iRun] && !isLowPtAndPfEle_[iRun] && isLowPtEle_[iRun])
	    objectPassed = getCandLepX_LowPtEle(selectedLowPtGsfTracks[iTM].at(j), lepton2TT,
						candLep2Dxy, candLep2Dz,candLep2DxyS, candLep2DzS,
						true, iTM, iRun,
						lowPtGsfTracksHandle, mvaSeedsHandle, gsfLinkLTHandle, gsfLinkPCHandle);
	  
	  else if(!isLepEle_[iRun])
	    objectPassed = getCandLepX_PFmuon(selectedMuons[iTM].at(j), candLep2, lepton2TT,
					      candLep2Dxy, candLep2Dz,candLep2DxyS, candLep2DzS,
					      true, iTM, iRun, PV);
	}
	else if((isLep2PFC || isLep2LT) && !isLepEle_[iRun]){
	  if(isLep2PFC) 
	    objectPassed = getCandX_PackedCand(selectedPackedCandidates[iTM].at(j-leptonNumber), candLep2, lepton2TT,
					       candLep2Dxy, candLep2Dz,candLep2DxyS, candLep2DzS,
					       false, true, false, iTM, iRun);
	  
	  else
	    objectPassed = getCandX_LostTracks(selectedLostTracks[iTM].at(j-leptonNumber-pfCandNumber), candLep2, lepton2TT,
					       candLep2Dxy, candLep2Dz,candLep2DxyS, candLep2DzS,
					       false, true, false, iTM, iRun);

	  //	  if(debug) std::cout << " candLep2->pt() = " << candLep2->pt() << std::endl;

	}
	else{ std::cout << " lepton 2 problem " << std::endl; }      

	if(!objectPassed) continue;
	if(debug) std::cout << " passed lepton 2 " << std::endl;

	float maxl1l2_dxyS = std::max(candLep2DxyS, candLep1DxyS);

	if(debug){
	  std::cout << " 1_vz =  " << vzL1 << " 2_vz = " << vzL2 << std::endl; 
	  //study first then in case if > 1 continue;
	}
	
	bool passedDiLepton = false;
	
	double LepLepLSBS = -1.;
	double LepLepLSBSErr = -1.;
	double LepLepVtx_Chi2 = -1.;
	double LepLepVtx_CL = -1.;
	
	math::XYZVector refitLepLep_2trks;
	RefCountedKinematicVertex refitVertexLepLep_2trks;

	if(debug) std::cout << " before save2TrkRefit_ (= " << save2TrkRefit_ << ")" << std::endl;
	if(save2TrkRefit_){
	  passedDiLepton = LepLepVertexRefitting(lepton1TT,
						 lepton2TT,
						 refitVertexLepLep_2trks,
						 refitLepLep_2trks, iRun);
	  
	  if (passedDiLepton){
	    pair<double,double> LepLepLS = computeLS(refitVertexLepLep_2trks, beamSpot);
	    LepLepLSBS = LepLepLS.first;
	    LepLepLSBSErr = LepLepLS.second;
	    
	    LepLepVtx_Chi2 = (double)refitVertexLepLep_2trks->chiSquared();
	    LepLepVtx_CL = TMath::Prob((double)refitVertexLepLep_2trks->chiSquared(),
				       int(rint(refitVertexLepLep_2trks->degreesOfFreedom())));
	  }
	}
	if(!isLep2PFL && !isLowPtEle_[iRun] && !isLowPtAndPfEle_[iRun] && save2TrkRefit_ &&
	   IPPV_llRefitVtx_min_ != -1 && LepLepLSBS/LepLepLSBSErr < IPPV_llRefitVtx_min_) continue;

	if(save2TrkRefit_ && !passedDiLepton) continue;
	const math::XYZPoint &leplepRefitVertex = (save2TrkRefit_ && passedDiLepton) ? math::XYZPoint(refitVertexLepLep_2trks->position()) : math::XYZPoint(0.,0.,0.);

	if(debug) std::cout << " end save2TrkRefit_ " << std::endl;

	//Kaon
	for (unsigned int k = 0; k < (pfCandNumber+lostChHadrTrackNumber); ++k) {
	  if(debug) std::cout << " i = " << i << " leptonNumber = " << leptonNumber << " j = " << j << " k = " << k << std::endl; 

	  if( ( (isLep1PFC || isLep1LT) && (i-leptonNumber) == k) ||
	      ( (isLep2PFC || isLep2LT) && (j-leptonNumber) == k) ) continue;
	  bool isKPFCand = k<pfCandNumber;

	  const reco::Candidate* kaon;
	  reco::TransientTrack kaonTT;
	  float candKaonDxy = -1;
          float candKaonDz = -1;
          float candKaonDxyS = -1;
          float candKaonDzS = -1;
	  float kaon_dxyFromRefitllVtx = -1;

	  if(isKPFCand)
	    objectPassed = getCandX_PackedCand(selectedPackedCandidates[iTM].at(k), kaon, kaonTT,
					       candKaonDxy, candKaonDz, candKaonDxyS, candKaonDzS,
					       false, false, true, iTM, iRun, (save2TrkRefit_) ? leplepRefitVertex : math::XYZPoint(0,0,0),
					       kaon_dxyFromRefitllVtx);

	  else
	    objectPassed = getCandX_LostTracks(selectedLostTracks[iTM].at(k-pfCandNumber),  kaon, kaonTT,
                                               candKaonDxy, candKaonDz, candKaonDxyS, candKaonDzS,
                                               false, false, true, iTM, iRun, (save2TrkRefit_) ? leplepRefitVertex : math::XYZPoint(0,0,0),
					       kaon_dxyFromRefitllVtx);

	  if(!objectPassed) continue;

	  //keep these here for the moment
	  float maxl1l2k_dxyS = std::max(maxl1l2_dxyS, candKaonDxyS);
	  if(!isLowPtEle_[iRun] && !isLowPtAndPfEle_[iRun] && !isLep2PFL && kll_dxyPV_min_ != -1 && std::abs(maxl1l2k_dxyS) < kll_dxyPV_min_) continue;

	  if(!isLowPtEle_[iRun] && !isLowPtAndPfEle_[iRun] && !isLep2PFL && lepKaon_dz_max_ > -1. &&
	     std::max(std::abs(vzL2 - vzK), std::abs(vzL1 - vzK)) > lepKaon_dz_max_ ) continue;

	  if(debug) std::cout << " passed kaon " << std::endl;	  

	  if(debug){
	    std::cout << " kaon.vz() =  " << vzK << std::endl;
	    //if > 1.5 continue
	  }
	  
	  pair<double,double> DCA = computeDCA(kaonTT, beamSpot);
	  double DCABS = DCA.first;
	  double DCABSErr = DCA.second;
	  
	  if(fabs(DCABS/DCABSErr)<DCASigMinKaon_) continue;
	  

	  //for Kst
	  const reco::Candidate* pion;
	  double candPionDxy = -1;
	  double candPionDz = -1;
	  double candPionDxyS = -1;
	  double candPionDzS = -1;

	  int candPionIndex = -1;

	  bool isPionPFCand = false;
	  double DCABS_pion = -1;
	  double DCABSErr_pion = -1;
	  
	  TLorentzVector refitKstTL;
	  RefCountedKinematicParticle refitKst;
	  RefCountedKinematicVertex refitVertexKst;
	  double KstLSBS = -1;
	  double KstLSBSErr = -1;
	  double KstVtx_Chi2 = -1; 
	  double KstVtx_CL = -1;
	  //double Kst_mass_err;
	  
	  TLorentzVector BToKstllTL;
	  RefCountedKinematicVertex refitVertexBToKstLepLep;
	  math::XYZVector refitBToKstLepLep;
	  math::XYZVector refitKPi;
	  double BToKstLSBS = -1;
	  double BToKstLSBSErr = -1;
	  double BToKstLepLepVtx_Chi2 = -1;
	  double BToKstLepLepVtx_CL = -1;
	  double BToKstcosAlpha = -1;
	  //for K 
	  TLorentzVector BToKllTL;
	  RefCountedKinematicVertex refitVertexBToKLepLep;
	  math::XYZVector refitBToKLepLep;
	  double BToKLSBS = -1;
	  double BToKLSBSErr = -1;
	  double BToKLepLepVtx_Chi2 = -1;
	  double BToKLepLepVtx_CL = -1;
	  double BToKcosAlpha = -1;
	  //
	  
	  if(debug) std::cout << " before pion " << std::endl;	  
	  if(isChKst_[iRun]){
	    for (unsigned int l = 0; l < (pfCandNumber+lostChHadrTrackNumber); ++l) {

	      candPionIndex = l;

	      if( k == l ||
		  ( (isLep1PFC || isLep1LT) && (i-leptonNumber) == l) ||
		  ( (isLep2PFC || isLep2LT) && (j-leptonNumber) == l) ) continue;
	      
	      isPionPFCand = l<pfCandNumber;

	      reco::TransientTrack pionTT;

	      if(isPionPFCand)
		objectPassed = getCandX_PackedCand(selectedPackedCandidates[iTM].at(l), pion, pionTT,
						   candPionDxy, candPionDz, candPionDxyS, candPionDzS,
						   false, false, false, iTM, iRun);

	      else
		objectPassed = getCandX_LostTracks(selectedLostTracks[iTM].at(l-pfCandNumber), pion, pionTT,
						   candPionDxy, candPionDz, candPionDxyS, candPionDzS,
						   false, false, false, iTM, iRun);

	      if(!objectPassed) continue;

	      //for the moment keep these here
	      if(!isLowPtEle_[iRun] && !isLowPtAndPfEle_[iRun] && !isLep2PFL && lepPion_dz_max_ > -1. &&
		 std::max(std::abs(vzL2 - vzP), std::abs(vzL1 - vzP)) > lepPion_dz_max_ ) continue;
	      if(kaonPion_dz_max_ > -1. && std::abs(vzK - vzP) > kaonPion_dz_max_) continue;

	      if(debug){
		std::cout << " pion.vz() =  " << vzP << std::endl;
		//if > 1.5 continue
	      }

	      if(debug) std::cout << " passed pion " << std::endl;
	      
	      pair<double,double> DCA_pion = computeDCA(pionTT, beamSpot);
	      DCABS_pion = DCA_pion.first;
	      DCABSErr_pion = DCA_pion.second;
	      
	      if(fabs(DCABS_pion/DCABSErr_pion)<DCASigMinPion_) continue;
	      	     
	      bool passed = KstVertexRefitting(kaonTT, pionTT,
					       refitVertexKst,
					       refitKst);
	      
	      if(!passed) continue;
	      
	      pair<double,double> KstLS = computeLS(refitVertexKst,beamSpot);
	      KstLSBS = KstLS.first;
	      KstLSBSErr = KstLS.second;
	      
	      KstVtx_Chi2 = (double)refitVertexKst->chiSquared();
	      KstVtx_CL = TMath::Prob((double)refitVertexKst->chiSquared(),
				      int(rint(refitVertexKst->degreesOfFreedom())));
	      
	      //Kst_mass_err = sqrt(refitKst->currentState().kinematicParametersError().matrix()(6,6));

	      passed = BToKstLepLepVertexRefitting(lepton1TT, lepton2TT, refitKst,
						   refitVertexBToKstLepLep,
						   refitKPi, iRun);
	      
	      if (!passed) continue;

	      pair<double,double> BToKstEELS = computeLS(refitVertexBToKstLepLep,beamSpot);
	      BToKstLSBS = BToKstEELS.first;
	      BToKstLSBSErr = BToKstEELS.second;

	      BToKstLepLepVtx_Chi2 = (double)refitVertexBToKstLepLep->chiSquared();
	      BToKstLepLepVtx_CL = TMath::Prob((double)refitVertexBToKstLepLep->chiSquared(),
					       int(rint(refitVertexBToKstLepLep->degreesOfFreedom())));

	      math::XYZVector l1(lepton1.Px(), lepton1.Py(), lepton1.Pz());
	      math::XYZVector l2(lepton2.Px(), lepton2.Py(), lepton2.Pz());
	      refitBToKstLepLep = l1+l2+refitKPi; 
	    
	      BToKstcosAlpha = computeCosAlpha(refitBToKstLepLep,refitVertexBToKstLepLep,beamSpot);
	      
	      if(BToKstLepLepVtx_CL < vtxCL_min_){
		if(debug)std::cout << " bad BToKstLepLepVtx_CL " << std::endl;
		continue;}
	      
	      refitKstTL.SetPtEtaPhiM(sqrt(refitKPi.perp2()), refitKPi.eta(),
				      refitKPi.phi(), refitKst->currentState().mass());
	      BToKstllTL = lepton1+lepton2+refitKstTL;
	      double massKstll = (BToKstllTL).Mag();
	      if(massKstll > Bmass_Kst_max_ || massKstll < Bmass_Kst_min_){
		if(debug)std::cout << " bad massKstll " << std::endl;
		continue;}
	    }
	  }
	  else{
	    if(debug) std::cout << " no pion B case " << std::endl;	  

	    bool passed = BToKLepLepVertexRefitting(lepton1TT, lepton2TT, kaonTT, 
						    refitVertexBToKLepLep, iRun);
	    
	    if (!passed) {
	      // if(save2TrkRefit_ && passedDiLepton) std::cout << " kaonDXY = " << kaon.dxy(leplepRefitVertex)
	      // 						     << " kaonDZ = " << kaon.dz(leplepRefitVertex) << std::endl;
	      continue;
	    }
	    if(debug) std::cout << " >>> BToKLepLepVertexRefitting ok " << std::endl;
	    pair<double,double> BToKLepLepLS = computeLS(refitVertexBToKLepLep,beamSpot);
	    BToKLSBS = BToKLepLepLS.first;
	    BToKLSBSErr = BToKLepLepLS.second;
	    
	    BToKLepLepVtx_Chi2 = (double)refitVertexBToKLepLep->chiSquared();
	    BToKLepLepVtx_CL = TMath::Prob((double)refitVertexBToKLepLep->chiSquared(),
					   int(rint(refitVertexBToKLepLep->degreesOfFreedom())));
	    
	    math::XYZVector l1(lepton1.Px(), lepton1.Py(), lepton1.Pz());
	    math::XYZVector l2(lepton2.Px(), lepton2.Py(), lepton2.Pz());
	    math::XYZVector lk(kaonTL.Px(), kaonTL.Py(), kaonTL.Pz());
	    refitBToKLepLep = l1+l2+lk; 

	    BToKcosAlpha = computeCosAlpha(refitBToKLepLep,refitVertexBToKLepLep,beamSpot);
	    
	    if(debug) std::cout << " BToKLepLepVtx_CL = " << BToKLepLepVtx_CL << std::endl;
	    if(BToKLepLepVtx_CL < vtxCL_min_){ 
	      if(debug)std::cout << " bad BToKLepLepVtx_CL " << std::endl; 
	      continue;}
	    
	    BToKllTL = lepton1+lepton2+kaonTL;
	    double massKll = (BToKllTL).Mag();        
	    if(massKll > Bmass_Kst_max_ || massKll < Bmass_Kst_min_) {
	      if(debug)std::cout << " bad massKll " << std::endl; 
	      continue;}
	  }
	
	  if(debug) std::cout << " result[iRun]->size() = " << result[iRun]->size()
			      << " worstTag_val_idx[iRun].first = " << worstTag_val_idx[iRun].first
			      << " isLep2PFL = " << isLep2PFL << std::endl;

	  if((nSelectedTriplets_ != -1 && result[iRun]->size() >= ((unsigned int)nSelectedTriplets_) ) &&
	     ((isChKst_[iRun] == true && (isLep2PFL ? (BToKstLepLepVtx_CL+2) : BToKstLepLepVtx_CL) < worstTag_val_idx[iRun].first) ||
	      (isChKst_[iRun] == false && (isLep2PFL ? (BToKLepLepVtx_CL+2) : BToKLepLepVtx_CL) < worstTag_val_idx[iRun].first)) ) continue;

	  if(debug) std::cout << " now filling " << std::endl;
	  pat::CompositeCandidate BToKstLLCand;

	  BToKstLLCand.addUserInt("configuration", (int)iRun);
	  BToKstLLCand.addUserInt("isEleCh", (int)isLepEle_[iRun]);
	  BToKstLLCand.addUserInt("isKstCh", (int)isChKst_[iRun]);
	  BToKstLLCand.addUserInt("isLowPtEle", (int)isLowPtEle_[iRun]);
	  BToKstLLCand.addUserInt("isLowPtAndPfEle", (int)isLowPtAndPfEle_[iRun]);
	  BToKstLLCand.addUserInt("isUseLostLeadMuonTracks", (int)useLostLeadMuonTracks_[iRun]);

	  BToKstLLCand.addUserInt("lep1_index", isLep1PFL ? (isLepEle_[iRun] ? selectedElectrons[iTM].at(i) : selectedMuons[iTM].at(i) ) : -1);
	  BToKstLLCand.addUserInt("lep1_lowPt_index", (isLep1LowPt && isLowPtAndPfEle_[iRun]) ? (i - selectedElectrons.size()) : ((isLep1LowPt && isLowPtEle_[iRun]) ? i : -1));

	  BToKstLLCand.addUserInt("lep2_index", isLep2PFL ? (isLepEle_[iRun] ? selectedElectrons[iTM].at(j) : selectedMuons[iTM].at(j) ) : -1);
	  BToKstLLCand.addUserInt("lep2_lowPt_index", (isLep2LowPt && isLowPtAndPfEle_[iRun]) ? (j - selectedElectrons.size()) : ((isLep2LowPt && isLowPtEle_[iRun]) ? j : -1));

	  BToKstLLCand.addUserInt("lep1_pfCand_index", isLep1PFC ? i-leptonNumber : -1);
	  BToKstLLCand.addUserInt("lep2_pfCand_index", isLep2PFC ? j-leptonNumber : -1);

	  BToKstLLCand.addUserInt("lep1_lostTrack_index", isLep1LT ? i-leptonNumber-pfCandNumber : -1);
	  BToKstLLCand.addUserInt("lep2_lostTrack_index", isLep2LT ? j-leptonNumber-pfCandNumber : -1);

	  BToKstLLCand.addUserInt("kaon_index", isKPFCand ? k : -1);
	  BToKstLLCand.addUserInt("kaon_lostTrack_index", isKPFCand ? -1 : k-pfCandNumber);

	  BToKstLLCand.addUserInt("lep1_isLowPt", (int)isLep1LowPt);
	  BToKstLLCand.addUserInt("lep2_isLowPt", (int)isLep2LowPt);
	  BToKstLLCand.addUserInt("lep1_isPFLep", (int)isLep1PFL);
	  BToKstLLCand.addUserInt("lep2_isPFLep", (int)isLep2PFL);
	  BToKstLLCand.addUserInt("lep1_isPFCand", (int)isLep1PFC);
	  BToKstLLCand.addUserInt("lep2_isPFCand", (int)isLep2PFC);
	  BToKstLLCand.addUserInt("lep1_isLT", (int)isLep1LT);
	  BToKstLLCand.addUserInt("lep2_isLT", (int)isLep2LT);
	  BToKstLLCand.addUserInt("kaon_isPFCand", (int)isKPFCand);
	  
	  if(debug) std::cout << " post BDT " << std::endl;

	  if(isChKst_[iRun]){
	    BToKstLLCand.addUserInt("pion_index", isPionPFCand ? candPionIndex : -1);
	    BToKstLLCand.addUserInt("pion_lostTrack_index", isPionPFCand ? -1 : candPionIndex-pfCandNumber);
	    BToKstLLCand.addUserInt("pion_isPFCand", (int)isPionPFCand);
	  }
	  else{
	    BToKstLLCand.addUserInt("pion_index", -1);
	    BToKstLLCand.addUserInt("pion_lostTrack_index", -1);
	    BToKstLLCand.addUserInt("pion_isPFCand", -1);
	  }
	
	  BToKstLLCand.addUserInt("muTrg_index", iTM);

	  BToKstLLCand.addUserFloat("kaon_dxy_wrtllVtx", kaon_dxyFromRefitllVtx);

	  BToKstLLCand.addUserInt("fitLepLep", (int)save2TrkRefit_);
	  BToKstLLCand.addUserInt("fitLepLepPassed", (int)passedDiLepton);
	  BToKstLLCand.addUserFloat("ll_mass", (lepton1+lepton2).Mag());	
	  BToKstLLCand.addUserFloat("ll_pt", (lepton1+lepton2).Pt());
	  BToKstLLCand.addUserFloat("ll_eta", (lepton1+lepton2).Eta());
	  BToKstLLCand.addUserFloat("ll_phi", (lepton1+lepton2).Phi());
	  BToKstLLCand.addUserFloat("ll_Lxy", (passedDiLepton)? (float) LepLepLSBS/LepLepLSBSErr : -1.);
	  BToKstLLCand.addUserFloat("ll_ctxy", (passedDiLepton)? (float) LepLepLSBS/(lepton1+lepton2).Pt() : -1.);
	  BToKstLLCand.addUserFloat("ll_Chi2_vtx", (passedDiLepton)? (float) LepLepVtx_Chi2 : -1.);
	  BToKstLLCand.addUserFloat("ll_CL_vtx", (passedDiLepton)? (float) LepLepVtx_CL : -1.);
	  BToKstLLCand.addUserFloat("maxl1l2_dxyS", maxl1l2_dxyS);

	  //if Kst
	  BToKstLLCand.addUserFloat("Kst_mass", isChKst_[iRun]? refitKstTL.Mag() : -1);
	  BToKstLLCand.addUserFloat("Kst_pt", isChKst_[iRun]? refitKstTL.Pt() : -1);
	  BToKstLLCand.addUserFloat("Kst_eta", isChKst_[iRun]? refitKstTL.Eta() : -1);
	  BToKstLLCand.addUserFloat("Kst_phi", isChKst_[iRun]? refitKstTL.Phi() : -1);
	  BToKstLLCand.addUserFloat("Kst_Lxy", isChKst_[iRun]? (float) KstLSBS/KstLSBSErr : -1);
	  BToKstLLCand.addUserFloat("Kst_ctxy", isChKst_[iRun]? (float) KstLSBS/BToKllTL.Pt() : -1);
	  BToKstLLCand.addUserFloat("Kst_Chi2_vtx", isChKst_[iRun]? (float) KstVtx_Chi2 : -1);
	  BToKstLLCand.addUserFloat("Kst_CL_vtx", isChKst_[iRun]? (float) KstVtx_CL : -1);
	
	  //B mass and co
	  BToKstLLCand.addUserFloat("B_pt", isChKst_[iRun]? BToKstllTL.Pt() : BToKllTL.Pt());
	  BToKstLLCand.addUserFloat("B_eta", isChKst_[iRun]? BToKstllTL.Eta() : BToKllTL.Eta());
	  BToKstLLCand.addUserFloat("B_phi", isChKst_[iRun]? BToKstllTL.Phi() : BToKllTL.Phi());
	  BToKstLLCand.addUserFloat("B_mass", isChKst_[iRun]? BToKstllTL.Mag() : BToKllTL.Mag());
	  
	  BToKstLLCand.addUserFloat("B_Lxy", isChKst_[iRun]? (float)BToKstLSBS/BToKstLSBSErr : (float) BToKLSBS/BToKLSBSErr);
	  BToKstLLCand.addUserFloat("B_ctxy", isChKst_[iRun]? (float)BToKstLSBS/BToKstllTL.Pt() :(float) BToKLSBS/BToKllTL.Pt());
	  BToKstLLCand.addUserFloat("B_Chi2_vtx", isChKst_[iRun]? (float) BToKstLepLepVtx_Chi2 : (float) BToKLepLepVtx_Chi2);
	  BToKstLLCand.addUserFloat("B_CL_vtx", isChKst_[iRun]? (float) BToKstLepLepVtx_CL : (float) BToKLepLepVtx_CL);
	  BToKstLLCand.addUserFloat("B_cosAlpha", isChKst_[iRun]? (float) BToKstcosAlpha : (float) BToKcosAlpha);
	  BToKstLLCand.addUserFloat("maxl1l2k_dxyS", maxl1l2k_dxyS );

	  if(result[iRun]->size() < ((unsigned int)nSelectedTriplets_ )|| nSelectedTriplets_ == -1){
	    result[iRun]->push_back(BToKstLLCand);
	    if(isChKst_[iRun]) resultTag[iRun].push_back(isLep2PFL ? (BToKstLepLepVtx_CL+2) : BToKstLepLepVtx_CL);
	    else resultTag[iRun].push_back(isLep2PFL ? (BToKLepLepVtx_CL+2) : BToKLepLepVtx_CL);

	    if(worstTag_val_idx[iRun].first == 0 || resultTag[iRun][resultTag[iRun].size()-1] < worstTag_val_idx[iRun].first){
	      worstTag_val_idx[iRun].first = resultTag[iRun][resultTag[iRun].size()-1];
	      worstTag_val_idx[iRun].second = resultTag[iRun].size()-1;
	    }
	  }
	  else{
	    if(debug) std::cout << " worst CL = " << worstTag_val_idx[iRun].first << " with idx = " << worstTag_val_idx[iRun].second << std::endl;
	    
	    resultTag[iRun][worstTag_val_idx[iRun].second] = isChKst_[iRun] ? (isLep2PFL ? (BToKstLepLepVtx_CL+2) : BToKstLepLepVtx_CL) :
	      (isLep2PFL ? (BToKLepLepVtx_CL+2) : BToKLepLepVtx_CL);
	    result[iRun]->at(worstTag_val_idx[iRun].second) = BToKstLLCand;

	    if(debug){
	      if(isChKst_[iRun]) std::cout << " post push_back result[iRun]->size() = " << result[iRun]->size()
					   << " BToKstLepLepVtx_CL = " << BToKstLepLepVtx_CL << std::endl;
	      else std::cout << " post push_back result[iRun]->size() = " << result[iRun]->size() << " BToKLepLepVtx_CL = " << BToKLepLepVtx_CL << std::endl;
	    }

	    if(debug) std::cout << " now looking for new worst " << std::endl;
	    float dummyVal = 10.;
	    for(unsigned int ij=0; ij<resultTag[iRun].size(); ++ij){
	      if(resultTag[iRun][ij] < dummyVal){
		dummyVal = resultTag[iRun][ij];
		worstTag_val_idx[iRun].first = dummyVal;
		worstTag_val_idx[iRun].second = ij;
	      }
	      if(debug) std::cout << ij << " resultTag[iRun][ij] = " << resultTag[iRun][ij] << " dummyVal = " << dummyVal
				  << " worstTag_val_idx[iRun].first = " << worstTag_val_idx[iRun].first
				  << " worstTag_val_idx[iRun].second = " << worstTag_val_idx[iRun].second << std::endl;
	    }
	  }


	  // if(save2TrkRefit_ && passedDiLepton) std::cout << " kaonDXY = " << kaon.dxy(leplepRefitVertex)
	  // 						 << " kaonDZ = " << kaon.dz(leplepRefitVertex) << std::endl;

	}//loop over k 
	
      }//subleading lepton
    }// loop over leading lepton
  }//triggering muons

  }//loop over configuration runs

  unsigned resultSize = 0;
  for(unsigned int iRun=0; iRun<isLepEle_.size(); ++iRun)
    resultSize += result[iRun]->size();
  tripletResults->resize(resultSize);
  for(unsigned int iRun=0; iRun<isLepEle_.size(); ++iRun){
    auto startPoint = (iRun == 0)? tripletResults->begin() : tripletResults->begin() + result[iRun-1]->size();
    std::copy(result[iRun]->begin(), result[iRun]->end(), startPoint);

    if(debug)std::cout << "iRun = " << iRun << " resultSize = " << resultSize << " tripletResults->size() = " << tripletResults->size() << std::endl;
  }

  if(debug) std::cout << " resultSize = " << resultSize << " tripletResults->size() = " << tripletResults->size()
		      << " result[0]->size() = " << result[0]->size() << " result[1]->size() = " << result[1]->size() << std::endl;

  iEvent.put(std::move(tripletResults));
  iEvent.put(std::move(resultPF), "PFCand");
  iEvent.put(std::move(resultLT), "LostTracks");
  iEvent.put(std::move(resultGsf), "LowPtGsfTrack");
}



bool BToKstllProducer::LepLepVertexRefitting(const reco::TransientTrack &lep1TT,
					     const reco::TransientTrack &lep2TT,
					     RefCountedKinematicVertex &refitVertex,
					     math::XYZVector &refitLepLep, unsigned int iRun){
    
  const reco::TransientTrack l1TT = lep1TT;
  const reco::TransientTrack l2TT = lep2TT;
    
  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;
  //KinematicConstrainedVertexFitter PartVtxFitter;
    
  std::vector<RefCountedKinematicParticle> Particles;
  double chi = 0.;
  double ndf = 0.;
  Particles.push_back(partFactory.particle(l1TT, lep1Mass_[iRun],chi,ndf, lep1MassErr_[iRun]));
  Particles.push_back(partFactory.particle(l2TT, lep2Mass_[iRun],chi,ndf, lep2MassErr_[iRun]));
  RefCountedKinematicTree leplepVertexFitTree = PartVtxFitter.fit(Particles);
    
  if ( !leplepVertexFitTree->isValid()) return false;
    
  leplepVertexFitTree->movePointerToTheTop();
  refitVertex = leplepVertexFitTree->currentDecayVertex();

  RefCountedKinematicParticle refitParticle = leplepVertexFitTree->currentParticle();
    
  if ( !refitVertex->vertexIsValid()) return false;

  refitLepLep = refitParticle->refittedTransientTrack().track().momentum();
   
  return true;    
}


bool BToKstllProducer::KstVertexRefitting(const reco::TransientTrack &kaonTT,
					  const reco::TransientTrack &pionTT,
					  RefCountedKinematicVertex &refitVertex,
                                          RefCountedKinematicParticle &refitKst){

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;
  //KinematicConstrainedVertexFitter PartVtxFitter;

  std::vector<RefCountedKinematicParticle> KstParticles;
  double chi = 0.;
  double ndf = 0.;
  KstParticles.push_back(partFactory.particle(kaonTT,KaonMass_,chi,ndf,KaonMassErr_));
  KstParticles.push_back(partFactory.particle(pionTT,PionMass_,chi,ndf,PionMassErr_));
  RefCountedKinematicTree KstVertexFitTree = PartVtxFitter.fit(KstParticles);

  if ( !KstVertexFitTree->isValid()) return false;

  KstVertexFitTree->movePointerToTheTop();
  refitVertex = KstVertexFitTree->currentDecayVertex();
  refitKst = KstVertexFitTree->currentParticle();

  if ( !refitVertex->vertexIsValid()) return false;

  return true;
}

bool BToKstllProducer::BToKstLepLepVertexRefitting(const reco::TransientTrack &lep1TT,
						   const reco::TransientTrack &lep2TT,
						   const RefCountedKinematicParticle &refitKst,
						   RefCountedKinematicVertex &refitVertex,
						   math::XYZVector &refitKPi, unsigned int iRun){

  const reco::TransientTrack l1TT = lep1TT;
  const reco::TransientTrack l2TT = lep2TT;
  const reco::TransientTrack KPiTT = refitKst->refittedTransientTrack();


  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;
  //KinematicConstrainedVertexFitter PartVtxFitter;
  
  float Kst_mass = refitKst->currentState().mass();
  float Kst_mass_err = sqrt(refitKst->currentState().kinematicParametersError().matrix()(6,6));
  if(KstMassConstraint_ > 0){
    Kst_mass = KstMassConstraint_;
    Kst_mass_err = KstMassErr_;
  }
  
  std::vector<RefCountedKinematicParticle> Particles;
  double chi = 0.;
  double ndf = 0.;
  Particles.push_back(partFactory.particle(l1TT,lep1Mass_[iRun],chi,ndf,lep1MassErr_[iRun]));
  Particles.push_back(partFactory.particle(l2TT,lep2Mass_[iRun],chi,ndf,lep2MassErr_[iRun]));
  Particles.push_back(partFactory.particle(KPiTT,Kst_mass,chi,ndf,Kst_mass_err));

  RefCountedKinematicTree VertexFitTree = PartVtxFitter.fit(Particles);

  if ( !VertexFitTree->isValid()) return false;

  VertexFitTree->movePointerToTheTop();
  refitVertex = VertexFitTree->currentDecayVertex();
  if ( !refitVertex->vertexIsValid()) return false;

  refitKPi = refitKst->refittedTransientTrack().track().momentum();

  return true;
}



bool BToKstllProducer::BToKLepLepVertexRefitting(const reco::TransientTrack &lep1TT,
						 const reco::TransientTrack &lep2TT,
						 const reco::TransientTrack &kaonTT,
						 RefCountedKinematicVertex &refitVertex,
						 unsigned int iRun){

  const reco::TransientTrack l1TT = lep1TT;
  const reco::TransientTrack l2TT = lep2TT;
  const reco::TransientTrack kTT = kaonTT;

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;
  //KinematicConstrainedVertexFitter PartVtxFitter;

  std::vector<RefCountedKinematicParticle> Particles;
  double chi = 0.;
  double ndf = 0.;
  Particles.push_back(partFactory.particle(l1TT,lep1Mass_[iRun],chi,ndf,lep1MassErr_[iRun]));
  Particles.push_back(partFactory.particle(l2TT,lep2Mass_[iRun],chi,ndf,lep2MassErr_[iRun]));
  Particles.push_back(partFactory.particle(kTT,KaonMass_,chi,ndf,KaonMassErr_));

  RefCountedKinematicTree VertexFitTree = PartVtxFitter.fit(Particles);

  if ( !VertexFitTree->isValid()) return false;

  VertexFitTree->movePointerToTheTop();
  refitVertex = VertexFitTree->currentDecayVertex();
  if ( !refitVertex->vertexIsValid()) return false;

  return true;
}

pair<double,double> BToKstllProducer::computeLS(RefCountedKinematicVertex refitVertex,
                                                reco::BeamSpot beamSpot){

  TVector v(2);
  v[0] = refitVertex->position().x()-beamSpot.position().x();
  v[1] = refitVertex->position().y()-beamSpot.position().y();

  TMatrix errVtx(2,2);
  errVtx(0,0) = refitVertex->error().cxx();
  errVtx(0,1) = refitVertex->error().matrix()(0,1);
  errVtx(1,0) = errVtx(0,1);
  errVtx(1,1) = refitVertex->error().cyy();

  TMatrix errBS(2,2);
  errBS(0,0) = beamSpot.covariance()(0,0);
  errBS(0,1) = beamSpot.covariance()(0,1);
  errBS(1,0) = beamSpot.covariance()(1,0);
  errBS(1,1) = beamSpot.covariance()(1,1);

  double LSBS = sqrt(v.Norm2Sqr());
  double LSBSErr = sqrt( v*(errVtx*v) + v*(errBS*v) ) / LSBS;

  pair<double,double> LS = make_pair(LSBS,LSBSErr);

  return LS;
}


pair<double,double> BToKstllProducer::computeLS(RefCountedKinematicVertex refitVertex,
                                                reco::Vertex pvertex){

  TVector v(2);
  v[0] = refitVertex->position().x()-pvertex.position().x();
  v[1] = refitVertex->position().y()-pvertex.position().y();

  TMatrix errVtx(2,2);
  errVtx(0,0) = refitVertex->error().cxx();
  errVtx(0,1) = refitVertex->error().matrix()(0,1);
  errVtx(1,0) = errVtx(0,1);
  errVtx(1,1) = refitVertex->error().cyy();

  TMatrix errBS(2,2);
  errBS(0,0) = pvertex.covariance()(0,0);
  errBS(0,1) = pvertex.covariance()(0,1);
  errBS(1,0) = pvertex.covariance()(1,0);
  errBS(1,1) = pvertex.covariance()(1,1);
    
  double LSBS = sqrt(v.Norm2Sqr());
  double LSBSErr = sqrt( v*(errVtx*v) + v*(errBS*v) ) / LSBS;
    
  pair<double,double> LS = make_pair(LSBS,LSBSErr);

  return LS;
}


double BToKstllProducer::computeCosAlpha(math::XYZVector& refitBToMLepLep,
					 RefCountedKinematicVertex refitVertex,
					 reco::BeamSpot beamSpot){
  
  TVector v(2);
  v[0] = refitVertex->position().x()-beamSpot.position().x();
  v[1] = refitVertex->position().y()-beamSpot.position().y();

  TVector w(2);
  w[0] = refitBToMLepLep.x();
  w[1] = refitBToMLepLep.y();

  double cosAlpha = v*w/sqrt(v.Norm2Sqr()*w.Norm2Sqr());
  return cosAlpha;
}


pair<double,double> BToKstllProducer::computeDCA(const reco::TransientTrack &hadronTT,
                                                 reco::BeamSpot beamSpot){

  TrajectoryStateClosestToPoint theDCAXBS = hadronTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );  
  
  double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
  double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
  pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
  return DCA;
}




DEFINE_FWK_MODULE(BToKstllProducer); 
