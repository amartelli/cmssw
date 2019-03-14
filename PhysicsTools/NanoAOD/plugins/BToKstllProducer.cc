#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Common/interface/RefToPtr.h"
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
    
  virtual void produce(edm::Event&, const edm::EventSetup&);
    
  bool LepLepVertexRefitting(const reco::TransientTrack &lep1TT,
			     const reco::TransientTrack &lep2TT,
                             RefCountedKinematicVertex &refitVertex,
                             math::XYZVector &refitLepLep);
    
  bool KstVertexRefitting(const reco::TransientTrack &kaonTT,
			  const reco::TransientTrack &pionTT,
			  RefCountedKinematicVertex &refitVertex,
			  RefCountedKinematicParticle &refitKst);			 

  bool BToKLepLepVertexRefitting(const reco::TransientTrack &lep1TT,
				 const reco::TransientTrack &lep2TT,
				 const reco::TransientTrack &kaonTT,
                                 RefCountedKinematicVertex &refitVertex);

  bool BToKstLepLepVertexRefitting(const reco::TransientTrack &lep1TT,
				   const reco::TransientTrack &lep2TT,
				   const RefCountedKinematicParticle &refitKst,
				   RefCountedKinematicVertex &refitVertex,
				   math::XYZVector &refitKPi);
    
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
  edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> PFCandSrc_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> lostSubLeadLepTrackSrc_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> lostChHadrTrackSrc_;
  
  int nSelectedTriplets_;

  bool isLepEle_;
  bool isChKst_;

  double ptMinLeadLep_;
  double etaMaxLeadLep_;
  
  double ptMinSubLeadLep_;
  double etaMaxSubLeadLep_;

  double ptMinKaon_;
  double etaMaxKaon_;
  double DCASigMinKaon_;

  double ptMinPion_;
  double etaMaxPion_;
  double DCASigMinPion_;

  bool diLepCharge_;
  bool KstCharge_;
  double JPsiMassConstraint_;
  double KstMassConstraint_;
  bool save2TrkRefit_;

  bool useLostSubLeadLepTracks_;
  bool useLostChHadrTracks_;

  double vtxCL_min_;
  double Bmass_min_;
  double Bmass_max_;
  double Bmass_Kst_min_;
  double Bmass_Kst_max_;

  double diLepton_dz_max_;
  double lepKaon_dz_max_;
  double lepPion_dz_max_;
  double kaonPion_dz_max_;
  double kaonRefitllVertex_dxy_max_;
  double kll_dxyPV_min_;
  double IPPV_llRefitVtx_min_;

  float ElectronMass_ = 0.5109989e-3;
  float ElectronMassErr_ = 3.1*1e-12;
  float MuonMass_ = 0.10565837;
  float MuonMassErr_ = 3.5*1e-9;
  float KaonMass_ = 0.493677;
  float KaonMassErr_ = 1.6e-5;
  float PionMass_ = 0.139570;
  float PionMassErr_ = 3.5e-7;
  
  float lep1Mass_; 
  float lep2Mass_; 
  float lep1MassErr_; 
  float lep2MassErr_; 


  //float JPsiMass_ = 3.096916;  //Configurable parameter
  float JPsiMassErr_ = 0.011;
  
  //float KstMass_ = 0.89176;  //Configurable parameter  
  float KstMassErr_ = 0.25e-3;

  bool debug;
};



BToKstllProducer::BToKstllProducer(const edm::ParameterSet &iConfig):
  beamSpotSrc_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
  vertexSrc_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ),
  electronSrc_( consumes<std::vector<pat::Electron>> ( iConfig.getParameter<edm::InputTag>( "electronCollection" ) ) ),
  muonSrc_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
  PFCandSrc_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
  lostSubLeadLepTrackSrc_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "lostSubLeadLepTrackCollection" ) ) ),
  lostChHadrTrackSrc_( consumes<edm::View<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "lostChHadrTrackCollection" ) ) ),
  nSelectedTriplets_( iConfig.getParameter<int>( "nSelectedTriplets" ) ),
  isLepEle_( iConfig.getParameter<bool>( "isLeptonElectron" ) ),
  isChKst_( iConfig.getParameter<bool>( "isChannelKst" ) ),
  ptMinLeadLep_( iConfig.getParameter<double>( (isLepEle_ == true) ? "LeadEleMinPt" : "LeadMuonMinPt" ) ),
  etaMaxLeadLep_( iConfig.getParameter<double>( (isLepEle_ == true) ? "LeadEleMaxEta" : "LeadMuonMaxEta" ) ),
  ptMinSubLeadLep_( iConfig.getParameter<double>( (isLepEle_ == true) ? "SubLeadEleMinPt" : "SubLeadMuonMinPt" ) ),
  etaMaxSubLeadLep_( iConfig.getParameter<double>( (isLepEle_ == true) ? "SubLeadEleMaxEta" : "SubLeadMuonMaxEta" ) ),
  ptMinKaon_( iConfig.getParameter<double>( "KaonMinPt" ) ),
  etaMaxKaon_( iConfig.getParameter<double>( "KaonMaxEta" ) ),
  DCASigMinKaon_( iConfig.getParameter<double>( "KaonMinDCASig" ) ),
  ptMinPion_( iConfig.getParameter<double>( "PionMinPt" ) ),
  etaMaxPion_( iConfig.getParameter<double>( "PionMaxEta" ) ),
  DCASigMinPion_( iConfig.getParameter<double>( "PionMinDCASig" ) ),
  diLepCharge_( iConfig.getParameter<bool>( "DiLeptonChargeCheck" ) ),
  KstCharge_( iConfig.getParameter<bool>( "KstarChargeCheck" ) ),

  JPsiMassConstraint_( iConfig.getParameter<double>( "JPsiMassConstraint" ) ),
  KstMassConstraint_( iConfig.getParameter<double>( "KstMassConstraint" ) ),

  save2TrkRefit_( iConfig.getParameter<bool>( "save2TrackRefit" ) ),
//  save4TrkRefit_( iConfig.getParameter<bool>( "save4TrackRefit" ) ),
  useLostSubLeadLepTracks_( iConfig.getParameter<bool>( "useLostSubLeadLepTracks" ) ),
  useLostChHadrTracks_( iConfig.getParameter<bool>( "useLostChHadrTracks" ) ),
  vtxCL_min_( iConfig.getParameter<double>( "vtxCL_min" ) ),
  Bmass_min_( iConfig.getParameter<double>( "Bmass_min" ) ),
  Bmass_max_( iConfig.getParameter<double>( "Bmass_max" ) ),
  Bmass_Kst_min_( iConfig.getParameter<double>( "Bmass_Kst_min" ) ),
  Bmass_Kst_max_( iConfig.getParameter<double>( "Bmass_Kst_max" ) ),
  diLepton_dz_max_( iConfig.getParameter<double>( "diLepton_dz_max" ) ),
  lepKaon_dz_max_( iConfig.getParameter<double>( "lepKaon_dz_max" ) ),
  lepPion_dz_max_( iConfig.getParameter<double>( "lepPion_dz_max" ) ),
  kaonPion_dz_max_( iConfig.getParameter<double>( "kaonPion_dz_max" ) ),
  kaonRefitllVertex_dxy_max_( iConfig.getParameter<double>( "kaonRefitllVertex_dxy_max" ) ),
  kll_dxyPV_min_( iConfig.getParameter<double>( "kll_dxyPV_min" ) ),
  IPPV_llRefitVtx_min_( iConfig.getParameter<double>( "IPPV_llRefitVtx_min" ) )
{
  lep1Mass_ = (isLepEle_) ? ElectronMass_ : MuonMass_;
  lep2Mass_ = lep1Mass_;
  lep1MassErr_ = (isLepEle_) ? ElectronMassErr_ : MuonMassErr_;
  lep2MassErr_ = lep1MassErr_;

  produces<pat::CompositeCandidateCollection>();
  //typedef std::vector<CompositeCandidate> CompositeCandidateCollection

  debug = false;
}


void BToKstllProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
  edm::ESHandle<MagneticField> bFieldHandle;
  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  edm::Handle<reco::VertexCollection> vertexHandle;
    
  iEvent.getByToken(beamSpotSrc_, beamSpotHandle);
    
  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("BToKstllProducer") << "No beam spot available from EventSetup" ;
  }
  
  reco::BeamSpot beamSpot = *beamSpotHandle;
  
  iEvent.getByToken(vertexSrc_, vertexHandle);
  const reco::Vertex & PV = vertexHandle->front();
  
  edm::Handle<std::vector<pat::Electron>> electronHandle;  
  edm::Handle<std::vector<pat::Muon>> muonHandle;
  edm::Handle<edm::View<pat::PackedCandidate>> pfCandHandle;
  edm::Handle<edm::View<pat::PackedCandidate>> lostSubLeadLepTrackHandle;
  edm::Handle<edm::View<pat::PackedCandidate>> lostChHadrTrackHandle;
  
  if(isLepEle_) iEvent.getByToken(electronSrc_, electronHandle);
  else iEvent.getByToken(muonSrc_, muonHandle);
  iEvent.getByToken(PFCandSrc_, pfCandHandle);
  if(useLostSubLeadLepTracks_) iEvent.getByToken(lostSubLeadLepTrackSrc_, lostSubLeadLepTrackHandle);
  if(useLostChHadrTracks_) iEvent.getByToken(lostChHadrTrackSrc_, lostChHadrTrackHandle);
  
  unsigned int leptonNumber = (isLepEle_ == true) ? electronHandle->size() : muonHandle->size();   
  unsigned int pfCandNumber = pfCandHandle->size();
  unsigned int subLeadLeptonTrackNumber = useLostSubLeadLepTracks_ ? (leptonNumber + pfCandNumber + lostSubLeadLepTrackHandle->size()) : leptonNumber;
  unsigned int lostChHadrTrackNumber = useLostChHadrTracks_ ? lostChHadrTrackHandle->size() : 0;
    
  // Output collection
  std::unique_ptr<pat::CompositeCandidateCollection> result( new pat::CompositeCandidateCollection );
  std::vector<float> resultTag;
  std::pair<float, int> worstTag_val_idx;

  std::map<std::pair<edm::ProductID,int>, unsigned int> checkLeptonsDuplicate;
  //std::vector<std::pair<edm::ProductID,int>> removableL2;

  if(debug){
    std::cout << " leptonNumber = " << leptonNumber << " pfCandNumber = " << pfCandNumber
	      << " subLeadLeptonTrackNumber = " << subLeadLeptonTrackNumber << " lostChHadrTrackNumber = " << lostChHadrTrackNumber << std::endl; 
    std::cout << " isLepEle_ = " << isLepEle_ << " isChKst_ = " << isChKst_ << std::endl; 
  }

  if(leptonNumber>0){

    // loop on leading lepton from pat::Lepton 
    // subleading from pat::Lepton or PFcandidate+LT
    // and then kaon (+pion) for K or Kstar
    // to build triplets
    for (unsigned int i = 0; i < leptonNumber; ++i) {
      checkLeptonsDuplicate.clear();

      float candLep1Dxy;
      float candLep1Dz;      
      float candLep1DxyS;
      float candLep1DzS;
      const reco::Candidate* candLep1;
      TLorentzVector lepton1;
      reco::TransientTrack lepton1TT; 
      int lepton1Charge;
      float lepton1VZ;
      if(isLepEle_){
	const pat::Electron & ele1 = (*electronHandle)[i];
	candLep1 = &ele1;

	if(debug)
	  std::cout << " ele1.sourceCandidatePtr(i).isNonnull() = " << ele1.sourceCandidatePtr(i).isNonnull()
		    << " ele1.sourceCandidatePtr(i).isAvailable() = " << ele1.sourceCandidatePtr(i).isAvailable() << std::endl;

	//just save ref a PFCand
	for(unsigned int ic = 0; ic < ele1.numberOfSourceCandidatePtrs(); ++ic){
	  reco::CandidatePtr dummyCandLep1 = ele1.sourceCandidatePtr(ic);
	  if(!ele1.sourceCandidatePtr(ic).isNonnull() || !ele1.sourceCandidatePtr(ic).isAvailable()) continue;
	  checkLeptonsDuplicate[std::pair<edm::ProductID,int>(dummyCandLep1.id(), dummyCandLep1.key())] = i;
	}

	//could implement ele ID criteria here

	lepton1TT = theTTBuilder->build(ele1.gsfTrack()); // it is build from the reco::Track - use buildfromGSF otherwise
	candLep1Dxy = ele1.dB(pat::Electron::PV2D);
	candLep1Dz = ele1.dB(pat::Electron::PVDZ);
	candLep1DxyS = candLep1Dxy / ele1.edB(pat::Electron::PV2D);
	candLep1DzS = candLep1Dz / ele1.edB(pat::Electron::PVDZ);
      }
      else{
	const pat::Muon & muon1 = (*muonHandle)[i];
	candLep1 = &muon1;

	if(debug) std::cout << " muon1.sourceCandidatePtr(i).isNonnull() = " << muon1.sourceCandidatePtr(i).isNonnull() << std::endl;

	//just save ref a PFCand
	for(unsigned int ic = 0; ic < muon1.numberOfSourceCandidatePtrs(); ++ic){
	  if(!muon1.sourceCandidatePtr(ic).isNonnull() || !muon1.sourceCandidatePtr(ic).isAvailable()) continue;
	  reco::CandidatePtr dummyCandLep1 = muon1.sourceCandidatePtr(ic);
	  checkLeptonsDuplicate[std::pair<edm::ProductID,int>(dummyCandLep1.id(), dummyCandLep1.key())] = i;
	}

	//have muon softID on leading muon !!
	if(!(muon1.isLooseMuon() && muon1.isSoftMuon(PV))) continue;
	lepton1TT = theTTBuilder->build(muon1.bestTrack()); // it is build from the reco::Track - generalTracks
	candLep1Dxy = muon1.dB(pat::Muon::PV2D);
	candLep1Dz = muon1.dB(pat::Muon::PVDZ);
	candLep1DxyS = candLep1Dxy/muon1.edB(pat::Muon::PV2D);
	candLep1DzS = candLep1Dz/muon1.edB(pat::Muon::PVDZ);
      }

      lepton1.SetPtEtaPhiM(candLep1->pt(), candLep1->eta(), candLep1->phi(), (isLepEle_) ? ElectronMass_ : MuonMass_);
      lepton1Charge = candLep1->charge();
      lepton1VZ = candLep1->vz();
      if(lepton1.Pt()<ptMinLeadLep_ || abs(lepton1.Eta()) > etaMaxLeadLep_) continue;
      
      if(debug) std::cout << " passed lepton 1 " << std::endl;

      //removableL2.clear();

      for (unsigned int j = i+1; j < subLeadLeptonTrackNumber; ++j) {
	bool isLep2PFL = (j < leptonNumber || useLostSubLeadLepTracks_ == false) ? true : false;
	bool isLep2PFC = (isLep2PFL == false && j < (pfCandNumber+leptonNumber)) ? true : false;
	bool isLep2LT = (j >= (pfCandNumber+leptonNumber)) ? true : false;

	if(debug) std::cout << " i = " << i << " j = " << j << " leptonNumber = " << leptonNumber <<std::endl;

	if(debug) std::cout << " isLep2PFL = " << isLep2PFL << " isLep2PFC = " << isLep2PFC << " isLep2LT = " << isLep2LT << std::endl;

	//cross clean also wrt 2nd lepton
	/*
	if(removableL2.size() != 0){
	  if(debug) std::cout << " pre clear checkLeptonsDuplicate.size = " << checkLeptonsDuplicate.size() << " removableL2 size = " << removableL2.size() << std::endl;

	  for(auto iLD : removableL2)
	    checkLeptonsDuplicate.erase(iLD);
	  removableL2.clear();

	  if(debug)
	    std::cout << " post clear checkLeptonsDuplicate.size = " << checkLeptonsDuplicate.size() << " removableL2 size = " << removableL2.size() << std::endl;
	}
	*/

	float candLep2Dxy;
	float candLep2Dz;
	float candLep2DxyS;
	float candLep2DzS;
	const reco::Candidate* candLep2;
	TLorentzVector lepton2;
	reco::TransientTrack lepton2TT; 
	int lepton2Charge;      
	float lepton2VZ;


	if(isLep2PFL){
	  if(i == j) continue;
	  if(isLepEle_){
	    const pat::Electron & ele2 = (*electronHandle)[j];
	    candLep2 = &ele2;

	    //just save ref a PFCand
	    for(unsigned int ic = 0; ic < ele2.numberOfSourceCandidatePtrs(); ++ic){
	      if(!ele2.sourceCandidatePtr(ic).isNonnull() || !ele2.sourceCandidatePtr(ic).isAvailable()) continue;
	      reco::CandidatePtr dummyCandLep2 = ele2.sourceCandidatePtr(ic);
	      std::pair<edm::ProductID,int> keyVal(dummyCandLep2.id(), dummyCandLep2.key());
	      checkLeptonsDuplicate[keyVal] = j;
	      //removableL2.push_back(keyVal);
	    }

	    lepton2TT = theTTBuilder->build(ele2.gsfTrack()); // it is build from the reco::Track - use buildfromGSF otherwise
	    candLep2Dxy = ele2.dB(pat::Electron::PV2D);
	    candLep2Dz = ele2.dB(pat::Electron::PVDZ);
	    candLep2DxyS = candLep2Dxy/ele2.edB(pat::Electron::PV2D);
	    candLep2DzS = candLep2Dz/ele2.edB(pat::Electron::PVDZ);
	  }
	  else{
	    const pat::Muon & muon2 = (*muonHandle)[j];
	    candLep2 = &muon2;

	    //just save ref a PFCand
	    for(unsigned int ic = 0; ic < muon2.numberOfSourceCandidatePtrs(); ++ic){
	      if(!muon2.sourceCandidatePtr(ic).isNonnull() || !muon2.sourceCandidatePtr(ic).isAvailable()) continue;
	      reco::CandidatePtr dummyCandLep2 = muon2.sourceCandidatePtr(ic);
	      std::pair<edm::ProductID,int> keyVal(dummyCandLep2.id(), dummyCandLep2.key());
              checkLeptonsDuplicate[keyVal] = j;
	      //removableL2.push_back(keyVal);
	    }

	    //could implement muon ID criteria here

	    lepton2TT = theTTBuilder->build(muon2.bestTrack()); // it is build from the reco::Track - generalTracks
	    candLep2Dxy = muon2.dB(pat::Muon::PV2D);
	    candLep2Dz = muon2.dB(pat::Muon::PVDZ);
	    candLep2DxyS = candLep2Dxy/muon2.edB(pat::Muon::PV2D);
	    candLep2DzS = candLep2Dz/muon2.edB(pat::Muon::PVDZ);
	  }
	}
	else if(isLep2PFC || isLep2LT){

	  edm::ProductID pair1 = isLep2PFC ? pfCandHandle.id() : lostSubLeadLepTrackHandle.id();
	  int pair2 = isLep2PFC ? (j-leptonNumber) : (j-leptonNumber-pfCandNumber);
	  std::pair<edm::ProductID, int> dummyPair(pair1, pair2);
	  std::map<std::pair<edm::ProductID,int>, unsigned int>::iterator checkLD_it = checkLeptonsDuplicate.find(dummyPair);
	  if(checkLD_it != checkLeptonsDuplicate.end()) { 
	    if(debug) std::cout << " found duplicate L2" << std::endl; 
	    continue;
	  }
	  if(debug) std::cout << " NON duplicate L2 " << std::endl;

	  if(isLepEle_){
	    const pat::PackedCandidate & ele2 = isLep2PFC ? (*pfCandHandle)[j-leptonNumber] : (*lostSubLeadLepTrackHandle)[j-leptonNumber-pfCandNumber];
	    candLep2 = &ele2;

	    //could implement ele ID criteria here
	    if(isLep2LT && debug) std::cout << " >>> ele2.hasTrackDetails() = " << ele2.hasTrackDetails() << " ele2.pdgId() = " << ele2.pdgId()  << std::endl;
	    if(!ele2.hasTrackDetails()) continue;
	    //exclude neutral should be safe do not ask too much ID
	    if(abs(ele2.pdgId()) == 0) continue;
	    //FIXME
	    if(ele2.pt() > 3.) continue;

	    lepton2TT = theTTBuilder->build(ele2.bestTrack());
	    candLep2Dxy = ele2.dxy();
	    candLep2Dz = ele2.dz();
	    candLep2DxyS = candLep2Dxy/ele2.dxyError();
	    candLep2DzS = candLep2Dz/ele2.dzError();
	  }
	  else{
	    const pat::PackedCandidate & muon2 = isLep2PFC ? (*pfCandHandle)[j-leptonNumber] : (*lostSubLeadLepTrackHandle)[j-leptonNumber-pfCandNumber];
	    candLep2 = &muon2;

	    if(debug)	    std::cout << " muon2 taken " << std::endl;
	    //could implement muon ID criteria here
	    if(isLep2LT && debug) std::cout << " >>> muon2.hasTrackDetails() = " << muon2.hasTrackDetails() << " muon2.pdgId() = " << muon2.pdgId()  << std::endl;
	    if(!muon2.hasTrackDetails()) continue;
	    //exclude neutral should be safe do not ask too much ID
	    if(abs(muon2.pdgId()) == 0) continue;   
	    //muonID for pT > 3 in EB => recover with pfCand
	    //FIXME
	    //if(std::abs(muon2.eta()) > 1.497 || muon2.pt() > 3.) continue;
	    if(muon2.pt() > 3.) continue;
	    if(debug) std::cout << " muon2 ok for candLep2 " << std::endl;

	    lepton2TT = theTTBuilder->build(muon2.bestTrack()); 
	    candLep2Dxy = muon2.dxy();
	    candLep2Dz = muon2.dz();
	    candLep2DxyS = candLep2Dxy/muon2.dxyError();
	    candLep2DzS = candLep2Dz/muon2.dzError();
	  }
	}
	else{
	  std::cout << " ERROR: not assigned subleading lepton " << std::endl;
	  return;
	}
	lepton2.SetPtEtaPhiM(candLep2->pt(), candLep2->eta(), candLep2->phi(), (isLepEle_) ? ElectronMass_ : MuonMass_);
	lepton2Charge = candLep2->charge();
	lepton2VZ = candLep2->vz();

	if(lepton1.Pt() < lepton2.Pt()) continue; //Lepton 1 is always saved as the leading one
	if(lepton2.Pt() < ptMinSubLeadLep_ || abs(lepton2.Eta()) > etaMaxSubLeadLep_) continue;
	if(debug && lepton1Charge != lepton2Charge) std::cout << " lepton1Charge = " << lepton1Charge << " lepton2Charge = " << lepton2Charge << std::endl;
	if(diLepCharge_ && lepton1Charge*lepton2Charge > 0) continue;
	// lepton1 and lepton2 belong to different collections need to check they are different
	if(!isLep2PFL && deltaR(lepton1.Eta(), lepton1.Phi(), lepton2.Eta(), lepton2.Phi()) < 0.01) continue;
	if(debug) std::cout << " passed lepton 2 " << std::endl;

	if(diLepton_dz_max_ > -1. && std::abs(lepton2VZ - lepton1VZ) > diLepton_dz_max_) continue;

	float maxl1l2_dxyS = std::max(candLep2DxyS, candLep1DxyS);

	if(debug){
	  std::cout << " 1_vz =  " << lepton1VZ << " 2_vz = " << lepton2VZ << std::endl; 
	  //study first then in case if > 1 continue;
	}
	
	bool passedDiLepton = false;
	
	double LepLepLSBS = -1.;
	double LepLepLSBSErr = -1.;
	double LepLepVtx_Chi2 = -1.;
	double LepLepVtx_CL = -1.;
	
	math::XYZVector refitLepLep_2trks;
	RefCountedKinematicVertex refitVertexLepLep_2trks;

	if(debug) std::cout << " before save2TrkRefit_ " << std::endl;
	if(save2TrkRefit_){
	  passedDiLepton = LepLepVertexRefitting(lepton1TT,
						 lepton2TT,
						 refitVertexLepLep_2trks,
						 refitLepLep_2trks);
	  
	  if (passedDiLepton){
	    pair<double,double> LepLepLS = computeLS(refitVertexLepLep_2trks, PV);
	    LepLepLSBS = LepLepLS.first;
	    LepLepLSBSErr = LepLepLS.second;
	    
	    LepLepVtx_Chi2 = (double)refitVertexLepLep_2trks->chiSquared();
	    LepLepVtx_CL = TMath::Prob((double)refitVertexLepLep_2trks->chiSquared(),
				       int(rint(refitVertexLepLep_2trks->degreesOfFreedom())));
	  }
	}
	if(save2TrkRefit_ && IPPV_llRefitVtx_min_ != -1 && LepLepLSBS/LepLepLSBSErr < IPPV_llRefitVtx_min_) continue;

	if(save2TrkRefit_ && !passedDiLepton) continue;
	const math::XYZPoint &leplepRefitVertex = (save2TrkRefit_ && passedDiLepton) ? math::XYZPoint(refitVertexLepLep_2trks->position()) : math::XYZPoint(0.,0.,0.);

	if(debug) std::cout << " end save2TrkRefit_ " << std::endl;

	//Kaon
	for (unsigned int k = 0; k < (pfCandNumber+lostChHadrTrackNumber); ++k) {
	  if(debug) std::cout << " i = " << i << " leptonNumber = " << leptonNumber << " j = " << j << " k = " << k << std::endl; 

	  if(!isLep2PFL && (j-leptonNumber) == k) continue;
	  bool isKPFCand = k<pfCandNumber;
	  const pat::PackedCandidate & kaon = isKPFCand ? (*pfCandHandle)[k] : (*lostChHadrTrackHandle)[k-pfCandNumber];
	  
	  edm::ProductID pair1k = isKPFCand ? pfCandHandle.id() : lostChHadrTrackHandle.id();
          int pair2k = isKPFCand ? k : (k-pfCandNumber);
	  std::pair<edm::ProductID, int> dummyPairk(pair1k, pair2k);
	  std::map<std::pair<edm::ProductID,int>, unsigned int>::iterator checkLD_itk = checkLeptonsDuplicate.find(dummyPairk);
          if(checkLD_itk != checkLeptonsDuplicate.end()){ 
	    if(debug) std::cout << " found duplicate k " << std::endl; 
	    continue;
	  }
	  if(debug) std::cout << " NON duplicate K " << std::endl;

	  if(!kaon.hasTrackDetails()) continue;
	  if(abs(kaon.pdgId()) != 211) // && abs(kaon.pdgId()) != 11 && abs(kaon.pdgId()) != 13) continue; //Charged hadrons
	    if(kaon.pt()<ptMinKaon_ || abs(kaon.eta())>etaMaxKaon_) continue;
	  if(deltaR(lepton1.Eta(), lepton1.Phi(), kaon.eta(), kaon.phi()) < 0.01 || 
	     deltaR(lepton2.Eta(), lepton2.Phi(), kaon.eta(), kaon.phi()) < 0.01 ) continue;

	  float kaon_dxyFromRefitllVtx = save2TrkRefit_ ? kaon.dxy(leplepRefitVertex) : -1;
	  if(save2TrkRefit_ && kaonRefitllVertex_dxy_max_ != -1 &&
	     std::abs(kaon_dxyFromRefitllVtx) > kaonRefitllVertex_dxy_max_) continue;

	  float kaon_dxyS = kaon.dxy()/kaon.dxyError();
	  float maxl1l2k_dxyS = std::max(maxl1l2_dxyS, kaon_dxyS);
	  if(kll_dxyPV_min_ != -1 && std::abs(maxl1l2k_dxyS) < kll_dxyPV_min_) continue;

	  if(debug) std::cout << " passed kaon " << std::endl;	  

	  if(lepKaon_dz_max_ > -1. && 
	     std::max(std::abs(lepton2VZ - kaon.vz()), std::abs(lepton1VZ - kaon.vz())) > lepKaon_dz_max_ ) continue;

	  if(debug){
	    std::cout << " kaon.vz() =  " << kaon.vz() << std::endl;
	    //if > 1.5 continue
	  }
	  
	  TLorentzVector kaonTL;
	  kaonTL.SetPtEtaPhiM(kaon.pt(), kaon.eta(), kaon.phi(), KaonMass_); 
	  const reco::TransientTrack kaonTT = theTTBuilder->build(kaon.bestTrack());
	  pair<double,double> DCA = computeDCA(kaonTT, beamSpot);
	  double DCABS = DCA.first;
	  double DCABSErr = DCA.second;
	  
	  if(fabs(DCABS/DCABSErr)<DCASigMinKaon_) continue;
	  
	  //for Kst
	  int candPionL = -1;
	  double candPionDxy = -1;
	  double candPionDz = -1;
	  double candPionDxyS = -1;
	  double candPionDzS = -1;
	  const reco::Candidate* candPion;
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
	  if(isChKst_){
	    for (unsigned int l = 0; l < (pfCandNumber+lostChHadrTrackNumber); ++l) {
	      candPionL = l;
	      if(!isLep2PFL && (k == l || (j-leptonNumber) == l)) continue;
	      
	      isPionPFCand = l<pfCandNumber;

	      edm::ProductID pair1p = isPionPFCand ? pfCandHandle.id() : lostChHadrTrackHandle.id();
	      int pair2p = isPionPFCand ? l : (l-pfCandNumber);
	      std::pair<edm::ProductID, int> dummyPairp(pair1p, pair2p);
	      std::map<std::pair<edm::ProductID,int>, unsigned int>::iterator checkLD_itp = checkLeptonsDuplicate.find(dummyPairp);
	      if(checkLD_itp != checkLeptonsDuplicate.end()) { 
		if(debug) std::cout << " found duplicate pi " << std::endl; 
		continue;
	      }
	      if(debug) std::cout << " NON duplicate Pi " << std::endl;

	      const pat::PackedCandidate & pion = isPionPFCand ? (*pfCandHandle)[l] : (*lostChHadrTrackHandle)[l-pfCandNumber];
	      candPion = &pion;
	      candPionDxy = pion.dxy();
	      candPionDz = pion.dz();
	      candPionDxyS = pion.dxyError();
	      candPionDzS = pion.dzError();

	      if(!pion.hasTrackDetails()) continue;
	      if(abs(pion.pdgId())!=211) continue; //Charged hadrons
	      if(pion.pt() < ptMinPion_ || abs(pion.eta()) > etaMaxPion_) continue;
	      if(KstCharge_ && kaon.charge()*pion.charge()>0) continue;
	      if(deltaR(lepton1.Eta(), lepton1.Phi(), pion.eta(), pion.phi()) < 0.01 || 
		 deltaR(lepton2.Eta(), lepton2.Phi(), pion.eta(), pion.phi()) < 0.01 ||
		 deltaR(kaon, pion) < 0.01) continue;
	      
	      if(lepPion_dz_max_ > -1. &&
		 std::max(std::abs(lepton2VZ - pion.vz()), std::abs(lepton1VZ - pion.vz())) > lepPion_dz_max_ ) continue;
	      if(kaonPion_dz_max_ > -1. && std::abs(kaon.vz() - pion.vz()) > kaonPion_dz_max_) continue;

	      if(debug){
		std::cout << " pion.vz() =  " << pion.vz() << std::endl;
		//if > 1.5 continue
	      }

	      if(debug) std::cout << " passed pion " << std::endl;
	      
	      const reco::TransientTrack pionTT = theTTBuilder->build(pion.bestTrack());
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
						   refitKPi);
	      
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
						    refitVertexBToKLepLep);
	    
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
	    
	    if(BToKLepLepVtx_CL < vtxCL_min_){ 
	      if(debug)std::cout << " bad BToKLepLepVtx_CL " << std::endl; 
	      continue;}
	    
	    BToKllTL = lepton1+lepton2+kaonTL;
	    double massKll = (BToKllTL).Mag();        
	    if(massKll > Bmass_Kst_max_ || massKll < Bmass_Kst_min_) {
	      if(debug)std::cout << " bad massKll " << std::endl; 
	      continue;}
	  }
	
	  if(debug) std::cout << " result->size() = " << result->size()
			      << " worstTag_val_idx.first = " << worstTag_val_idx.first
			      << " isLep2PFL = " << isLep2PFL << std::endl;

	  if((nSelectedTriplets_ != -1 && result->size() >= (unsigned int)nSelectedTriplets_) &&
	     ((isChKst_ == true && (isLep2PFL ? (BToKstLepLepVtx_CL+2) : BToKstLepLepVtx_CL) < worstTag_val_idx.first) ||
	      (isChKst_ == false && (isLep2PFL ? (BToKLepLepVtx_CL+2) : BToKLepLepVtx_CL) < worstTag_val_idx.first)) ) continue;

	  if(debug) std::cout << " now filling " << std::endl;
	  pat::CompositeCandidate BToKstLLCand;
	  BToKstLLCand.addDaughter( *(candLep1) , "lep1");
	  BToKstLLCand.addDaughter( *(candLep2) , "lep2");
	  BToKstLLCand.addDaughter( kaon, "kaon");
	  
	  BToKstLLCand.addUserInt("isEleCh", (int)isLepEle_);
	  BToKstLLCand.addUserInt("isKstCh", (int)isChKst_);
	  
	  BToKstLLCand.addUserInt("lep1_index", i);
	  BToKstLLCand.addUserInt("lep2_index", isLep2PFL ? j : -1);
	  BToKstLLCand.addUserInt("lep2_pfCand_index", isLep2PFC ? j-leptonNumber : -1);
	  BToKstLLCand.addUserInt("lep2_lostTrack_index", isLep2LT ? j-leptonNumber-pfCandNumber : -1);
	  BToKstLLCand.addUserInt("kaon_index", isKPFCand ? k : -1);
	  BToKstLLCand.addUserInt("kaon_lostTrack_index", isKPFCand ? -1 : k-pfCandNumber);

	  BToKstLLCand.addUserInt("lep2_isPFLep", (int)isLep2PFL);
	  BToKstLLCand.addUserInt("lep2_isPFCand", (int)isLep2PFC);
	  BToKstLLCand.addUserInt("kaon_isPFCand", (int)isKPFCand);
	  
	  if(isChKst_){
	    BToKstLLCand.addDaughter( *(candPion), "pion");	      
	    BToKstLLCand.addUserInt("pion_index", isPionPFCand ? candPionL : -1);
	    BToKstLLCand.addUserInt("pion_lostTrack_index", isPionPFCand ? -1 : candPionL-pfCandNumber);
	    BToKstLLCand.addUserInt("pion_isPFCand", (int)isPionPFCand);
	  }
	  else{
	    BToKstLLCand.addUserInt("pion_index", -1);
	    BToKstLLCand.addUserInt("pion_lostTrack_index", -1);
	    BToKstLLCand.addUserInt("pion_isPFCand", -1);
	  }
	  
	  BToKstLLCand.addUserFloat("lep1_pt",     lepton1.Pt());
	  BToKstLLCand.addUserFloat("lep1_eta",    lepton1.Eta());
	  BToKstLLCand.addUserFloat("lep1_phi",    lepton1.Phi());
	  BToKstLLCand.addUserInt("lep1_charge",   candLep1->charge());
	  BToKstLLCand.addUserFloat("lep1_dxy",    candLep1Dxy);
	  BToKstLLCand.addUserFloat("lep1_dxyS",   candLep1DxyS);
	  BToKstLLCand.addUserFloat("lep1_dz",     candLep1Dz);
	  BToKstLLCand.addUserFloat("lep1_dzS",    candLep1DzS);
	  BToKstLLCand.addUserFloat("lep1_vz",     candLep1->vz());

	  BToKstLLCand.addUserFloat("lep2_pt",     lepton2.Pt());
	  BToKstLLCand.addUserFloat("lep2_eta",    lepton2.Eta());
	  BToKstLLCand.addUserFloat("lep2_phi",    lepton2.Phi());
	  BToKstLLCand.addUserInt("lep2_charge",   candLep2->charge());
	  BToKstLLCand.addUserFloat("lep2_dxy",    candLep2Dxy);
	  BToKstLLCand.addUserFloat("lep2_dxyS",   candLep2DxyS);
	  BToKstLLCand.addUserFloat("lep2_dz",     candLep2Dz);
	  BToKstLLCand.addUserFloat("lep2_dzS",    candLep2DzS);
	  BToKstLLCand.addUserFloat("lep2_vz",     candLep2->vz());

	  BToKstLLCand.addUserFloat("kaon_pt",    kaonTL.Pt());
	  BToKstLLCand.addUserFloat("kaon_eta",   kaonTL.Eta());
	  BToKstLLCand.addUserFloat("kaon_phi",   kaonTL.Phi());
	  BToKstLLCand.addUserInt("kaon_charge",  kaon.charge());
	  BToKstLLCand.addUserFloat("kaon_DCASig", DCABS/DCABSErr);
	  // BToKstLLCand.addUserFloat("kaon_dxy",    kaon.dxy());
	  BToKstLLCand.addUserFloat("kaon_dxyS",   kaon_dxyS);
	  // BToKstLLCand.addUserFloat("kaon_dz",    kaon.dz());
	  BToKstLLCand.addUserFloat("kaon_dzS",   kaon.dz()/kaon.dzError());
	  BToKstLLCand.addUserFloat("kaon_vz",    kaon.vz());
	  BToKstLLCand.addUserFloat("kaon_dxy_wrtllVtx", kaon_dxyFromRefitllVtx);

	  //if Kst
	  BToKstLLCand.addUserFloat("pion_pt", isChKst_ ? candPion->pt() : -1);
	  BToKstLLCand.addUserFloat("pion_eta", isChKst_ ? candPion->eta() : -1);
	  BToKstLLCand.addUserFloat("pion_phi", isChKst_ ? candPion->phi() : -1);
	  BToKstLLCand.addUserInt("pion_charge", isChKst_ ? candPion->charge() : -1);
	  BToKstLLCand.addUserFloat("pion_DCASig", isChKst_ ? DCABS_pion/DCABSErr_pion : -1);
	  BToKstLLCand.addUserFloat("pion_dxy", isChKst_ ? candPionDxy : -1);
	  BToKstLLCand.addUserFloat("pion_dxyS", isChKst_ ? candPionDxyS : -1);
	  BToKstLLCand.addUserFloat("pion_dz", isChKst_ ? candPionDz : -1);
	  BToKstLLCand.addUserFloat("pion_dzS", isChKst_ ? candPionDzS : -1);
	  BToKstLLCand.addUserFloat("pion_vz", isChKst_ ? candPion->vz() : -1);

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
	  BToKstLLCand.addUserFloat("Kst_mass", isChKst_? refitKstTL.Mag() : -1);
	  BToKstLLCand.addUserFloat("Kst_pt", isChKst_? refitKstTL.Pt() : -1);
	  BToKstLLCand.addUserFloat("Kst_eta", isChKst_? refitKstTL.Eta() : -1);
	  BToKstLLCand.addUserFloat("Kst_phi", isChKst_? refitKstTL.Phi() : -1);
	  BToKstLLCand.addUserFloat("Kst_Lxy", isChKst_? (float) KstLSBS/KstLSBSErr : -1);
	  BToKstLLCand.addUserFloat("Kst_ctxy", isChKst_? (float) KstLSBS/BToKllTL.Pt() : -1);
	  BToKstLLCand.addUserFloat("Kst_Chi2_vtx", isChKst_? (float) KstVtx_Chi2 : -1);
	  BToKstLLCand.addUserFloat("Kst_CL_vtx", isChKst_? (float) KstVtx_CL : -1);
	
	  //B mass and co
	  BToKstLLCand.addUserFloat("B_pt", isChKst_? BToKstllTL.Pt() : BToKllTL.Pt());
	  BToKstLLCand.addUserFloat("B_eta", isChKst_? BToKstllTL.Eta() : BToKllTL.Eta());
	  BToKstLLCand.addUserFloat("B_phi", isChKst_? BToKstllTL.Phi() : BToKllTL.Phi());
	  BToKstLLCand.addUserFloat("B_mass", isChKst_? BToKstllTL.Mag() : BToKllTL.Mag());
	  
	  BToKstLLCand.addUserFloat("B_Lxy", isChKst_? (float)BToKstLSBS/BToKstLSBSErr : (float) BToKLSBS/BToKLSBSErr);
	  BToKstLLCand.addUserFloat("B_ctxy", isChKst_? (float)BToKstLSBS/BToKstllTL.Pt() :(float) BToKLSBS/BToKllTL.Pt());
	  BToKstLLCand.addUserFloat("B_Chi2_vtx", isChKst_? (float) BToKstLepLepVtx_Chi2 : (float) BToKLepLepVtx_Chi2);
	  BToKstLLCand.addUserFloat("B_CL_vtx", isChKst_? (float) BToKstLepLepVtx_CL : (float) BToKLepLepVtx_CL);
	  BToKstLLCand.addUserFloat("B_cosAlpha", isChKst_? (float) BToKstcosAlpha : (float) BToKcosAlpha);
	  BToKstLLCand.addUserFloat("maxl1l2k_dxyS", maxl1l2k_dxyS );

	  if(result->size() < (unsigned int)nSelectedTriplets_ || nSelectedTriplets_ == -1){
	    result->push_back(BToKstLLCand);
	    if(isChKst_) resultTag.push_back(isLep2PFL ? (BToKstLepLepVtx_CL+2) : BToKstLepLepVtx_CL);
	    else resultTag.push_back(isLep2PFL ? (BToKLepLepVtx_CL+2) : BToKLepLepVtx_CL);
	  }
	  else{
	    float dummyVal = 10.;
	    for(unsigned int ij=0; ij<resultTag.size(); ++ij){
	      if(resultTag[ij] < dummyVal){
		dummyVal = resultTag[ij];
		worstTag_val_idx.first = dummyVal;
		worstTag_val_idx.second = ij;
	      }
	    }
	    if(debug) std::cout << " worst CL = " << worstTag_val_idx.first << " with idx = " << worstTag_val_idx.second << std::endl;
	    resultTag[worstTag_val_idx.second] = isChKst_ ? (isLep2PFL ? (BToKstLepLepVtx_CL+2) : BToKstLepLepVtx_CL) :
	      (isLep2PFL ? (BToKLepLepVtx_CL+2) : BToKLepLepVtx_CL);
	    result->at(worstTag_val_idx.second) = BToKstLLCand;
	  }

	  if(debug){
	    if(isChKst_) std::cout << " post push_back result->size() = " << result->size() << " BToKstLepLepVtx_CL = " << BToKstLepLepVtx_CL << std::endl;
	    else std::cout << " post push_back result->size() = " << result->size() << " BToKLepLepVtx_CL = " << BToKLepLepVtx_CL << std::endl;
	  }

	  // if(save2TrkRefit_ && passedDiLepton) std::cout << " kaonDXY = " << kaon.dxy(leplepRefitVertex)
	  // 						 << " kaonDZ = " << kaon.dz(leplepRefitVertex) << std::endl;

	}//loop over k 
	
      }//subleading lepton
    }// loop over leading lepton
    
  }//if leptons > 1

  iEvent.put(std::move(result));
}



bool BToKstllProducer::LepLepVertexRefitting(const reco::TransientTrack &lep1TT,
					     const reco::TransientTrack &lep2TT,
					     RefCountedKinematicVertex &refitVertex,
					     math::XYZVector &refitLepLep){
    
  const reco::TransientTrack l1TT = lep1TT;
  const reco::TransientTrack l2TT = lep2TT;
    
  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;
  //KinematicConstrainedVertexFitter PartVtxFitter;
    
  std::vector<RefCountedKinematicParticle> Particles;
  double chi = 0.;
  double ndf = 0.;
  Particles.push_back(partFactory.particle(l1TT, lep1Mass_,chi,ndf, lep1MassErr_));
  Particles.push_back(partFactory.particle(l2TT, lep2Mass_,chi,ndf, lep2MassErr_));
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
						   math::XYZVector &refitKPi){

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
  Particles.push_back(partFactory.particle(l1TT,lep1Mass_,chi,ndf,lep1MassErr_));
  Particles.push_back(partFactory.particle(l2TT,lep2Mass_,chi,ndf,lep2MassErr_));
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
						 RefCountedKinematicVertex &refitVertex){

  const reco::TransientTrack l1TT = lep1TT;
  const reco::TransientTrack l2TT = lep2TT;
  const reco::TransientTrack kTT = kaonTT;

  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;
  //KinematicConstrainedVertexFitter PartVtxFitter;

  std::vector<RefCountedKinematicParticle> Particles;
  double chi = 0.;
  double ndf = 0.;
  Particles.push_back(partFactory.particle(l1TT,lep1Mass_,chi,ndf,lep1MassErr_));
  Particles.push_back(partFactory.particle(l2TT,lep2Mass_,chi,ndf,lep2MassErr_));
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
