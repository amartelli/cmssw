// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"



#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"



#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"



#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

//#include "RecoNtuples/HGCalHitCalibration/interface/AEvent.h"
//#include "RecoNtuples/HGCalHitCalibration/interface/AObData.h"

#include <string>
#include <map>

class HGCalHitCalibration : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
//
// constructors and destructor
//
  explicit HGCalHitCalibration(const edm::ParameterSet&);
  ~HGCalHitCalibration();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

 // ----------member data ---------------------------

  edm::EDGetTokenT<HGCRecHitCollection> _recHitsEE;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsFH;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsBH;
  edm::EDGetTokenT<std::vector<TrackingVertex> > _vtx;
  edm::EDGetTokenT<std::vector<TrackingParticle> > _part;
  edm::EDGetTokenT<std::vector<SimCluster> > _simClusters;
  edm::EDGetTokenT<std::vector<CaloParticle> > _caloParticles;

  std::string                detector;
  int                        algo;
  HGCalDepthPreClusterer     pre;
  bool                       rawRecHits;
  hgcal::RecHitTools         recHitTools;


  TH1F* h_Vtx_x;
  TH1F* h_Vtx_y;
  TH1F* h_Vtx_z;

  TH1F* h_Vtx_dvx;
  TH1F* h_Vtx_dvy;
  TH1F* h_Vtx_dvz;

  TH2F* h_EoP_CPene_vs_frWithTime;
  TH2F* h_maxHitWithTime_vs_layer;
  TH2F* h_fracmaxHitWithTime_vs_layer;
  TH2F* h_maxHit_vs_layer;

  TH1F* h_frWithTime_100;
  TH1F* h_frWithTime_200; 
  TH1F* h_frWithTime_300;

  TH1F* h_EoP_CPene_100_calib_fraction;
  TH1F* h_EoP_CPene_100_calib;

  TH1F* h_EoP_CPene_200_calib_fraction;
  TH1F* h_EoP_CPene_200_calib;

  TH1F* h_EoP_CPene_300_calib_fraction;
  TH1F* h_EoP_CPene_300_calib;


  TH1F* LayerOccupancy;

  std::vector<float> Energy_layer_calib;
  std::vector<float> Energy_layer_calib_fraction;

  std::vector<int> nHits_layer;
  std::vector<int> nHitsWithTime_layer;
};



HGCalHitCalibration::HGCalHitCalibration(const edm::ParameterSet& iConfig) :
  detector(iConfig.getParameter<std::string >("detector")),
  rawRecHits(iConfig.getParameter<bool>("rawRecHits"))
{
  //now do what ever initialization is needed
  usesResource("TFileService");

  if(detector=="all") {
    _recHitsEE = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"));
    _recHitsFH = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits"));
    _recHitsBH = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEBRecHits"));
    algo = 1;
  }else if(detector=="EM") {
    _recHitsEE = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCEERecHits"));
    algo = 2;
  }else if(detector=="HAD") {
    _recHitsFH = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEFRecHits"));
    _recHitsBH = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit","HGCHEBRecHits"));
    algo = 3;
  }
  _vtx = consumes<std::vector<TrackingVertex> >(edm::InputTag("mix","MergedTrackTruth"));
  _part = consumes<std::vector<TrackingParticle> >(edm::InputTag("mix","MergedTrackTruth"));
  _simClusters = consumes<std::vector<SimCluster> >(edm::InputTag("mix","MergedCaloTruth"));
  _caloParticles = consumes<std::vector<CaloParticle> >(edm::InputTag("mix","MergedCaloTruth"));

  /*
  agpc = new AGenPartCollection();
  arhc = new ARecHitCollection();
  if (rawRecHits)
    arhc_raw = new ARecHitCollection();
  ascc = new ASimClusterCollection();
  acpc = new ACaloParticleCollection();
  event = new AEvent();
  */

  edm::Service<TFileService> fs;
  h_EoP_CPene_vs_frWithTime = fs->make<TH2F>("h_EoP_CPene_vs_frWithTime", "",  1000, 0., 1., 1000, -0.5, 2.5);

  h_maxHitWithTime_vs_layer = fs->make<TH2F>("h_maxHitWithTime_vs_layer", "",  60, 0., 60., 1000., 0., 500.);
  h_fracmaxHitWithTime_vs_layer = fs->make<TH2F>("h_fracmaxHitWithTime_vs_layer", "",   60, 0., 60., 1000, 0., 1.);
  h_maxHit_vs_layer = fs->make<TH2F>("h_maxHit_vs_layer", "",  60, 0., 60., 1000., 0., 500.);

  h_frWithTime_100 = fs->make<TH1F>("h_frWithTime_100", "", 1000, 0., 1.);
  h_frWithTime_200 = fs->make<TH1F>("h_frWithTime_100", "", 1000, 0., 1.);
  h_frWithTime_300 = fs->make<TH1F>("h_frWithTime_100", "", 1000, 0., 1.);

  h_EoP_CPene_100_calib_fraction = fs->make<TH1F>("h_EoP_CPene_100_calib_fraction", "", 1000, -0.5, 2.5);
  h_EoP_CPene_100_calib = fs->make<TH1F>("h_EoP_CPene_100_calib", "", 1000, -0.5, 2.5);

  h_EoP_CPene_200_calib_fraction = fs->make<TH1F>("h_EoP_CPene_200_calib_fraction", "", 1000, -0.5, 2.5);
  h_EoP_CPene_200_calib = fs->make<TH1F>("h_EoP_CPene_200_calib", "", 1000, -0.5, 2.5);

  h_EoP_CPene_300_calib_fraction = fs->make<TH1F>("h_EoP_CPene_300_calib_fraction", "", 1000, -0.5, 2.5);
  h_EoP_CPene_300_calib = fs->make<TH1F>("h_EoP_CPene_300_calib", "", 1000, -0.5, 2.5);

  h_Vtx_x = fs->make<TH1F>("h_Vtx_x", "", 1000, -10., 10.);
  h_Vtx_y = fs->make<TH1F>("h_Vtx_y", "", 1000, -10., 10.);
  h_Vtx_z = fs->make<TH1F>("h_Vtx_z", "", 1000, -10., 10.);

  h_Vtx_dvx = fs->make<TH1F>("h_Vtx_dvx", "", 1000, -10., 10.);
  h_Vtx_dvy = fs->make<TH1F>("h_Vtx_dvy", "", 1000, -10., 10.);
  h_Vtx_dvz = fs->make<TH1F>("h_Vtx_dvz", "", 1000, -10., 10.);

  LayerOccupancy = fs->make<TH1F>("LayerOccupancy", "", 60, 0., 60.);
}

HGCalHitCalibration::~HGCalHitCalibration()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


void
HGCalHitCalibration::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  /*
  agpc->clear();
  arhc->clear();
  arhc_raw->clear();
  acdc->clear();
  amcc->clear();
  ascc->clear();
  apfcc->clear();
  acpc->clear();
  */

  Energy_layer_calib.clear();
  Energy_layer_calib_fraction.clear();
  nHits_layer.clear();
  nHitsWithTime_layer.clear();
  for(unsigned int ij=0; ij<60; ++ij){
    Energy_layer_calib.push_back(0.);
    Energy_layer_calib_fraction.push_back(0.);
    nHits_layer.push_back(0.);
    nHitsWithTime_layer.push_back(0.);
  }

  recHitTools.getEventSetup(iSetup);

  //  int npart = 0;
  /*
  int nhit  = 0;
  int nhit_raw = 0;
  int nsimclus = 0;
  int ncalopart = 0;
  */

  Handle<HGCRecHitCollection> recHitHandleEE;
  Handle<HGCRecHitCollection> recHitHandleFH;
  Handle<HGCRecHitCollection> recHitHandleBH;

  Handle<std::vector<TrackingVertex> > vtxHandle;
  Handle<std::vector<TrackingParticle> > partHandle;
  iEvent.getByToken(_vtx,vtxHandle);
  iEvent.getByToken(_part,partHandle);
  const std::vector<TrackingVertex>& vtxs = *vtxHandle;
  const std::vector<TrackingParticle>& part = *partHandle;

  Handle<std::vector<SimCluster> > simClusterHandle;
  Handle<std::vector<CaloParticle> > caloParticleHandle;

  iEvent.getByToken(_simClusters, simClusterHandle);
  iEvent.getByToken(_caloParticles, caloParticleHandle);

  const std::vector<SimCluster>& simClusters = *simClusterHandle;
  const std::vector<CaloParticle>& caloParticles = *caloParticleHandle;

  float vx = 0.;
  float vy = 0.;
  float vz = 0.;
  if(vtxs.size()!=0){
    vx = vtxs[0].position().x();
    vy = vtxs[0].position().y();
    vz = vtxs[0].position().z();
  }

  h_Vtx_x->Fill(vx);
  h_Vtx_y->Fill(vy);
  h_Vtx_z->Fill(vz);

  // TODO: should fall back to beam spot if no vertex
  //  npart = part.size();
  for(unsigned int i=0;i<part.size();++i){
    if(part[i].parentVertex()->nGenVertices()>0){
      float dvx=0.;
      float dvy=0.;
      float dvz=0.;
      if(part[i].decayVertices().size()==1){
	 dvx=part[i].decayVertices()[0]->position().x();
	 dvy=part[i].decayVertices()[0]->position().y();
	 dvz=part[i].decayVertices()[0]->position().z();

	 h_Vtx_dvx->Fill(dvx);
	 h_Vtx_dvy->Fill(dvy);
	 h_Vtx_dvz->Fill(dvz);
      }
      //      agpc->push_back(AGenPart(part[i].eta(),part[i].phi(),part[i].pt(),part[i].energy(),dvx,dvy,dvz,part[i].pdgId()));
    }
  }
  //make a map detid-rechit
  std::map<DetId,const HGCRecHit*> hitmap;
  switch(algo){
  case 1:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      iEvent.getByToken(_recHitsFH,recHitHandleFH);
      iEvent.getByToken(_recHitsBH,recHitHandleBH);
      const auto& rechitsEE = *recHitHandleEE;
      const auto& rechitsFH = *recHitHandleFH;
      const auto& rechitsBH = *recHitHandleBH;
      for(unsigned int i = 0; i < rechitsEE.size(); ++i){
	hitmap[rechitsEE[i].detid()] = &rechitsEE[i];
      }
      for(unsigned int i = 0; i < rechitsFH.size(); ++i){
	hitmap[rechitsFH[i].detid()] = &rechitsFH[i];
      }
      for(unsigned int i = 0; i < rechitsBH.size(); ++i){
	hitmap[rechitsBH[i].detid()] = &rechitsBH[i];
      }
      break;
    }
  case 2:
    {
      iEvent.getByToken(_recHitsEE,recHitHandleEE);
      const HGCRecHitCollection& rechitsEE = *recHitHandleEE;
      for(unsigned int i = 0; i < rechitsEE.size(); i++){
	hitmap[rechitsEE[i].detid()] = &rechitsEE[i];
      }
      break;
    }
  case 3:
    {
      iEvent.getByToken(_recHitsFH,recHitHandleFH);
      iEvent.getByToken(_recHitsBH,recHitHandleBH);
      const auto& rechitsFH = *recHitHandleFH;
      const auto& rechitsBH = *recHitHandleBH;
      for(unsigned int i = 0; i < rechitsFH.size(); i++){
	hitmap[rechitsFH[i].detid()] = &rechitsFH[i];
      }
      for(unsigned int i = 0; i < rechitsBH.size(); i++){
	hitmap[rechitsBH[i].detid()] = &rechitsBH[i];
      }
      break;
    }
  default:
    break;
  }

  // loop over caloParticles
  for (std::vector<CaloParticle>::const_iterator it_caloPart = caloParticles.begin(); it_caloPart != caloParticles.end(); ++it_caloPart) {
    //    if(it_caloPart->eta() < 0) continue;
    //    ++ncalopart;
    const SimClusterRefVector simClusterRefVector = it_caloPart->simClusters();
    //    std::vector<uint32_t> simClusterIndex;


    Energy_layer_calib.clear();
    Energy_layer_calib_fraction.clear();
    nHits_layer.clear();
    nHitsWithTime_layer.clear();
    for(unsigned int ij=0; ij<60; ++ij) {
      Energy_layer_calib.push_back(0.);
      Energy_layer_calib_fraction.push_back(0.);
      nHits_layer.push_back(0.);
      nHitsWithTime_layer.push_back(0.);
    }

    int totHits = 0;
    int totHitsWithTime = 0;

    //    int seedCalo = -1; // 0 == EE 1 == EH 2 == BH
    int seedDet = 0;
    float seedEnergy = 0.;
    //    float seedLayer = 0.;

    int simClusterCount = 0;
    //    float hitEnergy = 0;
    //    float hitFraction = 0;

    for (CaloParticle::sc_iterator it_sc = simClusterRefVector.begin(); it_sc != simClusterRefVector.end(); ++it_sc) {
      const SimCluster simCluster = (*(*it_sc));
      ++simClusterCount;
      //      std::cout << ">>> simCluster.energy() = " << simCluster.energy() << std::endl;
      //      float hitSumEnergy = 0;
      const std::vector<std::pair<uint32_t,float> > hits_and_fractions = simCluster.hits_and_fractions();

      //loop over hits      
      for (std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) {
	unsigned int hitlayer = recHitTools.getLayerWithOffset(it_haf->first);
	DetId hitid = (it_haf->first); 
	
	//	bool isTValid = false;
	float hitTime = 0.;

        // dump raw RecHits and match
	if (rawRecHits) {
	  //	  if (algo < 3) {
	  if (hitid.det() == DetId::Forward && hitid.subdetId() == HGCEE) {
	    const HGCRecHitCollection& rechitsEE = *recHitHandleEE;
	    //	    std::cout << " >>> in EE " << std::endl;
	    // loop over EE RecHits
	    for (HGCRecHitCollection::const_iterator it_hit = rechitsEE.begin(); it_hit < rechitsEE.end(); ++it_hit) {
	      // int clusterIndex = -1;
	      // int flags = 0x0;
	      
	      const HGCalDetId detid = it_hit->detid();
	      unsigned int layer = recHitTools.getLayerWithOffset(detid);
	      if(detid == hitid){
		//	std::cout << " >>> trovato " << std::endl;
		if(layer != hitlayer) std::cout << " >>> problem " << std::endl;
		Energy_layer_calib[layer] += it_hit->energy();
		Energy_layer_calib_fraction[layer] += it_hit->energy()*it_haf->second;
		hitTime = it_hit->time();
		//		isTValid = it_hit->isTimeValid();
		//		hitEnergy = it_hit->energy();
		//		hitSumEnergy += it_hit->energy();
		//		hitFraction = it_haf->second;

		LayerOccupancy->Fill(layer);
		if(seedEnergy < it_hit->energy()){
		  seedEnergy = it_hit->energy();
		  seedDet = recHitTools.getSiThickness(detid);
		  //		  seedCalo = 0;
		  //		  seedLayer = layer;
		}
		break;
	      }
	    }
	  }
	  //     if (algo != 2) {
	  if (hitid.det() == DetId::Forward && hitid.subdetId() == HGCHEF) {
	    const HGCRecHitCollection& rechitsFH = *recHitHandleFH;
	    //	    std::cout << " >>> in EH " << std::endl;
	    // loop over HEF RecHits
	    for (HGCRecHitCollection::const_iterator it_hit = rechitsFH.begin(); it_hit < rechitsFH.end(); ++it_hit) {
	      // int clusterIndex = -1;
	      // int flags = 0x0;
	      
	      const HGCalDetId detid = it_hit->detid();
	      unsigned int layer = recHitTools.getLayerWithOffset(detid);

	      if(detid == hitid){
		//	std::cout << " >>> trovato " << std::endl;
		if(layer != hitlayer) std::cout << " >>> problem " << std::endl;
		Energy_layer_calib[layer] += it_hit->energy();
		Energy_layer_calib_fraction[layer] += it_hit->energy()*it_haf->second;
		hitTime = it_hit->time();
		//		isTValid = it_hit->isTimeValid();

		//		hitSumEnergy += it_hit->energy();
		//		hitEnergy = it_hit->energy();
		//		hitFraction = it_haf->second;

		LayerOccupancy->Fill(layer);
		if(seedEnergy < it_hit->energy()){
		  seedEnergy = it_hit->energy();
		  seedDet = recHitTools.getSiThickness(detid);
		  //		  seedCalo = 1;
		  //		  seedLayer = layer;
		}
		break;
	      }
	    }
	  }
	  //	  std::cout << " >>> layer = " << hitlayer << " hitEnergy = " << hitEnergy << " fraction = " << hitFraction << std::endl;
	  // for(unsigned int iC=0; iC<Energy_layer_calib.size(); ++iC){
	  //   std::cout << " >>> Energy_layer_calib.at(iC)" << Energy_layer_calib.at(iC) << std::endl;
	  // }


	  
	  if (hitid.det() == DetId::Forward && hitid.subdetId() == HGCHEB) {
	    const HGCRecHitCollection& rechitsBH = *recHitHandleBH;
	    //	    std::cout << " >>> in BH " << std::endl;
	    // loop over BH RecHits
	    for (HGCRecHitCollection::const_iterator it_hit = rechitsBH.begin(); it_hit < rechitsBH.end(); ++it_hit) {
	      const HcalDetId detid = it_hit->detid();
	      unsigned int layer = recHitTools.getLayerWithOffset(detid);

	      if(detid == hitid){
		//		std::cout << " >>> trovato " << std::endl;
		if(layer != hitlayer) std::cout << " >>> problem " << std::endl;
		Energy_layer_calib[layer] += it_hit->energy();
		Energy_layer_calib_fraction[layer] += it_hit->energy()*it_haf->second;
		hitTime = it_hit->time();
		//		isTValid = it_hit->isTimeValid();

		//		hitSumEnergy += it_hit->energy();
		//		hitEnergy = it_hit->energy();
		//		hitFraction = it_haf->second;
		//		LayerOccupancy->Fill(layer);
		if(seedEnergy < it_hit->energy()){
		  seedEnergy = it_hit->energy();
		  seedDet = recHitTools.getSiThickness(detid);
		  //		  seedCalo = 2;
		  //		  seedLayer = layer;
		}
		break;
	      }
	    }
	  }
	  
	  ++totHits;
	  if(hitTime != -1){
	    ++totHitsWithTime;
	    nHitsWithTime_layer[hitlayer] += 1;
	  }
	  nHits_layer[hitlayer] += 1;


	  //	  std::cout << " isTimeValid() = " << isTValid << " hitTime = " << hitTime << std::endl; 
	}//end recHits
      }// end simHit
      //      std::cout << " >>> hitSumEnergy = " << hitSumEnergy << std::endl; 

    }//end simCluster


    float sumCalibRecHitCalib = 0;
    float sumCalibRecHitCalib_fraction = 0;
    // std::cout << " calib size  = " << Energy_layer_calib.size() 
    // 	      << " calib fraction size  = " << Energy_layer_calib_fraction.size() 
    // 	      << " fraction size  = " << Energy_layer_fraction.size() << std::endl;

    for(unsigned int iL=0; iL<Energy_layer_calib.size(); ++iL){
      sumCalibRecHitCalib += Energy_layer_calib[iL];
      sumCalibRecHitCalib_fraction += Energy_layer_calib_fraction[iL];
    }

    //    std::cout << "seedDet = " << seedDet << std::endl;
    
    if(seedDet == 100){
      h_EoP_CPene_100_calib_fraction->Fill(sumCalibRecHitCalib_fraction / it_caloPart->energy());
      h_EoP_CPene_100_calib->Fill(sumCalibRecHitCalib / it_caloPart->energy());

      //      std::cout << ">> sumCalibRecHitCalib = " << sumCalibRecHitCalib << std::endl;
      //      std::cout << ">> it_caloPart->energy() = " << it_caloPart->energy() << std::endl;
    }
    if(seedDet == 200){
      h_EoP_CPene_200_calib_fraction->Fill(sumCalibRecHitCalib_fraction / it_caloPart->energy());
      h_EoP_CPene_200_calib->Fill(sumCalibRecHitCalib / it_caloPart->energy());

      //      std::cout << ">> sumCalibRecHitCalib = " << sumCalibRecHitCalib << std::endl;
      //      std::cout << ">> it_caloPart->energy() = " << it_caloPart->energy() << std::endl;
    }
    if(seedDet == 300){
      h_EoP_CPene_300_calib_fraction->Fill(sumCalibRecHitCalib_fraction / it_caloPart->energy());
      h_EoP_CPene_300_calib->Fill(sumCalibRecHitCalib / it_caloPart->energy());

      //      std::cout << ">> sumCalibRecHitCalib = " << sumCalibRecHitCalib << std::endl;
      //      std::cout << ">> it_caloPart->energy() = " << it_caloPart->energy() << std::endl;
    }


    float frWT = 1.*totHitsWithTime/totHits;
    h_EoP_CPene_vs_frWithTime->Fill(frWT, sumCalibRecHitCalib_fraction/it_caloPart->energy());
    //    std::cout << " >>> sumCalibRecHitCalib_fraction = " << sumCalibRecHitCalib_fraction << std::endl;
    if(seedDet == 100) h_frWithTime_100->Fill(frWT);
    if(seedDet == 200) h_frWithTime_200->Fill(frWT);  
    if(seedDet == 300) h_frWithTime_300->Fill(frWT);

    for(unsigned int ij=0; ij<nHitsWithTime_layer.size(); ++ij){
      h_maxHitWithTime_vs_layer->Fill(ij, nHitsWithTime_layer[ij]);
      h_maxHit_vs_layer->Fill(ij, nHits_layer[ij]);
      h_fracmaxHitWithTime_vs_layer->Fill(ij, 1.*nHitsWithTime_layer[ij]/nHits_layer[ij]);
    }

    //    std::cout << " >>> simClusterCount = " << simClusterCount << std::endl;
  }//end caloparticle


}

void
HGCalHitCalibration::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
HGCalHitCalibration::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCalHitCalibration::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalHitCalibration);
