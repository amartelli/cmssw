#ifndef __newpf_PFClusterProducer_H__
#define __newpf_PFClusterProducer_H__

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoParticleFlow/PFClusterProducer/interface/RecHitTopologicalCleanerBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/SeedFinderBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/InitialClusteringStepBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterBuilderBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFCPositionCalculatorBase.h"
#include "RecoParticleFlow/PFClusterProducer/interface/PFClusterEnergyCorrectorBase.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include <memory>

#include "RecoLocalCalo/HGCalRecHitDump/interface/HGCalImagingAlgo.h"
#include "RecoLocalCalo/HGCalRecHitDump/interface/HGCalDepthPreClusterer.h"
#include "RecoLocalCalo/HGCalRecHitDump/interface/HGCalMultiCluster.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

/*
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
*/

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

class HGCalClusterTestProducer : public edm::stream::EDProducer<> {
 public:
  HGCalClusterTestProducer(const edm::ParameterSet&);
  ~HGCalClusterTestProducer() { }

  virtual void produce(edm::Event&, const edm::EventSetup&);

 private:

  edm::EDGetTokenT<HGCRecHitCollection> hits_ee_token;
  edm::EDGetTokenT<HGCRecHitCollection> hits_ef_token;
  edm::EDGetTokenT<HGCRecHitCollection> hits_eb_token;
  edm::EDGetTokenT<edm::View<reco::Candidate> > genParticlesToken;

  reco::CaloCluster::AlgoId algoId;

  HGCalImagingAlgo *algo;
  bool doSharing;
  std::string detector;

  HGCalImagingAlgo::VerbosityLevel verbosity;
};

DEFINE_FWK_MODULE(HGCalClusterTestProducer);

HGCalClusterTestProducer::HGCalClusterTestProducer(const edm::ParameterSet &ps) :
  algoId(reco::CaloCluster::undefined),
  algo(0),doSharing(ps.getParameter<bool>("doSharing")),
  detector(ps.getParameter<std::string >("detector")),              //one of EE, EF BH or "both" or "all"
  verbosity((HGCalImagingAlgo::VerbosityLevel)ps.getUntrackedParameter<unsigned int>("verbosity",3)){
  double ecut = ps.getParameter<double>("ecut");
  double delta_c = ps.getParameter<double>("deltac");
  double kappa = ps.getParameter<double>("kappa");

  if(detector=="all"){
    hits_ee_token = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit:HGCEERecHits"));
    hits_ef_token = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit:HGCHEFRecHits"));
    hits_eb_token = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit:HGCHEBRecHits"));
    algoId = reco::CaloCluster::hgcal_mixed;
  }if(detector=="both"){
    hits_ee_token = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit:HGCEERecHits"));
    hits_ef_token = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit:HGCHEFRecHits"));
    algoId = reco::CaloCluster::hgcal_mixed;
  }else if(detector=="EE"){
    hits_ee_token = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit:HGCEERecHits"));
    algoId = reco::CaloCluster::hgcal_em;
  }else if(detector=="EF"){
    hits_ef_token = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit:HGCHEFRecHits"));
    algoId = reco::CaloCluster::hgcal_had;
  }else if(detector=="EB"){
    hits_eb_token = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit:HGCHEBRecHits"));
    algoId = reco::CaloCluster::hgcal_had;
  }
  if(doSharing){
    double showerSigma =  ps.getParameter<double>("showerSigma");
    algo = new HGCalImagingAlgo(delta_c, kappa, ecut, showerSigma, 0, algoId, verbosity);
  }else{
    algo = new HGCalImagingAlgo(delta_c, kappa, ecut, 0, algoId, verbosity);
  }

  genParticlesToken = consumes<edm::View<reco::Candidate> >( edm::InputTag("genParticles") );

  std::cout << "Constructing HGCalClusterTestProducer" << std::endl;

  produces<std::vector<reco::BasicCluster> >();
  produces<std::vector<reco::BasicCluster> >("sharing");
}

void HGCalClusterTestProducer::produce(edm::Event& evt,
				       const edm::EventSetup& es) {
  edm::ESHandle<HGCalGeometry> ee_geom;
  es.get<IdealGeometryRecord>().get("HGCalEESensitive",ee_geom);

  edm::ESHandle<HGCalGeometry> ef_geom;
  es.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",ef_geom);

  edm::ESHandle<HGCalGeometry> eb_geom;
  es.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive",eb_geom);

  edm::Handle<HGCRecHitCollection> ee_hits;
  edm::Handle<HGCRecHitCollection> ef_hits;
  edm::Handle<HGCRecHitCollection> eb_hits;

  edm::Handle<edm::View<reco::Candidate> > genParticlesH;
  evt.getByToken(genParticlesToken, genParticlesH);
  const auto& genParticles = *genParticlesH;

  std::unique_ptr<std::vector<reco::BasicCluster> > clusters( new std::vector<reco::BasicCluster> ),
    clusters_sharing( new std::vector<reco::BasicCluster> );

  algo->reset();
  //  switch(algoId){
  switch(algoId){
  case reco::CaloCluster::hgcal_em:
    evt.getByToken(hits_ee_token,ee_hits);
    algo->setGeometry(ee_geom.product());
    algo->populate(*ee_hits);
    break;
  case  reco::CaloCluster::hgcal_had:
    if(detector=="EF"){
      evt.getByToken(hits_ef_token,ef_hits);
      algo->setGeometry(ef_geom.product());
      algo->populate(*ef_hits);
      break;
    }
    if(detector=="EB"){
      evt.getByToken(hits_eb_token,ef_hits);
      algo->setGeometry(eb_geom.product());
      algo->populate(*eb_hits);
      break;
    }
  case reco::CaloCluster::hgcal_mixed:
    if(detector=="both"){
      evt.getByToken(hits_ee_token,ee_hits);
      algo->setGeometry(ee_geom.product());
      algo->populate(*ee_hits);
      evt.getByToken(hits_ef_token,ef_hits);
      algo->setGeometry(ef_geom.product());
      algo->populate(*ef_hits);
      break;
    }
    if(detector=="all"){
      evt.getByToken(hits_ee_token,ee_hits);
      algo->setGeometry(ee_geom.product());
      algo->populate(*ee_hits);
      evt.getByToken(hits_ef_token,ef_hits);
      algo->setGeometry(ef_geom.product());
      algo->populate(*ef_hits);
      evt.getByToken(hits_eb_token,eb_hits);
      algo->setGeometry(eb_geom.product());
      algo->populate(*eb_hits);
      break;
    }
  default:
    break;
  }
  algo->makeClusters();
  *clusters = algo->getClusters(false);
  if(doSharing)
    *clusters_sharing = algo->getClusters(true);

  std::cout << "Density based cluster size: " << clusters->size() << std::endl;
  if(doSharing)
    std::cout << "Sharing clusters size     : " << clusters_sharing->size() << std::endl;

  for(unsigned i = 0; i < genParticles.size(); ++i ) {
    std::cout << genParticles[i].pt() << ' ' << genParticles[i].eta() << ' ' << genParticles[i].phi() << std::endl;
  }

  evt.put(std::move(clusters));
  evt.put(std::move(clusters_sharing),"sharing");
}

#endif
