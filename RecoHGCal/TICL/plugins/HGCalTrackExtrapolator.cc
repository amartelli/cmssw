// Author: Felice Pantaleo - felice.pantaleo@cern.ch
// Date: 09/2018
// Copyright CERN

// user include files
#include <vector>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

#include "RecoHGCal/TICL/interface/PatternRecognitionAlgoBase.h"
#include "RecoHGCal/TICL/interface/Trackster.h"
#include "RecoHGCal/TICL/plugins/HGCalTrackExtrapolator.h"

#include "PatternRecognitionbyCA.h"
#include "PatternRecognitionbyMultiClusters.h"

#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

DEFINE_FWK_MODULE(HGCalTrackExtrapolator);

HGCalTrackExtrapolator::HGCalTrackExtrapolator(const edm::ParameterSet& ps)
  : myAlgo_(std::make_unique<PatternRecognitionbyCA>(ps)),
    cutTk_(ps.getParameter<std::string>("cutTk")),
    propName_(ps.getParameter<std::string>("propagator")) {
  trks_token_ = consumes<reco::TrackCollection>(
      ps.getParameter<edm::InputTag>("tracker_tracks"));
  clusters_token_ = consumes<std::vector<reco::CaloCluster>>(
      ps.getParameter<edm::InputTag>("hgcal_layerclusters"));
  filtered_layerclusters_mask_token_ = consumes<std::vector<std::pair<unsigned int, float>>>(
      ps.getParameter<edm::InputTag>("filtered_layerclusters_mask"));
  original_layerclusters_mask_token_ =
      consumes<std::vector<float>>(ps.getParameter<edm::InputTag>("original_layerclusters_mask"));
  produces<std::vector<Trackster>>("TrackstersByCA");
  produces<std::vector<float>>();  // Mask to be applied at the next iteration

}

void HGCalTrackExtrapolator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracker_tracks", edm::InputTag("generalTracks"));
  desc.add<std::string>("cutTk", "1.48 < abs(eta) < 3.0 && pt > 0.5 && p > 1 && quality(\"highPurity\") && hitPattern().numberOfLostHits(\"MISSING_OUTER_HITS\") < 10");
  desc.add<std::string>("propagator", "PropagatorWithMaterial");
  desc.add<edm::InputTag>("hgcal_layerclusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<edm::InputTag>("filtered_layerclusters_mask",
                          edm::InputTag("FilteredLayerClusters", "iterationLabelGoesHere"));
  desc.add<edm::InputTag>("original_layerclusters_mask",
                          edm::InputTag("hgcalLayerClusters", "InitialLayerClustersMask"));
  desc.add<int>("algo_verbosity", 0);
  desc.add<double>("min_cos_theta", 0.915);
  desc.add<double>("min_cos_pointing", -1.);
  desc.add<int>("missing_layers", 0);
  desc.add<int>("min_clusters_per_ntuplet", 10);
  descriptions.add("HGCalTrackExtrapolator", desc);
}


void HGCalTrackExtrapolator::buildFirstLayers(){

  const CaloSubdetectorGeometry *subGeom = geom_->getSubdetectorGeometry(DetId::HGCalEE, ForwardSubdetector::ForwardEmpty);

  std::vector<float>  rmax(2, 0), rmin(2, 9e9), zDisk(2, 0);
  const std::vector<DetId> & ids = subGeom->getValidDetIds();

  for (auto & i : ids) {
    int layer =  HGCSiliconDetId(i).layer();
    if(layer != 1) continue;

    const GlobalPoint & pos = geom_->getPosition(i);
    float z = pos.z();
    int side = int(z > 0);
    zDisk[side] = z;

    float rho = pos.perp();
    if (rho > rmax[side]) rmax[side] = rho;
    if (rho < rmin[side]) rmin[side] = rho;
  }

  for (int iSide=0; iSide<2; ++iSide){
    firstDisk[iSide] = new GeomDet(Disk::build(Disk::PositionType(0,0,zDisk[iSide]), Disk::RotationType(), 
					       SimpleDiskBounds(rmin[iSide], rmax[iSide], zDisk[iSide]-0.5, zDisk[iSide]+0.5)).get() );
  }
}


bool HGCalTrackExtrapolator::propagateToFirstLayers(const reco::TrackRef &tk, std::vector<float>& impactP){

  FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(*tk, bfield_.product());

  //PropagationDirection direction = alongMomentum;
  const Propagator &prop = (*propagator_);
  int iSide = int(tk->eta() > 0);
  TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk[iSide]->surface());

  if (!tsos.isValid()){
    //printf(" Destination zside %+1d, z = %+8.2f\n", iSide, firstDisk[iSide]->toGlobal(LocalPoint(0,0,0)).z());
    //printf("         --> propagation failed.\n");    
    return false;
  }

  GlobalPoint gp = tsos.globalPosition();
  //    printf( " Prop point: global eta %+5.2f phi %+5.2f  x = %+8.2f y = %+8.2f z = %+8.2f rho = %8.2f\n", 

  const math::XYZVector dirRefProp = math::XYZVector(tsos.globalMomentum().x(), tsos.globalMomentum().y(), tsos.globalMomentum().z());

  impactP[0] = dirRefProp.eta();
  impactP[1] = dirRefProp.phi();
  impactP[2] = gp.x();
  impactP[3] = gp.y();
  impactP[4] = gp.z();

  return true;
}


void HGCalTrackExtrapolator::produce(edm::Event& evt, const edm::EventSetup& es) {

  //  std::cout << " \n in HGCalTrackExtrapolator " << std::endl;
  auto result = std::make_unique<std::vector<Trackster>>();
  auto output_mask = std::make_unique<std::vector<float>>();

  edm::Handle<reco::TrackCollection> tracks_h;
  edm::Handle<std::vector<reco::CaloCluster>> cluster_h;
  edm::Handle<std::vector<std::pair<unsigned int, float>>> filtered_layerclusters_mask_h;
  edm::Handle<std::vector<float>> original_layerclusters_mask_h;

  evt.getByToken(trks_token_, tracks_h);
  evt.getByToken(clusters_token_, cluster_h);
  evt.getByToken(filtered_layerclusters_mask_token_, filtered_layerclusters_mask_h);
  evt.getByToken(original_layerclusters_mask_token_, original_layerclusters_mask_h);

  const auto& layerClusters = *cluster_h;
  const auto& inputClusterMask = *filtered_layerclusters_mask_h;
  std::unique_ptr<std::vector<std::pair<unsigned int, float>>> filteredLayerClusters;

  //all this block should be done once per run not per event
  //to propagate tracks to 1st layer
  es.get<CaloGeometryRecord>().get(caloGeom_);
  geom_ = caloGeom_.product();
  es.get<TrackingComponentsRecord>().get(propName_, propagator_);
  es.get<IdealMagneticFieldRecord>().get(bfield_);
  buildFirstLayers();

  //should create a better container also keeping track of the track idx
  std::vector<std::vector<float>> pointRefDir;
  unsigned int itrack = 0;
  if(tracks_h->size() == 0) std::cout << " empty trks " << std::endl;
  for (const reco::Track &tk : *tracks_h) { ++itrack;
    if (!cutTk_(tk)){
      continue;
    }

    std::vector<float> locDir(5,0.);
    bool keepTrack = propagateToFirstLayers(reco::TrackRef(tracks_h,itrack-1), locDir);

    //std::cout << " propagated:   eta = " << locDir[0] << " phi= " << locDir[1] << std::endl;
    if(keepTrack) pointRefDir.emplace_back(locDir);
  }

  for (int iSide =0; iSide < 2; ++iSide)
    delete firstDisk[iSide];

  if(pointRefDir.size() != 0)  myAlgo_->makeTrackstersSeeded(evt, es, layerClusters, inputClusterMask, *result, pointRefDir);
  //  else {std::cout << " non seeding track - continue " << std::endl;}

  // Now update the global mask and put it into the event
  output_mask->reserve(original_layerclusters_mask_h->size());
  // Copy over the previous state
  std::copy(std::begin(*original_layerclusters_mask_h), std::end(*original_layerclusters_mask_h),
	    std::back_inserter(*output_mask));
  
  
  evt.put(std::move(result), "TrackstersByCA");
  evt.put(std::move(output_mask));
  
}
