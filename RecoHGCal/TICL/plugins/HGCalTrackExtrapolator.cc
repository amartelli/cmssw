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
#include "DataFormats/TICL/interface/Trackster.h"
#include "RecoHGCal/TICL/interface/PatternRecognitionAlgoBase.h"
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
  : myAlgo_(std::make_unique<ticl::PatternRecognitionbyCA>(ps)),
    cutTk_(ps.getParameter<std::string>("cutTk")),
    propName_(ps.getParameter<std::string>("propagator")) {
  trks_token_ = consumes<reco::TrackCollection>(ps.getParameter<edm::InputTag>("tracker_tracks"));
  clusters_token_ = consumes<std::vector<reco::CaloCluster>>(ps.getParameter<edm::InputTag>("layer_clusters"));
  filtered_layerclusters_mask_token_ = consumes<std::vector<float>>(ps.getParameter<edm::InputTag>("filtered_mask"));
  original_layerclusters_mask_token_ = consumes<std::vector<float>>(ps.getParameter<edm::InputTag>("original_mask"));
  layer_clusters_tiles_token_ = consumes<ticl::TICLLayerTiles>(ps.getParameter<edm::InputTag>("layer_clusters_tiles"));
  produces<std::vector<Trackster>>();
  produces<std::vector<float>>();  // Mask to be applied at the next iteration

  detectorName_ = "HGCalEESensitive";
}

void HGCalTrackExtrapolator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracker_tracks", edm::InputTag("generalTracks"));
  desc.add<std::string>("cutTk", "1.48 < abs(eta) < 3.0 && pt > 0.5 && p > 1 && quality(\"highPurity\") && hitPattern().numberOfLostHits(\"MISSING_OUTER_HITS\") < 10");
  desc.add<std::string>("propagator", "PropagatorWithMaterial");
  desc.add<edm::InputTag>("layer_clusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<edm::InputTag>("filtered_mask", edm::InputTag("FilteredLayerClusters", "iterationLabelGoesHere"));
  desc.add<edm::InputTag>("original_mask", edm::InputTag("hgcalLayerClusters", "InitialLayerClustersMask"));
  desc.add<edm::InputTag>("layer_clusters_tiles", edm::InputTag("TICLLayerTileProducer"));
  desc.add<int>("algo_verbosity", 0);
  desc.add<double>("min_cos_theta", 0.915);
  desc.add<double>("min_cos_pointing", -1.);
  desc.add<int>("missing_layers", 0);
  desc.add<int>("min_clusters_per_ntuplet", 10);
  descriptions.add("HGCalTrackExtrapolator", desc);
}


void HGCalTrackExtrapolator::beginRun(edm::Run const& iEvent, edm::EventSetup const& es){

  edm::ESHandle<HGCalDDDConstants> hdc;
  es.get<IdealGeometryRecord>().get(detectorName_, hdc);
  hgcons_ = hdc.product();

  buildFirstLayers();

  es.get<IdealMagneticFieldRecord>().get(bfield_);
  es.get<TrackingComponentsRecord>().get(propName_, propagator_);
}


void HGCalTrackExtrapolator::endRun(edm::Run const&, edm::EventSetup const&){

  for (int iSide = 0; iSide < 2; ++iSide)
    delete firstDisk_[iSide];
}


void HGCalTrackExtrapolator::buildFirstLayers(){

  float zVal = hgcons_->waferZ(1, true);
  std::pair<double, double> rMinMax = hgcons_->rangeR(zVal, true);
  //  std::cout << " z = " << zVal << " min = " << rMinMax.first << " max = " << rMinMax.second << std::endl;

  for (int iSide=0; iSide<2; ++iSide){
    float zSide = (iSide == 0) ? (-1.*zVal) : zVal;
    firstDisk_[iSide] = new GeomDet(Disk::build(Disk::PositionType(0,0,zSide), Disk::RotationType(),
						SimpleDiskBounds(rMinMax.first, rMinMax.second, zSide-0.5, zSide+0.5)).get() );
  }
}


bool HGCalTrackExtrapolator::propagateToFirstLayers(const reco::TrackRef &tk, SeedingRegion& impactP){

  FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(*tk, bfield_.product());

  //PropagationDirection direction = alongMomentum;
  const Propagator &prop = (*propagator_);
  int iSide = int(tk->eta() > 0);

  TrajectoryStateOnSurface tsos = prop.propagate(fts, firstDisk_[iSide]->surface());
  if (!tsos.isValid()){
    //printf(" Destination zside %+1d, z = %+8.2f\n", iSide, firstDisk[iSide]->toGlobal(LocalPoint(0,0,0)).z());
    //printf("         --> propagation failed.\n");
    return false;
  }

  impactP.point = GlobalPoint(tsos.globalPosition());
  impactP.direction = GlobalVector(tsos.globalMomentum().x(), tsos.globalMomentum().y(), tsos.globalMomentum().z());
  impactP.zSide = iSide;

  return true;
}


void HGCalTrackExtrapolator::produce(edm::Event& evt, const edm::EventSetup& es) {

  auto result = std::make_unique<std::vector<Trackster>>();
  auto output_mask = std::make_unique<std::vector<float>>();

  edm::Handle<reco::TrackCollection> tracks_h;
  edm::Handle<std::vector<reco::CaloCluster>> cluster_h;
  edm::Handle<std::vector<float>> filtered_layerclusters_mask_h;
  edm::Handle<std::vector<float>> original_layerclusters_mask_h;
  edm::Handle<ticl::TICLLayerTiles> layer_clusters_tiles_h;

  evt.getByToken(trks_token_, tracks_h);
  evt.getByToken(clusters_token_, cluster_h);
  evt.getByToken(filtered_layerclusters_mask_token_, filtered_layerclusters_mask_h);
  evt.getByToken(original_layerclusters_mask_token_, original_layerclusters_mask_h);
  evt.getByToken(layer_clusters_tiles_token_, layer_clusters_tiles_h);

  const auto& layerClusters = *cluster_h;
  const auto& inputClusterMask = *filtered_layerclusters_mask_h;
  const auto& layer_clusters_tiles = *layer_clusters_tiles_h;


  //should create a better container also keeping track of the track idx
  std::vector<SeedingRegion> pointRefDir;
  int itrack = -1;
  if(tracks_h->size() == 0) std::cout << " empty trks " << std::endl;
  for (const reco::Track &tk : *tracks_h) {
    ++itrack;
    if (!cutTk_(tk)){
      continue;
    }

    SeedingRegion point;
    bool keepTrack = propagateToFirstLayers(reco::TrackRef(tracks_h,itrack), point);

    if(keepTrack){
      //FIXME RA                                                                
      int firstLayer = (point.zSide > 0) ? 52 : 0;
      point.etaBin = layer_clusters_tiles[firstLayer].getEtaBin(point.point.eta());
      point.phiBin = layer_clusters_tiles[firstLayer].getPhiBin(point.point.phi());
      point.index = itrack;

      point.etaStart =  point.etaBin - 2;
      point.etaLength = 5;
      point.phiStart = point.phiBin - 2;
      point.phiLength = 5;
      pointRefDir.emplace_back(point);
    }
  }

  if(pointRefDir.size() != 0)
    myAlgo_->makeTracksters(evt, es, layerClusters, inputClusterMask, layer_clusters_tiles, *result, pointRefDir);


  // Now update the global mask and put it into the event
  output_mask->reserve(original_layerclusters_mask_h->size());
  // Copy over the previous state
  std::copy(std::begin(*original_layerclusters_mask_h),
            std::end(*original_layerclusters_mask_h),
            std::back_inserter(*output_mask));

  //RA maybe not
  /*
  // Mask the used elements, accordingly
  for (auto const& trackster : *result) {
    for (auto const v : trackster.vertices) {
      // TODO(rovere): for the moment we mask the layer cluster completely. In
      // the future, properly compute the fraction of usage.
      (*output_mask)[v] = 0.;
    }
  }
  */
  evt.put(std::move(result));
  evt.put(std::move(output_mask));
}
