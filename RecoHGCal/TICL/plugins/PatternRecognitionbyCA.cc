// Author: Felice Pantaleo, Marco Rovere - felice.pantaleo@cern.ch, marco.rovere@cern.ch
// Date: 11/2018
#include <algorithm>
#include <set>
#include <vector>
#include <cmath>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/TICL/interface/TICLLayerTile.h"
#include "PatternRecognitionbyCA.h"
#include "HGCGraph.h"

using namespace ticl;

PatternRecognitionbyCA::PatternRecognitionbyCA(const edm::ParameterSet &conf) : PatternRecognitionAlgoBase(conf) {
  theGraph_ = std::make_unique<HGCGraph>();
  min_cos_theta_ = (float)conf.getParameter<double>("min_cos_theta");
  min_cos_pointing_ = (float)conf.getParameter<double>("min_cos_pointing");
  missing_layers_ = conf.getParameter<int>("missing_layers");
  min_clusters_per_ntuplet_ = conf.getParameter<int>("min_clusters_per_ntuplet");
}

PatternRecognitionbyCA::~PatternRecognitionbyCA(){};


std::vector<SeedingRegion> PatternRecognitionbyCA::createSeedingRegions(std::vector<PropagationSeedingPoint>& points,
									const ticl::TICLLayerTiles &tiles){
  std::vector<SeedingRegion> regions;
  for(auto point : points){
    // 0 for z<0 and rhtools_.lastLayerFH() for z > 0
    //take values from ticl::constants::nLayers
    //assuming that coherence is checked elsewhere
    int firstLayer = (point.p.z() > 0) ? rhtools_.lastLayerFH() : 0;
    auto etaBin = tiles[firstLayer].getEtaBin(point.p.eta());
    auto phiBin = tiles[firstLayer].getPhiBin(point.p.phi());
    regions.emplace_back(etaBin - 2, 5, phiBin - 2, 5, int(point.p.z() > 0));
  }
  return regions;
}

void PatternRecognitionbyCA::makeTracksters(const edm::Event &ev, const edm::EventSetup &es,
                                            const std::vector<reco::CaloCluster> &layerClusters,
                                            const std::vector<float> &mask,
                                            const ticl::TICLLayerTiles &tiles,
                                            std::vector<Trackster> &result) {

  // for unseeded iterations create 2 global seeding regions
  // spanning all eta and phi bins for both endcaps
  std::vector<PropagationSeedingPoint> pointV;
  std::vector<SeedingRegion> globalRegion;

  PropagationSeedingPoint point;
  point.p = GlobalPoint(0,0,0);
  point.v = GlobalVector(0,0,0);
  point.index = -1;

  pointV.emplace_back(point);
  globalRegion.emplace_back(0, ticl::constants::nEtaBins, 0, ticl::constants::nPhiBins, 0);

  pointV.emplace_back(point);
  globalRegion.emplace_back(0, ticl::constants::nEtaBins, 0, ticl::constants::nPhiBins, 1);

  makeTracksters(ev, es, layerClusters, mask, tiles, result, pointV, globalRegion);
  return;
}

void PatternRecognitionbyCA::makeTracksters(const edm::Event &ev, const edm::EventSetup &es,
                                            const std::vector<reco::CaloCluster> &layerClusters,
					    const std::vector<float> &mask,
                                            const ticl::TICLLayerTiles &tiles,
                                            std::vector<Trackster> &result, std::vector<PropagationSeedingPoint>& points) {

  std::vector<SeedingRegion> regions = createSeedingRegions(points, tiles);

  makeTracksters(ev, es, layerClusters, mask, tiles, result, points, regions);
  return;
}


void PatternRecognitionbyCA::makeTracksters(const edm::Event &ev, const edm::EventSetup &es,
                                            const std::vector<reco::CaloCluster> &layerClusters,
                                            const std::vector<float> &mask,
                                            const ticl::TICLLayerTiles &tiles,
					    std::vector<Trackster> &result, std::vector<PropagationSeedingPoint>& points,
					    std::vector<SeedingRegion>& regions) {
  rhtools_.getEventSetup(es);

  theGraph_->setVerbosity(algo_verbosity_);
  theGraph_->clear();
  if (algo_verbosity_ > None) {
    LogDebug("HGCPatterRecoByCA") << "Making Tracksters with CA" << std::endl;
  }
  std::vector<HGCDoublet::HGCntuplet> foundNtuplets;
  std::vector<uint8_t> layer_cluster_usage(layerClusters.size(), 0);
  theGraph_->makeAndConnectDoublets(tiles,
                                    ticl::constants::nEtaBins,
                                    ticl::constants::nPhiBins,
                                    layerClusters,
                                    mask,
                                    2,
                                    2,
                                    min_cos_theta_,
                                    min_cos_pointing_,
                                    missing_layers_,
                                    rhtools_.lastLayerFH(),
				    points, regions);

  //RAFIX
  //change passing a vector of index
  //also will return a vector of reco index of the tracksters
  theGraph_->findNtuplets(foundNtuplets, min_clusters_per_ntuplet_);


  //#ifdef FP_DEBUG
  const auto &doublets = theGraph_->getAllDoublets();
  int tracksterId = 0;
  for (auto const &ntuplet : foundNtuplets) {
    std::set<unsigned int> effective_cluster_idx;
    for (auto const &doublet : ntuplet) {
      auto innerCluster = doublets[doublet].innerClusterId();
      auto outerCluster = doublets[doublet].outerClusterId();
      effective_cluster_idx.insert(innerCluster);
      effective_cluster_idx.insert(outerCluster);
      if (algo_verbosity_ > Advanced) {
        LogDebug("HGCPatterRecoByCA") << "New doublet " << doublet << " for trackster: " << result.size() << " InnerCl "
                                      << innerCluster << " " << layerClusters[innerCluster].x() << " "
                                      << layerClusters[innerCluster].y() << " " << layerClusters[innerCluster].z()
                                      << " OuterCl " << outerCluster << " " << layerClusters[outerCluster].x() << " "
                                      << layerClusters[outerCluster].y() << " " << layerClusters[outerCluster].z()
                                      << " " << tracksterId << std::endl;
      }
    }
    for (auto const i : effective_cluster_idx) {
      layer_cluster_usage[i]++;
      LogDebug("HGCPatterRecoByCA") << "LayerID: " << i << " count: " << (int)layer_cluster_usage[i] << std::endl;
    }
    // Put back indices, in the form of a Trackster, into the results vector
    Trackster tmp;
    tmp.vertices.reserve(effective_cluster_idx.size());
    tmp.vertex_multiplicity.resize(effective_cluster_idx.size(), 0);
    std::copy(std::begin(effective_cluster_idx), std::end(effective_cluster_idx), std::back_inserter(tmp.vertices));
    result.push_back(tmp);
    tracksterId++;
  }

  for (auto &trackster : result) {
    for (size_t i = 0; i < trackster.vertices.size(); ++i) {
      assert(i < trackster.vertex_multiplicity.size());
      trackster.vertex_multiplicity[i] = layer_cluster_usage[trackster.vertices[i]];
      LogDebug("HGCPatterRecoByCA") << "LayerID: " << trackster.vertices[i]
                                    << " count: " << (int)trackster.vertex_multiplicity[i] << std::endl;
    }
  }
}
