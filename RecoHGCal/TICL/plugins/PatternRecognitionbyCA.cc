// Author: Felice Pantaleo, Marco Rovere - felice.pantaleo@cern.ch, marco.rovere@cern.ch
// Date: 11/2018
// Copyright CERN
#include <algorithm>
#include <set>
#include <vector>
#include <cmath>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PatternRecognitionbyCA.h"

void PatternRecognitionbyCA::fillHistogram(
    const std::vector<reco::CaloCluster> &layerClusters,
    const std::vector<std::pair<unsigned int, float>> &mask) {
  if (algo_verbosity_ > None) {
    LogDebug("HGCPatterRecoByCA") << "filling eta/phi histogram per Layer" << std::endl;
  }
  for (auto &m : mask) {
    auto lcId = m.first;
    const auto &lc = layerClusters[lcId];
    // getting the layer Id from the detId of the first hit of the layerCluster
    const auto firstHitDetId = lc.hitsAndFractions()[0].first;
    int layer = rhtools_.getLayerWithOffset(firstHitDetId) +
      rhtools_.lastLayerFH() * ((rhtools_.zside(firstHitDetId) + 1) >> 1) - 1;
    assert(layer >= 0);
    auto etaBin = getEtaBin(lc.eta());
    auto phiBin = getPhiBin(lc.phi());
    histogram_[layer][globalBin(etaBin, phiBin)].push_back(lcId);
    if (algo_verbosity_ > Advanced) {
      LogDebug("HGCPatterRecoByCA")
	<< "Adding layerClusterId: " << lcId << " into bin [eta,phi]: [ " << etaBin << ", "
	<< phiBin << "] for layer: " << layer << std::endl;
    }
  }
}

std::vector<SeedingRegion> PatternRecognitionbyCA::createSeedingRegions(std::vector<PropagationSeedingPoint>& points){
  std::vector<SeedingRegion> regions;
  for(auto point : points){
    auto etaBin = getEtaBin(point.p.eta());
    auto phiBin = getPhiBin(point.p.phi());
    regions.emplace_back(etaBin - 2, 5, phiBin - 2, 5, int(point.p.z() > 0));
  }
  return regions;
}

void PatternRecognitionbyCA::makeTracksters(const edm::Event &ev, const edm::EventSetup &es,
                                            const std::vector<reco::CaloCluster> &layerClusters,
                                            const std::vector<std::pair<unsigned int, float>> &mask,
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
  globalRegion.emplace_back(0, int(nEtaBins_), 0, int(nPhiBins_), 0);

  pointV.emplace_back(point);
  globalRegion.emplace_back(0, int(nEtaBins_), 0, int(nPhiBins_), 1);

  makeTracksters(ev, es, layerClusters, mask, result, pointV, globalRegion);
  return;
}

void PatternRecognitionbyCA::makeTracksters(const edm::Event &ev, const edm::EventSetup &es,
                                            const std::vector<reco::CaloCluster> &layerClusters,
                                            const std::vector<std::pair<unsigned int, float>> &mask,
                                            std::vector<Trackster> &result, std::vector<PropagationSeedingPoint>& points) {

  std::vector<SeedingRegion> regions = createSeedingRegions(points);

  makeTracksters(ev, es, layerClusters, mask, result, points, regions);
  return;
}


void PatternRecognitionbyCA::makeTracksters(const edm::Event &ev, const edm::EventSetup &es,
                                            const std::vector<reco::CaloCluster> &layerClusters,
                                            const std::vector<std::pair<unsigned int, float>> &mask,
					    std::vector<Trackster> &result, std::vector<PropagationSeedingPoint>& points,
					    std::vector<SeedingRegion>& regions) {
  rhtools_.getEventSetup(es);

  clearHistogram();
  theGraph_.setVerbosity(algo_verbosity_);
  theGraph_.clear();
  if (algo_verbosity_ > None) {
    LogDebug("HGCPatterRecoByCA") << "Making Tracksters with CA" << std::endl;
  }
  std::vector<HGCDoublet::HGCntuplet> foundNtuplets;
  fillHistogram(layerClusters, mask);

  theGraph_.makeAndConnectDoublets(histogram_, nEtaBins_, nPhiBins_, layerClusters,
                                   2, 2, min_cos_theta_, min_cos_pointing_, missing_layers_,
                                   rhtools_.lastLayerFH(), points, regions);

  //RAFIX
  //change passing a vector of index
  //also will return a vector of reco index of the tracksters
  theGraph_.findNtuplets(foundNtuplets, min_clusters_per_ntuplet_);
  //#ifdef FP_DEBUG
  const auto &doublets = theGraph_.getAllDoublets();
  int tracksterId = 0;
  for (auto &ntuplet : foundNtuplets) {
    std::set<unsigned int> effective_cluster_idx;
    for (auto &doublet : ntuplet) {
      auto innerCluster = doublets[doublet].innerClusterId();
      auto outerCluster = doublets[doublet].outerClusterId();
      effective_cluster_idx.insert(innerCluster);
      effective_cluster_idx.insert(outerCluster);
      if (algo_verbosity_ > Advanced) {
        LogDebug("HGCPatterRecoByCA")
	  << "New doublet " << doublet << " for trackster: " << result.size() << " InnerCl "
	  << innerCluster << " " << layerClusters[innerCluster].x() << " "
	  << layerClusters[innerCluster].y() << " " << layerClusters[innerCluster].z()
	  << " OuterCl " << outerCluster << " " << layerClusters[outerCluster].x() << " "
	  << layerClusters[outerCluster].y() << " " << layerClusters[outerCluster].z() << " "
	  << tracksterId << std::endl;
      }
    }
    // Put back indices, in the form of a Trackster, into the results vector
    Trackster tmp;
    tmp.vertices.reserve(effective_cluster_idx.size());
    std::copy(std::begin(effective_cluster_idx), std::end(effective_cluster_idx),
              std::back_inserter(tmp.vertices));
    result.push_back(tmp);
    tracksterId++;
  }
  //#endif
}
