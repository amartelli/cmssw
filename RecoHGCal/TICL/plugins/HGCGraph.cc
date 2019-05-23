// Author: Felice Pantaleo - felice.pantaleo@cern.ch
// Date: 11/2018
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/TICL/interface/Common.h"
#include "PatternRecognitionbyCA.h"
#include "HGCDoublet.h"
#include "HGCGraph.h"

void HGCGraph::makeAndConnectDoublets(const ticl::TICLLayerTiles &histo,
                                      int nEtaBins,
                                      int nPhiBins,
                                      const std::vector<reco::CaloCluster> &layerClusters,
                                      const std::vector<float> &mask,
                                      int deltaIEta, int deltaIPhi, float minCosTheta,
                                      float minCosPointing, int missing_layers, int maxNumberOfLayers,
				      std::vector<PropagationSeedingPoint>& points, std::vector<SeedingRegion>& regions) {

  isOuterClusterOfDoublets_.clear();
  isOuterClusterOfDoublets_.resize(layerClusters.size());
  allDoublets_.clear();
  theRootDoublets_.clear();
  int indexRegion = -1;
  for(auto seedR : regions){
    ++indexRegion;

    for (int il = 0; il < maxNumberOfLayers - 1; ++il) {

      for (int outer_layer = 0; outer_layer < std::min(1 + missing_layers, maxNumberOfLayers - 1 - il); ++outer_layer) {
        int currentInnerLayerId = il + maxNumberOfLayers * seedR.zSide;
        int currentOuterLayerId = currentInnerLayerId + 1 + outer_layer;
        auto const &outerLayerHisto = histo[currentOuterLayerId];
        auto const &innerLayerHisto = histo[currentInnerLayerId];

	int etaS = seedR.etaStart;
	int phiS = seedR.phiStart;

	for (int oeta = std::max(0, etaS); oeta < std::min(nEtaBins, etaS + seedR.etaLength); ++oeta) {
          auto offset = oeta * nPhiBins;
	  for (int ophi_it = phiS; ophi_it < phiS + seedR.phiLength; ++ophi_it) {
	    int ophi = (( ophi_it % nPhiBins + nPhiBins) % nPhiBins);
            for (auto outerClusterId : outerLayerHisto[offset + ophi]) {
              // Skip masked clusters
              if (mask[outerClusterId] == 0.)
                continue;
              const auto etaRangeMin = std::max(0, oeta - deltaIEta);
              const auto etaRangeMax = std::min(oeta + deltaIEta, nEtaBins);

              for (int ieta = etaRangeMin; ieta < etaRangeMax; ++ieta) {
                // wrap phi bin
                for (int phiRange = 0; phiRange < 2 * deltaIPhi + 1; ++phiRange) {
                  // The first wrapping is to take into account the
                  // cases in which we would have to seach in
                  // negative bins. The second wrap is mandatory to
                  // account for all other cases, since we add in
                  // between a full nPhiBins slot.
                  auto iphi = ((ophi + phiRange - deltaIPhi) % nPhiBins + nPhiBins) % nPhiBins;
                  for (auto innerClusterId : innerLayerHisto[ieta * nPhiBins + iphi]) {
                    // Skip masked clusters
                    if (mask[innerClusterId] == 0.)
                      continue;
                    auto doubletId = allDoublets_.size();
                    allDoublets_.emplace_back(innerClusterId, outerClusterId, doubletId,
                                              &layerClusters, points[indexRegion].index);
                    if (verbosity_ > Advanced) {
                      LogDebug("HGCGraph")
                          << "Creating doubletsId: " << doubletId << " layerLink in-out: [" << currentInnerLayerId
                          << ", " << currentOuterLayerId << "] clusterLink in-out: [" << innerClusterId << ", "
                          << outerClusterId << "]" << std::endl;
                    }
                    isOuterClusterOfDoublets_[outerClusterId].push_back(doubletId);
                    auto &neigDoublets = isOuterClusterOfDoublets_[innerClusterId];
                    auto &thisDoublet = allDoublets_[doubletId];
                    if (verbosity_ > Expert) {
                      LogDebug("HGCGraph")
                          << "Checking compatibility of doubletId: " << doubletId
                          << " with all possible inners doublets link by the innerClusterId: " << innerClusterId
                          << std::endl;
                    }
                    bool isRootDoublet = thisDoublet.checkCompatibilityAndTag(
	                allDoublets_, neigDoublets, points[indexRegion].v,
			minCosTheta, minCosPointing, verbosity_ > Advanced);
                    if (isRootDoublet)
                      theRootDoublets_.push_back(doubletId);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  // #ifdef FP_DEBUG
  if (verbosity_ > None) {
    LogDebug("HGCGraph") << "number of Root doublets " << theRootDoublets_.size() << " over a total number of doublets "
                         << allDoublets_.size() << std::endl;
  }
  // #endif
}


//RA: as argument also give vector with tracks indices
// return vector of reco tracksters = as done now
// return also vector of tracks indices of reco tracksters
void HGCGraph::findNtuplets(std::vector<HGCDoublet::HGCntuplet> &foundNtuplets,
                            const unsigned int minClustersPerNtuplet) {
  HGCDoublet::HGCntuplet tmpNtuplet;
  tmpNtuplet.reserve(minClustersPerNtuplet);
  for (auto rootDoublet : theRootDoublets_) {
    tmpNtuplet.clear();
    allDoublets_[rootDoublet].findNtuplets(allDoublets_, tmpNtuplet);
    if (tmpNtuplet.size() > minClustersPerNtuplet) {
      foundNtuplets.push_back(tmpNtuplet);
    }
  }
}
