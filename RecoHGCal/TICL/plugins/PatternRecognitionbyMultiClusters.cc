#include "PatternRecognitionbyMultiClusters.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

void ticl::PatternRecognitionbyMultiClusters::makeTracksters(const edm::Event& ev,
                                                             const edm::EventSetup& es,
                                                             const std::vector<reco::CaloCluster>& layerClusters,
                                                             const std::vector<float>& mask,
                                                             const ticl::TICLLayerTiles& tiles,
                                                             const std::vector<ticl::TICLSeedingRegion>& regions,
                                                             std::vector<Trackster>& result) {
  LogDebug("HGCPatterRecoTrackster") << "making Tracksters" << std::endl;
}
