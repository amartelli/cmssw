import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import *
from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import *

# patch particle flow clusters for HGC into local reco sequence 
# (for now until global reco is going with some sort of clustering)
from RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGC_cff import *
from RecoParticleFlow.PFClusterProducer.particleFlowClusterHGC_cfi import *
from RecoLocalCalo.HGCalRecProducers.hgcalLayerClusters_cfi import hgcalLayerClusters
from RecoLocalCalo.HGCalRecProducers.hgcal3DClusters_cfi import hgcal3DClusters

hgcalLocalRecoSequence = cms.Sequence( HGCalUncalibRecHit+
                                       HGCalRecHit+
                                       hgcalLayerClusters+
                                       hgcal3DClusters+
                                       particleFlowRecHitHGCSeq+
                                       particleFlowClusterHGCal+
                                       particleFlowClusterHGCalFromMultiCl )
