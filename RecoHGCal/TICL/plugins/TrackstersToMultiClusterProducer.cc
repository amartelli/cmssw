#ifndef __RecoHGCal_TICL_TracksterToMultiClusterProducer_H__
#define __RecoHGCal_TICL_TracksterToMultiClusterProducer_H__

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
#include "RecoHGCal/TICL/interface/Trackster.h"

#include <string>
#include <vector>
#include <array>
#include <algorithm>

class TrackstersToMultiClusterProducer : public edm::stream::EDProducer<> {
 public:
  TrackstersToMultiClusterProducer(const edm::ParameterSet&);
  ~TrackstersToMultiClusterProducer() override { }
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  void produce(edm::Event&, const edm::EventSetup&) override;

 private:
  std::string label_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> layer_clusters_token_;
  edm::EDGetTokenT<std::vector<Trackster> > tracksters_token_;
};

DEFINE_FWK_MODULE(TrackstersToMultiClusterProducer);

TrackstersToMultiClusterProducer::TrackstersToMultiClusterProducer(const edm::ParameterSet &ps)
  : label_(ps.getParameter<std::string>("label")) {
    layer_clusters_token_ =
      consumes<std::vector<reco::CaloCluster>>(
          ps.getParameter<edm::InputTag>("LayerClusters"));
    tracksters_token_ =
      consumes<std::vector<Trackster>>(
          ps.getParameter<edm::InputTag>("Tracksters"));
    produces<std::vector<reco::HGCalMultiCluster> >(label_);
}


void TrackstersToMultiClusterProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // hgcalMultiClusters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("Tracksters", edm::InputTag("Tracksters", "TrackstersByCA"));
  desc.add<edm::InputTag>("LayerClusters", edm::InputTag("hgcalLayerClusters"));
  desc.add<std::string>("label", "MultiClustersFromTracksterByCA");
  desc.addUntracked<unsigned int>("verbosity", 3);
  descriptions.add("TrackstersToMultiCluster", desc);
}

void TrackstersToMultiClusterProducer::produce(edm::Event& evt, const edm::EventSetup& es) {

  auto multiclusters = std::make_unique<std::vector<reco::HGCalMultiCluster>>();
  edm::Handle<std::vector<Trackster> > tracksterHandle;
  evt.getByToken(tracksters_token_, tracksterHandle);

  edm::Handle<std::vector<reco::CaloCluster>> layer_clustersHandle;
  evt.getByToken(layer_clusters_token_, layer_clustersHandle);

  auto const & tracksters = *tracksterHandle;
  auto const & layerClusters = *layer_clustersHandle;

  edm::PtrVector<reco::BasicCluster> clusterPtrs;
  for( unsigned i = 0; i < layerClusters.size(); ++i ) {
    edm::Ptr<reco::BasicCluster> ptr(layer_clustersHandle,i);
    clusterPtrs.push_back(ptr);
  }

  std::for_each(std::begin(tracksters), std::end(tracksters),
    [&](auto const & trackster) {
      std::array<double, 3> baricenter{{0., 0., 0.}};
      double total_weight = 0.;
      reco::HGCalMultiCluster temp;
      std::for_each(std::begin(trackster.vertices), std::end(trackster.vertices),
        [&](unsigned int idx) {
          temp.push_back(clusterPtrs[idx]);
          auto weight = clusterPtrs[idx]->energy();
          total_weight += weight;
          baricenter[0] += clusterPtrs[idx]->x()*weight;
          baricenter[1] += clusterPtrs[idx]->y()*weight;
          baricenter[2] += clusterPtrs[idx]->z()*weight;
        });
      std::transform(std::begin(baricenter),
        std::end(baricenter), std::begin(baricenter),
        [&total_weight](double val)->double {return val/total_weight;}
      );
      temp.setEnergy(total_weight);
      temp.setPosition(math::XYZPoint(baricenter[0], baricenter[1], baricenter[2]));
      temp.setAlgoId(reco::CaloCluster::hgcal_em);
      multiclusters->push_back(temp);
    }
  );

  evt.put(std::move(multiclusters), label_);
}

#endif
