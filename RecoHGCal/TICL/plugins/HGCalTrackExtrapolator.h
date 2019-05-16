#ifndef __RecoHGCal_TICL_HGCalTrackExtrapolator_H__
#define __RecoHGCal_TICL_HGCalTrackExtrapolator_H__
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "RecoHGCal/TICL/interface/PatternRecognitionAlgoBase.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

class HGCalTrackExtrapolator : public edm::stream::EDProducer<> {
 public:
  HGCalTrackExtrapolator(const edm::ParameterSet &);
  ~HGCalTrackExtrapolator() override {}
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  void produce(edm::Event &, const edm::EventSetup &) override;

  void buildFirstLayers();
  bool propagateToFirstLayers(const reco::TrackRef &tk, std::vector<float>& impactP);

 private:
  edm::EDGetTokenT<reco::TrackCollection> trks_token_;
  edm::EDGetTokenT<std::vector<reco::CaloCluster>> clusters_token_;
  edm::EDGetTokenT<std::vector<std::pair<unsigned int, float>>> filtered_layerclusters_mask_token_;
  edm::EDGetTokenT<std::vector<float>> original_layerclusters_mask_token_;

  std::unique_ptr<PatternRecognitionAlgoBase> myAlgo_;

  const StringCutObjectSelector<reco::Track> cutTk_;
  const std::string propName_;
  const CaloGeometry *geom_;

  GeomDet* firstDisk[2];
  edm::ESHandle<CaloGeometry> caloGeom_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<MagneticField> bfield_;
};

#endif
