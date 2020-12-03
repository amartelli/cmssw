#ifndef RECOHGCAL_TICL_TRACKSTERSTOF_H
#define RECOHGCAL_TICL_TRACKSTERSTOF_H

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "FWCore/Framework/interface/Event.h"

#include <vector>

typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>> XYZPointF;

namespace ticl {
  class TrackstersTOF{
  public:

    TrackstersTOF();

    void setParameters(edm::ConsumesCollector& sumes);

    void initialize(const edm::EventSetup &es);

    float correctTOFforTracksters(edm::Event &,
				  //std::vector<Trackster> &,
				  const Trackster &,
				  const std::vector<reco::Track> &,
				  const std::vector<reco::Vertex> &,
				  const GlobalPoint &,				  
				  const GlobalPoint &,
				  const GlobalPoint &, const float &, 
				  const float &, const float &,
				  const XYZPointF &, const float &);
    
    float correctTOFforTracksters(edm::Event &,
				  const Trackster &,
				  const reco::Track &,
				  //const size_t trkIdx,
				  const std::vector<reco::Vertex> &,
				  const GlobalPoint &,
				  const XYZPointF &, const float &);

  private:

    float GetPropagatedTrackLength(const GlobalPoint &startingPoint,
				   const GlobalPoint &endPoint,
				   const reco::Track &track);
    

    /*edm::EDGetTokenT<reco::TrackCollection> tracks_token_; */
    /* edm::EDGetTokenT<reco::VertexCollection> vertices_token_; */
    
    edm::EDGetTokenT<std::vector<CaloParticle> > _caloParticles;

    edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfield_token_;
    edm::ESHandle<MagneticField> bfield_;
  };

} // namespace ticl
#endif
