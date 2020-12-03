#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackstersTOF.h"
#include "TPrincipal.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "DataFormats/Math/interface/GeantUnits.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <iostream>
#include <set>

#include <Eigen/Core>
#include <Eigen/Dense>

class ConsumesCollector;
using namespace ticl;

TrackstersTOF::TrackstersTOF() {
};

void TrackstersTOF::setParameters(edm::ConsumesCollector &sumes){
  bfield_token_ = sumes.esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>();
  caloParticles_token_ = sumes.consumes<std::vector<CaloParticle> >(edm::InputTag("mix","MergedCaloTruth"));
}

void TrackstersTOF::initialize(const edm::EventSetup &es){
  bfield_ = es.getHandle(bfield_token_);
}

float TrackstersTOF::correctTOFforTracksters(edm::Event &evt,
					     const Trackster &trackster,
					     const std::vector<reco::Track> &tracks,
					     const std::vector<reco::Vertex> &vertex,
					     const GlobalPoint &showerStart,
					     const GlobalPoint &showerStart_LC,
					     const GlobalPoint &showerSeed_LC, const float &seedEnergy_LC, 
					     const float &seedTime_LC, const float &seedTimeE_LC,
					     const XYZPointF &vtxPos, const float & vtxT){

    edm::Handle<std::vector<CaloParticle> > caloParticleHandle;
    evt.getByToken(caloParticles_token_, caloParticleHandle);
    const std::vector<CaloParticle>& caloParticles = *caloParticleHandle;
    
    std::cout << " caloParticles.size() = " << caloParticles.size() << std::endl;
    for(std::vector<CaloParticle>::const_iterator it_caloPart = caloParticles.begin(); it_caloPart != caloParticles.end(); ++it_caloPart){
    std::cout << "   " << it_caloPart->simClusters().size() << std::endl;
    if(it_caloPart->simClusters().size() > 2) continue;
    
    if(trackster.barycenter().eta() * it_caloPart->eta() < 0) continue;

    if(trackster.timeError() == -1) return -99;
  
    //get the vertex from the track (PCA to BS?)
    //get trackster position in HGCAL: would get the barycenter as the trackster time should be closer to the core time
    
    GlobalPoint hgcPoint(trackster.barycenter().x(), trackster.barycenter().y(), trackster.barycenter().z());
    // std::cout << " position = " << hgcPoint << " showerStart = " << showerStart 
    // 	      << " diff = " << sqrt(pow(hgcPoint.x() - showerStart.x(), 2) + pow(hgcPoint.y() - showerStart.y(), 2) + pow(hgcPoint.z() - showerStart.z(), 2)) << std::endl; 



    float min_dz = 0.;

    //gen vtx position
    GlobalPoint vtxPoint(vtxPos.x(), vtxPos.y(), vtxPos.z());

    float trkL0 = hgcPoint.mag();
    float trkL0st = showerStart.mag();
    float trkL0st_LC = showerStart_LC.mag();
    float trkL0seed_LC = showerSeed_LC.mag();

    float trkL = 0.;
    float trkL_showerStart = 0.;
    float trkL_showerStart_LC = 0.;
    float trkL_showerSeed_LC = 0.;
    constexpr double c_cm_ns = geant_units::operators::convertMmToCm(CLHEP::c_light);  // [mm/ns] -> [cm/ns]                                                                           
    constexpr double c_inv = 1.0 / c_cm_ns;
    //std::cout << " diff wrt origin = " << trkL0 - trkL0st << " tof diff = " << (trkL0 - trkL0st) * c_inv << std::endl;

    float beta_pi = 0.;
    int charge = -1;
    float trackP = 1.;
    float correctedT_st = 0.;
    float correctedT = 0.;
    if(trackster.seedID() != edm::ProductID()) {
      auto trackIdx = trackster.seedIndex();
      auto const &track = tracks[trackIdx];
      trackP = track.p();

      trkL = GetPropagatedTrackLength(vtxPoint, hgcPoint, track);
      trkL_showerStart = GetPropagatedTrackLength(vtxPoint, showerStart, track);
      trkL_showerStart_LC = GetPropagatedTrackLength(vtxPoint, showerStart_LC, track);
      trkL_showerSeed_LC = GetPropagatedTrackLength(vtxPoint, showerSeed_LC, track);

      constexpr double m_pi = 0.13957018;
      constexpr double m_pi_inv2 = 1.0 / m_pi / m_pi;
      constexpr double m_k = 0.493677;
      constexpr double m_k_inv2 = 1.0 / m_k / m_k;
      constexpr double m_p = 0.9382720813;
      constexpr double m_p_inv2 = 1.0 / m_p / m_p;

      float gammasq_pi = 1. + track.p2() * m_pi_inv2;
      beta_pi = std::sqrt(1. - 1. / gammasq_pi);
      float dt_pi = (trkL0 - trkL / beta_pi) * c_inv;
      correctedT = trackster.time() +  dt_pi;
      float dt_pi_st = (trkL0st - trkL_showerStart / beta_pi) * c_inv;
      correctedT_st = trackster.time() +  dt_pi_st;
      
      float gammasq_k = 1. + track.p2() * m_k_inv2;
      float beta_k = std::sqrt(1. - 1. / gammasq_k);
      float dt_k = (trkL0 - trkL / beta_k) * c_inv;
      float correctedT_k = trackster.time() +  dt_k;
      charge = 1;
    }//charged
    else{
      trkL = sqrt(pow(vtxPoint.x() - hgcPoint.x(), 2) + pow(vtxPoint.y() - hgcPoint.y(), 2) + pow(vtxPoint.z() - hgcPoint.z(), 2));
      trkL_showerStart = sqrt(pow(vtxPoint.x() - showerStart.x(), 2) + pow(vtxPoint.y() - showerStart.y(), 2) + pow(vtxPoint.z() - showerStart.z(), 2));
      trkL_showerStart_LC = sqrt(pow(vtxPoint.x() - showerStart_LC.x(), 2) + pow(vtxPoint.y() - showerStart_LC.y(), 2) + pow(vtxPoint.z() - showerStart_LC.z(), 2));
      trkL_showerSeed_LC = sqrt(pow(vtxPoint.x() - showerSeed_LC.x(), 2) + pow(vtxPoint.y() - showerSeed_LC.y(), 2) + pow(vtxPoint.z() - showerSeed_LC.z(), 2));

      beta_pi = 1.;
      charge = 0;      
      float dt_photon = (trkL0 - trkL) * c_inv;
      correctedT = trackster.time() +  dt_photon;
      float dt_photon_st = (trkL0st - trkL_showerStart) * c_inv ;//+  (trkL0 - trkL0st) * c_inv;
      correctedT_st = trackster.time() +  dt_photon_st;
      float correctedT_k = 0.;
    }
    // std::cout << " error = " << trackster.timeError() << " *20 = " << 20. * trackster.timeError() << " *15 = " << 15. * trackster.timeError() << std::endl;  
    // std::cout << " trkL = " << trkL << " trkL0 = " << trkL0 << " trkL0st = " << trkL0st << " trkL_showerStart = " << trkL_showerStart << std::endl;

    std::cout << "charge " << charge << " t " << trackster.time() << " correctedT " << correctedT << " correctedT_st " << correctedT_st
	      << " tError " << trackster.timeError() << " raw_pt " << trackster.raw_pt() << " trackp " << trackP
	      << " energy " << trackster.raw_energy() << " vtxT " << vtxT << " vtxPos " << vtxPos.x() << " " << vtxPos.y() << " " << vtxPos.z()
	      << " seedE " << seedEnergy_LC << " seedT " << seedTime_LC << " seedTError " << seedTimeE_LC
	      << " trkL0 " << trkL0 << " trkL " << trkL << " trkL0st " << trkL0st << " trkLst " << trkL_showerStart 
	      << " trkL0st_LC " << trkL0st_LC << " trkLst_LC " << trkL_showerStart_LC << " trkL0seed_LC " << trkL0seed_LC << " trkL_showerSeed_LC " << trkL_showerSeed_LC
	      << " showerStart " << showerStart.x() << " " << showerStart.y() << " " << showerStart.z()
	      << " showerStart_LC " << showerStart_LC.x() << " " << showerStart_LC.y() << " " << showerStart_LC.z()
	      << " showerSeed_LC " << showerSeed_LC.x() << " " << showerStart_LC.y() << " " << showerStart_LC.z()
	      << " hgcPoint " << hgcPoint.x() << " " << hgcPoint.y() << " " << hgcPoint.z() 
	      << " c_inv " << c_inv << " beta_pi " << beta_pi << " nSC " << it_caloPart->simClusters().size()
	      << " pull_st " << (correctedT_st - vtxT) / (10. * trackster.timeError()) 
	      << " pull_stBias " << (correctedT_st - vtxT + 0.01) / (10. * trackster.timeError())
	      << " pull " << (correctedT - vtxT) / (10. * trackster.timeError()) 
	      << " pull_Bias " << (correctedT - vtxT + 0.01) / (10. * trackster.timeError()) 
	      << " CPxyz " << it_caloPart->momentum().x() << " " << it_caloPart->momentum().y() << " " << it_caloPart->momentum().z()  << std::endl;


    //  std::cout << " pullPID = " << correctedT - correctedT_k << std::endl;          

    }//caloparticles

    return -1.;
}

float TrackstersTOF::correctTOFforTracksters(edm::Event &evt, const Trackster &trackster,
					     const reco::Track &track,
					     //const std::vector<reco::CaloCluster> &layerClusters,
					     const std::vector<reco::Vertex> &vertex,
					     const GlobalPoint &showerStart,
					     const XYZPointF &vtxPos, const float & vtxT){

  //caloparticles just to x-check while developing
  edm::Handle<std::vector<CaloParticle> > caloParticleHandle;
  evt.getByToken(caloParticles_token_, caloParticleHandle);
  const std::vector<CaloParticle>& caloParticles = *caloParticleHandle;
  
  std::cout << " caloParticles.size() = " << caloParticles.size() << std::endl;
  for(std::vector<CaloParticle>::const_iterator it_caloPart = caloParticles.begin(); it_caloPart != caloParticles.end(); ++it_caloPart){
    std::cout << "   " << it_caloPart->simClusters().size() << std::endl;
    if(it_caloPart->simClusters().size() > 2) continue;

    if(trackster.barycenter().z() * it_caloPart->momentum().z() < 0) continue;
    if(trackster.timeError() == -1) return -99;

    //get trackster position in HGCAL: would get the barycenter as the trackster time should be closer to the core time
    GlobalPoint hgcPoint(trackster.barycenter().x(), trackster.barycenter().y(), trackster.barycenter().z());
    std::cout << " position = " << hgcPoint << " showerStart = " << showerStart 
	      << " diff = " << sqrt(pow(hgcPoint.x() - showerStart.x(), 2) + pow(hgcPoint.y() - showerStart.y(), 2) + pow(hgcPoint.z() - showerStart.z(), 2)) << std::endl; 


    //  reco::TrackRef tkR(reco::TrackRef(&tracks, trkIdx)); 
    // std::cout << "\n  in TrackstersTOF trackster pt = " << trackster.raw_pt() 
    // 	    
    // 	    << " track pt = " << track.pt() << std::endl;

    //get the vertex from the track: the closest vtx in dz? tobe checked
    float min_dz = 99.;
    int bestVtx_Idx = -1;
    
    for(unsigned i=0; i<vertex.size(); ++i){
      auto const vtx = vertex[i];
      //std::cout << " pos = " << vtx.position() << " time = " << vtx.t() << std::endl;
      float current_dz = std::abs(track.dz(vtx.position()));
      //std::cout << " dz = " << current_dz << std::endl;
      if(current_dz < min_dz){
	min_dz = current_dz;
	bestVtx_Idx = i;
	// for(auto iT = iV.tracks_begin(); iT != iV.tracks_end(); ++iT){
	//   std::cout << " pt tracks = " << (*iT)->pt() << std::endl;
	// }
      }
    }
    
    //  std::cout << " best vtx = " << bestVtx_Idx << std::endl;
    if(bestVtx_Idx == -1) return -99;
    
    //  GlobalPoint vtxPoint(vertex[bestVtx_Idx].position().x(), vertex[bestVtx_Idx].position().y(), vertex[bestVtx_Idx].position().z());
    GlobalPoint vtxPoint(vtxPos.x(), vtxPos.y(), vtxPos.z());
    
    float trkL0 = hgcPoint.mag();
    float trkL0st = showerStart.mag();
    float trkL = GetPropagatedTrackLength(vtxPoint, hgcPoint, track);
    float trkL_showerStart = GetPropagatedTrackLength(vtxPoint, showerStart, track);
    //  std::cout << " trkL = " << trkL << " trkL0 = " << trkL0 << std::endl;
    
    constexpr double m_pi = 0.13957018;
    constexpr double m_pi_inv2 = 1.0 / m_pi / m_pi;
    constexpr double m_k = 0.493677;
    constexpr double m_k_inv2 = 1.0 / m_k / m_k;
    constexpr double m_p = 0.9382720813;
    constexpr double m_p_inv2 = 1.0 / m_p / m_p;
    constexpr double c_cm_ns = geant_units::operators::convertMmToCm(CLHEP::c_light);  // [mm/ns] -> [cm/ns]
    constexpr double c_inv = 1.0 / c_cm_ns;
    
    float gammasq_pi = 1. + track.p2() * m_pi_inv2;
    float beta_pi = std::sqrt(1. - 1. / gammasq_pi);
    float dt_pi = (trkL0 - trkL / beta_pi) * c_inv;
    float correctedT = trackster.time() +  dt_pi;
    
    float dt_pi_st = (trkL0st - trkL_showerStart / beta_pi) * c_inv;
    float correctedT_st = trackster.time() +  dt_pi_st;
    //  std::cout << " original t = " << trackster.time() << " correzione0 = " << trkL0*c_inv << " correzioneL = " << trkL / beta_pi * c_inv << " c_inv = " << c_inv << std::endl;
    //std::cout << " >>> corrected time = " << correctedT << " genT = " << vtxT << " pull = " << (correctedT - vtxT) / trackster.timeError() << " min_dz = " << min_dz << std::endl;
    
    float gammasq_k = 1. + track.p2() * m_k_inv2;
    float beta_k = std::sqrt(1. - 1. / gammasq_k);
    float dt_k = (trkL0 - trkL / beta_k) * c_inv;
    float correctedT_k = trackster.time() +  dt_k;
    //  std::cout << " pullPID = " << correctedT - correctedT_k << std::endl; 
    
    std::cout << " minDz " << min_dz << " pullO " << (trackster.time() - vtxT) / trackster.timeError() 
	      << " pullGhgc " << (correctedT - vtxT) / trackster.timeError() 
	      << " pullGenst " << (correctedT_st - vtxT) / trackster.timeError() 
	      << " pullPID " << (correctedT - correctedT_k) / trackster.timeError() 
	      << " pullGenk " << (correctedT_k - vtxT) / trackster.timeError() 
	      << std::endl;
    
    
    //std::cout << " track time = " << track.time() << " recoVtx t = " << vertex[bestVtx_Idx].t() << " track error t = " << track.timeError() << std::endl;
    //  std::cout << " pullR = " << (correctedT - vertex[bestVtx_Idx].t) / trackster.timeError() << std::endl;
    
  }
  return -99.;
  
  
}



float TrackstersTOF::GetPropagatedTrackLength(const GlobalPoint &startingPoint, 
					      const GlobalPoint &endPoint, 
					      const reco::Track &track){
  // if charged compute the track length exploiting HelixPropagator
  GlobalVector startingMomentum(track.px(), track.py(), track.pz());
  auto magField = bfield_.product();
  FreeTrajectoryState trajectory(startingPoint, startingMomentum, track.charge(), magField);
  SteppingHelixPropagator propagator(magField);
  float propagatedTrackL_ = propagator.propagateWithPath(trajectory, endPoint).second;            
  
  return propagatedTrackL_;
}
