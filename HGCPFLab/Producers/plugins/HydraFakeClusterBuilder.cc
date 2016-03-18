#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HGCPFLab/DataFormats/interface/HighRapidityDevRecoAssociation.h"
#include "HGCPFLab/DataFormats/interface/HydraWrapper.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
//#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

// We probably don't need all of these
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/FlatTrd.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "SimDataFormats/CaloTest/interface/HGCalTestNumbering.h"

#include "TH2.h"
using namespace std;
using namespace edm;
using namespace reco;

class HydraFakeClusterBuilder : public EDProducer
{
public:
    HydraFakeClusterBuilder( const ParameterSet & );
    
private:
    void produce( Event &, const EventSetup & ) override;

    EDGetTokenT<View<Hydra> > tokenHydra_;
    unique_ptr<HydraWrapper> hydraObj;

    bool useGenParticles_; // if true use SimTracks corresponding to GenParticles
                           // if false use SimTracks at calorimeter face
    
    bool splitRecHits_; // if true, split the rechits in proportion to the contributing simhits
                        // if false, assign each rechit fully to the simtrack that contributed most

    bool debugPrint_;
    double minDebugEnergy_;

    unsigned int GetHGCLayer(const DetId& detid, const ForwardSubdetector& subdet) const;
    
    TH1F * h_ene_track;
	TH1F * h_pt_track;
	TH1F * h_eta_track;
	TH1F * h_phi_track;
	TH1F * h_X_track;
	TH1F * h_Y_track;
    TH1F * h_layer;
    TH1F * h_det;
    TH2D * h_cel_vs_X;
    TH2D * h_cel_vs_Y;
    TH2D * h_cel_vs_Z;
    TH2D * h_cel_vs_XY;
    

    //    CalibHGC m_calibEE, m_calibHEF, m_calibHEB;
    //    bool calibInitialized;
    // https://github.com/lgray/cmssw/blob/topic_cheated_reco/RecoParticleFlow/PandoraTranslator/plugins/PandoraCMSPFCandProducer.cc#L390-L477
};

HydraFakeClusterBuilder::HydraFakeClusterBuilder( const ParameterSet &iConfig ) :
    tokenHydra_( consumes<View<Hydra> >( iConfig.getParameter<InputTag> ( "HydraTag" ) ) ),
    useGenParticles_( iConfig.getParameter<bool>( "UseGenParticles" ) ),
    splitRecHits_( iConfig.getParameter<bool>( "SplitRecHits" ) ),
    debugPrint_( iConfig.getUntrackedParameter<bool>( "DebugPrint" , true ) ),
    minDebugEnergy_( iConfig.getUntrackedParameter<double>( "MinDebugEnergy", 10. ) )
    //    m_calibEE(ForwardSubdetector::HGCEE,"EE",debugPrint), m_calibHEF(ForwardSubdetector::HGCHEF,"HEF",debugPrint), m_calibHEB(ForwardSubdetector::HGCHEB,"HEB",debugPrint), calibInitialized(false),
{
    edm::Service<TFileService> fs; 
    produces<reco::PFClusterCollection>();
    h_ene_track = fs->make<TH1F>("h_ene_track","Energy of the track associated to cluster",1000,0,1000);
    h_pt_track = fs->make<TH1F>("h_pt_track","Pt of the track associated to cluster",150,0,150);
    h_eta_track = fs->make<TH1F>("h_eta_track","Eta of the track associated to cluster",500,-4,4);
    h_phi_track = fs->make<TH1F>("h_phi_track","Phi of the track associated to cluster",500,-3.5,3.5);
    h_X_track   = fs->make<TH1F>("h_X_track","X of the track associated to cluster",300,-150,150);
    h_Y_track   = fs->make<TH1F>("h_Y_track","Y of the track associated to cluster",300,-150,150);
	h_layer   = fs->make<TH1F>("h_layer","h_layer",50,0,50);
    h_det   = fs->make<TH1F>("h_det","h_det",15,0,15);
    h_cel_vs_X   = fs->make<TH2D>("h_cel_vs_X","h_cel_vs_X",300,-150,150,500,0,500);
    h_cel_vs_Y   = fs->make<TH2D>("h_cel_vs_Y","h_cel_vs_Y",300,-150,150,500,0,500);
    h_cel_vs_Z   = fs->make<TH2D>("h_cel_vs_Z","h_cel_vs_Z",820,-410,410,500,0,500);
    h_cel_vs_XY	 = fs->make<TH2D>("h_cel_vs_XY","h_cel_vs_XY",300,-150,150,300,-150,150);

}

void HydraFakeClusterBuilder::produce( Event &iEvent, const EventSetup & )
{
    std::cout << " HydraFakeClusterBuilder::produce useGenParticles=" << useGenParticles_ << std::endl;
    //    if ( splitRecHits_ ) {
    //        throw cms::Exception("NotImplemented")
    //            << "splitRecHits option not available just yet";
    //    }
    Handle<View<Hydra> > HydraHandle;
    iEvent.getByToken(tokenHydra_, HydraHandle);
    assert ( HydraHandle->size() == 1 );
    hydraObj.reset( new HydraWrapper( HydraHandle->ptrAt(0)) );

    auto_ptr<reco::PFClusterCollection> pfClusterCol ( new reco::PFClusterCollection );

    std::vector<Index_t> tracksToBecomeClusters;

    for ( unsigned i = 0 ; i < hydraObj->simTrackSize() ; i++ ) {
        edm::Ptr<SimTrack> currentTrack = hydraObj->simTrack( i );
        if (useGenParticles_) {
            if (!currentTrack->noGenpart() ) {
                tracksToBecomeClusters.push_back(i);
            }
        } else {
            unsigned nHitsImmediate = hydraObj->simHitsFromSimTrack( i, false ).size();
            unsigned nHitsAllDescendants = hydraObj->simHitsFromSimTrack( i, true ).size();
            if ( nHitsImmediate > 0 && nHitsAllDescendants == nHitsImmediate ) {
                tracksToBecomeClusters.push_back(i);
            }
        }
    }

    if (debugPrint_) {
        std::cout << "There are " << tracksToBecomeClusters.size() << " tracks to become clusters: " << std::endl;
        for ( auto i : tracksToBecomeClusters ) {
            edm::Ptr<SimTrack> currentTrack = hydraObj->simTrack( i );
            auto p = currentTrack->momentum();
            if (p.E() > minDebugEnergy_) {
                std::cout << "   Track for cheated cluster " << i << " has energy eta geantId pdgId " << p.E() << " " << p.Eta() 
                          << " " << currentTrack->trackId() << " " << currentTrack->type() << std::endl;
            }
        }
    }

    for ( auto i : tracksToBecomeClusters ) {
        reco::PFCluster temp;
        edm::Ptr<SimTrack> currentTrack = hydraObj->simTrack( i );
        auto simHits = hydraObj->simHitsFromSimTrack( i, true );
        for ( unsigned j = 0 ; j < hydraObj->recHitSize(); j++ ) {
            edm::Ptr<reco::PFRecHit> currentRecHit = hydraObj->recHit( j );
            std::cout << " hasSimHitFromRecHit( j ) = "  << hydraObj->hasSimHitFromRecHit( j ) << std::endl;
            if ( hydraObj->hasSimHitFromRecHit( j ) ) {
                if ( !splitRecHits_ ) {
                    auto result = hydraObj->simHitAndFractionFromRecHit( j );
                    
                    //Returns an iterator to the first element in the range [first,last) 
                    //that compares equal to val. If no such element is found, the function returns last.
                    auto match = std::find(simHits.begin(), simHits.end(), result.first );
                    
                    if ( match != simHits.end() ) {
						std::cout << " match true! "<< std::endl;
                        if ( true ) {
                            std::cout << "   Track " << i << "  simHit "<<(result.first)->energy() << " matches to rechit " << j << " with fraction " << result.second;
                            std::cout << "     e(track) = " << currentTrack->momentum().E() << " e(hit)=" << currentRecHit->energy();
                            std::cout << std::endl;
                        }
                        edm::Ref<reco::PFRecHitCollection> currentRecHitRef = hydraObj->recHitRef( j );
                        if ( false && debugPrint_ ) std::cout << "    energy of currentRecHitRef: " << currentRecHitRef->energy() << std::endl;
                        std::cout << "       Unsplit: we use a fraction of 1, but it's really " << result.second << std::endl;
                        reco::PFRecHitFraction rhf( currentRecHitRef, 1. ); // unsplit case
                        temp.addRecHitFraction( rhf );
                        
						int layer = 999;
                        int det = 999;
						
                        if(true) std::cout <<  "recHit Id = " << currentRecHitRef->detId() << std::endl;
						DetId hitid(currentRecHitRef->detId());
						layer = GetHGCLayer( hitid, (ForwardSubdetector)hitid.subdetId() );
                        det = hitid.subdetId();
                        
						if( det == 3 ){ 
                        	layer = GetHGCLayer( hitid, (ForwardSubdetector)hitid.subdetId() );
                        }else if( det == 4 ){
                        	layer = 28 + GetHGCLayer( hitid, (ForwardSubdetector)hitid.subdetId() );
                        }else if( det == 5 ){
                        	layer = 28 + 12 + GetHGCLayer( hitid, (ForwardSubdetector)hitid.subdetId() );
                        }
						std::cout << "layer: " <<(unsigned int) ((HGCalDetId)(hitid)).layer(); //<< " waffer: " <<(unsigned int) ((HGCalDetId)(hitid)).sector();
                        //std::cout << " cel type: " << (unsigned int) ((HGCalDetId)(hitid)).subsector()
                        std::cout << " cell: " << (unsigned int) ((HGCalDetId)(hitid)).cell() << std::endl;
						if(debugPrint_) std::cout << "recHit Id = "<< currentRecHitRef->detId() << "  subdet = " << det << "  layer = " << layer << std::endl;
                        
                        h_layer->Fill(layer);
                        h_det->Fill(det);
                        h_cel_vs_X->Fill(currentRecHitRef->position().x(), (unsigned int) ((HGCalDetId)(hitid)).cell());
						h_cel_vs_Y->Fill(currentRecHitRef->position().y(), (unsigned int) ((HGCalDetId)(hitid)).cell());
                        h_cel_vs_Z->Fill(currentRecHitRef->position().z(), (unsigned int) ((HGCalDetId)(hitid)).cell());
						h_cel_vs_XY->Fill(currentRecHitRef->position().x(), currentRecHitRef->position().y(), (unsigned int) ((HGCalDetId)(hitid)).cell() ); 
                    }
                } else {
                    auto result_vec = hydraObj->simHitsAndFractionsFromRecHit( j );
                    for ( auto result : result_vec ) {
                        auto match = std::find(simHits.begin(), simHits.end(), result.first );
                        if ( match != simHits.end() ) {
                            if ( false && debugPrint_ ) {
                                std::cout << "   Track " << i << " matches to rechit " << j << " with fraction " << result.second;
                                std::cout << "     e(track) = " << currentTrack->momentum().E() << " e(hit)=" << currentRecHit->energy();
                                std::cout << std::endl;
                            }
                            edm::Ref<reco::PFRecHitCollection> currentRecHitRef = hydraObj->recHitRef( j );
                            if ( false && debugPrint_ ) std::cout << "    energy of currentRecHitRef: " << currentRecHitRef->energy() << std::endl;
                            std::cout << "     Split: fraction=" << result.second << std::endl;
                            reco::PFRecHitFraction rhf( currentRecHitRef, result.second ); // split case
                            temp.addRecHitFraction( rhf );
                        }
                    }
                }
            }
        }

        // Cluster position calculation is not log weighted yet!
        float cl_energy = 0.;
        float cl_x =0., cl_y =0., cl_z =0.;
        int n_h = 0;
        for( const reco::PFRecHitFraction& rhf : temp.recHitFractions() ) {
            const reco::PFRecHitRef& refhit = rhf.recHitRef();
            const double rh_fraction = rhf.fraction();
            const double rh_rawenergy = refhit->energy();
            const double rh_energy = rh_rawenergy * rh_fraction;   
            cl_energy += rh_energy;
            auto cl_pos = refhit->position();
            cl_x += cl_pos.x();
            cl_y += cl_pos.y();
            cl_z += cl_pos.z();
            n_h++;
        }
        if (n_h > 0) {
            temp.setEnergy(cl_energy);
            temp.setPosition(math::XYZPoint(cl_x/n_h,cl_y/n_h,cl_z/n_h));
            if ( currentTrack->momentum().E() > minDebugEnergy_ ) {
                if (debugPrint_) std::cout << "Track " << i << " energy=" << currentTrack->momentum().E() << " fake cluster energy=" << temp.energy() << std::endl;
                if (debugPrint_) std::cout << "Track " << i << " eta=" << currentTrack->momentum().Eta() << " fake cluster eta=" << temp.eta() << std::endl;
                if (debugPrint_) std::cout << "Track " << i << " phi=" << currentTrack->momentum().Phi() << " fake cluster phi=" << temp.phi() << std::endl;
            }
            pfClusterCol->push_back( temp );
            h_ene_track->Fill(currentTrack->momentum().E());
            h_pt_track->Fill(currentTrack->momentum().Pt());            
            h_eta_track->Fill(currentTrack->momentum().Eta());
            h_phi_track->Fill(currentTrack->momentum().Phi());
            h_X_track->Fill(currentTrack->trackerSurfacePosition().X());
            h_Y_track->Fill(currentTrack->trackerSurfacePosition().Y());
            
        } else {
            if ( false ) std::cout << "   cluster has no hits, skipping" << std::endl;
        }
    }
        

    iEvent.put ( pfClusterCol );
}


unsigned int HydraFakeClusterBuilder::GetHGCLayer(const DetId& detid, const ForwardSubdetector& subdet) const {
    unsigned int layer = 0;
	
    layer = (unsigned int) ((HGCalDetId)(detid)).layer() ;
	
    return layer;
}


DEFINE_FWK_MODULE( HydraFakeClusterBuilder );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
