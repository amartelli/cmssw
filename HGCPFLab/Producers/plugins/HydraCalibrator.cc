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
#include "DataFormats/HGCRecHit/interface/HGCUncalibratedRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

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
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"

#include "HGCPFLab/Producers/interface/MipCalibrator.h"

#include <map>
#include <string>
#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
using namespace std;
using namespace edm;
using namespace reco;

class HydraCalibrator : public EDProducer
{
public:
    HydraCalibrator( const ParameterSet & );
    
private:
    void produce( Event &, const EventSetup & ) override;

    EDGetTokenT<View<Hydra> > tokenHydra_;
    unique_ptr<HydraWrapper> hydraObj;
    std::vector<edm::InputTag> uncalibRHit_;

    edm::ESHandle<HGCalGeometry> hgceeGeoHandle_;
    edm::ESHandle<HGCalGeometry> hgchefGeoHandle_;
    edm::ESHandle<HGCalGeometry> hgchebGeoHandle_;   

    //const HGCalDDDConstants* hgcons_;


    bool useGenParticles_; // if true use SimTracks corresponding to GenParticles
                           // if false use SimTracks at calorimeter face
    
    bool splitRecHits_; // if true, split the rechits in proportion to the contributing simhits
                        // if false, assign each rechit fully to the simtrack that contributed most

    bool debugPrint_;
    double minDebugEnergy_;

    unsigned int GetHGCLayer(const DetId& detid, const ForwardSubdetector& subdet) const;

    MipCalibrator mipCalib_;
    
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


    TProfile* h_amplitude_vsLayer;
    TProfile* h_energy_vsLayer;
    TProfile* h_absorberEnergy_vsLayer;
    TProfile* h_Mip_vsLayer_abs;
    TProfile* h_Mip_vsLayer_absP1;
    TProfile* h_Mip_vsLayer_act;

    TH2F* h_energy_vs_amplitude;
    TProfile* tp_energy_vs_amplitude;
    
    TH2F* h_energy_vs_recHE;
    TProfile* tp_energy_vs_recHE;

    TH1F* h_EoP_hitTrack;
    TProfile* tp_EoP_hitTrack_vsP;

    TH2F* h_thickness_vs_R;
    TProfile2D* h_thickness_xyMap_R50;
    TProfile2D* h_thickness_xyMap_R80;
    TProfile2D* h_thickness_xyMap_R110;

    TProfile2D* h_R_layer;

    std::map<int, float> nMIP_layer;
    std::map<int, float> Energy_layer_act;
    std::map<int, float> Energy_layer_abs;

    //    CalibHGC m_calibEE, m_calibHEF, m_calibHEB;
    //    bool calibInitialized;
    // https://github.com/lgray/cmssw/blob/topic_cheated_reco/RecoParticleFlow/PandoraTranslator/plugins/PandoraCMSPFCandProducer.cc#L390-L477
};

HydraCalibrator::HydraCalibrator( const ParameterSet &iConfig ) :
    tokenHydra_( consumes<View<Hydra> >( iConfig.getParameter<InputTag> ( "HydraTag" ) ) ),
    useGenParticles_( iConfig.getParameter<bool>( "UseGenParticles" ) ),
    splitRecHits_( iConfig.getParameter<bool>( "SplitRecHits" ) ),
    debugPrint_( iConfig.getUntrackedParameter<bool>( "DebugPrint" , true ) ),
    minDebugEnergy_( iConfig.getUntrackedParameter<double>( "MinDebugEnergy", 10. ) )
    //    m_calibEE(ForwardSubdetector::HGCEE,"EE",debugPrint), m_calibHEF(ForwardSubdetector::HGCHEF,"HEF",debugPrint), m_calibHEB(ForwardSubdetector::HGCHEB,"HEB",debugPrint), calibInitialized(false),
{
    uncalibRHit_ = iConfig.getParameter<std::vector<InputTag> >("HGCalUncalibRecHitCollection");
    for( const auto& tag : uncalibRHit_ ) {
        consumes<HGCUncalibratedRecHitCollection>(tag);
    }


    //    mipCalib_(iConfig);

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

    h_amplitude_vsLayer = fs->make<TProfile>("h_amplitude_vsLayer", "", 50, 0., 50.);
    h_energy_vsLayer = fs->make<TProfile>("h_energy_vsLayer", "", 50, 0., 50.);
    h_absorberEnergy_vsLayer = fs->make<TProfile>("h_absorberEnergy_vsLayer", "", 50, 0., 50.);

    h_Mip_vsLayer_abs = fs->make<TProfile>("h_Mip_vsLayer_abs", "", 50, 0., 50.);
    h_Mip_vsLayer_absP1 = fs->make<TProfile>("h_Mip_vsLayer_absP1", "", 50, 0., 50.);
    h_Mip_vsLayer_act = fs->make<TProfile>("h_Mip_vsLayer_act", "", 50, 0., 50.);

    h_energy_vs_amplitude = fs->make<TH2F>("h_energy_vs_amplitude", "", 1000, 0., 2., 1000, 0., 2.);
    tp_energy_vs_amplitude = fs->make<TProfile>("tp_energy_vs_amplitude", "", 1000, 0., 2.);

    h_energy_vs_recHE = fs->make<TH2F>("h_energy_vs_recHE", "", 1000, 0., 2., 1000, 0., 2.);
    tp_energy_vs_recHE = fs->make<TProfile>("tp_energy_vs_recHE", "", 1000, 0., 2.);

    h_EoP_hitTrack = fs->make<TH1F>("h_EoP_hitTrack", "", 500, -1., 2.);
    tp_EoP_hitTrack_vsP = fs->make<TProfile>("h_EoP_hitTrack", "", 200, 0., 100.);

    h_thickness_vs_R = fs->make<TH2F>("h_thickness_vs_R", "", 300, 0., 300., 5, 0., 5.);
    h_thickness_xyMap_R50 = fs->make<TProfile2D>("h_thickness_xyMap_R50", "", 400, -200., 200., 400, -200., 200.);
    h_thickness_xyMap_R80 = fs->make<TProfile2D>("h_thickness_xyMap_R80", "", 400, -200., 200., 400, -200., 200.);
    h_thickness_xyMap_R110 = fs->make<TProfile2D>("h_thickness_xyMap_R110", "", 400, -200., 200., 400, -200., 200.);

    h_R_layer = fs->make<TProfile2D>("h_R_layer", "", 60, 0., 60., 200, 0., 200.);


    nMIP_layer.clear();
    Energy_layer_act.clear();
    Energy_layer_abs.clear();

    // //convecrt in GeV
    // for(unsigned int iC=0; iC<dEdX_abs.size(); ++iC) dEdX_abs[iC] = dEdX_abs[iC]/1000.;

}

void HydraCalibrator::produce( Event &iEvent, const EventSetup& iSetup)
{

    iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",hgceeGeoHandle_) ;
    iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",hgchefGeoHandle_) ;
    iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive",hgchebGeoHandle_) ;

    const HGCalGeometry* geom = nullptr;

    // const HGCalTopology& topo = geom->topology();
    // const HGCalDDDConstants& dddConst = topo.dddConstants();



    // edm::ESHandle<HGCalDDDConstants> pHGNDC;
    // iSetup.get<IdealGeometryRecord>().get(nameSense_,pHGNDC);
    // const HGCalDDDConstants hgdc(*pHGNDC);


    //    std::cout << " HydraCalibrator::produce useGenParticles=" << useGenParticles_ << std::endl;
    //    if ( splitRecHits_ ) {
    //        throw cms::Exception("NotImplemented")
    //            << "splitRecHits option not available just yet";
    //    }
    Handle<View<Hydra> > HydraHandle;
    iEvent.getByToken(tokenHydra_, HydraHandle);
    assert ( HydraHandle->size() == 1 );
    hydraObj.reset( new HydraWrapper( HydraHandle->ptrAt(0)) );

    vector<Handle<HGCUncalibratedRecHitCollection>> uncalibRecHits;
    for( const auto& tag : uncalibRHit_ ) {
        std::cout << tag << std::endl;
        uncalibRecHits.emplace_back( Handle<HGCUncalibratedRecHitCollection>() );
        iEvent.getByLabel(tag, uncalibRecHits.back());

    }

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
    /*
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
    */




    for ( auto i : tracksToBecomeClusters ) {
        reco::PFCluster temp;
        edm::Ptr<SimTrack> currentTrack = hydraObj->simTrack( i );

        nMIP_layer.clear();
        Energy_layer_act.clear();
        Energy_layer_abs.clear();

        auto simHits = hydraObj->simHitsFromSimTrack( i, true );
        for ( unsigned j = 0 ; j < hydraObj->recHitSize(); j++ ) {
            edm::Ptr<reco::PFRecHit> currentRecHit = hydraObj->recHit( j );
            //            std::cout << " hasSimHitFromRecHit( j ) = "  << hydraObj->hasSimHitFromRecHit( j ) << std::endl;
            if ( hydraObj->hasSimHitFromRecHit( j ) ) {
                if ( !splitRecHits_ ) {
                    auto result = hydraObj->simHitAndFractionFromRecHit( j );
                    
                    //Returns an iterator to the first element in the range [first,last) 
                    //that compares equal to val. If no such element is found, the function returns last.
                    auto match = std::find(simHits.begin(), simHits.end(), result.first );
                    
                    if ( match != simHits.end() ) {
                        //						std::cout << " match true! "<< std::endl;
                        /*
                        if ( true ) {
                            std::cout << "   Track " << i << "  simHit "<<(result.first)->energy() << " matches to rechit " << j << " with fraction " << result.second;
                            std::cout << "     e(track) = " << currentTrack->momentum().E() << " e(hit)=" << currentRecHit->energy();
                            std::cout << std::endl;
                        }
                        */
                        edm::Ref<reco::PFRecHitCollection> currentRecHitRef = hydraObj->recHitRef( j );
                        if ( false && debugPrint_ ) std::cout << "    energy of currentRecHitRef: " << currentRecHitRef->energy() << std::endl;
                        //                        std::cout << "       Unsplit: we use a fraction of 1, but it's really " << result.second << std::endl;
                        reco::PFRecHitFraction rhf( currentRecHitRef, 1. ); // unsplit case
                        temp.addRecHitFraction( rhf );
                        

						int layer = 999;
                        int det = 999;
                        std::string detBlock = "";

                        //                        if(true) std::cout <<  "recHit Id = " << currentRecHitRef->detId() << std::endl;
						DetId hitid(currentRecHitRef->detId());
						layer = GetHGCLayer( hitid, (ForwardSubdetector)hitid.subdetId() );
                        det = hitid.subdetId();
                        //                        double cellType = ((HGCalDetId)(currentRecHitRef->detId())).waferType();

                        float rh_x = currentRecHitRef->position().x();
                        float rh_y = currentRecHitRef->position().y();
                        float rh_R2 = pow(rh_x,2)+pow(rh_y,2); 

                        /*
                        if(det != 5){
                            // check thickness and choose 100 200 300
                            if(cellType == 1) detBlock = "100";
                            else{
                                if(rh_R2 > 130.*130.) detBlock = "300";
                                else detBlock = "200";
                            }
                        }
                        */

						if( det == 3 ){ 
                        	layer = GetHGCLayer( hitid, (ForwardSubdetector)hitid.subdetId() );
                            geom = hgceeGeoHandle_.product();     
                        }else if( det == 4 ){
                        	layer = 28 + GetHGCLayer( hitid, (ForwardSubdetector)hitid.subdetId() );
                            geom = hgchefGeoHandle_.product();
                        }else if( det == 5 ){
                        	layer = 28 + 12 + GetHGCLayer( hitid, (ForwardSubdetector)hitid.subdetId() );
                            detBlock = "BH";
                            geom = hgchebGeoHandle_.product();
                        }

                        const HGCalTopology& topo = geom->topology();
                        const HGCalDDDConstants& dddConst = topo.dddConstants();
                        //int sec = (unsigned int) ((HGCalDetId)(hitid)).sector();
                        int cellThickness = dddConst.waferTypeL((unsigned int) ((HGCalDetId)(hitid)).wafer() );
                        //                        std::cout << " >>>>>> cellThickness = " << cellThickness << std::endl;
                        if(cellThickness == 1) detBlock = "100";
                        if(cellThickness == 2) detBlock = "200";
                        if(cellThickness == 3) detBlock = "300";

                        h_thickness_vs_R->Fill(sqrt(rh_R2), cellThickness);
                        h_R_layer->Fill(layer, sqrt(rh_R2), cellThickness);
                        if(sqrt(rh_R2) < 50) h_thickness_xyMap_R50->Fill(rh_x, rh_y, cellThickness);
                        if(sqrt(rh_R2) > 70 && sqrt(rh_R2) < 90) h_thickness_xyMap_R80->Fill(rh_x, rh_y, cellThickness);
                        if(sqrt(rh_R2) > 110) h_thickness_xyMap_R110->Fill(rh_x, rh_y, cellThickness);

                        //						std::cout << "layer: " <<(unsigned int) ((HGCalDetId)(hitid)).layer(); //<< " waffer: " <<(unsigned int) ((HGCalDetId)(hitid)).sector();
                        //std::cout << " cel type: " << (unsigned int) ((HGCalDetId)(hitid)).subsector()
                        //                        std::cout << " cell: " << (unsigned int) ((HGCalDetId)(hitid)).cell() << std::endl;
                        //						if(debugPrint_) std::cout << "recHit Id = "<< currentRecHitRef->detId() << "  subdet = " << det << "  layer = " << layer << std::endl;


                        //                        bool amplitudeSet = false;
                        double amplitude = -10.;

                        //match recHit with uncalib
                        for(unsigned k = 0; k < uncalibRecHits.size(); ++k){
                            //    std::cout << " k "<< k << std::endl;
                            for(HGCUncalibratedRecHitCollection::const_iterator itr = uncalibRecHits[k]->begin(); itr != uncalibRecHits[k]->end(); ++itr){
                                DetId uncDetId( DetId(itr->id()));
                                if(uncDetId == hitid){
                                    amplitude = double(itr->amplitude());
                                    //          amplitudeSet = true;
                                    break;
                                }
                            }
                        }

                        if(amplitude != -10.){
                            //dEdX_act[detBlock] =  dedX in Si of correct thickness
                            nMIP_layer[layer] += amplitude * mipCalib_.GetAct_fC2MIP(detBlock);

                            float energyRecHit = amplitude * mipCalib_.GetAct_fC2MIP(detBlock) * mipCalib_.GetAct_MIP2keV(detBlock) * 1.e-06;  // in keV??? => * 1.e-6
                            Energy_layer_act[layer] += energyRecHit;

                            //                            h_Mip_vsLayer_act->Fill(nL, nMIP_layer[nL]);


                            h_energy_vs_amplitude->Fill(amplitude, energyRecHit);
                            tp_energy_vs_amplitude->Fill(amplitude, energyRecHit);
                            
                            h_energy_vs_recHE->Fill(currentRecHitRef->energy(), energyRecHit);
                            tp_energy_vs_recHE->Fill(currentRecHitRef->energy(), energyRecHit);
                        
                            h_amplitude_vsLayer->Fill(layer, amplitude);
                            h_energy_vsLayer->Fill(layer, energyRecHit);

                            //                            h_Mip_vsLayer_abs->Fill(layer, nMIP_layer[layer]);
                        }
                        else std::cout << " amplitude == 10 PROBLEM!!! " << std::endl;



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
                            // if ( false && debugPrint_ ) {
                            //     std::cout << "   Track " << i << " matches to rechit " << j << " with fraction " << result.second;
                            //     std::cout << "     e(track) = " << currentTrack->momentum().E() << " e(hit)=" << currentRecHit->energy();
                            //     std::cout << std::endl;
                            // }
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
                // if (debugPrint_) std::cout << "Track " << i << " energy=" << currentTrack->momentum().E() << " fake cluster energy=" << temp.energy() << std::endl;
                // if (debugPrint_) std::cout << "Track " << i << " eta=" << currentTrack->momentum().Eta() << " fake cluster eta=" << temp.eta() << std::endl;
                // if (debugPrint_) std::cout << "Track " << i << " phi=" << currentTrack->momentum().Phi() << " fake cluster phi=" << temp.phi() << std::endl;
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
    

        float sumCalibRecHit = 0.;   
        for(unsigned int nL=1; nL<(1+nMIP_layer.size()); ++nL){
            //            nMIP_layer[nL] = nMIP_layer[nL];
            if(nL == 1) {
                Energy_layer_abs[nL] = nMIP_layer[nL] * mipCalib_.GetAbs_MIP2keV(nL) * 1.e-03;   // in MeV
                h_Mip_vsLayer_abs->Fill(nL, nMIP_layer[nL]);        
                h_Mip_vsLayer_absP1->Fill(nL, nMIP_layer[nL+1]);        
            }
            else{
                Energy_layer_abs[nL] = (nMIP_layer[nL-1]+nMIP_layer[nL])/2. * mipCalib_.GetAbs_MIP2keV(nL) * 1.e-03; // in MeV
                h_Mip_vsLayer_abs->Fill(nL, (nMIP_layer[nL-1]+nMIP_layer[nL])/2.);        
                h_Mip_vsLayer_absP1->Fill(nL, (nMIP_layer[nL-1]+nMIP_layer[nL])/2.);        
            }
            sumCalibRecHit += Energy_layer_abs[nL] + Energy_layer_act[nL];
            h_absorberEnergy_vsLayer->Fill(nL, Energy_layer_abs[nL]);
                        
            //                h_Mip_vsLayer_abs->Fill(nL, nMIP_layer[nL]);
            // if(nL != nMIP_layer.size()+1)            h_Mip_vsLayer_absP1->Fill(nL, nMIP_layer[nL+1]);
            
            //std::cout << " nL = " << nL << " >>> nMIP_layer[layer] = " << nMIP_layer[nL] << std::endl;
            h_Mip_vsLayer_abs->Fill(nL, nMIP_layer[nL]);
            h_Mip_vsLayer_act->Fill(nL, nMIP_layer[nL]);
        }        
        h_EoP_hitTrack->Fill(sumCalibRecHit / currentTrack->momentum().E());
        tp_EoP_hitTrack_vsP->Fill(currentTrack->momentum().E(), sumCalibRecHit / currentTrack->momentum().E());
        
        
    }
    
    
    iEvent.put ( pfClusterCol );
    


}


unsigned int HydraCalibrator::GetHGCLayer(const DetId& detid, const ForwardSubdetector& subdet) const {
    unsigned int layer = 0;
	
    layer = (unsigned int) ((HGCalDetId)(detid)).layer() ;
	
    return layer;
}


DEFINE_FWK_MODULE( HydraCalibrator );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
