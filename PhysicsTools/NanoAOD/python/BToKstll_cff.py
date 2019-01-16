import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import *

BToKstll=cms.EDProducer("BToKstllProducer",
                        beamSpot = cms.InputTag("offlineBeamSpot"),
                        vertexCollection=cms.InputTag("offlineSlimmedPrimaryVertices"),
                        electronCollection = cms.InputTag("linkedObjects","electrons"), #same collection as in NanoAOD
                        muonCollection = cms.InputTag("linkedObjects","muons"), #same collection as in NanoAOD
                        PFCandCollection = cms.InputTag("packedPFCandidates"),
                        lostSubLeadLepTrackCollection = cms.InputTag("lostTracks"),
                        lostChHadrTrackCollection = cms.InputTag("lostTracks"),

                        nSelectedTriplets = cms.int32(50),  #50
                        isLeptonElectron = cms.bool(False),
                        isChannelKst = cms.bool(False),

                        #case electron
                        LeadEleMinPt = cms.double(1.),
                        LeadEleMaxEta = cms.double(2.4),
                        SubLeadEleMinPt = cms.double(1.),
                        SubLeadEleMaxEta = cms.double(2.4),
                        #case muon
                        LeadMuonMinPt = cms.double(1.),
                        LeadMuonMaxEta = cms.double(2.4),
                        SubLeadMuonMinPt = cms.double(1.),
                        SubLeadMuonMaxEta = cms.double(2.4),

                        KaonMinPt = cms.double(1.),
                        KaonMaxEta = cms.double(2.4),
                        KaonMinDCASig = cms.double(-1.),
                        PionMinPt = cms.double(1.),
                        PionMaxEta = cms.double(2.4),
                        PionMinDCASig = cms.double(-1.),
                        
                        ## following for lepton + lepton + track
                        #diLepton_dz_max = cms.double(-1),
                        #lepKaon_dz_max = cms.double(-1),
                        #lepPion_dz_max = cms.double(-1),
                        #kaonPion_dz_max = cms.double(-1),
                        #kaonRefitllVertex_dxy_max = cms.double(-1), # > 0.2 ?
                        #kll_dxyPV_min = cms.double(-1),  #<1 ?
                        #IPPV_llRefitVtx_min = cms.double(-1),  #<1 ?
                        ###

                        ## following for lepton + track + track
                        diLepton_dz_max = cms.double(2.),
                        lepKaon_dz_max = cms.double(2.),
                        lepPion_dz_max = cms.double(2.),
                        kaonPion_dz_max = cms.double(2.),
                        kaonRefitllVertex_dxy_max = cms.double(0.02), # > 0.2 ?
                        kll_dxyPV_min = cms.double(0.02),  #<1 ?
                        IPPV_llRefitVtx_min = cms.double(1.),  #<1 ?
                        ###

                        DiLeptonChargeCheck = cms.bool(True),
                        KstarChargeCheck = cms.bool(True),
                        JPsiMassConstraint = cms.double(-1), #2-trk refitting uses measured di-ele mass
                        KstMassConstraint = cms.double(0.89176), #2-trk refitting uses nominal K*(892) mass
                        save2TrackRefit = cms.bool(True),
                        #save4TrackRefit = cms.bool(False),
                        useLostSubLeadLepTracks = cms.bool(False),
                        useLostChHadrTracks = cms.bool(False),
                        vtxCL_min = cms.double(1.e-3),  #e-3
                        Bmass_min = cms.double(2.),     # 2.
                        Bmass_max = cms.double(8.),     #8.
                        Bmass_Kst_min = cms.double(2.),
                        Bmass_Kst_max = cms.double(8.)
                      )

BToKstllTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
                             src=cms.InputTag("BToKstll"),
                             cut=cms.string(""),
                             name=cms.string("BToKstll"),
                             doc=cms.string("BToKstll Variable"),
                             singleton=cms.bool(False),
                             extension=cms.bool(False),
                             variables=cms.PSet(
                                 isEleCh=Var("userInt('isEleCh')", int, doc="electrons final state"),
                                 isKstCh=Var("userInt('isKstCh')", int, doc="Kstart final state"),

                                 lep1_index=Var("userInt('lep1_index')", int,doc="index of leading lepton in pfLepton collection"),
                                 lep2_index=Var("userInt('lep2_index')", int,doc="index of subleading lepton in pfLepton collection"),
                                 lep2_pfCand_index=Var("userInt('lep2_pfCand_index')", int,doc="index of subleading lepton in pfCandidate collection"),
                                 lep2_lostTrack_index=Var("userInt('lep2_lostTrack_index')", int,doc="index of subleading lepton in LostTrack collection"),
                                 kaon_index=Var("userInt('kaon_index')", int,doc="PFCand index of corresponding kaon"),
                                 kaon_lostTrack_index=Var("userInt('kaon_lostTrack_index')", int,doc="LostTrack index of corresponding kaon"),
                                 lep2_isPFLep=Var("userInt('lep2_isPFLep')", int,doc="flag is lepton2 from PFLepton collection"),
                                 lep2_isPFCand=Var("userInt('lep2_isPFCand')", int,doc="flag is lepton2 from PFCandidate collection"),
                                 kaon_isPFCand=Var("userInt('kaon_isPFCand')", int,doc="flag is kaon from PFCandidate collection"),
                                 pion_index=Var("userInt('pion_index')", int,doc="PFCand index of corresponding pion"),
                                 pion_lostTrack_index=Var("userInt('pion_lostTrack_index')", int,doc="LostTrack index of corresponding pion"),
                                 pion_isPFCand=Var("userInt('pion_isPFCand')", int,doc="flag is pion from PFCand"),
                                 
                                 lep1_pt=Var("userFloat('lep1_pt')", float,doc="pt of leading lepton"),
                                 lep1_eta=Var("userFloat('lep1_eta')", float,doc="eta of leading lepton"),
                                 lep1_phi=Var("userFloat('lep1_phi')", float,doc="phi of leading lepton"),
                                 lep1_charge=Var("userInt('lep1_charge')", int,doc="charge of leading lepton"),
                                 lep1_dxy=Var("userFloat('lep1_dxy')", float,doc="dxy of leading lepton (with sign) wrt first PV, in cm"),
                                 lep1_dxyS=Var("userFloat('lep1_dxyS')", float,doc="dxy/err of leading lepton (with sign) wrt first PV, in cm"),
                                 lep1_dz=Var("userFloat('lep1_dz')", float,doc="dz of leading lepton (with sign) wrt first PV, in cm"),
                                 lep1_dzS=Var("userFloat('lep1_dzS')", float,doc="dz/err of leading lepton (with sign) wrt first PV, in cm"),
                                 lep1_vz=Var("userFloat('lep1_vz')", float,doc="z coordinate of vertex position for lep1"),

                                 lep2_pt=Var("userFloat('lep2_pt')", float,doc="pt of subleading lepton"),
                                 lep2_eta=Var("userFloat('lep2_eta')", float,doc="eta of subleading lepton"),
                                 lep2_phi=Var("userFloat('lep2_phi')", float,doc="phi of subleading lepton"),
                                 lep2_charge=Var("userInt('lep2_charge')", int,doc="charge of subleading lepton"),
                                 lep2_dxy=Var("userFloat('lep2_dxy')", float,doc="dxy of subleading lepton (with sign) wrt first PV, in cm"),
                                 lep2_dxyS=Var("userFloat('lep2_dxyS')", float,doc="dxy/err of subleading lepton (with sign) wrt first PV, in cm"),
                                 lep2_dz=Var("userFloat('lep2_dz')", float,doc="dz of subleading lepton (with sign) wrt first PV, in cm"),
                                 lep2_dzS=Var("userFloat('lep2_dzS')", float,doc="dz/err of subleading lepton (with sign) wrt first PV, in cm"),
                                 lep2_vz=Var("userFloat('lep2_vz')", float,doc="z coordinate of vertex position for lep2"),

                                 kaon_pt=Var("userFloat('kaon_pt')", float,doc="pt of kaon"),
                                 kaon_eta=Var("userFloat('kaon_eta')", float,doc="eta of kaon"),
                                 kaon_phi=Var("userFloat('kaon_phi')", float,doc="phi of kaon"),
                                 kaon_charge=Var("userInt('kaon_charge')", int,doc="charge of kaon"),
                                 kaon_DCASig=Var("userFloat('kaon_DCASig')", float,doc="significance of xy-distance of closest approach kaon-beamspot"),
                                 kaon_dxy=Var("daughter('kaon').dxy()", float,doc="dxy of kaon"),
                                 kaon_dxyS=Var("userFloat('kaon_dxyS')", float,doc="dxy/err of kaon"),
                                 kaon_dz=Var("daughter('kaon').dz()", float,doc="dz of kaon"),
                                 kaon_dzS=Var("userFloat('kaon_dzS')", float,doc="dz/err of kaon"),
                                 kaon_vz=Var("userFloat('kaon_vz')", float,doc="z coordinate of vertex position for kaon"),
                                 kaon_dxy_wrtllVtx=Var("userFloat('kaon_dxy_wrtllVtx')", float,doc="impact parameter of kaon wrt ll refit vtx"),

                                 pion_pt=Var("userFloat('pion_pt')", float,doc="pt of pion"),
                                 pion_eta=Var("userFloat('pion_eta')", float,doc="eta of pion"),
                                 pion_phi=Var("userFloat('pion_phi')", float,doc="phi of pion"),
                                 pion_charge=Var("userInt('pion_charge')", int,doc="charge of pion"),
                                 pion_DCASig=Var("userFloat('pion_DCASig')", float,doc="significance of xy-distance of closest approach pion-beamspot"),
                                 pion_dxy=Var("userFloat('pion_dxy')", float,doc="dxy of pion"),
                                 pion_dxyS=Var("userFloat('pion_dxyS')", float,doc="dxy/err of pion"),
                                 pion_dz=Var("userFloat('pion_dz')", float,doc="dz of pion"),
                                 pion_dzS=Var("userFloat('pion_dzS')", float,doc="dz/err of pion"),
                                 pion_vz=Var("userFloat('pion_vz')", float,doc="z coordinate of vertex position for pion"),

                                 fitLepLep=Var("userInt('fitLepLep')", int, doc="flag 1 if lepton-lepton fit is performed"),
                                 fitLepLepPassed=Var("userInt('fitLepLep')", int, doc="flag 1 if lepton-lepton fit is ok"),
                                 ll_mass=Var("userFloat('ll_mass')", float,doc="dilepton mass"),
                                 ll_pt=Var("userFloat('ll_mass')", float,doc="dilepton pt"),
                                 ll_eta=Var("userFloat('ll_mass')", float,doc="dilepton eta"),
                                 ll_phi=Var("userFloat('ll_mass')", float,doc="dilepton phi"),
                                 ll_Lxy=Var("userFloat('ll_mass')", float,doc="significance of dilepton vertex-beamspot xy-separation"),
                                 ll_ctxy=Var("userFloat('ll_mass')", float,doc="dielectron vertex-beamspot xy-separation/pt"),
                                 ll_Chi2_vtx=Var("userFloat('ll_Chi2_vtx')", float,doc="dilepton vertex chi2"),
                                 ll_CL_vtx=Var("userFloat('ll_Chi2_vtx')", float,doc="dilepton vertex chi2 vertex probability"),
                                 maxl1l2_dxyS=Var("userFloat('maxl1l2_dxyS')", float,doc="max l1 l2 transverse displacement wrt PV"),

                                 Kst_mass=Var("userFloat('Kst_mass')", float,doc="K* mass (refitted)"),
                                 Kst_pt=Var("userFloat('Kst_pt')", float,doc="K* pt (refitted)"),
                                 Kst_eta=Var("userFloat('Kst_eta')", float,doc="K* eta (refitted)"),
                                 Kst_phi=Var("userFloat('Kst_phi')", float,doc="K* phi (refitted)"),
                                 #Kst_mass_err=Var("userFloat('Kst_mass_err')", float,doc="error on K* mass"),
                                 Kst_Lxy=Var("userFloat('Kst_Lxy')", float,doc="significance of K* vertex-beamspot xy-separation"),
                                 Kst_ctxy=Var("userFloat('Kst_ctxy')", float,doc="K* vertex-beamspot xy-separation/pt"),
                                 Kst_Chi2_vtx=Var("userFloat('Kst_Chi2_vtx')", float,doc="K* vertex chi2"),
                                 Kst_CL_vtx=Var("userFloat('Kst_CL_vtx')", float,doc="K* chi2 vertex probability"),

                                 B_pt=Var("userFloat('B_pt')", float,doc="pt of B candidate"),
                                 B_eta=Var("userFloat('B_eta')", float,doc="eta of B candidate"),
                                 B_phi=Var("userFloat('B_phi')", float,doc="phi of B candidate"),
                                 B_mass=Var("userFloat('B_mass')", float,doc="mass of B candidate"),
                                 B_Lxy=Var("userFloat('B_Lxy')", float,doc="significance of B vertex-beamspot xy-separation"),
                                 B_ctxy=Var("userFloat('B_ctxy')", float,doc="B vertex-beamspot xy-separation/pt"),
                                 B_Chi2_vtx=Var("userFloat('B_Chi2_vtx')", float,doc="B vertex chi2"),
                                 B_CL_vtx=Var("userFloat('B_CL_vtx')", float,doc="B chi2 vertex probability"),
                                 B_cosAlpha=Var("userFloat('B_cosAlpha')", float,doc="cosine of angle between B momentum and vertex-beamspot separation"),
                                 maxl1l2k_dxyS=Var("userFloat('maxl1l2k_dxyS')", float,doc="max l1 l2 k transverse displacement wrt PV"),
                                )
                             )

BToKstllSequence=cms.Sequence(BToKstll)
BToKstllTables=cms.Sequence(BToKstllTable)
