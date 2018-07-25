import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import *

BToKstee=cms.EDProducer("BToKsteeProducer",
                        beamSpot=cms.InputTag("offlineBeamSpot"),
                        electronCollection=cms.InputTag("linkedObjects","electrons"), #same collection as in NanoAOD
                        PFCandCollection=cms.InputTag("packedPFCandidates"),
                        ElectronMinPt=cms.double(1.),
                        ElectronMaxEta=cms.double(2.4),
                        KaonMinPt=cms.double(1.),
                        KaonMaxEta=cms.double(2.4),
                        KaonMinDCASig=cms.double(3.3),
                        PionMinPt=cms.double(1.),
                        PionMaxEta=cms.double(2.4),
                        PionMinDCASig=cms.double(3.3),
                        DiElectronChargeCheck=cms.bool(False),
                        KstarChargeCheck=cms.bool(True),
                        JPsiMassConstraint=cms.double(-1), #2-trk refitting uses measured di-ele mass
                        KstMassConstraint=cms.double(0.89176), #2-trk refitting uses nominal K*(892) mass
                        save2TrackRefit=cms.bool(True),
                        save4TrackRefit=cms.bool(True)
                      )

BToKsteeTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
                             src=cms.InputTag("BToKstee"),
                             cut=cms.string(""),
                             name=cms.string("BToKstee"),
                             doc=cms.string("BToKstee Variable"),
                             singleton=cms.bool(False),
                             extension=cms.bool(False),
                             variables=cms.PSet(
                                ele1_pt=Var("userFloat('ele1_pt')", float,doc="pt of leading electron (refitted)"),
                                ele1_eta=Var("userFloat('ele1_eta')", float,doc="eta of leading electron (refitted)"),
                                ele1_phi=Var("userFloat('ele1_phi')", float,doc="phi of leading electron (refitted)"),
                                ele1_charge=Var("userFloat('ele1_charge')", int,doc="charge of leading electron"),
                                ele1_dxy=Var("daughter('ele1').dB('PV2D')", float,doc="dxy of leading electron (with sign) wrt first PV, in cm"),
                                ele1_dz=Var("daughter('ele1').dB('PVDZ')", float,doc="dz of leading electron (with sign) wrt first PV, in cm"),
                                ele2_pt=Var("userFloat('ele2_pt')", float,doc="pt of subleading electron (refitted)"),
                                ele2_eta=Var("userFloat('ele2_eta')", float,doc="eta of subleading electron (refitted)"),
                                ele2_phi=Var("userFloat('ele2_phi')", float,doc="phi of subleading electron (refitted)"),
                                ele2_charge=Var("userFloat('ele2_charge')", int,doc="charge of subleading electron"),
                                ele2_dxy=Var("daughter('ele2').dB('PV2D')", float,doc="dxy of subleading electron (with sign) wrt first PV, in cm"),
                                ele2_dz=Var("daughter('ele2').dB('PVDZ')", float,doc="dz of subleading electron (with sign) wrt first PV, in cm"),
                                kaon_pt=Var("userFloat('kaon_pt')", float,doc="pt of kaon (refitted)"),
                                kaon_eta=Var("userFloat('kaon_eta')", float,doc="eta of kaon (refitted)"),
                                kaon_phi=Var("userFloat('kaon_phi')", float,doc="phi of kaon (refitted)"),
                                kaon_charge=Var("userFloat('kaon_charge')", int,doc="charge of kaon"),
                                kaon_DCASig=Var("userFloat('kaon_DCASig')", float,doc="significance of xy-distance of closest approach kaon-beamspot"),
                                kaon_dxy=Var("daughter('kaon').dxy()", float,doc="dxy of kaon (not refitted)"),
                                kaon_dz=Var("daughter('kaon').dz()", float,doc="dz of kaon (not refitted)"),
                                pion_pt=Var("userFloat('pion_pt')", float,doc="pt of kaon (refitted)"),
                                pion_eta=Var("userFloat('pion_eta')", float,doc="eta of kaon (refitted)"),
                                pion_phi=Var("userFloat('pion_phi')", float,doc="phi of kaon (refitted)"),
                                pion_charge=Var("userFloat('pion_charge')", int,doc="charge of kaon"),
                                pion_DCASig=Var("userFloat('pion_DCASig')", float,doc="significance of xy-distance of closest approach kaon-beamspot"),
                                pion_dxy=Var("daughter('pion').dxy()", float,doc="dxy of kaon (not refitted)"),
                                pion_dz=Var("daughter('pion').dz()", float,doc="dz of kaon (not refitted)"),
                                ee_pt=Var("userFloat('ee_pt')", float,doc="dielectron pt (refitted)"),
                                ee_eta=Var("userFloat('ee_eta')", float,doc="dielectron eta (refitted)"),
                                ee_phi=Var("userFloat('ee_phi')", float,doc="dielectron phi (refitted)"),
                                ee_mass=Var("userFloat('ee_mass')", float,doc="dielectron mass (refitted)"),
                                ee_mass_err=Var("userFloat('ee_mass_err')", float,doc="error on dielectron mass"),
                                ee_Lxy=Var("userFloat('ee_Lxy')", float,doc="significance of dielectron vertex-beamspot xy-separation"),
                                ee_ctxy=Var("userFloat('ee_ctxy')", float,doc="dielectron vertex-beamspot xy-separation/pt"),
                                ee_CL_vtx=Var("userFloat('ee_CL_vtx')", float,doc="dielectron chi2 vertex probability"),
                                Kst_pt=Var("userFloat('Kst_pt')", float,doc="K* pt (refitted)"),
                                Kst_eta=Var("userFloat('Kst_eta')", float,doc="K* eta (refitted)"),
                                Kst_phi=Var("userFloat('Kst_phi')", float,doc="K* phi (refitted)"),
                                Kst_mass=Var("userFloat('Kst_mass')", float,doc="K* mass (refitted)"),
                                Kst_mass_err=Var("userFloat('Kst_mass_err')", float,doc="error on K* mass"),
                                Kst_Lxy=Var("userFloat('Kst_Lxy')", float,doc="significance of K* vertex-beamspot xy-separation"),
                                Kst_ctxy=Var("userFloat('Kst_ctxy')", float,doc="K* vertex-beamspot xy-separation/pt"),
                                Kst_CL_vtx=Var("userFloat('Kst_CL_vtx')", float,doc="K* chi2 vertex probability"),
                                pt=Var("userFloat('pt')", float,doc="pt of BToKstee candidate (3-trk refitted)"),
                                eta=Var("userFloat('eta')", float,doc="eta of BToKstee candidate (3-trk refitted)"),
                                phi=Var("userFloat('phi')", float,doc="phi of BToKstee candidate (3-trk refitted)"),
                                mass=Var("userFloat('mass')", float,doc="mass of BToKstee candidate (3-trk refitted)"),
                                mass_err=Var("userFloat('mass_err')", float,doc="error on mass of BToKstee candidate (3-trk refitted)"),
                                Lxy=Var("userFloat('Lxy')", float,doc="significance of BToKstee vertex-beamspot xy-separation (3-trk refitted)"),
                                ctxy=Var("userFloat('ctxy')", float,doc="BToKstee vertex-beamspot xy-separation/pt (3-trk refitted)"),
                                CL_vtx=Var("userFloat('CL_vtx')", float,doc="BToKstee chi2 vertex probability (3-trk refitted)"),
                                cosAlpha=Var("userFloat('cosAlpha')", float,doc="cosine of angle between BToKmumu momentum and vertex-beamspot separation (3-trk refitted)"),
                                pt_2trk=Var("userFloat('pt_2trk')", float,doc="pt of BToKstee candidate (2-trk refitted)"),
                                eta_2trk=Var("userFloat('eta_2trk')", float,doc="eta of BToKstee candidate (2-trk refitted)"),
                                phi_2trk=Var("userFloat('phi_2trk')", float,doc="phi of BToKstee candidate (2-trk refitted)"),
                                mass_2trk=Var("userFloat('mass_2trk')", float,doc="mass of BToKstee candidate (2-trk refitted)"),
                                mass_err_2trk=Var("userFloat('mass_err_2trk')", float,doc="error on mass of BToKstee candidate (2-trk refitted)"),
                                Lxy_2trk=Var("userFloat('Lxy_2trk')", float,doc="significance of BToKstee vertex-beamspot xy-separation (2-trk refitted)"),
                                ctxy_2trk=Var("userFloat('ctxy_2trk')", float,doc="BToKstee vertex-beamspot xy-separation/pt (2-trk refitted)"),
                                CL_vtx_2trk=Var("userFloat('CL_vtx_2trk')", float,doc="BToKstee chi2 vertex probability (2-trk refitted)"),
                                cosAlpha_2trk=Var("userFloat('cosAlpha_2trk')", float,doc="cosine of angle between BToKmumu momentum and vertex-beamspot separation (2-trk refitted)"),
                                pt_4trk=Var("userFloat('pt_4trk')", float,doc="pt of BToKstee candidate (4-trk refitted)"),
                                eta_4trk=Var("userFloat('eta_4trk')", float,doc="eta of BToKstee candidate (4-trk refitted)"),
                                phi_4trk=Var("userFloat('phi_4trk')", float,doc="phi of BToKstee candidate (4-trk refitted)"),
                                mass_4trk=Var("userFloat('mass_4trk')", float,doc="mass of BToKstee candidate (4-trk refitted)"),
                                mass_err_4trk=Var("userFloat('mass_err_4trk')", float,doc="error on mass of BToKstee candidate (4-trk refitted)"),
                                Lxy_4trk=Var("userFloat('Lxy_4trk')", float,doc="significance of BToKstee vertex-beamspot xy-separation (4-trk refitted)"),
                                ctxy_4trk=Var("userFloat('ctxy_4trk')", float,doc="BToKstee vertex-beamspot xy-separation/pt (4-trk refitted)"),
                                CL_vtx_4trk=Var("userFloat('CL_vtx_4trk')", float,doc="BToKstee chi2 vertex probability (4-trk refitted)"),
                                cosAlpha_4trk=Var("userFloat('cosAlpha_4trk')", float,doc="cosine of angle between BToKmumu momentum and vertex-beamspot separation (4-trk refitted)"),
                                )
                             )

BToKsteeSequence=cms.Sequence(BToKstee)
BToKsteeTables=cms.Sequence(BToKsteeTable)
