from PhysicsTools.NanoAOD.common_cff import *
import FWCore.ParameterSet.Config as cms


BToKpipi=cms.EDProducer("BToKpipiProducer",
                        beamSpot=cms.InputTag("offlineBeamSpot"),
                        PFCandCollection=cms.InputTag("packedPFCandidates"),
                        MinPt=cms.double(2.),
                        MinPtPiB=cms.double(2.),
                        MaxEta=cms.double(2.4),
                        MinDCASig=cms.double(2.),
                        D0MinMass=cms.double(0.),
                        D0MaxMass=cms.double(10.),
                        D0MinCLVtx=cms.double(0.),
                        D0ChargeCheck=cms.bool(True),
                        KaonChargeCheck=cms.bool(True),
                        BMinMass=cms.double(0.),
                        BMaxMass=cms.double(10.),
                        BMinCLVtx=cms.double(0.),
                        D0MassConstraint=cms.double(1.86484),
                        save3TrackRefit=cms.bool(True)
                        )

BToKpipiTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer",
                             src=cms.InputTag("BToKpipi"),
                             cut=cms.string(""),
                             name=cms.string("BToKpipi"),
                             doc=cms.string("BToKpipi Variable"),
                             singleton=cms.bool(False),
                             extension=cms.bool(False),
                             variables=cms.PSet(
                                kaon_pt=Var("userFloat('kaon_pt')", float,doc="pt of kaon (refitted)"),
                                kaon_eta=Var("userFloat('kaon_eta')", float,doc="eta of kaon (refitted)"),
                                kaon_phi=Var("userFloat('kaon_phi')", float,doc="phi of kaon (refitted)"),
                                kaon_charge=Var("userFloat('kaon_charge')", int,doc="charge of kaon"),
                                kaon_DCASig=Var("userFloat('kaon_DCASig')", float,doc="significance of xy-distance of closest approach kaon-beamspot"),
                                kaon_dxy=Var("daughter('kaon').dxy()", float,doc="dxy of kaon (not refitted)"),
                                kaon_dz=Var("daughter('kaon').dz()", float,doc="dz of kaon (not refitted)"),
                                piD0_pt=Var("userFloat('piD0_pt')", float,doc="pt of pion from D0->Kpi (refitted)"),
                                piD0_eta=Var("userFloat('piD0_eta')", float,doc="eta of pion from D0->Kpi (refitted)"),
                                piD0_phi=Var("userFloat('piD0_phi')", float,doc="phi of pion from D0->Kpi (refitted)"),
                                piD0_charge=Var("userFloat('piD0_charge')", int,doc="charge of pion from D0->Kpi"),
                                piD0_DCASig=Var("userFloat('piD0_DCASig')", float,doc="significance of xy-distance of closest approach piD0-beamspot"),
                                piD0_dxy=Var("daughter('piD0').dxy()", float,doc="dxy of kaon (not refitted)"),
                                piD0_dz=Var("daughter('piD0').dz()", float,doc="dz of kaon (not refitted)"),
                                piBu_pt=Var("userFloat('piBu_pt')", float,doc="pt of pion from Bu->D0pi (refitted)"),
                                piBu_eta=Var("userFloat('piBu_eta')", float,doc="eta of pion from Bu->D0pi (refitted)"),
                                piBu_phi=Var("userFloat('piBu_phi')", float,doc="phi of pion from Bu->D0pi (refitted)"),
                                piBu_charge=Var("userFloat('piBu_charge')", int,doc="charge of pion from Bu->D0pi"),
                                piBu_DCASig=Var("userFloat('piBu_DCASig')", float,doc="significance of xy-distance of closest approach piBu-beamspot"),
                                piBu_dxy=Var("daughter('piBu').dxy()", float,doc="dxy of piBu (not refitted)"),
                                piBu_dz=Var("daughter('piBu').dz()", float,doc="dz of piBu (not refitted)"),
                                Kpi_pt=Var("userFloat('Kpi_pt')", float,doc="D0->Kpi pt (refitted)"),
                                Kpi_eta=Var("userFloat('Kpi_eta')", float,doc="D0->Kpi eta (refitted)"),
                                Kpi_phi=Var("userFloat('Kpi_phi')", float,doc="D0->Kpi phi (refitted)"),
                                Kpi_mass=Var("userFloat('Kpi_mass')", float,doc="D0->Kpi mass (refitted)"),
                                Kpi_mass_err=Var("userFloat('Kpi_mass_err')", float,doc="error on D0->Kpi mass"),
                                Kpi_Lxy=Var("userFloat('Kpi_Lxy')", float,doc="significance of D0->Kpi vertex-beamspot xy-separation"),
                                Kpi_ctxy=Var("userFloat('Kpi_ctxy')", float,doc="D0->Kpi vertex-beamspot xy-separation/pt"),
                                Kpi_CL_vtx=Var("userFloat('Kpi_CL_vtx')", float,doc="D0->Kpi chi2 vertex probability"),
                                pt=Var("userFloat('pt')", float,doc="pt of BToKpipi candidate (2-trk refitted)"),
                                eta=Var("userFloat('eta')", float,doc="eta of BToKpipi candidate (2-trk refitted)"),
                                phi=Var("userFloat('phi')", float,doc="phi of BToKpipi candidate (2-trk refitted)"),
                                mass=Var("userFloat('mass')", float,doc="mass of BToKpipi candidate (2-trk refitted)"),
                                mass_err=Var("userFloat('mass_err')", float,doc="error on mass of BToKpipi candidate(2-trk refitted)"),
                                Lxy=Var("userFloat('Lxy')", float,doc="significance of BToKpipi vertex-beamspot xy-separation (2-trk refitted)"),
                                ctxy=Var("userFloat('ctxy')", float,doc="BToKpipi vertex-beamspot xy-separation/pt (2-trk refitted)"),
                                CL_vtx=Var("userFloat('CL_vtx')", float,doc="BToKpipi chi2 vertex probability (2-trk refitted)"),
                                cosAlpha=Var("userFloat('cosAlpha')", float,doc="cosine of angle between BToKpipi momentum and vertex-beamspot separation (2-trk refitted)"),
                                pt_3trk=Var("userFloat('pt_3trk')", float,doc="pt of BToKpipi candidate (3-trk refitted)"),
                                eta_3trk=Var("userFloat('eta_3trk')", float,doc="eta of BToKpipi candidate (3-trk refitted)"),
                                phi_3trk=Var("userFloat('phi_3trk')", float,doc="phi of BToKpipi candidate (3-trk refitted)"),
                                mass_3trk=Var("userFloat('mass_3trk')", float,doc="mass of BToKpipi candidate (3-trk refitted)"),
                                mass_err_3trk=Var("userFloat('mass_err_3trk')", float,doc="error on mass of BToKpipi candidate(3-trk refitted)"),
                                Lxy_3trk=Var("userFloat('Lxy_3trk')", float,doc="significance of BToKpipi vertex-beamspot xy-separation (3-trk refitted)"),
                                ctxy_3trk=Var("userFloat('ctxy_3trk')", float,doc="BToKpipi vertex-beamspot xy-separation/pt (3-trk refitted)"),
                                CL_vtx_3trk=Var("userFloat('CL_vtx_3trk')", float,doc="BToKpipi chi2 vertex probability (3-trk refitted)"),
                                cosAlpha_3trk=Var("userFloat('cosAlpha_3trk')", float,doc="cosine of angle between BToKpipi momentum and vertex-beamspot separation (3-trk refitted)"),
                                )
                             )

BToKpipiSequence=cms.Sequence(BToKpipi)
BToKpipiTables=cms.Sequence(BToKpipiTable)
