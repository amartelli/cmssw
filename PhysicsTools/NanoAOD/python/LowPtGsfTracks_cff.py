from PhysicsTools.NanoAOD.common_cff import *
import FWCore.ParameterSet.Config as cms


LowPtGsfTrackTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer",
                                    src=cms.InputTag("BToKstll:LowPtGsfTrack"),
                                    cut=cms.string(""),
                                    name=cms.string("LowPtGsfTrack"),
                                    doc=cms.string("LowPtGsfTrack Variable in pc format"),
                                    singleton=cms.bool(False),
                                    extension=cms.bool(False),
                                    variables=cms.PSet(pt = Var("userFloat('pt')",float,doc="pt Mode"),
                                                       eta = Var("userFloat('eta')",float,doc="eta Mode"),
                                                       phi = Var("userFloat('phi')",float,doc="phi Mode"),
                                                       charge = Var("userFloat('charge')",int,doc="charge Mode"),
                                                       seedBDT_unbiased = Var("userFloat('seedBDT_unbiased')", float,doc="seed BDT unbiased"),
                                                       seedBDT_ptbiased = Var("userFloat('seedBDT_ptbiased')", float,doc="seed BDT ptbiased"),
                                                       dz = Var("userFloat('dz')",float,doc="dz (with sign) wrt first PV, in cm"),
                                                       dxy = Var("userFloat('dxy')",float,doc="dxy (with sign) wrt first PV, in cm"),
                                                       dzS = Var("userFloat('dzS')", float, doc="dz/err (with sign) wrt first PV, in cm"),
                                                       dxyS = Var("userFloat('dxyS')", float, doc="dxy/err (with sign) wrt first PV, in cm"),
                                                       vz = Var("userFloat('vz')", float, doc="z coordinate of vertex position, in cm"),
                                                   ),
                                )



LowPtGsfTrkTables=cms.Sequence(LowPtGsfTrackTable)
