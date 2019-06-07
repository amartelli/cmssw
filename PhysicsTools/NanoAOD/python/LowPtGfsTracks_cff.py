from PhysicsTools.NanoAOD.common_cff import *
import FWCore.ParameterSet.Config as cms


LowPtGsfTrackTable=cms.EDProducer("SimpleGsfTrackFlatTableProducer",
                                  src=cms.InputTag("BToKstll:LowPtGsfTrack"),
                                  cut=cms.string(""),
                                  name=cms.string("LowPtGsfTrk"),
                                  doc=cms.string("LowPtGsfTrack Variable"),
                                  singleton=cms.bool(False),
                                  extension=cms.bool(False),
                                  variables=cms.PSet(pt = Var("ptMode()",float,doc="pt"),
                                                     eta = Var("etaMode()",float,doc="eta"),
                                                     phi = Var("phiMode()",float,doc="phi"),
                                                     charge = Var("chargeMode()",int,doc="charge"),
                                                 ),
                              )


LowPtGsfTrkTables=cms.Sequence(LowPtGsfTrackTable)
