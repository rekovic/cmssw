import FWCore.ParameterSet.Config as cms

l1tMuonJets = cms.EDProducer("L1TMuonJetProducer",
  
  MuonStubCollectionInputTag = cms.InputTag("simEmtfDigis"),

  L1TkMuStubCollectionInputTag = cms.InputTag("L1TkMuonStubS12"),
  
  L1EMTFTrackCollectionInputTag = cms.InputTag("simEmtfDigis"),
  L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),

  max_dR  = cms.double(0.9),
  max_dZ  = cms.double(1.0),

  debug   = cms.untracked.bool(True) 

)
