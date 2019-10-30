import FWCore.ParameterSet.Config as cms

l1tMuonJets = cms.EDProducer("L1TMuonJetProducer",
  
  MuonStubCollectionInputTag = cms.InputTag("simEmtfDigis"),

  L1TkMuStubCollectionInputTag = cms.InputTag("L1TkMuonStubS12"),
  
  L1EMTFTrackCollectionInputTag = cms.InputTag("simEmtfDigis"),
  L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),

  max_dR_3_TKMUSTUB  = cms.double(0.6),
  max_dZ_3_TKMUSTUB  = cms.double(0.6),

  max_dR_2_TKMUSTUB_1_MUSTUB  = cms.double(0.6),
  max_dZ_2_TKMUSTUB_1_MUSTUB  = cms.double(0.6),

  max_dR_1_TKMUSTUB_2_MUSTUB  = cms.double(0.7),
  max_dZ_1_TKMUSTUB_2_MUSTUB  = cms.double(0.85),

  stubStation = cms.uint32(1),
  stubSubsystem = cms.uint32(4),
  stubMinQuality = cms.uint32(3),
  debug   = cms.untracked.bool(True) 

)

l1tMuonJetsDebug = l1tMuonJets.clone()
l1tMuonJetsDebug.stubMinQuality = cms.uint32(1)
l1tMuonJetsDebug.debug   = cms.untracked.bool(True)
