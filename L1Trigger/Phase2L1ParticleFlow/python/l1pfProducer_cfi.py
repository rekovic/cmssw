import FWCore.ParameterSet.Config as cms

from math import sqrt

l1pfProducer = cms.EDProducer("L1TPFProducer",
     tracks = cms.InputTag('pfTracksFromL1Tracks'),
     muons = cms.InputTag('simGmtStage2Digis',),
     emClusters = cms.VInputTag(cms.InputTag('pfClustersFromHGC3DClustersEM'), cms.InputTag('pfClustersFromL1EGClusters')),
     hadClusters = cms.VInputTag(cms.InputTag('pfClustersFromCombinedCalo:calibrated')),
     emPtCut  = cms.double(0.5),
     hadPtCut = cms.double(1.0),
     trkPtCut    = cms.double(2.0),
     trkMinStubs = cms.uint32(4),
     trkMaxChi2  = cms.double(15),
     etaCharged  = cms.double(2.5),
     puppiDr     = cms.double(0.3),
     puppiPtCut  = cms.double(4.0), # for old algo
     puppiEtaCuts       = cms.vdouble(1.5, 2.5, 3.0, 5.5), 
     puppiPtCuts        = cms.vdouble(0.0, 3.0, 6.0, 8.0),
     puppiPtCutsPhotons = cms.vdouble(0.0, 3.0, 6.0, 8.0),
     vtxRes      = cms.double(0.333),
     vtxAlgo     = cms.string("TP"),
     vtxAdaptiveCut = cms.bool(True),
     linking = cms.PSet(
        algo = cms.string("PFAlgo3"),
        # track -> mu linking configurables
        trackMuDR    = cms.double(0.2), # accounts for poor resolution of standalone, and missing propagations
        trackMuMatch = cms.string("boxBestByPtRatio"), # also drBestByPtRatio
        # track -> em linking configurables
        trackEmDR   = cms.double(0.04), # 1 Ecal crystal size is 0.02, and ~2 cm in HGCal is ~0.007
        trackEmUseAlsoTrackSigma = cms.bool(True), # also use the track uncertainty for electron linking
        # em -> calo linking configurables
        emCaloDR    = cms.double(0.10),    # 1 Hcal tower size is ~0.09
        caloEmPtMinFrac = cms.double(0.5), # Calo object must have an EM Et at least half of that of the EM cluster to allow linking
        emCaloUseAlsoCaloSigma = cms.bool(True), # also use the track uncertainty for electron linking
        # track -> calo linking configurables
        trackCaloLinkMetric = cms.string("bestByDRPt"),
        #trackCaloLinkMetric = cms.string("bestByDR"),
        trackCaloDR = cms.double(0.12),
        trackCaloNSigmaLow  = cms.double(2.0),
        trackCaloNSigmaHigh = cms.double(sqrt(1.5)), # optimization gave 1.2, but sqrt(1.5) is easier to implement in hardware
        useTrackCaloSigma = cms.bool(True), # take the uncertainty on the calo cluster from the track, for linking purposes
        sumTkCaloErr2 = cms.bool(True), # add up track calo errors in quadrature instead of linearly
        rescaleTracks = cms.bool(False), # if tracks exceed the calo, rescale the track momenta
        useCaloTrkWeightedAverage = cms.bool(False), # do the weighted average of track & calo pTs if it's a 1-1 link 
        # how to deal with unlinked tracks
        maxInvisiblePt = cms.double(10.0), # max allowed pt of a track with no calo energy
        tightTrackMinStubs = cms.uint32(6),
        tightTrackMaxChi2  = cms.double(50),
        tightTrackMaxInvisiblePt = cms.double(20),
        # how to deal with neutrals
        ecalPriority  = cms.bool(True), # take first ecal energy when making neutrals
        # other features not turned on: reliniking of neutrals to track-matched calo clusters with track excess
        caloReLink  = cms.bool(False),
        caloReLinkDR = cms.double(0.3),
        caloReLinkThreshold = cms.double(0.5),
        # other features not turned on: matching too high pt tracks to calo but rescaling track pt (not implemented in PFAlgo3)
        rescaleUnmatchedTrack = cms.bool(False),
     ),
     debug = cms.untracked.int32(0),
)


