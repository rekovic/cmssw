import FWCore.ParameterSet.Config as cms 
from Configuration.StandardSequences.Eras import eras
from Configuration.ProcessModifiers.convertHGCalDigisSim_cff import convertHGCalDigisSim

# For old samples use the digi converter
#process = cms.Process('DIGI',eras.Phase2,convertHGCalDigisSim)
process = cms.Process('DIGI',eras.Phase2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)

# Input source
process.source = cms.Source("PoolSource",
       fileNames = cms.untracked.vstring('/store/relval/CMSSW_10_4_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/103X_upgrade2023_realistic_v2_2023D21noPU-v1/20000/F4344045-AEDE-4240-B7B1-27D2CF96C34E.root'),
       inputCommands=cms.untracked.vstring(
           'keep *',
           'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
           'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
           'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
           'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
           'drop l1tEMTFTrack2016s_simEmtfDigis__HLT',
           )
       )

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('SingleElectronPt10_cfi nevts:10'),
    name = cms.untracked.string('Applications')
)

# Output definition
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("ntuple.root")
    )

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# load HGCAL TPG simulation
process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')
from L1Trigger.L1THGCalUtilities.hgcalTriggerChain import HGCalTriggerChain
import L1Trigger.L1THGCalUtilities.vfe as vfe
import L1Trigger.L1THGCalUtilities.concentrator as concentrator
import L1Trigger.L1THGCalUtilities.clustering2d as clustering2d
import L1Trigger.L1THGCalUtilities.clustering3d as clustering3d

chain = HGCalTriggerChain()
chain.register_vfe("Floatingpoint7", lambda p : vfe.create_compression(p, 4, 3, True))
chain.register_concentrator("Supertriggercell", concentrator.create_supertriggercell)
chain.register_concentrator("Bestchoice", lambda p,i : concentrator.create_bestchoice(p,i, triggercells=12))
chain.register_backend1("Dummy", clustering2d.create_dummy)
chain.register_backend2("Histothreshold", clustering3d.create_histoThreshold)

chain.register_chain('Floatingpoint7', 'Supertriggercell', 'Dummy', 'Histothreshold')
chain.register_chain('Floatingpoint7', 'Bestchoice', 'Dummy', 'Histothreshold')

process = chain.create_sequences(process)

# Remove towers from sequence
process.hgcalTriggerPrimitives.remove(process.hgcalTowerMap)
process.hgcalTriggerPrimitives.remove(process.hgcalTower)

process.hgcl1tpg_step = cms.Path(process.hgcalTriggerPrimitives)
print(process.hgcl1tpg_step)

# load ntuplizer
#process.load('L1Trigger.L1THGCalUtilities.hgcalTriggerNtuples_cff')
#process.ntuple_multicluster_supertriggercell =  process.ntuple_multicluster.clone()
#process.ntuple_multicluster_supertriggercell.Multiclusters = cms.InputTag('hgcalStage2SuperTriggerCell:HGCalBackendLayer2Processor3DClustering')
#process.ntuple_multicluster_bestchoice =  process.ntuple_multicluster.clone()
#process.ntuple_multicluster_bestchoice.Multiclusters = cms.InputTag('hgcalStage2BestChoice:HGCalBackendLayer2Processor3DClustering')
#
#process.hgcalTriggerNtuplizerSuperTriggerCell = process.hgcalTriggerNtuplizer.clone()
#process.hgcalTriggerNtuplizerSuperTriggerCell.Ntuples = cms.VPSet(process.ntuple_multicluster_supertriggercell)
#
#process.hgcalTriggerNtuplizerBestChoice = process.hgcalTriggerNtuplizer.clone()
#process.hgcalTriggerNtuplizerBestChoice.Ntuples = cms.VPSet(process.ntuple_multicluster_bestchoice)
#
#ntuples = cms.Sequence(process.hgcalTriggerNtuplizerSuperTriggerCell*process.hgcalTriggerNtuplizerBestChoice)
#process.globalReplace("hgcalTriggerNtuples",ntuples)
#
#process.ntuple_step = cms.Path(process.hgcalTriggerNtuples)

# Schedule definition
process.schedule = cms.Schedule(process.hgcl1tpg_step)#, process.ntuple_step)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

