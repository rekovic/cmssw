import FWCore.ParameterSet.Config as cms

def custom_ntuples_V9(process):
    ntuples = process.hgcalTriggerNtuplizer.Ntuples
    for ntuple in ntuples:
        if ntuple.NtupleName=='HGCalTriggerNtupleHGCDigis' or \
           ntuple.NtupleName=='HGCalTriggerNtupleHGCTriggerCells':
            ntuple.bhSimHits = cms.InputTag('g4SimHits:HGCHitsHEback')
    return process



def create_ntuple(process, inputs,
        ntuple_list=[
            'event',
            'gen', 'genjet', 'gentau',
            'digis',
            'triggercells',
            'clusters', 'multicluster',
            'tower'
            ]
        ):
    vpset = []
    for ntuple in ntuple_list:
        pset = getattr(process, 'ntuple_'+ntuple).clone()
        if ntuple=='triggercells':
            pset.TriggerCells = cms.InputTag('{}:HGCalConcentratorProcessorSelection'.format(inputs[0]))
            pset.Multiclusters = cms.InputTag('{}:HGCalBackendLayer2Processor3DClustering'.format(inputs[2]))
        elif ntuple=='clusters':
            pset.Clusters = cms.InputTag('{}:HGCalBackendLayer1Processor2DClustering'.format(inputs[1]))
            pset.Multiclusters = cms.InputTag('{}:HGCalBackendLayer2Processor3DClustering'.format(inputs[2]))
        elif ntuple=='multiclusters':
            pset.Multiclusters = cms.InputTag('{}:HGCalBackendLayer2Processor3DClustering'.format(inputs[2]))
        vpset.append(pset)
    ntuplizer = process.hgcalTriggerNtuplizer.clone()
    ntuplizer.Ntuples = cms.VPSet(*vpset)
    print(ntuplizer.Ntuples)
    return ntuplizer




