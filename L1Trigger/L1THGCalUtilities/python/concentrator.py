import FWCore.ParameterSet.Config as cms
import SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi as digiparam

def create_supertriggercell(process, inputs):
    producer = process.hgcalConcentratorProducer.clone() 
    producer.ProcessorParameters.Method = cms.string('superTriggerCellSelect')
    producer.InputTriggerCells = cms.InputTag('{}:HGCalVFEProcessorSums'.format(inputs))
    producer.InputTriggerSums = cms.InputTag('{}:HGCalVFEProcessorSums'.format(inputs))
    return producer


def create_bestchoice(process, inputs,
       triggercells=12
        ):
    adcSaturationBH_MIP = digiparam.hgchebackDigitizer.digiCfg.feCfg.adcSaturation_fC
    adcNbitsBH = digiparam.hgchebackDigitizer.digiCfg.feCfg.adcNbits
    producer = process.hgcalConcentratorProducer.clone() 
    producer.ProcessorParameters.Method = cms.string('bestChoiceSelect')
    producer.ProcessorParameters.NData = cms.uint32(triggercells)
    producer.ProcessorParameters.linLSB = cms.double(100./1024.)
    producer.ProcessorParameters.adcsaturationBH = adcSaturationBH_MIP
    producer.ProcessorParameters.adcnBitsBH = adcNbitsBH
    producer.InputTriggerCells = cms.InputTag('{}:HGCalVFEProcessorSums'.format(inputs))
    producer.InputTriggerSums = cms.InputTag('{}:HGCalVFEProcessorSums'.format(inputs))
    return producer
