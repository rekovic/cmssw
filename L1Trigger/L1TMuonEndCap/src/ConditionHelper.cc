#include "L1Trigger/L1TMuonEndCap/interface/ConditionHelper.h"

#include "CondFormats/L1TObjects/interface/L1TMuonEndCapParams.h"
#include "CondFormats/DataRecord/interface/L1TMuonEndcapParamsRcd.h"

#include "CondFormats/L1TObjects/interface/L1TMuonEndCapForest.h"
#include "CondFormats/DataRecord/interface/L1TMuonEndCapForestRcd.h"

#include "L1Trigger/L1TMuonEndCap/interface/PtAssignmentEngine.h"

#include "FWCore/Framework/interface/EventSetup.h"

ConditionHelper::ConditionHelper():
  params_cache_id_(0ULL), forest_cache_id_(0ULL) {
}

ConditionHelper::~ConditionHelper() {
}

void ConditionHelper::checkAndUpdateConditions(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Pull configuration from the EventSetup
  auto& params_setup = iSetup.get<L1TMuonEndcapParamsRcd>();
  if (params_setup.cacheIdentifier() != params_cache_id_) {
    params_setup.get(params_);

    // with the magic above you can use params_->fwVersion to change emulator's behavior
    // ...

    // reset cache id
    params_cache_id_ = params_setup.cacheIdentifier();
  }

  // Pull pt LUT from the EventSetup
  auto& forest_setup = iSetup.get<L1TMuonEndCapForestRcd>();
  if (forest_setup.cacheIdentifier() != forest_cache_id_) {
    forest_setup.get(forest_);

    // at this point we want to reload the newly pulled pT LUT
    // ...

    // reset cache id
    forest_cache_id_ = forest_setup.cacheIdentifier();
  }

  // Debug
  //std::cout << "Run number: " << iEvent.id().run() << " fw_version: " << get_fw_version()
  //    << " pt_lut_version: " << get_pt_lut_version() << " pc_lut_version: " << get_pc_lut_version()
  //    << std::endl;
}

unsigned int ConditionHelper::get_fw_version() const {
  return params_->firmwareVersion_;
}

unsigned int ConditionHelper::get_pt_lut_version() const {
  return params_->PtAssignVersion_;
}

unsigned int ConditionHelper::get_pc_lut_version() const {
  // Not yet implemented in O2O; default to coordinate LUTs from beginning of 2017
  return 2;
  // Requires change to CondFormats/L1TObjects/interface/L1TMuonEndCapParams.h
  // return params_->PrimConvVersion_;  
}
