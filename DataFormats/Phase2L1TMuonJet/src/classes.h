#include "Rtypes.h"

#include "DataFormats/Phase2L1TMuonJet/interface/MuonJet.h"
#include "DataFormats/Phase2L1TMuonJet/interface/MuonJetFwd.h"

namespace DataFormats_Phase2L1TMuonJet
{
  struct dictionary 
  {
    l1t::MuonJet l1muonjet;
    l1t::MuonJetCollection l1muonjetCollection;
    edm::Wrapper<l1t::MuonJetCollection> l1muonjetCWrapper;
  };
}
