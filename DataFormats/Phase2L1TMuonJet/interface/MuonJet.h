// Class to hold MuonJet objects for PhaseII Level1 Trigger
//
// Author:  Vladimir Rekovic  2019.10.05
#ifndef DataFormats_Phase2L1TMuonJet_MuonJet_H
#define DataFormats_Phase2L1TMuonJet_MuonJet_H

#include <vector>
#include <Math/VectorUtil.h>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/L1TMuon/interface/EMTFHit.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"

namespace l1t {

  class MuonJet : public  L1Candidate
  {

    public:
      typedef TTTrack< Ref_Phase2TrackerDigi_ >  L1TTTrackType;
      typedef std::vector< L1TTTrackType > L1TTTrackCollection;
      typedef edm::Ref<L1TTTrackCollection> L1TTTrackRef;
      typedef std::vector<edm::Ref<L1TTTrackCollection > > L1TTrackRefVector;

      // default constructor
      MuonJet(); 

      MuonJet(const vector<float> vpt, const vector<float> veta, const vector<float> vphi, const vector<float> zvtx); 
      MuonJet(const L1TkMuonParticle &, const L1TkMuonParticle &, const L1TkMuonParticle &);
      MuonJet(const L1TkMuonParticle &, const L1TkMuonParticle &, const EMTFHit &);

      virtual ~MuonJet() {}

      const vector<float> & getPt() const    { return pt_; }
      const vector<float> & getEta() const   { return eta_; }
      const vector<float> & getPhi() const   { return phi_; }
      const vector<float> & getZVtx() const  { return zvtx_; }

      void configure(float const& maxDeltaEta, float const& maxDeltaPhi, float const& maxDeltaR, float const& maxDeltaZ);
      void process();
      bool isValid();

      void print();

    private:

      int n_muons_;
      std::vector<float> pt_; //with charge sign
      std::vector<float> eta_;
      std::vector<float> phi_;
      std::vector<float> zvtx_;

      float maxDeltaPt_;
      float maxDeltaEta_;
      float maxDeltaPhi_;
      float maxDeltaR_;
      float maxDeltaZ_;
      float maxDeltaM_;

      float deltaPt_;
      float deltaEta_;
      float deltaPhi_;
      float deltaR_;
      float deltaZ_;
      float deltaM_;

      float mass_;
      int totalCharge_;
  

      EMTFHitRefVector stubs_;
      L1TTrackRefVector l1tracks_;

  };

}

#endif
