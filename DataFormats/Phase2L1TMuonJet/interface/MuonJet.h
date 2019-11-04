// Class to hold MuonJet objects for PhaseII Level1 Trigger
//
// Author:  Vladimir Rekovic  2019.10.05
#ifndef DataFormats_Phase2L1TMuonJet_MuonJet_H
#define DataFormats_Phase2L1TMuonJet_MuonJet_H

#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    
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

      enum MuonType
      {
        k_NONE = 0,
        k_MUSTUB = 1,
        k_TKMUSTUB =2 
       };

      enum MuonJetType
      {
        ONE_MUSTUB,
        ONE_TKMUSTUB,
        TWO_MUSTUB,
        TWO_TKMUSTUB,
        THREE_MUSTUB,
        ONE_TKMUSTUB_TWO_MUSTUB,
        TWO_TKMUSTUB_ONE_MUSTUB,
        THREE_TKMUSTUB
      };

      // default constructor
      MuonJet(); 

      MuonJet(const vector<float> vpt, const vector<float> veta, const vector<float> vphi, const vector<float> zvtx); 

      MuonJet(const L1TkMuonParticle &, const L1TkMuonParticle &, const L1TkMuonParticle &);
      MuonJet(const L1TkMuonParticle &, const L1TkMuonParticle &);
      MuonJet(const L1TkMuonParticle &, const L1TkMuonParticle &, const EMTFHit &);
      MuonJet(const L1TkMuonParticle & ,const L1TkMuonParticle &, const EMTFHitRef);
      MuonJet(const L1TkMuonParticle &, const EMTFHit &, const EMTFHit &);
      MuonJet(const L1TkMuonParticle &, const EMTFHitRef, const EMTFHitRef);
      MuonJet(const EMTFHit &, const EMTFHit &, const EMTFHit &);
      MuonJet(const EMTFHitRef, const EMTFHitRef,const EMTFHitRef);

      virtual ~MuonJet() {}

      const vector<float> & getPt() const    { return pt_; }
      const vector<float> & getEta() const   { return eta_; }
      const vector<float> & getPhi() const   { return phi_; }
      const vector<float> & getZVtx() const  { return zvtx_; }

      const vector<int> & getCharge()   const  { return charge_; }
      const vector<int> & getMuonType()   const  { return muonType_; }
      
      const vector<int> & getStubQuality()  const  { return stubQual_; }

      const vector<float> & getStubEta() const   { return stubEta_; }
      const vector<float> & getStubPhi() const   { return stubPhi_; }
      const vector<float> & getStubBend() const   { return stubBend_; }
      const std::vector<EMTFHitRef> & getStubRef() const   { return stubRef_; }
      const std::vector<EMTFHitRefVector> & getStubRefs() const   { return stubRefs_; }

      float getMass()         const { return mass_; }
      float getMass12()       const { return mass12_; }

      float getDeltaR()       const { return deltaR_; }
      float getDeltaEta()       const { return deltaEta_; }
      float getDeltaPhi()       const { return deltaPhi_; }
      float getDeltaPhiAtVtx()       const { return deltaPhiAtVtx_; }
      float getDeltaPt()       const { return deltaPt_; }
      float getDeltaZ()       const { return deltaZ_; }

      float getDeltaRMin()       const { return deltaRMin_; }
      float getDeltaPtMin()       const { return deltaPtMin_; }
      float getDeltaEtaMin()       const { return deltaEtaMin_; }
      float getDeltaPhiMin()       const { return deltaPhiMin_; }
      float getDeltaPhiAtVtxMin()       const { return deltaPhiAtVtxMin_; }
      float getDeltaZMin()       const { return deltaZMin_; }

      int   getTotalCharge()    const { return totalCharge_; }
      MuonJetType getType()   const { return type_; }
      float getPtFromBendingAngle(const EMTFHit & muStub);
      float getPhiAtVertex(const EMTFHit & muStub);
      int   getChargeFromBendingAngle(const EMTFHit & muStub);


      void configure(float const& maxDeltaEta, float const& maxDeltaPhi, float const& maxDeltaR, float const& maxDeltaZ);
      void process();
      void sortWithPt();
      bool isValid();
      bool areValidStubs();
      bool isValidTkMuStub(const EMTFHitRef);
      bool isValidMuStub(const EMTFHitRef);
      bool const areUniqueStubs ();
      bool const areUniqueTkMuStubs ();
      bool const isUniqueStubRefs ();
      void print();

    private:

      MuonJetType type_;

      int n_muons_;

      std::vector<float> pt_; //with charge sign
      std::vector<float> eta_;
      std::vector<float> phi_;
      std::vector<float> zvtx_;

      std::vector<int> charge_;

      std::vector<int> muonType_; // TKMUSTUB, MUSTUB

      std::vector<int> stubQual_;

      std::vector<float> stubEta_;
      std::vector<float> stubPhi_;
      std::vector<float> stubBend_;
      std::vector<EMTFHitRef> stubRef_; // this holds the principle stubRef
      std::vector<EMTFHitRefVector> stubRefs_; // this holds all the stubRefs associated with singl mu object

      float maxDeltaPt_;
      float maxDeltaEta_;
      float maxDeltaPhi_;
      float maxDeltaR_;
      float maxDeltaZ_;
      float maxDeltaM_;

      float deltaPt_;  // store max delta value
      float deltaEta_;  // store max delta value
      float deltaPhi_;  // store max delta value
      float deltaPhiAtVtx_;  // store max delta value
      float deltaR_;  // store max delta value
      float deltaZ_;  // store max delta value
      float deltaM_;  // store max delta value

      float deltaPtMin_;  // store min delta value
      float deltaEtaMin_;  // store min delta value
      float deltaPhiMin_;  // store min delta value
      float deltaPhiAtVtxMin_;  // store min delta value
      float deltaRMin_;  // store min delta value
      float deltaZMin_;  // store min delta value
      float deltaMMin_;  // store min delta value

      float mass_;
      float mass12_;

      int totalCharge_;
  

      EMTFHitRefVector stubs_;
      L1TTrackRefVector l1tracks_;

      template <typename T>
      void reorder( std::vector<T> & data, std::vector<std::size_t> const & order )
      {
         std::vector<T> tmp;         // create an empty vector
         tmp.reserve( data.size() ); // ensure memory and avoid moves in the vector
         for ( std::size_t i = 0; i < order.size(); ++i ) {
            tmp.push_back( data[order[i]] );
         }
         data.swap( tmp );          // swap vector contents
      }

      template <typename T>
      vector<size_t> sort_indexes(const vector<T> &v) {

        // initialize original index locations
        vector<size_t> idx(v.size());
        iota(idx.begin(), idx.end(), 0);

        // sort indexes based on comparing values in v
        sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

        return idx;

      }

  };

}

#endif
