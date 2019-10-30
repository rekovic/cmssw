#include "DataFormats/Phase2L1TMuonJet/interface/MuonJet.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include "TLorentzVector.h"

using namespace std;
using namespace l1t;

MuonJet::MuonJet() {

  pt_.reserve(3);
  eta_.reserve(3);
  phi_.reserve(3);
  zvtx_.reserve(3);
  charge_.reserve(3);
  muonType_.reserve(3);
  stubQual_.reserve(3);
  stubEta_.reserve(3);
  stubPhi_.reserve(3);
  stubBend_.reserve(3);
  stubRef_.reserve(3);

}

MuonJet::MuonJet(const vector<float> vpt, const vector<float> veta, const  vector<float> vphi, const  vector<float> vz): pt_(vpt), eta_(veta), phi_(vphi), zvtx_(vz)
{
}

MuonJet::MuonJet(const L1TkMuonParticle & tkMuStub_1, const L1TkMuonParticle & tkMuStub_2, const L1TkMuonParticle & tkMuStub_3) 
{

  type_ = THREE_TKMUSTUB;

  muonType_  .push_back(k_TKMUSTUB);
  muonType_  .push_back(k_TKMUSTUB);
  muonType_  .push_back(k_TKMUSTUB);

  pt_   .push_back(tkMuStub_1.pt());
  pt_   .push_back(tkMuStub_2.pt());
  pt_   .push_back(tkMuStub_3.pt());

  eta_  .push_back(tkMuStub_1.eta());
  eta_  .push_back(tkMuStub_2.eta());
  eta_  .push_back(tkMuStub_3.eta());

  phi_  .push_back(tkMuStub_1.phi());
  phi_  .push_back(tkMuStub_2.phi());
  phi_  .push_back(tkMuStub_3.phi());

  stubEta_  .push_back(tkMuStub_1.getMuStubRef()->Eta_sim());
  stubEta_  .push_back(tkMuStub_2.getMuStubRef()->Eta_sim());
  stubEta_  .push_back(tkMuStub_3.getMuStubRef()->Eta_sim());

  stubPhi_  .push_back(tkMuStub_1.getMuStubRef()->Phi_sim() * TMath::Pi()/180.);
  stubPhi_  .push_back(tkMuStub_2.getMuStubRef()->Phi_sim() * TMath::Pi()/180.);
  stubPhi_  .push_back(tkMuStub_3.getMuStubRef()->Phi_sim() * TMath::Pi()/180.);

  stubBend_  .push_back(tkMuStub_1.getMuStubRef()->Bend());
  stubBend_  .push_back(tkMuStub_2.getMuStubRef()->Bend());
  stubBend_  .push_back(tkMuStub_3.getMuStubRef()->Bend());

  stubRef_  .push_back(tkMuStub_1.getMuStubRef());
  stubRef_  .push_back(tkMuStub_2.getMuStubRef());
  stubRef_  .push_back(tkMuStub_3.getMuStubRef());

  stubRefs_  .push_back(tkMuStub_1.getMuStubRefs());
  stubRefs_  .push_back(tkMuStub_2.getMuStubRefs());
  stubRefs_  .push_back(tkMuStub_3.getMuStubRefs());

  zvtx_ .push_back(tkMuStub_1.getTrkzVtx());
  zvtx_ .push_back(tkMuStub_2.getTrkzVtx());
  zvtx_ .push_back(tkMuStub_3.getTrkzVtx());

  charge_ .push_back(tkMuStub_1.charge());
  charge_ .push_back(tkMuStub_2.charge());
  charge_ .push_back(tkMuStub_3.charge());

  stubQual_ .push_back(tkMuStub_1.getMuStubRef()->Quality());
  stubQual_ .push_back(tkMuStub_2.getMuStubRef()->Quality());
  stubQual_ .push_back(tkMuStub_3.getMuStubRef()->Quality());

}

MuonJet::MuonJet(const L1TkMuonParticle & tkMuStub_1, const L1TkMuonParticle & tkMuStub_2) 
{
  
  type_ = TWO_TKMUSTUB;
  
  muonType_  .push_back(k_TKMUSTUB);
  muonType_  .push_back(k_TKMUSTUB);
  muonType_  .push_back(k_NONE);

  pt_   .push_back(tkMuStub_1.pt());
  pt_   .push_back(tkMuStub_2.pt());
  pt_   .push_back(-99);

  eta_  .push_back(tkMuStub_1.eta());
  eta_  .push_back(tkMuStub_2.eta());
  eta_  .push_back(-99);

  phi_  .push_back(tkMuStub_1.phi());
  phi_  .push_back(tkMuStub_2.phi());
  phi_  .push_back(-99);

  stubEta_  .push_back(tkMuStub_1.getMuStubRef()->Eta_sim());
  stubEta_  .push_back(tkMuStub_2.getMuStubRef()->Eta_sim());
  stubEta_  .push_back(-99);

  stubPhi_  .push_back(tkMuStub_1.getMuStubRef()->Phi_sim() * TMath::Pi()/180.);
  stubPhi_  .push_back(tkMuStub_2.getMuStubRef()->Phi_sim() * TMath::Pi()/180.);
  stubPhi_  .push_back(-99);

  stubBend_  .push_back(tkMuStub_1.getMuStubRef()->Bend());
  stubBend_  .push_back(tkMuStub_2.getMuStubRef()->Bend());
  stubBend_  .push_back(-99);

  stubRef_  .push_back(tkMuStub_1.getMuStubRef());
  stubRef_  .push_back(tkMuStub_2.getMuStubRef());
  stubRef_  .push_back(tkMuStub_1.getMuStubRef());  // FIXME: temporary push

  stubRefs_  .push_back(tkMuStub_1.getMuStubRefs());
  stubRefs_  .push_back(tkMuStub_2.getMuStubRefs());
  stubRefs_  .push_back(tkMuStub_1.getMuStubRefs());  // FIXME: temporary push

  zvtx_ .push_back(tkMuStub_1.getTrkzVtx());
  zvtx_ .push_back(tkMuStub_2.getTrkzVtx());
  zvtx_ .push_back(-99);

  charge_ .push_back(tkMuStub_1.charge());
  charge_ .push_back(tkMuStub_2.charge());
  charge_ .push_back(-99); 

  stubQual_ .push_back(tkMuStub_1.getMuStubRef()->Quality());
  stubQual_ .push_back(tkMuStub_2.getMuStubRef()->Quality());
  stubQual_  .push_back(-99);

}

MuonJet::MuonJet(const L1TkMuonParticle & tkMuStub_1, const L1TkMuonParticle & tkMuStub_2, const EMTFHit & muStub_1) 
{
  
  type_ = TWO_TKMUSTUB_ONE_MUSTUB;
  
  muonType_  .push_back(k_TKMUSTUB);
  muonType_  .push_back(k_TKMUSTUB);
  muonType_  .push_back(k_MUSTUB);

  pt_   .push_back(tkMuStub_1.pt());
  pt_   .push_back(tkMuStub_2.pt());
  pt_   .push_back(getPtFromBendingAngle(muStub_1));

  eta_  .push_back(tkMuStub_1.eta());
  eta_  .push_back(tkMuStub_2.eta());
  eta_  .push_back(muStub_1.Eta_sim());

  phi_  .push_back(tkMuStub_1.phi());
  phi_  .push_back(tkMuStub_2.phi());
  //phi_  .push_back(muStub_1.Phi_sim() * TMath::Pi()/180.);
  phi_  .push_back(getPhiAtVertex(muStub_1));

  stubEta_  .push_back(tkMuStub_1.getMuStubRef()->Eta_sim());
  stubEta_  .push_back(tkMuStub_2.getMuStubRef()->Eta_sim());
  stubEta_  .push_back(muStub_1.Eta_sim());

  stubPhi_  .push_back(tkMuStub_1.getMuStubRef()->Phi_sim() * TMath::Pi()/180.);
  stubPhi_  .push_back(tkMuStub_2.getMuStubRef()->Phi_sim() * TMath::Pi()/180.);
  stubPhi_  .push_back(muStub_1.Phi_sim() * TMath::Pi()/180.);

  stubBend_  .push_back(tkMuStub_1.getMuStubRef()->Bend());
  stubBend_  .push_back(tkMuStub_2.getMuStubRef()->Bend());
  stubBend_  .push_back(muStub_1.Bend());

  stubRef_  .push_back(tkMuStub_1.getMuStubRef());
  stubRef_  .push_back(tkMuStub_2.getMuStubRef());
  stubRef_  .push_back(tkMuStub_1.getMuStubRef());  // FIXME: temporary push
  //stubRef_  .push_back(muStub_1);

  zvtx_ .push_back(tkMuStub_1.getTrkzVtx());
  zvtx_ .push_back(tkMuStub_2.getTrkzVtx());
  zvtx_ .push_back(-99);

  charge_ .push_back(tkMuStub_1.charge());
  charge_ .push_back(tkMuStub_2.charge());
  charge_ .push_back(getChargeFromBendingAngle(muStub_1)); 

  stubQual_ .push_back(tkMuStub_1.getMuStubRef()->Quality());
  stubQual_ .push_back(tkMuStub_2.getMuStubRef()->Quality());
  stubQual_ .push_back(muStub_1.Quality());

}


MuonJet::MuonJet(const L1TkMuonParticle & tkMuStub_1, const L1TkMuonParticle & tkMuStub_2, const EMTFHitRef muStub_1) 
{
  
  type_ = TWO_TKMUSTUB_ONE_MUSTUB;
  
  muonType_  .push_back(k_TKMUSTUB);
  muonType_  .push_back(k_TKMUSTUB);
  muonType_  .push_back(k_MUSTUB);

  pt_   .push_back(tkMuStub_1.pt());
  pt_   .push_back(tkMuStub_2.pt());
  pt_   .push_back(getPtFromBendingAngle(*muStub_1));

  eta_  .push_back(tkMuStub_1.eta());
  eta_  .push_back(tkMuStub_2.eta());
  eta_  .push_back(muStub_1->Eta_sim());

  phi_  .push_back(tkMuStub_1.phi());
  phi_  .push_back(tkMuStub_2.phi());
  //phi_  .push_back(muStub_1->Phi_sim() * TMath::Pi()/180.);
  phi_  .push_back(getPhiAtVertex(*muStub_1));

  stubEta_  .push_back(tkMuStub_1.getMuStubRef()->Eta_sim());
  stubEta_  .push_back(tkMuStub_2.getMuStubRef()->Eta_sim());
  stubEta_  .push_back(muStub_1->Eta_sim());

  stubPhi_  .push_back(tkMuStub_1.getMuStubRef()->Phi_sim() * TMath::Pi()/180.);
  stubPhi_  .push_back(tkMuStub_2.getMuStubRef()->Phi_sim() * TMath::Pi()/180.);
  stubPhi_  .push_back(muStub_1->Phi_sim() * TMath::Pi()/180.);

  stubBend_  .push_back(tkMuStub_1.getMuStubRef()->Bend());
  stubBend_  .push_back(tkMuStub_2.getMuStubRef()->Bend());
  stubBend_  .push_back(muStub_1->Bend());

  stubRef_  .push_back(tkMuStub_1.getMuStubRef());
  stubRef_  .push_back(tkMuStub_2.getMuStubRef());
  stubRef_  .push_back(muStub_1);

  stubRefs_  .push_back(tkMuStub_1.getMuStubRefs());
  stubRefs_  .push_back(tkMuStub_2.getMuStubRefs());
  EMTFHitRefVector singleRefVector(1,muStub_1);
  stubRefs_  .push_back(singleRefVector);

  zvtx_ .push_back(tkMuStub_1.getTrkzVtx());
  zvtx_ .push_back(tkMuStub_2.getTrkzVtx());
  zvtx_ .push_back(-99);

  charge_ .push_back(tkMuStub_1.charge());
  charge_ .push_back(tkMuStub_2.charge());
  charge_ .push_back(getChargeFromBendingAngle(*muStub_1)); 

  stubQual_ .push_back(tkMuStub_1.getMuStubRef()->Quality());
  stubQual_ .push_back(tkMuStub_2.getMuStubRef()->Quality());
  stubQual_ .push_back(muStub_1->Quality());

}

MuonJet::MuonJet(const L1TkMuonParticle & tkMuStub_1, const EMTFHit & muStub_1, const EMTFHit & muStub_2) 
{
  
  type_ = ONE_TKMUSTUB_TWO_MUSTUB;
  
  muonType_  .push_back(k_TKMUSTUB);
  muonType_  .push_back(k_MUSTUB);
  muonType_  .push_back(k_MUSTUB);

  pt_   .push_back(tkMuStub_1.pt());
  pt_   .push_back(getPtFromBendingAngle(muStub_1));
  pt_   .push_back(getPtFromBendingAngle(muStub_2));

  eta_  .push_back(tkMuStub_1.eta());
  eta_  .push_back(muStub_1.Eta_sim());
  eta_  .push_back(muStub_2.Eta_sim());

  phi_  .push_back(tkMuStub_1.phi());
  phi_  .push_back(muStub_1.Phi_sim() * TMath::Pi()/180.);
  phi_  .push_back(muStub_2.Phi_sim() * TMath::Pi()/180.);

  stubEta_  .push_back(tkMuStub_1.getMuStubRef()->Eta_sim());
  stubEta_  .push_back(muStub_1.Eta_sim());
  stubEta_  .push_back(muStub_2.Eta_sim());

  stubPhi_  .push_back(tkMuStub_1.getMuStubRef()->Phi_sim() * TMath::Pi()/180.);
  stubPhi_  .push_back(muStub_1.Phi_sim() * TMath::Pi()/180.);
  stubPhi_  .push_back(muStub_2.Phi_sim() * TMath::Pi()/180.);

  stubBend_  .push_back(tkMuStub_1.getMuStubRef()->Bend());
  stubBend_  .push_back(muStub_1.Bend());
  stubBend_  .push_back(muStub_2.Bend());

  stubRef_  .push_back(tkMuStub_1.getMuStubRef());
  stubRef_  .push_back(tkMuStub_1.getMuStubRef());  // FIXME: temporary push
  stubRef_  .push_back(tkMuStub_1.getMuStubRef());  // FIXME: temporary push
  //stubRef_  .push_back(muStub_1);
  //stubRef_  .push_back(muStub_2);

  zvtx_ .push_back(tkMuStub_1.getTrkzVtx());
  zvtx_ .push_back(-99);
  zvtx_ .push_back(-99);

  charge_ .push_back(tkMuStub_1.charge());
  charge_ .push_back(getChargeFromBendingAngle(muStub_1)); 
  charge_ .push_back(getChargeFromBendingAngle(muStub_2)); 

  stubQual_  .push_back(tkMuStub_1.getMuStubRef()->Quality());
  stubQual_  .push_back(muStub_1.Quality());
  stubQual_  .push_back(muStub_2.Quality());

}

MuonJet::MuonJet(const L1TkMuonParticle & tkMuStub_1, const EMTFHitRef muStub_1, const EMTFHitRef muStub_2) 
{
  
  type_ = ONE_TKMUSTUB_TWO_MUSTUB;
  
  muonType_  .push_back(k_TKMUSTUB);
  muonType_  .push_back(k_MUSTUB);
  muonType_  .push_back(k_MUSTUB);

  pt_   .push_back(tkMuStub_1.pt());
  pt_   .push_back(getPtFromBendingAngle(*muStub_1));
  pt_   .push_back(getPtFromBendingAngle(*muStub_2));

  eta_  .push_back(tkMuStub_1.eta());
  eta_  .push_back(muStub_1->Eta_sim());
  eta_  .push_back(muStub_2->Eta_sim());

  phi_  .push_back(tkMuStub_1.phi());
  phi_  .push_back(muStub_1->Phi_sim() * TMath::Pi()/180.);
  phi_  .push_back(muStub_2->Phi_sim() * TMath::Pi()/180.);

  stubEta_  .push_back(tkMuStub_1.getMuStubRef()->Eta_sim());
  stubEta_  .push_back(muStub_1->Eta_sim());
  stubEta_  .push_back(muStub_2->Eta_sim());

  stubPhi_  .push_back(tkMuStub_1.getMuStubRef()->Phi_sim() * TMath::Pi()/180.);
  stubPhi_  .push_back(muStub_1->Phi_sim() * TMath::Pi()/180.);
  stubPhi_  .push_back(muStub_2->Phi_sim() * TMath::Pi()/180.);

  stubBend_  .push_back(tkMuStub_1.getMuStubRef()->Bend());
  stubBend_  .push_back(muStub_1->Bend());
  stubBend_  .push_back(muStub_2->Bend());

  stubRef_  .push_back(tkMuStub_1.getMuStubRef());
  stubRef_  .push_back(muStub_1);
  stubRef_  .push_back(muStub_2);

  stubRefs_  .push_back(tkMuStub_1.getMuStubRefs());
  EMTFHitRefVector singleRefVector_1(1,muStub_1);
  stubRefs_  .push_back(singleRefVector_1);
  EMTFHitRefVector singleRefVector_2(1,muStub_2);
  stubRefs_  .push_back(singleRefVector_2);

  zvtx_ .push_back(tkMuStub_1.getTrkzVtx());
  zvtx_ .push_back(-99);
  zvtx_ .push_back(-99);

  charge_ .push_back(tkMuStub_1.charge());
  charge_ .push_back(getChargeFromBendingAngle(*muStub_1)); 
  charge_ .push_back(getChargeFromBendingAngle(*muStub_2)); 

  stubQual_  .push_back(tkMuStub_1.getMuStubRef()->Quality());
  stubQual_  .push_back(muStub_1->Quality());
  stubQual_  .push_back(muStub_2->Quality());

}

MuonJet::MuonJet(const EMTFHit & muStub_1, const EMTFHit & muStub_2, const EMTFHit & muStub_3) 
{
  
  type_ = THREE_MUSTUB;
  
  muonType_  .push_back(k_MUSTUB);
  muonType_  .push_back(k_MUSTUB);
  muonType_  .push_back(k_MUSTUB);

  pt_   .push_back(getPtFromBendingAngle(muStub_1));
  pt_   .push_back(getPtFromBendingAngle(muStub_2));
  pt_   .push_back(getPtFromBendingAngle(muStub_3));

  eta_  .push_back(muStub_1.Eta_sim());
  eta_  .push_back(muStub_2.Eta_sim());
  eta_  .push_back(muStub_3.Eta_sim());

  phi_  .push_back(muStub_1.Phi_sim() * TMath::Pi()/180.);
  phi_  .push_back(muStub_2.Phi_sim() * TMath::Pi()/180.);
  phi_  .push_back(muStub_3.Phi_sim() * TMath::Pi()/180.);

  stubEta_  .push_back(muStub_1.Eta_sim());
  stubEta_  .push_back(muStub_2.Eta_sim());
  stubEta_  .push_back(muStub_3.Eta_sim());

  stubPhi_  .push_back(muStub_1.Phi_sim() * TMath::Pi()/180.);
  stubPhi_  .push_back(muStub_2.Phi_sim() * TMath::Pi()/180.);
  stubPhi_  .push_back(muStub_3.Phi_sim() * TMath::Pi()/180.);

  stubBend_  .push_back(muStub_1.Bend());
  stubBend_  .push_back(muStub_2.Bend());
  stubBend_  .push_back(muStub_3.Bend());

  //stubRef_  .push_back(muStub_1);
  //stubRef_  .push_back(muStub_2);
  //stubRef_  .push_back(muStub_3);

  zvtx_ .push_back(-99);
  zvtx_ .push_back(-99);
  zvtx_ .push_back(-99);

  charge_ .push_back(getChargeFromBendingAngle(muStub_1)); 
  charge_ .push_back(getChargeFromBendingAngle(muStub_2)); 
  charge_ .push_back(getChargeFromBendingAngle(muStub_3)); 

  stubQual_  .push_back(muStub_1.Quality());
  stubQual_  .push_back(muStub_2.Quality());
  stubQual_  .push_back(muStub_3.Quality());

}

void MuonJet::configure(float const& maxDeltaEta, float const& maxDeltaPhi, float const& maxDeltaR, float const& maxDeltaZ) 
{

  maxDeltaEta_ = maxDeltaEta;
  maxDeltaPhi_ = maxDeltaPhi;
  maxDeltaR_ = maxDeltaR;
  maxDeltaZ_ = maxDeltaZ;

}

void  MuonJet::process() {

  // first sort according to pt
  //
  // /////////////////////////
  sortWithPt();

  // calculate delta eta,phi,R
  //
  // /////////////////////////
  const float MU_MASS = 0.105658;
  math::PtEtaPhiMLorentzVectorF LV[3];

  for (int i=0; i<3; i++) 
  {
    LV[i].SetCoordinates(pt_[i], eta_[i], phi_[i], MU_MASS);
  }

  std::vector<float> deltaR_v;
  deltaR_v.push_back(ROOT::Math::VectorUtil::DeltaR(LV[0],LV[1]));
  deltaR_v.push_back(ROOT::Math::VectorUtil::DeltaR(LV[0],LV[2]));
  deltaR_v.push_back(ROOT::Math::VectorUtil::DeltaR(LV[1],LV[2]));
  
  math::PtEtaPhiMLorentzVectorF LV3 = LV[0] + LV[1] + LV[2];

  math::PtEtaPhiMLorentzVectorF LV01  = LV[0] + LV[1];

 // in case of TWO_TKMUSTUB_ONE_MUSTUB, 
 // use the two leggs of muonType_ k_TKMUSTUB
 // to be used for mass12 calculation
  if(type_ == TWO_TKMUSTUB_ONE_MUSTUB) {

    if(muonType_[0] == k_MUSTUB) LV01  = LV[1] + LV[2];
    if(muonType_[1] == k_MUSTUB) LV01  = LV[0] + LV[2];
    if(muonType_[2] == k_MUSTUB) LV01  = LV[0] + LV[1];
 
 }

  deltaR_ = TMath::MaxElement(3,deltaR_v.data());

  // in case of TWO_TKMUSTUB deltaR is the one betwen those two objects
  if(type_ == TWO_TKMUSTUB) deltaR_ = TMath::MaxElement(1,deltaR_v.data());

  // calculate mass
  //
  // /////////////////////////
  mass_ = LV3.M();

  mass12_ = LV01.M();

  // calculate total charge
  // 
  // /////////////////////////
  totalCharge_ = charge_[0] + charge_[1] + charge_[2];

  // calculate delta Z
  //
  // /////////////////////////
  float deltaZ = -99;

  for (int i=0; i<3; i++) 
  {
    for (int j=i+1; j<3; j++) 
    { 

      // calculate delta z
      // only calculate if valid z
      if(zvtx_[i] == -99 || zvtx_[j] == -99) continue;

      float temp_deltaZ = abs(zvtx_[i]-zvtx_[j]);
      //cout << "    deltaZ = " << deltaZ << " temp_deltaZ = " << temp_deltaZ << endl;

      if (temp_deltaZ > deltaZ) 
        deltaZ = temp_deltaZ;
    }
  }

  deltaZ_ = deltaZ;

  // in case of ONE_TKMUSTUB_TWO_MUSTUB or THREE_MUSTUB deltaZ is set to -99;s
  if(type_ == ONE_TKMUSTUB_TWO_MUSTUB || type_ == THREE_MUSTUB) deltaZ_ = -99;


  return;
}

bool MuonJet::isValid()
{

  bool rc = true;

  if (deltaZ_ > maxDeltaZ_) return false;
  if (deltaR_ > maxDeltaR_) return false;
  if (abs(totalCharge_) > 1) return false;

  //  TkMuStub muste be 
  //  type 1 if eta (1.2 - 2.0)
  //  type 4 if eta (2.0 - 2.4)
  if(! areValidStubs() ) return false;

  // The principle stubs of TkMuStub object are unique, by construction of TkMuStub object
  // This is a sanity check;
  if(type_ == THREE_TKMUSTUB && ! areUniqueTkMuStubs() ) return false;

  // Check that stub of MUSTUB type leg are not among stubs associated with another leg
  if(type_ != THREE_TKMUSTUB && ! areUniqueStubs() ) return false;


  return rc;
}

bool MuonJet::areValidStubs()
{

  bool rc = false;
  bool isValidMuon_0 = false;
  bool isValidMuon_1 = false;
  bool isValidMuon_2 = false;

  cout << "areValidStubs():  muonType_[0] = " << muonType_[0] << " muonType_[1] = " << muonType_[1] << " muonType_[2] = " << muonType_[2] << endl;

  if(muonType_[0] == k_TKMUSTUB) {
    if(isValidTkMuStub(stubRef_[0])) isValidMuon_0 = true;
  }
  else if(muonType_[0] == k_MUSTUB) {
    if(isValidMuStub(stubRef_[0])) isValidMuon_0 = true;
  }
  else {
    isValidMuon_0 = false;
  }

  if(muonType_[1] == k_TKMUSTUB) {
    if(isValidTkMuStub(stubRef_[1])) isValidMuon_1 = true;
  }
  else if(muonType_[1] == k_MUSTUB) {
    if(isValidMuStub(stubRef_[1])) isValidMuon_1 = true;
  }
  else {
    isValidMuon_1 = false;
  }

  if(muonType_[2] == k_TKMUSTUB) {
    if(isValidTkMuStub(stubRef_[2])) isValidMuon_2 = true;
  }
  else if(muonType_[2] == k_MUSTUB) {
    if(isValidMuStub(stubRef_[2])) isValidMuon_2 = true;
  }
  else {
    isValidMuon_2 = true;
  }

  cout << "   Validation of muons (isValidMuon_0, isValidMuon_1, isValidMuon_2) = (" << isValidMuon_0 << "," << isValidMuon_1 << "," << isValidMuon_2 << ")" << endl;
  if(isValidMuon_0 && isValidMuon_1 && isValidMuon_2) rc = true;

  return rc;

}

bool MuonJet::isValidMuStub(const EMTFHitRef pStub)
{

  bool rc = true;
  float astubEta  = abs(pStub->Eta_sim());
  int stubType    = pStub->Subsystem();
  int stubStation = pStub->Station();
  int stubRing    = pStub->Ring();

  cout << "Validation of stub: (eta, phi, type, station, ring = ( " << pStub->Eta_sim() << "," << pStub->Phi_sim() * TMath::Pi()/180.  << "," << stubType << "," << stubStation << "," << stubRing << ")" << endl;

  // in eta (1.2 - 1.6) must be hit from ME1/2
  if(astubEta > 1.2 && astubEta < 1.6 && (stubStation != 1 || stubRing != 2 || stubType != EMTFHit::kCSC) ) rc = false;
  // in eta (1.6 - 2.0) must be hit from ME1/1
  if(astubEta > 1.6 && astubEta < 2.0 && (stubStation != 1 || stubRing != 1 || stubType != EMTFHit::kCSC) ) rc = false;
  // in eta (2.0 - 2.8) must be hit from ME0
  if(astubEta > 2.0 && astubEta < 2.8 && (stubStation != 1 || stubRing != 1 || stubType != EMTFHit::kME0) ) rc = false;

  return rc;

}

bool MuonJet::isValidTkMuStub(const EMTFHitRef pStub)
{

  bool rc = true;

  //float astubEta  = abs(pStub->Eta_sim());
  //int stubType    = pStub->Subsystem();

  //if(astubEta > 1.2 && astubEta < 2.0 && stubType != EMTFHit::kCSC && stubType != EMTFHit::kRPC) rc = false;
  //if(astubEta > 2.0 && astubEta < 2.8 && stubType == EMTFHit::kME0) rc = true;
  

  return rc;

}

bool const MuonJet::areUniqueTkMuStubs ()
{

  bool rc = true;

  if (abs(stubEta_[0] - stubEta_[1]) == 0) rc = false;
  if (abs(stubEta_[0] - stubEta_[2]) == 0) rc = false;
  if (abs(stubEta_[1] - stubEta_[2]) == 0) rc = false;

  if (abs(stubPhi_[0] - stubPhi_[1]) == 0) rc = false;
  if (abs(stubPhi_[0] - stubPhi_[2]) == 0) rc = false;
  if (abs(stubPhi_[1] - stubPhi_[2]) == 0) rc = false;

  if(rc == false) edm::LogWarning("MuonJet::areUniqueTkMuStubs") << "WARNING: Not passing sanity check: principle stubs of TkMuStub object within MuonJet are not unique. " << endl;
  return rc;

}


bool const MuonJet::areUniqueStubs ()
{

  bool rc = true;

  // remove duplicates from overlapped 
  if (abs(stubEta_[0] - stubEta_[1]) < 0.01)
    if (abs(stubPhi_[0] - stubPhi_[1]) < 0.02) rc = false;

  if (abs(stubEta_[0] - stubEta_[2]) < 0.01)
    if (abs(stubPhi_[0] - stubPhi_[2]) < 0.02) rc = false;

  if (abs(stubEta_[1] - stubEta_[2]) < 0.01)
    if (abs(stubPhi_[1] - stubPhi_[2]) < 0.02) rc = false;


  // FIXME:  Implement
  // for each MUSTUB type, check that MuStubRef is not common with any of MuStubRefs of other legs
  // The above checks with  difference in eta and phi might be OK for now.
   
  for (int i =0; i<3; i++) {
    
    if(muonType_[i] != k_MUSTUB) continue;

    const EMTFHitRef i_muonMainStub = stubRef_[i];

    for (int j =i+1; j<3; j++) {

      for (auto it = stubRefs_[j].begin(); it != stubRefs_[j].end(); ++it) {

        const EMTFHitRef j_muonAssociatedStub = *it;

        if(i_muonMainStub == j_muonAssociatedStub) cout << "NOT UNIQUE STUBS !!!!!! " << endl;

        float delta_eta = abs(i_muonMainStub->Eta_sim() - j_muonAssociatedStub->Eta_sim());
        float delta_phi = abs(i_muonMainStub->Phi_sim()* TMath::Pi()/180. - j_muonAssociatedStub->Phi_sim()* TMath::Pi()/180.);
        int delta_chamber = -99;
        if(i_muonMainStub->Subsystem() == j_muonAssociatedStub->Subsystem() &&  i_muonMainStub->Station() == j_muonAssociatedStub->Station()) 
          delta_chamber = abs( i_muonMainStub->Chamber() - i_muonMainStub->Chamber() );

        if (delta_eta < 0.02 && delta_phi < 0.02) {
          if (delta_chamber == 0) cout << "NOT UNIQUE STUBS DUE TO DELTA RAYS !!!!!! " << endl;
          if (delta_chamber == 1) cout << "NOT UNIQUE STUBS DUE TO OVERLAPS !!!!!! " << endl;
        }

        /*
        if ( delta_eta < 0.02 && delta_phi < 0.02) {
          
          cout << "NOT UNIQUE STUBS: TWO SUTB TOO CLOSE!!!!!! " << endl;
          cout << " subsystem = " << i_muonMainStub->Subsystem() <<  "station = " << i_muonMainStub->Station() << "chamber = " << i_muonMainStub->Chamber() << " neighbor = " << i_muonMainStub->Neighbor() << endl;
          cout << " subsystem = " << j_muonAssociatedStub->Subsystem() <<  "station = " << j_muonAssociatedStub->Station() << "chamber = " << j_muonAssociatedStub->Chamber() << " neighbor = " << j_muonAssociatedStub->Neighbor() << endl;

        }
        */

      } // end for it
     
   } // end for j

  } // end for i

  cout << "areUniqueStubs(): rc == " << rc << endl;
  return rc;

}


void MuonJet::print() 
{

  printf("%10s %5s = %1d \n","MuonJet:","type",type_);
  printf("%15s %8s %8s %8s %10s %10s %10s %10s %8s %8s %8s %8s %8s %8s %8s \n","muon","pt","eta","phi","stubEta","stubPhi", "stubRho","stubBend","z","charge","muType","station","subsyst","chamber","stubQual");

  for (int i=0; i<3; i++) 
  printf("%15d %8.2f %8.2f %8.2f %10.4f %10.4f %10.1f %10f %8.2f %8d %8d %8d %8d %8d %8d \n", i, pt_[i], eta_[i], phi_[i], stubEta_[i], stubPhi_[i], stubRef_[i]->Rho_sim(), stubBend_[i], zvtx_[i], charge_[i], muonType_[i], stubRef_[i]->Station(), stubRef_[i]->Subsystem(), stubRef_[i]->Chamber(), stubQual_[i]);

  printf("%15s = %2.2f , %5s = %2.2f\n", "maxDeltaR", maxDeltaR_, "deltaR" , deltaR_);
  printf("%15s = %2.2f , %5s = %2.2f\n", "maxDeltaZ", maxDeltaZ_, "deltaZ" , deltaZ_);
  printf("%15s = %2.2f \n", "mass", mass_);
  printf("%15s = %2.2f \n", "mass12", mass12_);
  printf("%15s = %3d \n", "totalCharge", totalCharge_);

  return;

}

void MuonJet::sortWithPt() 
{
  
  vector<std::size_t> indicesPtDescend;
  for (auto i: sort_indexes(pt_)) indicesPtDescend.push_back(i);

  reorder(pt_,indicesPtDescend);
  reorder(eta_,indicesPtDescend);
  reorder(phi_,indicesPtDescend);
  reorder(stubEta_,indicesPtDescend);
  reorder(stubPhi_,indicesPtDescend);
  reorder(stubBend_,indicesPtDescend);
  reorder(zvtx_,indicesPtDescend);
  reorder(charge_,indicesPtDescend);
  reorder(stubQual_,indicesPtDescend);
  reorder(muonType_,indicesPtDescend);
  reorder(stubRef_,indicesPtDescend);
  reorder(stubRefs_,indicesPtDescend);

  return;

}

float MuonJet::getPhiAtVertex(const EMTFHit & muStub) {

  float pt = getPtFromBendingAngle(muStub);
  float rho = muStub.Rho_sim() / 100.0; // in meters
  float eta = muStub.Eta_sim();
  float phi = muStub.Phi_sim() * TMath::Pi()/180.; // in radians
  int charge = getChargeFromBendingAngle(muStub);

  int signEtaCharge = eta * charge /abs(eta * charge);

  float delta_phi = 0.15 * 3.8 * rho * signEtaCharge / pt;

  float phi_vtx = phi + delta_phi;

  return phi_vtx;


}

float MuonJet::getPtFromBendingAngle(const EMTFHit & muStub) {

  float stubBend = muStub.Bend();
  float stubPt = -99;

  // convert bend angle phi to rad in ME0 
  // ////////////////////////////////////
  // ME bend is in units of 1/4 strip, 
  // there are 384 strips in 20 deg chamber.
  // so unit of bending is quivalent to 0.00363 rad = (20/384 *4 *3.14/180)
  if (muStub.Subsystem() == 4) { // ME0

    stubBend *= 3.63*2.0/1000.; // in rad 
    stubPt = 1.0 / stubBend; // in GeV

  } // end if ME0

  // convert bend angle to mrad in CSC 
  // Encoded bend in 6 bits, which corresponds to 1/32-strip unit
  // strip unit = 0.00296 rad for ME1/1 (from Muon TDR) 
  // ////////////////////////////////
  if (muStub.Subsystem() == 1) { // CSC
      
    float stripSize = -1; // in mrad
    
    if(muStub.Station()== 1) {
      if( muStub.Ring() == 1) stripSize =  2.96; // ME1/1
      if( muStub.Ring() == 2) stripSize =  2.33; // ME1/2
      if( muStub.Ring() == 3) stripSize =  2.16; // ME1/3
    }

    if(muStub.Station()== 2) {
      if( muStub.Ring() == 1) stripSize =  4.65; // ME2/1
      if( muStub.Ring() == 2) stripSize =  2.33; // ME2/1
    }

    if(muStub.Station()== 3) {
      if( muStub.Ring() == 1) stripSize =  4.65; // ME3/1
      if( muStub.Ring() == 2) stripSize =  2.33; // ME3/1
    }
    
    if(muStub.Station()== 4) {
      if(muStub.Ring() == 1) stripSize =  4.65; // ME4/1
      if(muStub.Ring() == 2) stripSize =  2.33; // ME4/1
    }
      
    stubBend *= stripSize*32.0/4.0/1000.; // in rad
    stubPt = 1.0 / stubBend; // in GeV

  } // end if CSC

  return abs(stubPt);

}


/*
float MuonJet::getPtFromBendingAngle(const EMTFHitRef muStub) {

  float stubBend = muStub->Bend();
  float stubPt = -99;

  // convert bend angle phi to rad in ME0 
  // ////////////////////////////////////
  // ME bend is in units of 1/4 strip, 
  // there are 384 strips in 20 deg chamber.
  // so unit of bending is quivalent to 0.00363 rad = (20/384 *4 *3.14/180)
  if (muStub->Subsystem() == 4) { // ME0

    stubBend *= 3.63/1000.; // in rad 
    stubPt = 1.0 / stubBend; // in GeV

  } // end if ME0

  // convert bend angle to mrad in CSC 
  // Encoded bend in 6 bits, which corresponds to 1/32-strip unit
  // strip unit = 0.00296 rad for ME1/1 (from Muon TDR) 
  // ////////////////////////////////
  if (muStub->Subsystem() == 1) { // CSC
      
    float stripSize = -1; // in mrad
    
    if(muStub->Station()== 1) {
      if( muStub->Ring() == 1) stripSize =  2.96; // ME1/1
      if( muStub->Ring() == 2) stripSize =  2.33; // ME1/2
      if( muStub->Ring() == 3) stripSize =  2.16; // ME1/3
    }

    if(muStub->Station()== 2) {
      if( muStub->Ring() == 1) stripSize =  4.65; // ME2/1
      if( muStub->Ring() == 2) stripSize =  2.33; // ME2/1
    }

    if(muStub->Station()== 3) {
      if( muStub->Ring() == 1) stripSize =  4.65; // ME3/1
      if( muStub->Ring() == 2) stripSize =  2.33; // ME3/1
    }
    
    if(muStub->Station()== 4) {
      if(muStub->Ring() == 1) stripSize =  4.65; // ME4/1
      if(muStub->Ring() == 2) stripSize =  2.33; // ME4/1
    }
      
    stubBend *= stripSize*32.0/1000.; // in rad
    stubPt = 1.0 / stubBend; // in GeV

  } // end if CSC

  return abs(stubPt);

}
*/

int MuonJet::getChargeFromBendingAngle(const EMTFHit & muStub) {

  int charge = -99;
  
  float stubBend = muStub.Bend();

  if(stubBend > 0) charge = -1;
  if(stubBend < 0) charge = 1;

  return charge;

}

