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
  qual_.reserve(3);

}

MuonJet::MuonJet(const vector<float> vpt, const vector<float> veta, const  vector<float> vphi, const  vector<float> vz): pt_(vpt), eta_(veta), phi_(vphi), zvtx_(vz)
{
}

MuonJet::MuonJet(const L1TkMuonParticle & tkMuStub_1, const L1TkMuonParticle & tkMuStub_2, const L1TkMuonParticle & tkMuStub_3) 
{

  type_ = THREE_TKMUSTUB;

  pt_   .push_back(tkMuStub_1.pt());
  pt_   .push_back(tkMuStub_2.pt());
  pt_   .push_back(tkMuStub_3.pt());

  eta_  .push_back(tkMuStub_1.eta());
  eta_  .push_back(tkMuStub_2.eta());
  eta_  .push_back(tkMuStub_3.eta());

  phi_  .push_back(tkMuStub_1.phi());
  phi_  .push_back(tkMuStub_2.phi());
  phi_  .push_back(tkMuStub_3.phi());

  zvtx_ .push_back(tkMuStub_1.getTrkzVtx());
  zvtx_ .push_back(tkMuStub_2.getTrkzVtx());
  zvtx_ .push_back(tkMuStub_3.getTrkzVtx());

  qual_ .push_back(99);
  qual_ .push_back(99);
  qual_ .push_back(99);

}

MuonJet::MuonJet(const L1TkMuonParticle & tkMuStub_1, const L1TkMuonParticle & tkMuStub_2, const EMTFHit & muStub_1) 
{
  
  type_ = TWO_TKMUSTUB_ONE_MUSTUB;

  pt_   .push_back(tkMuStub_1.pt());
  pt_   .push_back(tkMuStub_2.pt());
  pt_   .push_back(-99);

  eta_  .push_back(tkMuStub_1.eta());
  eta_  .push_back(tkMuStub_2.eta());
  eta_  .push_back(muStub_1.Eta_sim());

  phi_  .push_back(tkMuStub_1.phi());
  phi_  .push_back(tkMuStub_2.phi());
  phi_  .push_back(muStub_1.Phi_sim() * TMath::Pi()/180.);

  zvtx_ .push_back(tkMuStub_1.getTrkzVtx());
  zvtx_ .push_back(tkMuStub_2.getTrkzVtx());
  zvtx_ .push_back(-99);

  qual_  .push_back(99);
  qual_  .push_back(99);
  qual_  .push_back(muStub_1.Quality());

}

void MuonJet::configure(float const& maxDeltaEta, float const& maxDeltaPhi, float const& maxDeltaR, float const& maxDeltaZ) 
{

  maxDeltaEta_ = maxDeltaEta;
  maxDeltaPhi_ = maxDeltaPhi;
  maxDeltaR_ = maxDeltaR;
  maxDeltaZ_ = maxDeltaZ;

}

void  MuonJet::process() {

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
  
  math::PtEtaPhiMLorentzVectorF LV01  = LV[0] + LV[1];
  math::PtEtaPhiMLorentzVectorF LV012 = LV01 + LV[2];

  deltaR_ = TMath::MaxElement(3,deltaR_v.data());

  mass_ = LV012.M();

  totalCharge_ = -99;

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

  // calculate delta eta,phi,R


  return;
}

bool MuonJet::isValid()
{

  bool rc = true;

  if (deltaZ_ > maxDeltaZ_) rc = false;
  if (deltaR_ > maxDeltaR_) rc = false;

  return rc;
}

void MuonJet::print() 
{

  printf("%10s %5s = %1d \n","MuonJet:","type",type_);
  printf("%15s %7s %7s %7s %7s %7s \n","muon","pt","eta","phi","z","qual");

  for (int i=0; i<3; i++) 
  printf("%15d %7.2f %7.2f %7.2f %7.2f %7d \n", i, pt_[i], eta_[i], phi_[i], zvtx_[i], qual_[i]);

  printf("%15s = %2.2f , %5s = %2.2f\n", "maxDeltaR", maxDeltaR_, "deltaR" , deltaR_);
  printf("%15s = %2.2f , %5s = %2.2f\n", "maxDeltaZ", maxDeltaZ_, "deltaZ" , deltaZ_);
  printf("%15s = %2.2f \n", "mass", mass_);
  printf("%15s = %3d \n", "totalCharge", totalCharge_);

  return;

}
