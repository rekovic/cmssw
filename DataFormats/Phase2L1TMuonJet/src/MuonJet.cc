#include "DataFormats/Phase2L1TMuonJet/interface/MuonJet.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace l1t;

MuonJet::MuonJet() {

  pt_.reserve(3);
  eta_.reserve(3);
  phi_.reserve(3);

}

MuonJet::MuonJet(const vector<float> vpt, const vector<float> veta, const  vector<float> vphi, const  vector<float> vz): pt_(vpt), eta_(veta), phi_(vphi), zvtx_(vz)
{
}

MuonJet::MuonJet(const L1TkMuonParticle & tkMuStub_1, const L1TkMuonParticle & tkMuStub_2, const L1TkMuonParticle & tkMuStub_3) 
{

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

}
