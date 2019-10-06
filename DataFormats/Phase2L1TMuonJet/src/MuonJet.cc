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

MuonJet::MuonJet(const vector<float> vpt, const vector<float> veta, const  vector<float> vphi): pt_(vpt), eta_(veta), phi_(vphi)
{
}
