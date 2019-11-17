// input: L1TkTracks and Muon Stubs 
//
// 
//
// Author: Vladimir Rekovic
// Version: 1   Date: 10.Feb.2019
// Descritiption: 
//      Match L1Tracks to EndCap muon stubs.
//      Stubs are taken to be inputs to EMTF, EMTFHitCollection
//      Use DynamicWindow matching algorithm which assumes muon 
//      detector object (here Stubs) have coordinates at 2nd station.  
//      Match the two and produce a collection of L1TkMuonParticle

// user include files
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"
#include "L1Trigger/L1TTrackMatch/interface/L1TkMuCorrDynamicWindows.h"
#include "L1Trigger/L1TMuonEndCap/interface/Common.h"

// system include files
#include <memory>
#include <string>

using namespace l1t;

class L1TkMuonStubProducer : public edm::EDProducer {
public:

  typedef TTTrack< Ref_Phase2TrackerDigi_ >  L1TTTrackType;
  typedef std::vector< L1TTTrackType > L1TTTrackCollectionType;

  struct PropState { //something simple, imagine it's hardware emulation
    PropState() :
      pt(-99),  eta(-99), phi(-99),
      sigmaPt(-99),  sigmaEta(-99), sigmaPhi(-99),
      valid(false) {}
    float pt;
    float eta;
    float phi;
    float sigmaPt;
    float sigmaEta;
    float sigmaPhi;
    bool valid;

  };

  enum AlgoType {
    kTP = 1,
    kDynamicWindows = 2
  };

  explicit L1TkMuonStubProducer(const edm::ParameterSet&);
  ~L1TkMuonStubProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  int findStubRefIndex( const edm::Handle<EMTFHitCollection>&, const EMTFHit & ) const;

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  // algo for endcap regions using dynamic windows for making the match trk + muStub
  void runOnMuonHitCollection(const edm::Handle<EMTFHitCollection>&,
                          const edm::Handle<L1TTTrackCollectionType>&,
                          L1TkMuonParticleCollection& tkMuons) const;

  void getExtendedEtaME0Muons(const edm::Handle<EMTFHitCollection>& , L1TkMuonParticleCollection&, EMTFHitCollection & ) const;

  void cleanStubs(const EMTFHitCollection &, EMTFHitCollection &) const;
  void printStubs(const EMTFHitCollection & , int) const;
  void printStub(const EMTFHit & ) const;

  // int emtfMatchAlgoVersion_ ;         
  AlgoType emtfMatchAlgoVersion_ ;         

  std::unique_ptr<L1TkMuCorrDynamicWindows> dwcorr_;
  bool requireBX0_;
  bool requirePhaseII_;
  int  mu_stub_station_;

  const edm::EDGetTokenT< EMTFTrackCollection >          emtfTCToken; // the track collection, directly from the EMTF and not formatted by GT
  const edm::EDGetTokenT< EMTFHitCollection >            emtfHCToken; // the hit collection, directly from the EMTF which stored the input Hits
  const edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > trackToken;
} ;


L1TkMuonStubProducer::L1TkMuonStubProducer(const edm::ParameterSet& iConfig) :
  emtfTCToken(consumes< EMTFTrackCollection >         (iConfig.getParameter<edm::InputTag>("L1EMTFTrackCollectionInputTag"))),
  emtfHCToken(consumes< EMTFHitCollection >           (iConfig.getParameter<edm::InputTag>("L1EMTFHitCollectionInputTag"))),
  trackToken (consumes< std::vector<TTTrack< Ref_Phase2TrackerDigi_> > > (iConfig.getParameter<edm::InputTag>("L1TrackInputTag")))
{
   
   // configuration of the EMTF algorithm type
   std::string emtfMatchAlgoVersionString = iConfig.getParameter<std::string>("emtfMatchAlgoVersion");
   std::transform(emtfMatchAlgoVersionString.begin(), emtfMatchAlgoVersionString.end(), emtfMatchAlgoVersionString.begin(), ::tolower); // make lowercase

   mu_stub_station_ = iConfig.getParameter<int>("mu_stub_station");
   requireBX0_ = iConfig.getParameter<bool>("require_BX0");
   requirePhaseII_ = iConfig.getParameter<bool>("require_PhaseII");

   if (emtfMatchAlgoVersionString == "dynamicwindows")
      emtfMatchAlgoVersion_ = kDynamicWindows;
   else
    throw cms::Exception("TkMuAlgoConfig") << "the ID of the EMTF algo matcher passed is invalid\n";
   

   produces<L1TkMuonParticleCollection>();
   produces<L1TkMuonParticleCollection>("ME0Ext");
   produces<EMTFHitCollection>("Cleaned");

   // initializations
   if (emtfMatchAlgoVersion_ == kDynamicWindows)
   {
      // FIXME: to merge eventually into an unique file with bith phi and theta boundaries
      std::string fIn_bounds_name = iConfig.getParameter<edm::FileInPath>("emtfcorr_boundaries").fullPath();
      std::string fIn_theta_name  = iConfig.getParameter<edm::FileInPath>("emtfcorr_theta_windows").fullPath();
      std::string fIn_phi_name    = iConfig.getParameter<edm::FileInPath>("emtfcorr_phi_windows").fullPath();
      std::string fIn_S1_theta_name  = iConfig.getParameter<edm::FileInPath>("emtfcorr_S1_theta_windows").fullPath();
      std::string fIn_S1_phi_name    = iConfig.getParameter<edm::FileInPath>("emtfcorr_S1_phi_windows").fullPath();
      auto bounds = L1TkMuCorrDynamicWindows::prepare_corr_bounds(fIn_bounds_name.c_str(), "h_dphi_l");
      TFile* fIn_theta = TFile::Open (fIn_theta_name.c_str());
      TFile* fIn_phi   = TFile::Open (fIn_phi_name.c_str());
      TFile* fIn_S1_theta = TFile::Open (fIn_theta_name.c_str());
      TFile* fIn_S1_phi   = TFile::Open (fIn_phi_name.c_str());
      dwcorr_ = std::unique_ptr<L1TkMuCorrDynamicWindows> (new L1TkMuCorrDynamicWindows(bounds, fIn_theta, fIn_phi, fIn_S1_theta, fIn_S1_phi));

      // files can be closed since the correlator code clones the TF1s
      fIn_theta->Close();
      fIn_phi->Close();

      // FIXME: more initialisation using the parameters passed from the cfg
      dwcorr_->set_safety_factor  (iConfig.getParameter<double>("final_window_factor"));
      dwcorr_->set_sf_initialrelax(iConfig.getParameter<double>("initial_window_factor"));
      
      dwcorr_->set_relaxation_pattern(
        iConfig.getParameter<double>("pt_start_relax"),
        iConfig.getParameter<double>("pt_end_relax")
        );
      
      dwcorr_->set_do_relax_factor(iConfig.getParameter<bool>("do_relax_factors"));

      //
      dwcorr_ -> set_n_trk_par       (iConfig.getParameter<int>("n_trk_par"));
      dwcorr_ -> set_min_trk_p       (iConfig.getParameter<double>("min_trk_p"));
      dwcorr_ -> set_max_trk_aeta    (iConfig.getParameter<double>("max_trk_aeta"));
      dwcorr_ -> set_max_trk_chi2    (iConfig.getParameter<double>("max_trk_chi2"));
      dwcorr_ -> set_min_trk_nstubs  (iConfig.getParameter<int>("min_trk_nstubs"));
      dwcorr_ -> set_do_trk_qual_presel(true);
   }
}

L1TkMuonStubProducer::~L1TkMuonStubProducer() {
}

// ------------ method called to produce the data  ------------
void
L1TkMuonStubProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // the L1Mu objects
  edm::Handle<EMTFTrackCollection> l1emtfTCH;
  edm::Handle<EMTFHitCollection> l1emtfHCH;

  iEvent.getByToken(emtfTCToken, l1emtfTCH);
  iEvent.getByToken(emtfHCToken, l1emtfHCH);

  // the L1Tracks
  edm::Handle<L1TTTrackCollectionType> l1tksH;
  iEvent.getByToken(trackToken, l1tksH);

  L1TkMuonParticleCollection  oc_endcap_tkmuonStub;
  L1TkMuonParticleCollection  oc_me0Extended_tkmuonStub;
  EMTFHitCollection           oc_cleanedStub;

  // process each of the MTF collections separately! -- we don't want to filter the muons
  //if (emtfMatchAlgoVersion_ == kDynamicWindows) 
    runOnMuonHitCollection(l1emtfHCH, l1tksH, oc_endcap_tkmuonStub);
  //else
    //throw cms::Exception("TkMuAlgoConfig") << "trying to run an invalid algorithm (this should never happen)\n";

  // now combine all trk muons into a single output collection!
  std::unique_ptr<L1TkMuonParticleCollection> oc_tkmuon(new L1TkMuonParticleCollection());
  for (const auto& p : {oc_endcap_tkmuonStub}){
    oc_tkmuon->insert(oc_tkmuon->end(), p.begin(), p.end());
  }

   getExtendedEtaME0Muons(l1emtfHCH,oc_me0Extended_tkmuonStub, oc_cleanedStub);

  // now combine all trk muons into a single output collection!
  std::unique_ptr<L1TkMuonParticleCollection> oc_me0Ext(new L1TkMuonParticleCollection());
  for (const auto& p : {oc_me0Extended_tkmuonStub}){
    oc_me0Ext->insert(oc_me0Ext->end(), p.begin(), p.end());
  }

  // now put cleaned stubs in a output collection!
  std::unique_ptr<EMTFHitCollection> oc_stub(new EMTFHitCollection());
  for (const auto& p : {oc_cleanedStub}){
    oc_stub->insert(oc_stub->end(), p.begin(), p.end());
  }

   
  // put the new track+muon objects in the event!
  iEvent.put( std::move(oc_tkmuon));
  // put the new me0  objects in the event!
  iEvent.put( std::move(oc_me0Ext),"ME0Ext");
  // put the clean stubs in the event!
  iEvent.put( std::move(oc_stub),"Cleaned");
};


void
L1TkMuonStubProducer::runOnMuonHitCollection(const edm::Handle<EMTFHitCollection>& muonStubH,
                                     const edm::Handle<L1TTTrackCollectionType>& l1tksH,
                                     L1TkMuonParticleCollection& tkMuons) const

{
  // collection of stubs from the event
  const EMTFHitCollection& l1muStubs = (*muonStubH.product());

  // collection for cleaned stubs
  EMTFHitCollection cleanedStubs;

  // reserve to size of original coolection
  cleanedStubs.reserve(l1muStubs.size());
  
  // fill collection of cleaned stubs
  cleanStubs(l1muStubs, cleanedStubs);
  
  const L1TTTrackCollectionType& l1trks = (*l1tksH.product());
  std::vector< std::vector<int> > matchedStubsPerTrack;
  auto corr_muStub_idxs = dwcorr_->find_match_stub(cleanedStubs, l1trks, matchedStubsPerTrack, mu_stub_station_, requireBX0_, requirePhaseII_);
  // it's a vector with as many entries as the L1TT vector.
  // >= 0 : the idx in the muStub vector of matched muStubs
  // < 0: no match
  
  // sanity check
  if (corr_muStub_idxs.size() != l1trks.size())
    throw cms::Exception("TkMuAlgoOutput") << "the size of tkmu indices does not match the size of input trk collection\n";
  
  for (uint il1ttrack = 0; il1ttrack < corr_muStub_idxs.size(); ++il1ttrack)
  {
    int muStub_idx = corr_muStub_idxs.at(il1ttrack);
    if (muStub_idx < 0)
      continue;

    std::vector<int> tkMatchedStubs = matchedStubsPerTrack.at(il1ttrack);
    // cout << "      TkMuStub: i_Track = " << il1ttrack << " closest_matched_stub = " << muStub_idx << " , ( all matched stubs = ";
    // for (auto ssm = tkMatchedStubs.begin(); ssm != tkMatchedStubs.end(); ++ssm)
    //  cout << *ssm << ' ';
    // cout << " ) " << endl;

    const L1TTTrackType& matchTk = l1trks[il1ttrack];
    const auto& p3 = matchTk.getMomentum(dwcorr_->get_n_trk_par());
    const auto& tkv3 = matchTk.getPOCA(dwcorr_->get_n_trk_par());
    const auto& curve = matchTk.getRInv( dwcorr_->get_n_trk_par());
    float p4e = sqrt(0.105658369*0.105658369 + p3.mag2() );
    math::XYZTLorentzVector l1tkp4(p3.x(), p3.y(), p3.z(), p4e);

    edm::Ptr< L1TTTrackType > l1tkPtr(l1tksH, il1ttrack);
    
    float trkisol = -999; // FIXME: now doing as in the TP algo

    const EMTFHit & muStub = cleanedStubs[muStub_idx];
    int i_ref = findStubRefIndex(muonStubH, muStub);

    edm::Ref< EMTFHitCollection > stubRef(muonStubH,i_ref);

    L1TkMuonParticle l1tkmuStub(l1tkp4, stubRef, l1tkPtr, trkisol);

    l1tkmuStub.setTrkzVtx( (float)tkv3.z() );
    
    // Set curvature
    l1tkmuStub.setTrackCurvature(curve);
    //Set charge
    int charge =  (curve > 0) ? 1 : -1;
    l1tkmuStub.setCharge(charge);

    // add references to stubs from station 1 and 2 which are matched to this track 
    for (auto ssm = tkMatchedStubs.begin(); ssm != tkMatchedStubs.end(); ++ssm) {

      if(*ssm == -1) continue; 

      // find index in the original list
      const EMTFHit & muStub = cleanedStubs[*ssm];
      int i_ref = findStubRefIndex(muonStubH, muStub);

      // store refs
      edm::Ref< EMTFHitCollection > matchedStubRef(muonStubH,i_ref);
      l1tkmuStub.addMuonStub(matchedStubRef);

    }
    
    tkMuons.push_back(l1tkmuStub);

  }


  return;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TkMuonStubProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void 
L1TkMuonStubProducer::cleanStubs(const EMTFHitCollection &  muStubs, EMTFHitCollection & cleanedStubs) const {

    // if empty collection don't do anything
    if(muStubs.size() == 0) return;
    
    // copy the first stub in the new collection
    const EMTFHit & muStub = muStubs[0];
    cleanedStubs.push_back(muStub);

    for(uint ms = 1; ms < muStubs.size(); ms++) {

      const EMTFHit & muStub = muStubs[ms];

      int n_duplicate = 0;

      for(uint i = 0; i <cleanedStubs.size(); i++) {

        const EMTFHit & cStub = cleanedStubs[i];
        int dSubsystem = cStub.Subsystem() - muStub.Subsystem();
        int dStation = cStub.Station() - muStub.Station();
        int dChamber = cStub.Chamber() - muStub.Chamber();
        int dBend = cStub.Bend() - muStub.Bend();
        float aDeltaPhi = abs(cStub.Phi_sim() * TMath::Pi()/180. - muStub.Phi_sim() * TMath::Pi()/180.);
        //cout << "   aDeltaPhi = " << aDeltaPhi << endl;

        // duplicate stubs defined as having same phi
        if(aDeltaPhi < 0.0001 && dSubsystem == 0 && dStation == 0 && dChamber <= 1 && dBend == 0) {

          //cout << "Found duplicates: " 
            //<< " (eta1,phi1) = (" << cStub.Eta_sim() << "," << cStub.Phi_sim() << ")" 
            //<< " (eta2,phi2) = (" << muStub.Eta_sim() << "," << muStub.Phi_sim() << ")" 
            //<< endl;
          n_duplicate++;

        } // end if

      } // end for cleaned

      // did not find any duplicate stubs and the stub is not Neighbor
      if(n_duplicate == 0 && muStub.Neighbor() == 0) {
        cleanedStubs.push_back(muStub);
      }

    } // end for muStubs

    /*
    printf("----------------- Cleaned stubs ---------------------------------------------------- \n");
    for(uint i = 0; i <cleanedStubs.size(); i++) {
        const EMTFHit & cStub = cleanedStubs[i];
        printStub(cStub);
    }
    printf("------------------------------------------------------------------------------- \n");
    */

}
void 
L1TkMuonStubProducer::printStubs(const EMTFHitCollection &  muStubs, int subsystem ) const {

    int subsystem_ = subsystem;
    int station_ = 1;
    int ring_ = 1;
    printf("----------------- stubs: subsystem = %1d  station = %1d   ring = %1d --------------- \n", subsystem_, station_, ring_ );
    for(uint ms = 0; ms < muStubs.size(); ms++) {

      const EMTFHit & muStub = muStubs[ms];
      if(muStub.Subsystem() == subsystem_ && muStub.Station() == station_ && muStub.Ring() == ring_) {
        printStub(muStub);
      }
    }

    ring_ = 2;
    printf("----------------- stubs: subsystem = %1d  station = %1d   ring = %1d --------------- \n", subsystem_, station_, ring_ );
    for(uint ms = 0; ms < muStubs.size(); ms++) {

      const EMTFHit & muStub = muStubs[ms];
      if(muStub.Subsystem() == subsystem_ && muStub.Station() == station_ && muStub.Ring() == ring_) {
        printStub(muStub);
      }
    }

    printf("------------------------------------------------------------------------------- \n");

}

void 
L1TkMuonStubProducer::printStub(const EMTFHit &  muStub ) const {

        printf("stub subsystem = %1d, station = %1d,  ring = %1d, chamber = %1d, strip = %1d, isNeighbor = %1d, quality = %1d, eta = %5.4f, phi = %5.4f, rho = %3.1f, bend = %3d, pattern = %3d, stubNum = %3d, BC0 = %1d, BX = %1d, valid = %1d \n", muStub.Subsystem(),  muStub.Station(), muStub.Ring(), muStub.Chamber(), muStub.Strip(), muStub.Neighbor(), muStub.Quality(), muStub.Eta_sim(), muStub.Phi_sim() * TMath::Pi()/180.,  muStub.Rho_sim(), muStub.Bend(), muStub.Pattern(), muStub.Stub_num(), muStub.BC0(), muStub.BX(), muStub.Valid());

}

int
L1TkMuonStubProducer::findStubRefIndex(const edm::Handle<EMTFHitCollection>& l1muonStubH, 
                                          const EMTFHit & muStub) const
{

  // dumb function.  it would be much more elegant to save the cleaned collection.
  int rc = -99;

  const EMTFHitCollection& l1muStubs = (*l1muonStubH.product());

  for(uint i = 0; i < l1muStubs.size(); i++) {

    const EMTFHit & eventStub = l1muStubs[i]; 

    if( eventStub.Eta_sim() == muStub.Eta_sim() &&
        eventStub.Phi_sim() == muStub.Phi_sim() &&
        eventStub.Rho_sim() == muStub.Rho_sim() ) return i;

  }

  return rc;

}

void
L1TkMuonStubProducer::getExtendedEtaME0Muons(const edm::Handle<EMTFHitCollection>& muonStubH, L1TkMuonParticleCollection& tkMuons, EMTFHitCollection & cleanedStubs) const

{
  const EMTFHitCollection& l1muStubs = (*muonStubH.product());

  // collection for cleaned stubs
  //EMTFHitCollection cleanedStubs;

  // reserve to size of original coolection
  cleanedStubs.reserve(l1muStubs.size());

  // fill collection of cleaned stubs
  cleanStubs(l1muStubs, cleanedStubs);

  //cout << "Numbe of cleaned stubs = " << cleanedStubs.size() << endl;

  for (auto muStub : cleanedStubs) {


    if(muStub.Subsystem() != EMTFHit::kME0) continue;

    
    float stubBend = muStub.Bend();
    stubBend *= 3.63/1000.; // in rad
    float pt = 0.25 * 1.0 / abs(stubBend); // in GeV

    float eta = muStub.Eta_sim();
    float phi = muStub.Phi_sim() * TMath::Pi()/180.;

    //cout << " Checkign the eta" << endl;
    if(abs(eta) < 2.4) continue;

    //cout << "Checking stub for pt = " << pt  << endl;
    if(pt < 5.0) continue;
    
    reco::Candidate::PolarLorentzVector muonLV(pt,eta,phi,0);
    L1TkMuonParticle muon(muonLV);

    int charge = 1;
    if(stubBend > 0) charge = -1;

    int quality = muStub.Quality();

    
    muon.setCharge(charge);
    muon.setQuality(quality);

    //cout << "Storing ME0Ext saving pt = " << pt << endl;
    // build a L1Candidate
    tkMuons.push_back( muon);

  }

  return;
}


//define this as a plug-in
DEFINE_FWK_MODULE(L1TkMuonStubProducer);



