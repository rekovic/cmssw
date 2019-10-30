// -*- C++ -*-
//
// Package:    L1Trigger/L1TMuonJet
// Class:      L1TMuonJetProducer
// 
/**\class L1TMuonJetProducer L1TMuonJetProducer.cc L1Trigger/L1TMuonJetProducer/plugins/L1TMuonJetProducer.cc

 Description: producer to fill MuonJet class with muon 3 muon stub information (pt, eta, phi)

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Vladimir Rekovic
//         Created:  Fri, 04 Oct 2019 16:13:59 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Phase2L1TMuonJet/interface/MuonJet.h"
#include "DataFormats/Phase2L1TMuonJet/interface/MuonJetFwd.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"
#include "L1Trigger/L1TMuonEndCap/interface/Common.h"
#include <vector>

#include "TMath.h"



using namespace l1t;
//
// class declaration
//

class L1TMuonJetProducer : public edm::stream::EDProducer<> {
   public:
      typedef TTTrack< Ref_Phase2TrackerDigi_ >  L1TTTrackType;
      typedef std::vector< L1TTTrackType > L1TTTrackCollectionType;
      explicit L1TMuonJetProducer(const edm::ParameterSet&);
      ~L1TMuonJetProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      void findMuonJets(const edm::Handle<EMTFHitCollection>&, const edm::Handle<L1TkMuonParticleCollection>&, MuonJetCollection& ) const;
      bool isAllowedMuonStub(const EMTFHit & ) const;
      void printStubs(const EMTFHitCollection & , int) const;
      void printStub(const EMTFHit & ) const;

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      const edm::EDGetTokenT< EMTFHitCollection >            muonStubCToken; // the hit collection, directly from the EMTF which stored the input Hits
      const edm::EDGetTokenT< L1TkMuonParticleCollection >   tkMuStubCToken; // the TkMuStub collection
      const edm::EDGetTokenT< EMTFTrackCollection >          emtfTCToken; // the track collection, directly from the EMTF and not formatted by GT
      const edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > trackToken;

      float max_dR_3_TKMUSTUB_;
      float max_dZ_3_TKMUSTUB_;

      float max_dR_2_TKMUSTUB_1_MUSTUB_;
      float max_dZ_2_TKMUSTUB_1_MUSTUB_;

      float max_dR_1_TKMUSTUB_2_MUSTUB_;
      float max_dZ_1_TKMUSTUB_2_MUSTUB_;

      int stubStation_;
      int stubSubsystem_;
      int stubMinQuality_;

      bool debug_;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
L1TMuonJetProducer::L1TMuonJetProducer(const edm::ParameterSet& iConfig) :
  muonStubCToken(consumes< EMTFHitCollection >           (iConfig.getParameter<edm::InputTag>("MuonStubCollectionInputTag"))),
  tkMuStubCToken(consumes< L1TkMuonParticleCollection>           (iConfig.getParameter<edm::InputTag>("L1TkMuStubCollectionInputTag"))),
  emtfTCToken(consumes< EMTFTrackCollection >         (iConfig.getParameter<edm::InputTag>("L1EMTFTrackCollectionInputTag"))),
  trackToken (consumes< std::vector<TTTrack< Ref_Phase2TrackerDigi_> > > (iConfig.getParameter<edm::InputTag>("L1TrackInputTag"))),
  max_dR_3_TKMUSTUB_(iConfig.getParameter<double>("max_dR_3_TKMUSTUB")) ,
  max_dZ_3_TKMUSTUB_(iConfig.getParameter<double>("max_dZ_3_TKMUSTUB")) ,
  max_dR_2_TKMUSTUB_1_MUSTUB_(iConfig.getParameter<double>("max_dR_2_TKMUSTUB_1_MUSTUB")) ,
  max_dZ_2_TKMUSTUB_1_MUSTUB_(iConfig.getParameter<double>("max_dZ_2_TKMUSTUB_1_MUSTUB")) ,
  max_dR_1_TKMUSTUB_2_MUSTUB_(iConfig.getParameter<double>("max_dR_1_TKMUSTUB_2_MUSTUB")) ,
  max_dZ_1_TKMUSTUB_2_MUSTUB_(iConfig.getParameter<double>("max_dZ_1_TKMUSTUB_2_MUSTUB")) ,
  stubStation_(iConfig.getParameter<unsigned int>("stubStation")),
  stubSubsystem_(iConfig.getParameter<unsigned int>("stubSubsystem")),
  stubMinQuality_(iConfig.getParameter<unsigned int>("stubMinQuality")),
  debug_(iConfig.getUntrackedParameter<bool>("debug", false))
{
   //register your products
   produces<MuonJetCollection>("3TkMuStub");
   produces<MuonJetCollection>("2TkMuStub");
   produces<MuonJetCollection>("2TkMuStub1MuStub");
   produces<MuonJetCollection>("1TkMuStub2MuStub");
   produces<MuonJetCollection>("3MuStub");
   
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
   //
   
  
}


L1TMuonJetProducer::~L1TMuonJetProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
L1TMuonJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   iEvent.put(std::make_unique<ExampleData2>(*pIn));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
  // the Mu stubs
  edm::Handle<EMTFHitCollection> muonStubCH;
  iEvent.getByToken(muonStubCToken, muonStubCH);

  // the L1TkMu objects
  edm::Handle<L1TkMuonParticleCollection> tkmuStubCH;
  iEvent.getByToken(tkMuStubCToken, tkmuStubCH);

  // the L1Mu objects
  //edm::Handle<EMTFTrackCollection> l1emtfTCH;
  //iEvent.getByToken(emtfTCToken, l1emtfTCH);

  // the L1Tracks
  //edm::Handle<L1TTTrackCollectionType> l1tksH;
  //iEvent.getByToken(trackToken, l1tksH);

  MuonJetCollection muonjets;

  findMuonJets(muonStubCH, tkmuStubCH, muonjets);
  

  auto out_3TkMuStub_muonjets          = std::make_unique<MuonJetCollection>(); 
  auto out_2TkMuStub_muonjets          = std::make_unique<MuonJetCollection>(); 
  auto out_2TkMuStub1Stub_muonjets     = std::make_unique<MuonJetCollection>(); 
  auto out_1TkMuStub2Stub_muonjets     = std::make_unique<MuonJetCollection>(); 
  auto out_3Stub_muonjets            = std::make_unique<MuonJetCollection>(); 


  for (const MuonJet& jet : muonjets) {

    // FIXME: make case switch
    
    if(jet.getType() == MuonJet::THREE_TKMUSTUB)              out_3TkMuStub_muonjets->push_back(jet);
    if(jet.getType() == MuonJet::TWO_TKMUSTUB)                out_2TkMuStub_muonjets->push_back(jet);
    if(jet.getType() == MuonJet::TWO_TKMUSTUB_ONE_MUSTUB)     out_2TkMuStub1Stub_muonjets->push_back(jet);
    if(jet.getType() == MuonJet::ONE_TKMUSTUB_TWO_MUSTUB)     out_1TkMuStub2Stub_muonjets->push_back(jet);
    if(jet.getType() == MuonJet::THREE_MUSTUB)                out_3Stub_muonjets->push_back(jet);

  }

  if(out_3TkMuStub_muonjets->size() != 0)       cout<< "At least one MuonJet::THREE_TKMUSTUB." << endl;
  if(out_2TkMuStub_muonjets->size() != 0)       cout<< "At least one MuonJet::TWO_TKMUSTUB." << endl;
  if(out_2TkMuStub1Stub_muonjets->size() != 0)  cout<< "At least one MuonJet::TWO_TKMUSTUB_ONE_MUSTUB." << endl;
  if(out_1TkMuStub2Stub_muonjets->size() != 0)  cout<< "At least one MuonJet::ONE_TKMUSTUB_TWO_MUSTUB." << endl;
  if(out_3Stub_muonjets->size() != 0)           cout<< "At least one MuonJet::THREE_MUSTUB." << endl;

  iEvent.put( std::move(out_3TkMuStub_muonjets),        "3TkMuStub");
  iEvent.put( std::move(out_2TkMuStub_muonjets),        "2TkMuStub");
  iEvent.put( std::move(out_2TkMuStub1Stub_muonjets),   "2TkMuStub1MuStub");
  iEvent.put( std::move(out_1TkMuStub2Stub_muonjets),   "1TkMuStub2MuStub");
  iEvent.put( std::move(out_3Stub_muonjets),            "3MuStub" );
 
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
L1TMuonJetProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
L1TMuonJetProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
L1TMuonJetProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
L1TMuonJetProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
L1TMuonJetProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
L1TMuonJetProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
void 
L1TMuonJetProducer::printStubs(const EMTFHitCollection &  muStubs, int subsystem ) const {

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
L1TMuonJetProducer::printStub(const EMTFHit &  muStub ) const {

        printf("stub subsystem = %1d, station = %1d,  ring = %1d, chamber = %1d, isNeighbor = %1d, quality = %1d, eta = %5.4f, phi = %5.4f, rho = %3.1f, bend = %3d, pattern = %3d, stubNum = %3d, BC0 = %1d, BX = %1d, valid = %1d \n", muStub.Subsystem(),  muStub.Station(), muStub.Ring(), muStub.Chamber(), muStub.Neighbor(), muStub.Quality(), muStub.Eta_sim(), muStub.Phi_sim() * TMath::Pi()/180.,  muStub.Rho_sim(), muStub.Bend(), muStub.Pattern(), muStub.Stub_num(), muStub.BC0(), muStub.BX(), muStub.Valid());

}

bool 
L1TMuonJetProducer::isAllowedMuonStub(const EMTFHit &  muStub) const {

  bool rc = true;

  int station   = muStub.Station();
  int subsystem = muStub.Subsystem();
  //int ring      = muStub.Ring();
  //int chamber   = muStub.Chamber();

  // FIXME:Selection of Stub station, ring, and subsystem is now done in MuonJet.
  // Should move it here, so can control it easier from configuration
  //
  // Station 
  if ( station != stubStation_) rc = false;

  // Subsystem --  only CSC and ME0 stubs
  //   if ( subsystem != stubSubsystem_) rc = false;
  if ( subsystem != EMTFHit::kCSC && subsystem != EMTFHit::kME0 ) rc = false;

  // only Quality 
  if ( muStub.Quality() < stubMinQuality_ ) rc = false;

  // don't want neighbors
  if ( muStub.Neighbor() ) { 
    if(debug_) cout << "   This is a neighbor! " << endl;
    rc = false;
  }

  // FIXME: 
  // There are 36 chambers of 10+ degrees each 
  // We have "front" and "rear" chambers.
  // Chambers are slightly larger tha 10 degrees, so they overlap,
  // and in this overlap a low pt muon can make hits in both 
  // "front" and "rear" overlaping parts of the chambers.
  // Don't want two hits from overlapping chambers
  // so, keep need to figure out how to clean them here
  // For now, do check if Delta phi or Delta eta of two
  // hits is less than 0.01, to identify two hits from 
  // overlaped parts of chambers.
  /*
  bool result = false;
  if (subsystem == TriggerPrimitive::kCSC) {
    bool isOverlapping = !(station == 1 && ring == 3);
    // not overlapping means back
    if(isOverlapping)
    {
      bool isEven = (chamber % 2 == 0);
      // odd chambers are bolted to the iron, which faces
      // forward in 1&2, backward in 3&4, so...
      result = (station < 3) ? isEven : !isEven;
    }
  }
  */

  return rc;

}

void
L1TMuonJetProducer::findMuonJets(const edm::Handle<EMTFHitCollection>& l1muonStubH,
                                     const edm::Handle<L1TkMuonParticleCollection>& l1tkmuonStubH,
                                     MuonJetCollection& out_muonJets) const
{
  const EMTFHitCollection& l1muStubs = (*l1muonStubH.product());
  const L1TkMuonParticleCollection& l1tkmuStubs = (*l1tkmuonStubH.product());

  if(debug_) printStubs(l1muStubs, EMTFHit::kCSC);

  // 1st loop over TkMuStubs
  for(uint i_tkms = 0; i_tkms < l1tkmuStubs.size(); i_tkms++) {

    const L1TkMuonParticle& tkmuStub = l1tkmuStubs[i_tkms]; 

    if(debug_) {
      cout << " tkmustub pt " << tkmuStub.pt() << " eta " << tkmuStub.eta() << " phi " << tkmuStub.phi() << " charge " << tkmuStub.charge()  << " z_vtx " << tkmuStub.getTrkzVtx() << "  stubRef (eta,phi,station,type,qual) = (" << tkmuStub.getMuStubRef()->Eta_sim() << "," << tkmuStub.getMuStubRef()->Phi_sim() * TMath::Pi()/180. << "," << tkmuStub.getMuStubRef()->Station() << "," << tkmuStub.getMuStubRef()->Subsystem() << "," << tkmuStub.getMuStubRef()->Quality() << ")"<< endl;
    }

    // 2nd loop - loop over muon stubs
    for(uint j_ms = 0; j_ms < l1muStubs.size(); j_ms++) {
      continue;

      // get 1st stub
      //const EMTFHit & muStub_1 = l1muStubs[j_ms];
      const EMTFHitRef muStub_1_Ref(l1muonStubH,j_ms);

      if (! isAllowedMuonStub(*muStub_1_Ref) ) continue;

      // 3rd loop - loop over muon stubs
      for(uint k_ms = j_ms+1; k_ms < l1muStubs.size(); k_ms++) {

        // get 2nd stub
        //const EMTFHit & muStub_2 = l1muStubs[k_ms];
        const EMTFHitRef muStub_2_Ref(l1muonStubH,k_ms);
        
        if (! isAllowedMuonStub(*muStub_2_Ref)) continue;

        MuonJet muonJet_1tms_2ms(l1tkmuStubs[i_tkms], muStub_1_Ref, muStub_2_Ref);

        muonJet_1tms_2ms.configure(max_dR_1_TKMUSTUB_2_MUSTUB_, max_dR_1_TKMUSTUB_2_MUSTUB_, max_dR_1_TKMUSTUB_2_MUSTUB_, max_dZ_1_TKMUSTUB_2_MUSTUB_);
        muonJet_1tms_2ms.process();

        // check if withn  max_dz, max_dR
        if (muonJet_1tms_2ms.isValid()) {

          cout << " pushing ONE_TKMUSTUB_TWO_MUSTUB" << endl;
          out_muonJets.push_back(muonJet_1tms_2ms);
          if(debug_) muonJet_1tms_2ms.print();
          
        } // end if valid

      } // end for k_ms muon stub

    } // end for j_ms muon stub
    

    // 2nd loop over - loop over TkMuStubs
    for (uint j_tkms = i_tkms + 1; j_tkms < l1tkmuStubs.size(); j_tkms++) {

      MuonJet muonJet_2tms(l1tkmuStubs[i_tkms],l1tkmuStubs[j_tkms]);
      muonJet_2tms.configure(max_dR_2_TKMUSTUB_1_MUSTUB_, max_dR_2_TKMUSTUB_1_MUSTUB_, max_dR_2_TKMUSTUB_1_MUSTUB_, max_dZ_2_TKMUSTUB_1_MUSTUB_);
      muonJet_2tms.process();

      // check if withn  max_dz, max_dR
      //if (muonJet_2tms.isValid()) 
      //{

        //cout << " pushing TWO_TKMUSTUB" << endl;
        out_muonJets.push_back(muonJet_2tms);
        //if(debug_) muonJet_2tms.print();
          
      //} // end if valid

      // 3rd loop - loop over TkMuStubs
      for (uint k_tkms = j_tkms + 1; k_tkms < l1tkmuStubs.size(); k_tkms++) {

        MuonJet muonJet_3tms(l1tkmuStubs[i_tkms],l1tkmuStubs[j_tkms],l1tkmuStubs[k_tkms]);
        muonJet_3tms.configure(max_dR_3_TKMUSTUB_, max_dR_3_TKMUSTUB_, max_dR_3_TKMUSTUB_, max_dZ_3_TKMUSTUB_);
        muonJet_3tms.process();

        // check if withn  max_dz, max_dR
        if (muonJet_3tms.isValid()) 
        {

          cout << " pushing THREE_TKMUSTUB" << endl;
          out_muonJets.push_back(muonJet_3tms);
          if(debug_) muonJet_3tms.print();
          
        } // end if valid
        
      } // end for k_tkms

        const L1TkMuonParticle& tkmuStub1 = l1tkmuStubs[i_tkms]; 
        const L1TkMuonParticle& tkmuStub2 = l1tkmuStubs[j_tkms]; 
        if(debug_) cout << " tkmuStub1 pt " << tkmuStub1.pt() << " eta " << tkmuStub1.eta() << " phi " << tkmuStub1.phi() << " charge " << tkmuStub1.charge()  << " z_vtx " << tkmuStub1.getTrkzVtx() << "  stubRef (eta,phi,station,type,qual) = (" << tkmuStub1.getMuStubRef()->Eta_sim() << "," << tkmuStub1.getMuStubRef()->Phi_sim() * TMath::Pi()/180. << "," << tkmuStub1.getMuStubRef()->Station() << "," << tkmuStub1.getMuStubRef()->Subsystem() << "," << tkmuStub1.getMuStubRef()->Quality() << ")"<< endl;
        if(debug_) cout << " tkmuStub2 pt " << tkmuStub2.pt() << " eta " << tkmuStub2.eta() << " phi " << tkmuStub2.phi() << " charge " << tkmuStub2.charge()  << " z_vtx " << tkmuStub2.getTrkzVtx() << "  stubRef (eta,phi,station,type,qual) = (" << tkmuStub2.getMuStubRef()->Eta_sim() << "," << tkmuStub2.getMuStubRef()->Phi_sim() * TMath::Pi()/180. << "," << tkmuStub2.getMuStubRef()->Station() << "," << tkmuStub2.getMuStubRef()->Subsystem() << "," << tkmuStub2.getMuStubRef()->Quality() << ")"<< endl;

      // 3rd loop - loop over muon stubs
      for(uint k_ms = 0; k_ms < l1muStubs.size(); k_ms++) {

        //const EMTFHit & muStub_3 = l1muStubs[k_ms];
        const EMTFHitRef muStub_3_Ref(l1muonStubH,k_ms);

        if(debug_) cout << "Trying stub " << k_ms << endl;
        if (! isAllowedMuonStub(*muStub_3_Ref)) continue;

        if(debug_) cout <<  "       stubRef (eta,phi,station,type,qual) = (" << muStub_3_Ref->Eta_sim() << "," << muStub_3_Ref->Phi_sim() * TMath::Pi()/180. << "," << muStub_3_Ref->Station() << "," << muStub_3_Ref->Subsystem() << "," << muStub_3_Ref->Quality() << ")"<< endl;
        MuonJet muonJet_2tms_1ms(l1tkmuStubs[i_tkms],l1tkmuStubs[j_tkms], muStub_3_Ref);

        if(debug_) cout << "Print before processing ::::::::::::::::::" << endl;
        if(debug_) muonJet_2tms_1ms.print();
        

        muonJet_2tms_1ms.configure(max_dR_2_TKMUSTUB_1_MUSTUB_, max_dR_2_TKMUSTUB_1_MUSTUB_, max_dR_2_TKMUSTUB_1_MUSTUB_, max_dZ_2_TKMUSTUB_1_MUSTUB_);
        muonJet_2tms_1ms.process();
        if(debug_) cout << "Print after processing ................. " << endl;
        if(debug_) muonJet_2tms_1ms.print();

        // check if withn  max_dz, max_dR
        if (muonJet_2tms_1ms.isValid()) {

          cout << " pushing TWO_TKMUSTUB_ONE_MUSTUB" << "(i_tkms,j_tkms,k_ms) = " << i_tkms << "," << j_tkms << "," << k_ms << endl;
          out_muonJets.push_back(muonJet_2tms_1ms);
          if(debug_) muonJet_2tms_1ms.print();
          
        } // end if valid

      } // end for k_ms muon stub

      continue;

    } // end for j_tkms

  } // end for i_tkms


  return;

}
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TMuonJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TMuonJetProducer);
