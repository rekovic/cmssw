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

      float max_dR_;
      float max_dZ_;
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
  max_dR_(iConfig.getParameter<double>("max_dR")) ,
  max_dZ_(iConfig.getParameter<double>("max_dZ")) ,
  debug_(iConfig.getUntrackedParameter<bool>("debug", false))
{
   //register your products
   produces<MuonJetCollection>("3Stub");
   produces<MuonJetCollection>("1TkMuStub2Stub");
   produces<MuonJetCollection>("2TkMuStub1Stub");
   produces<MuonJetCollection>("3TkMuStub");
   
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
  

  auto out_3Stub_muonjets              = std::make_unique<MuonJetCollection>(); 

  auto out_1TkMuStub2Stub_muonjets     = std::make_unique<MuonJetCollection>(); 
  auto out_2TkMuStub1Stub_muonjets     = std::make_unique<MuonJetCollection>(); 
  auto out_3TkMuStub_muonjets          = std::make_unique<MuonJetCollection>(); 


  for (const auto& p : {muonjets}){
    out_3TkMuStub_muonjets->insert(out_3TkMuStub_muonjets->end(), p.begin(), p.end());
  }

  iEvent.put( std::move(out_3Stub_muonjets),            "3Stub" );
  iEvent.put( std::move(out_1TkMuStub2Stub_muonjets),   "1TkMuStub2Stub");
  iEvent.put( std::move(out_2TkMuStub1Stub_muonjets),   "2TkMuStub1Stub");
  iEvent.put( std::move(out_3TkMuStub_muonjets),        "3TkMuStub");
 
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
L1TMuonJetProducer::findMuonJets(const edm::Handle<EMTFHitCollection>& l1muonStubH,
                                     const edm::Handle<L1TkMuonParticleCollection>& l1tkmuonStubH,
                                     MuonJetCollection& muonJets) const
{
  const EMTFHitCollection& l1muStubs = (*l1muonStubH.product());
  const L1TkMuonParticleCollection& l1tkmuStubs = (*l1tkmuonStubH.product());

  // 1st loop over TkMuStubs
  for(uint i_tkms = 0; i_tkms < l1tkmuStubs.size(); i_tkms++) {

    const L1TkMuonParticle& tkmuStub = l1tkmuStubs[i_tkms]; 

    if(debug_) {
      cout << " tkmustub pt " << tkmuStub.pt() << " eta " << tkmuStub.eta() << " phi " << tkmuStub.phi() << " charge " << tkmuStub.charge()  << " z_vtx " << tkmuStub.getTrkzVtx() << endl;
    }
    

    // 2nd loop over TkMuStubs
    for (uint j_tkms = i_tkms + 1; j_tkms < l1tkmuStubs.size(); j_tkms++) {

      // 3rd loop over TkMuStubs
      for (uint k_tkms = j_tkms + 1; k_tkms < l1tkmuStubs.size(); k_tkms++) {

        MuonJet muonJet(l1tkmuStubs[i_tkms],l1tkmuStubs[j_tkms],l1tkmuStubs[k_tkms]);

        /*
        */

        muonJet.configure(max_dR_, max_dR_, max_dR_, max_dZ_);
        muonJet.process();
        //if(debug_) muonJet.print(); 

        // check if Valid() - max_dz, max_dR
        if (muonJet.isValid()) {

          muonJets.push_back(muonJet);
          if(debug_) muonJet.print();
          
        } // end if valid
        
      } // end for k_tkms


      // loop over muon stubs
      for(uint k_ms = 0; k_ms < l1muStubs.size(); k_ms++) {

        const EMTFHit & muStub = l1muStubs[k_ms];
        
        // only 1st station
        if (muStub.Station() != 1) continue;

        // only CSC
        if (muStub.Subsystem() != 1) continue;

        MuonJet muonJet_2TkMuStub_1MuStub(l1tkmuStubs[i_tkms],l1tkmuStubs[j_tkms],l1muStubs[k_ms]);

        muonJet_2TkMuStub_1MuStub.configure(max_dR_, max_dR_, max_dR_, max_dZ_);
        muonJet_2TkMuStub_1MuStub.process();

        // check if Valid() - max_dz, max_dR
        if (muonJet_2TkMuStub_1MuStub.isValid()) {

          muonJets.push_back(muonJet_2TkMuStub_1MuStub);
          if(debug_) muonJet_2TkMuStub_1MuStub.print();
          
        } // end if valid

      } // end for k_ms muon stub

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
