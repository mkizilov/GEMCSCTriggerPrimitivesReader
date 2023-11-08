#include <assert.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

//user include files below
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/GEMDigi/interface/GEMPadDigiClusterCollection.h"
#include "DataFormats/GEMDigi/interface/GEMPadDigiCluster.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/CSCDigi/interface/CSCALCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCALCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigi.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;

class GEMCSCTriggerPrimitivesReader : public edm::one::EDAnalyzer<> {
public:
  explicit GEMCSCTriggerPrimitivesReader(const edm::ParameterSet&);
  ~GEMCSCTriggerPrimitivesReader(){};
private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  // virtual void beginJob() ;
  // virtual void endJob() ;
  int eventsAnalyzed;
    
  // Run number, Event number
  int RUN_;
  int Event_;

  edm::EDGetTokenT<CSCALCTDigiCollection> alcts_d_token_;
  edm::EDGetTokenT<CSCCLCTDigiCollection> clcts_d_token_;

  edm::EDGetTokenT<CSCALCTDigiCollection> alcts_e_token_;
  edm::EDGetTokenT<CSCCLCTDigiCollection> clcts_e_token_;

  // Producer's labels
  std::string lctProducerData_;
  std::string lctProducerEmul_;

  bool debug;
};


GEMCSCTriggerPrimitivesReader::GEMCSCTriggerPrimitivesReader(const edm::ParameterSet& iConfig) : eventsAnalyzed(0) {
    lctProducerData_ = iConfig.getUntrackedParameter<string>("CSCLCTProducerData", "cscunpacker");
    lctProducerEmul_ = iConfig.getUntrackedParameter<string>("CSCLCTProducerEmul", "cscTriggerPrimitiveDigis");
    debug = iConfig.getParameter<bool>("debug");

    //consumes Data
    alcts_d_token_ = consumes<CSCALCTDigiCollection>(edm::InputTag(lctProducerData_, "MuonCSCALCTDigi"));
    clcts_d_token_ = consumes<CSCCLCTDigiCollection>(edm::InputTag(lctProducerData_, "MuonCSCCLCTDigi"));
    // lcts_tmb_d_token_ = consumes<CSCCorrelatedLCTDigiCollection>(edm::InputTag(lctProducerData_, "MuonCSCCorrelatedLCTDigi"));

    //consumes Emul
    alcts_e_token_ = consumes<CSCALCTDigiCollection>(edm::InputTag(lctProducerEmul_));
    clcts_e_token_ = consumes<CSCCLCTDigiCollection>(edm::InputTag(lctProducerEmul_));
    // lcts_tmb_e_token_ = consumes<CSCCorrelatedLCTDigiCollection>(edm::InputTag(lctProducerEmul_));
  }

void
GEMCSCTriggerPrimitivesReader::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    ++eventsAnalyzed;

    RUN_ = iEvent.id().run();
    Event_ = iEvent.id().event();
    cout<<"RUN: "<<RUN_<<" Event: "<<Event_<<endl;

    //get Data
    edm::Handle<CSCALCTDigiCollection> alcts_data;
    edm::Handle<CSCCLCTDigiCollection> clcts_data;
    iEvent.getByToken(alcts_d_token_, alcts_data);
    iEvent.getByToken(clcts_d_token_, clcts_data);
    if (!alcts_data.isValid()) {
    edm::LogWarning("L1CSCTPEmulatorWrongInput")
        << "+++ Warning: Collection of ALCTs with label MuonCSCALCTDigi"
        << " requested, but not found in the event... Skipping the rest +++\n";
    return;
    }
    if (!clcts_data.isValid()) {
      edm::LogWarning("L1CSCTPEmulatorWrongInput")
          << "+++ Warning: Collection of CLCTs with label MuonCSCCLCTDigi"
          << " requested, but not found in the event... Skipping the rest +++\n";
      return;
    }

    //Print Data clct
    cout << "Data CLCTs" << clcts_data;
  


    //get Emul
    edm::Handle<CSCALCTDigiCollection> alcts_emul;
    edm::Handle<CSCCLCTDigiCollection> clcts_emul;
    iEvent.getByToken(alcts_e_token_, alcts_emul);
    iEvent.getByToken(clcts_e_token_, clcts_emul);
    if (!alcts_emul.isValid()) {
      edm::LogWarning("L1CSCTPEmulatorWrongInput")
          << "+++ Warning: Collection of emulated ALCTs"
          << " requested, but not found in the event... Skipping the rest +++\n";
      return;
    }
    if (!clcts_emul.isValid()) {
      edm::LogWarning("L1CSCTPEmulatorWrongInput")
          << "+++ Warning: Collection of emulated CLCTs"
          << " requested, but not found in the event... Skipping the rest +++\n";
      return;
    }

}

DEFINE_FWK_MODULE(GEMCSCTriggerPrimitivesReader);