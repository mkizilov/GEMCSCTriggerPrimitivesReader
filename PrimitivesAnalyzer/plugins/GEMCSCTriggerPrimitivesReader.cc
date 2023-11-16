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

using namespace edm;

class GEMCSCTriggerPrimitivesReader : public edm::one::EDAnalyzer<> {
public:
  explicit GEMCSCTriggerPrimitivesReader(const edm::ParameterSet&);
  ~GEMCSCTriggerPrimitivesReader(){};
private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  TTree* bookTTree_ALCT();
  TTree* bookTTree_CLCT();
  TTree* bookTTree_LCT();
  void SaveALCTs(const CSCALCTDigiCollection* alcts, bool is_data, bool is_emul);
  void SaveCLCTs(const CSCCLCTDigiCollection* clcts, bool is_data, bool is_emul);
  void SaveLCTs(const CSCCorrelatedLCTDigiCollection* lcts, bool is_data, bool is_emul);

  int eventsAnalyzed;

  // TTree variables
  edm::Service<TFileService> fs;
  TTree* tree;
  TTree* ALCT_tree;
  TTree* CLCT_tree;
  TTree* LCT_tree;
  TTree* t;
  bool t_is_data;
  bool t_is_emul;
  int t_endcap;
  int t_station;
  int t_ring;
  int t_chamber;
  int t_RUN;
  int t_Event;
  int t_eventsAnalyzed;
  bool t_isValid;
  int t_quality;
  int t_keyWG;
  int t_strip;
  int t_quadStrip;
  int t_eightStrip;
  int t_stripType;
  int t_bend;
  int t_slope;
  int t_bx;
  int t_pattern;
  int t_run3pattern;


    
  // Run number, Event number
  int RUN_;
  int Event_;


  edm::EDGetTokenT<CSCALCTDigiCollection> alcts_d_token_;
  edm::EDGetTokenT<CSCCLCTDigiCollection> clcts_d_token_;
  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> lcts_tmb_d_token_;

  edm::EDGetTokenT<CSCALCTDigiCollection> alcts_e_token_;
  edm::EDGetTokenT<CSCCLCTDigiCollection> clcts_e_token_;
  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> lcts_tmb_e_token_;

  // // Producer's labels
  std::string lctProducerData_;
  std::string lctProducerEmul_;

  bool debug;
};


GEMCSCTriggerPrimitivesReader::GEMCSCTriggerPrimitivesReader(const edm::ParameterSet& iConfig) : eventsAnalyzed(0) {
    lctProducerData_ = iConfig.getUntrackedParameter<std::string>("CSCLCTProducerData", "muonCSCDigis");
    lctProducerEmul_ = iConfig.getUntrackedParameter<std::string>("CSCLCTProducerEmul", "cscTriggerPrimitiveDigis");
    std::cout<<"lctProducerData: "<<lctProducerData_<<std::endl;
    std::cout<<"lctProducerEmul: "<<lctProducerEmul_<<std::endl;
    debug = iConfig.getParameter<bool>("debug");

    //consumes Data
    alcts_d_token_ = consumes<CSCALCTDigiCollection>(edm::InputTag(lctProducerData_, "MuonCSCALCTDigi"));
    clcts_d_token_ = consumes<CSCCLCTDigiCollection>(edm::InputTag(lctProducerData_, "MuonCSCCLCTDigi"));
    lcts_tmb_d_token_ = consumes<CSCCorrelatedLCTDigiCollection>(edm::InputTag(lctProducerData_, "MuonCSCCorrelatedLCTDigi"));

    //consumes Emul
    alcts_e_token_ = consumes<CSCALCTDigiCollection>(edm::InputTag(lctProducerEmul_));
    clcts_e_token_ = consumes<CSCCLCTDigiCollection>(edm::InputTag(lctProducerEmul_));
    lcts_tmb_e_token_ = consumes<CSCCorrelatedLCTDigiCollection>(edm::InputTag(lctProducerEmul_));

    // //TTree varialbes
    // RUN_ = 0;
    // Event_ = 0;
    // t_is_data = false;
    // t_is_emul = false;
    // t_endcap = 0;
    // t_station = 0;
    // t_ring = 0;
    // t_chamber = 0;
    //book TTree
    ALCT_tree = bookTTree_ALCT();
    CLCT_tree = bookTTree_CLCT();
    LCT_tree = bookTTree_LCT();

  }

void
GEMCSCTriggerPrimitivesReader::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    ++eventsAnalyzed;

    RUN_ = iEvent.id().run();
    Event_ = iEvent.id().event();
    if(debug)std::cout<<"RUN: "<<RUN_<<" Event: "<<Event_<<std::endl;

    //get Data
    edm::Handle<CSCALCTDigiCollection> alcts_data;
    edm::Handle<CSCCLCTDigiCollection> clcts_data;
    edm::Handle<CSCCorrelatedLCTDigiCollection> lcts_tmb_data;
    iEvent.getByToken(alcts_d_token_, alcts_data);
    iEvent.getByToken(clcts_d_token_, clcts_data);
    iEvent.getByToken(lcts_tmb_d_token_, lcts_tmb_data);


    //get Emul
    edm::Handle<CSCALCTDigiCollection> alcts_emul;
    edm::Handle<CSCCLCTDigiCollection> clcts_emul;
    edm::Handle<CSCCorrelatedLCTDigiCollection> lcts_tmb_emul;
    iEvent.getByToken(alcts_e_token_, alcts_emul);
    iEvent.getByToken(clcts_e_token_, clcts_emul);
    iEvent.getByToken(lcts_tmb_e_token_, lcts_tmb_emul);
    
        // Save ALCTs
    if (alcts_data.isValid()) {
      bool is_data = true;
      bool is_emul = false;
      SaveALCTs(alcts_data.product(), is_data, is_emul);
    }
    else if(debug) std::cout<<"Data ALCTs not found"<<std::endl;

    if (alcts_emul.isValid()) {
      bool is_data = false;
      bool is_emul = true;
      SaveALCTs(alcts_emul.product(), is_data, is_emul);
    }
    else if(debug) std::cout<<"Emul ALCTs not found"<<std::endl;
    // Save CLCTs
    if (clcts_data.isValid()) {
      bool is_data = true;
      bool is_emul = false;
      SaveCLCTs(clcts_data.product(), is_data, is_emul);
    }
    else if(debug) std::cout<<"Data CLCTs not found"<<std::endl;
    if (clcts_emul.isValid()) {
      bool is_data = false;
      bool is_emul = true;
      SaveCLCTs(clcts_emul.product(), is_data, is_emul);
    }
    else if(debug) std::cout<<"Emul CLCTs not found"<<std::endl;
    // Save LCTs
    if (lcts_tmb_data.isValid()) {
      bool is_data = true;
      bool is_emul = false;
      SaveLCTs(lcts_tmb_data.product(), is_data, is_emul);
    }
    else if(debug) std::cout<<"Data LCTs not found"<<std::endl;
    if (lcts_tmb_emul.isValid()) {
      bool is_data = false;
      bool is_emul = true;
      SaveLCTs(lcts_tmb_emul.product(), is_data, is_emul);
    }
    else if(debug) std::cout<<"Emul LCTs not found"<<std::endl;

}






void GEMCSCTriggerPrimitivesReader::SaveALCTs(const CSCALCTDigiCollection* alcts, bool is_data, bool is_emul)
{
  for (int endc = 1; endc <= 2; endc++) {
    for (int stat = 1; stat <= 4; stat++) {
      for (int ring = 1; ring <= 4; ring++) {
        for (int cham = 1; cham <= 36; cham++) {
          // Calculate DetId.  0th layer means whole chamber.
          CSCDetId detid(endc, stat, ring, cham, 0);
          std::vector<CSCALCTDigi> alctV;
          const auto& range = alcts->get(detid);
          for (auto digiIt = range.first; digiIt != range.second; digiIt++) {
            if ((*digiIt).isValid()) {
              std::cout << "ALCT "  <<(*digiIt) << " data:emul "<< is_data <<":"<< is_emul << std::endl;
              alctV.push_back(*digiIt);
              //Fill TTree
              t_RUN = RUN_;
              t_Event = Event_;
              t_eventsAnalyzed = eventsAnalyzed;
              t_is_data = is_data;
              t_is_emul = is_emul;
              t_endcap = endc;
              t_station = stat;
              t_ring = ring;
              t_chamber = cham;
              ALCT_tree->Fill();
            }
          }
        }
      }
    }
  }
}

void GEMCSCTriggerPrimitivesReader::SaveCLCTs(const CSCCLCTDigiCollection* clcts, bool is_data, bool is_emul)
{
  for (int endc = 1; endc <= 2; endc++) {
    for (int stat = 1; stat <= 4; stat++) {
      for (int ring = 1; ring <= 4; ring++) {
        for (int cham = 1; cham <= 36; cham++) {
          // Calculate DetId.  0th layer means whole chamber.
          CSCDetId detid(endc, stat, ring, cham, 0);
          std::vector<CSCCLCTDigi> clctV;
          const auto& range = clcts->get(detid);
          for (auto digiIt = range.first; digiIt != range.second; digiIt++) {
            if ((*digiIt).isValid()) {
              if (debug) std::cout << "CLCT " << (*digiIt) <<" data:emul "<< is_data <<":"<< is_emul << std::endl;
              clctV.push_back(*digiIt);
              //Fill TTree
              t_RUN = RUN_;
              t_Event = Event_;
              t_eventsAnalyzed = eventsAnalyzed;
              t_is_data = is_data;
              t_is_emul = is_emul;
              t_endcap = endc;
              t_station = stat;
              t_ring = ring;
              t_chamber = cham;

              CLCT_tree->Fill();
            }
          }
        }
      }
    }
  }
}

void GEMCSCTriggerPrimitivesReader::SaveLCTs(const CSCCorrelatedLCTDigiCollection* lcts, bool is_data, bool is_emul)
{
  for (int endc = 1; endc <= 2; endc++) {
    for (int stat = 1; stat <= 4; stat++) {
      for (int ring = 1; ring <= 4; ring++) {
        for (int cham = 1; cham <= 36; cham++) {
          // Calculate DetId.  0th layer means whole chamber.
          CSCDetId detid(endc, stat, ring, cham, 0);
          std::vector<CSCCorrelatedLCTDigi> lctV;
          const auto& range = lcts->get(detid);
          for (auto digiIt = range.first; digiIt != range.second; digiIt++) {
            if ((*digiIt).isValid()) {
              if (debug) std::cout << "LCT " << (*digiIt) <<" data:emul "<< is_data <<":"<< is_emul << std::endl;
              lctV.push_back(*digiIt);
              //Fill TTree
              t_RUN = RUN_;
              t_Event = Event_;
              t_eventsAnalyzed = eventsAnalyzed;
              t_is_data = is_data;
              t_is_emul = is_emul;
              t_endcap = endc;
              t_station = stat;
              t_ring = ring;
              t_chamber = cham;
              t_isValid = (*digiIt).isValid();
              t_quality = (*digiIt).getQuality();
              t_keyWG = (*digiIt).getKeyWG();
              t_strip = (*digiIt).getStrip();
              t_quadStrip = (*digiIt).getStrip(4);
              t_eightStrip = (*digiIt).getStrip(8);
              t_stripType = (*digiIt).getStripType();
              t_bend = (*digiIt).getBend();
              t_slope = (*digiIt).getSlope();
              t_bx = (*digiIt).getBX();
              t_pattern = (*digiIt).getPattern();
              t_run3pattern = (*digiIt).getRun3Pattern();


              LCT_tree->Fill();
            }
          }
        }
      }
    }
  }
}

TTree* GEMCSCTriggerPrimitivesReader::bookTTree_ALCT() {
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("ALCT_tree", "GEMCSCTriggerPrimitivesReader");
  t->Branch("RUN", &t_RUN, "RUN/I");
  t->Branch("Event", &t_Event, "Event/I");
  t->Branch("eventsAnalyzed", &t_eventsAnalyzed, "eventsAnalyzed/I");
  t->Branch("is_data", &t_is_data, "is_data/O");
  t->Branch("is_emul", &t_is_emul, "is_emul/O");
  t->Branch("endcap", &t_endcap, "endcap/I");
  t->Branch("station", &t_station, "station/I");
  t->Branch("ring", &t_ring, "ring/I");
  t->Branch("chamber", &t_chamber, "chamber/I");
  return t;
}
TTree* GEMCSCTriggerPrimitivesReader::bookTTree_CLCT(){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("CLCT_tree", "GEMCSCTriggerPrimitivesReader");
  t->Branch("RUN", &t_RUN, "RUN/I");
  t->Branch("Event", &t_Event, "Event/I");
  t->Branch("eventsAnalyzed", &t_eventsAnalyzed, "eventsAnalyzed/I");
  t->Branch("is_data", &t_is_data, "is_data/O");
  t->Branch("is_emul", &t_is_emul, "is_emul/O");
  t->Branch("endcap", &t_endcap, "endcap/I");
  t->Branch("station", &t_station, "station/I");
  t->Branch("ring", &t_ring, "ring/I");
  t->Branch("chamber", &t_chamber, "chamber/I");
  return t;

}
TTree* GEMCSCTriggerPrimitivesReader::bookTTree_LCT(){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("LCT_tree", "GEMCSCTriggerPrimitivesReader");
  t->Branch("RUN", &t_RUN, "RUN/I");
  t->Branch("Event", &t_Event, "Event/I");
  t->Branch("eventsAnalyzed", &t_eventsAnalyzed, "eventsAnalyzed/I");
  t->Branch("is_data", &t_is_data, "is_data/O");
  t->Branch("is_emul", &t_is_emul, "is_emul/O");
  t->Branch("endcap", &t_endcap, "endcap/I");
  t->Branch("station", &t_station, "station/I");
  t->Branch("ring", &t_ring, "ring/I");
  t->Branch("chamber", &t_chamber, "chamber/I");
  t->Branch("isValid", &t_isValid, "isValid/O");
  t->Branch("quality", &t_quality, "quality/I");
  t->Branch("keyWG", &t_keyWG, "keyWG/I");
  t->Branch("strip", &t_strip, "strip/I");
  t->Branch("quadStrip", &t_quadStrip, "quadStrip/I");
  t->Branch("eightStrip", &t_eightStrip, "eightStrip/I");
  t->Branch("stripType", &t_stripType, "stripType/I");
  t->Branch("bend", &t_bend, "bend/I");
  t->Branch("slope", &t_slope, "slope/I");
  t->Branch("bx", &t_bx, "bx/I");
  t->Branch("pattern", &t_pattern, "pattern/I");
  t->Branch("run3pattern", &t_run3pattern, "run3pattern/I");
  return t;
}

DEFINE_FWK_MODULE(GEMCSCTriggerPrimitivesReader);