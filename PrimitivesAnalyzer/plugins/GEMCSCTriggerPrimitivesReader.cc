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

  int t_RUN;
  int t_Event;
  int t_eventsAnalyzed;

  std::vector<bool> t_is_data;
  std::vector<bool> t_is_emul;
  std::vector<int> t_endcap;
  std::vector<int> t_station;
  std::vector<int> t_ring;
  std::vector<int> t_chamber;
  std::vector<bool> t_isValid;
  std::vector<int> t_quality;
  std::vector<int> t_keyWG;
  std::vector<int> t_strip;
  std::vector<int> t_quadStrip;
  std::vector<int> t_eightStrip;
  std::vector<int> t_stripType;
  std::vector<int> t_bend;
  std::vector<int> t_slope;
  std::vector<int> t_bx;
  std::vector<int> t_pattern;
  std::vector<int> t_run3pattern;
  std::vector<int> t_fractional_slope;
  std::vector<int> t_quartStripBit;
  std::vector<int> t_getEighthStripBit;
  std::vector<int> t_cfeb;
  std::vector<int> t_trknmb;
  std::vector<int> t_getKeyStrip;
  std::vector<int> t_fractional_strip;
  std::vector<int> t_full_bx;
  std::vector<int> t_comp_code;
  std::vector<int> t_collisionB;
  std::vector<int> t_accel;


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
    debug = iConfig.getParameter<bool>("debug");
    if(debug) std::cout<<"lctProducerData: "<<lctProducerData_<<std::endl;
    if(debug) std::cout<<"lctProducerEmul: "<<lctProducerEmul_<<std::endl;

    //consumes Data
    alcts_d_token_ = consumes<CSCALCTDigiCollection>(edm::InputTag(lctProducerData_, "MuonCSCALCTDigi"));
    clcts_d_token_ = consumes<CSCCLCTDigiCollection>(edm::InputTag(lctProducerData_, "MuonCSCCLCTDigi"));
    lcts_tmb_d_token_ = consumes<CSCCorrelatedLCTDigiCollection>(edm::InputTag(lctProducerData_, "MuonCSCCorrelatedLCTDigi"));

    //consumes Emul
    alcts_e_token_ = consumes<CSCALCTDigiCollection>(edm::InputTag(lctProducerEmul_));
    clcts_e_token_ = consumes<CSCCLCTDigiCollection>(edm::InputTag(lctProducerEmul_));
    lcts_tmb_e_token_ = consumes<CSCCorrelatedLCTDigiCollection>(edm::InputTag(lctProducerEmul_));

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
  t_RUN = RUN_;
  t_Event = Event_;
  t_eventsAnalyzed = eventsAnalyzed;
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
              if(debug) std::cout << "ALCT "  <<(*digiIt) << " data:emul "<< is_data <<":"<< is_emul << std::endl;
              alctV.push_back(*digiIt);
              //Fill TTree
              t_is_data.push_back(is_data);
              t_is_emul.push_back(is_emul);
              t_endcap.push_back(endc);
              t_station.push_back(stat);
              t_ring.push_back(ring);
              t_chamber.push_back(cham);
              t_isValid.push_back((*digiIt).isValid());
              t_quality.push_back((*digiIt).getQuality());
              t_accel.push_back((*digiIt).getAccelerator());
              t_collisionB.push_back((*digiIt).getCollisionB());
              t_keyWG.push_back((*digiIt).getKeyWG());
              t_bx.push_back((*digiIt).getBX());
              t_trknmb.push_back((*digiIt).getTrknmb());
            }
          }
        }
      }
    }
  }
  ALCT_tree->Fill();
  // Clear vectors
  t_is_data.clear();
  t_is_emul.clear();
  t_endcap.clear();
  t_station.clear();
  t_ring.clear();
  t_chamber.clear();
  t_isValid.clear();
  t_quality.clear();
  t_accel.clear();
  t_collisionB.clear();
  t_keyWG.clear();
  t_bx.clear();
  t_trknmb.clear();

}

void GEMCSCTriggerPrimitivesReader::SaveCLCTs(const CSCCLCTDigiCollection* clcts, bool is_data, bool is_emul)
{
  t_RUN = RUN_;
  t_Event = Event_;
  t_eventsAnalyzed = eventsAnalyzed;
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
              t_is_data.push_back(is_data);
              t_is_emul.push_back(is_emul);
              t_endcap.push_back(endc);
              t_station.push_back(stat);
              t_ring.push_back(ring);
              t_chamber.push_back(cham);  
              t_isValid.push_back((*digiIt).isValid());
              t_quality.push_back((*digiIt).getQuality());
              t_pattern.push_back((*digiIt).getPattern());
              t_run3pattern.push_back((*digiIt).getRun3Pattern());
              t_slope.push_back((*digiIt).getSlope());
              t_fractional_slope.push_back((*digiIt).getFractionalSlope());
              t_stripType.push_back((*digiIt).getStripType());
              t_bend.push_back((*digiIt).getBend());
              t_quartStripBit.push_back((*digiIt).getQuartStripBit());
              t_getEighthStripBit.push_back((*digiIt).getEighthStripBit());
              t_cfeb.push_back((*digiIt).getCFEB());
              t_trknmb.push_back((*digiIt).getTrknmb());
              t_getKeyStrip.push_back((*digiIt).getKeyStrip());
              t_fractional_strip.push_back((*digiIt).getFractionalStrip());
              t_full_bx.push_back((*digiIt).getFullBX());
              t_comp_code.push_back((*digiIt).getCompCode());
              t_bx.push_back((*digiIt).getBX());
              t_strip.push_back((*digiIt).getStrip());
            }
          }
        }
      }
    }
  }
  CLCT_tree->Fill();
  // Clear vectors
  t_is_data.clear();
  t_is_emul.clear();
  t_endcap.clear();
  t_station.clear();
  t_ring.clear();
  t_chamber.clear();
  t_isValid.clear();
  t_quality.clear();
  t_pattern.clear();
  t_run3pattern.clear();
  t_slope.clear();
  t_fractional_slope.clear();
  t_stripType.clear();
  t_bend.clear();
  t_quartStripBit.clear();
  t_getEighthStripBit.clear();
  t_cfeb.clear();
  t_trknmb.clear();
  t_getKeyStrip.clear();
  t_fractional_strip.clear();
  t_full_bx.clear();
  t_comp_code.clear();
  t_bx.clear();
  t_strip.clear();
}

void GEMCSCTriggerPrimitivesReader::SaveLCTs(const CSCCorrelatedLCTDigiCollection* lcts, bool is_data, bool is_emul)
{
  t_RUN = RUN_;
  t_Event = Event_;
  t_eventsAnalyzed = eventsAnalyzed;
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
              t_is_data.push_back(is_data);
              t_is_emul.push_back(is_emul);
              t_endcap.push_back(endc);
              t_station.push_back(stat);
              t_ring.push_back(ring);
              t_chamber.push_back(cham);
              t_isValid.push_back((*digiIt).isValid());
              t_quality.push_back((*digiIt).getQuality());
              t_keyWG.push_back((*digiIt).getKeyWG());
              t_strip.push_back((*digiIt).getStrip());
              t_quadStrip.push_back((*digiIt).getStrip(4));
              t_eightStrip.push_back((*digiIt).getStrip(8));
              t_stripType.push_back((*digiIt).getStripType());
              t_bend.push_back((*digiIt).getBend());
              t_slope.push_back((*digiIt).getSlope());
              t_bx.push_back((*digiIt).getBX());
              t_pattern.push_back((*digiIt).getPattern());
              t_run3pattern.push_back((*digiIt).getRun3Pattern());
              
            }
          }
        }
      }
    }
  }
  LCT_tree->Fill();
  // Clear vectors
  t_is_data.clear();
  t_is_emul.clear();
  t_endcap.clear();
  t_station.clear();
  t_ring.clear();
  t_chamber.clear();
  t_isValid.clear();
  t_quality.clear();
  t_keyWG.clear();
  t_strip.clear();
  t_quadStrip.clear();
  t_eightStrip.clear();
  t_stripType.clear();
  t_bend.clear();
  t_slope.clear();
  t_bx.clear();
  t_pattern.clear();
  t_run3pattern.clear();
}

TTree* GEMCSCTriggerPrimitivesReader::bookTTree_ALCT() {
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("ALCT_tree", "GEMCSCTriggerPrimitivesReader");
  t->Branch("RUN", &t_RUN, "RUN/I");
  t->Branch("Event", &t_Event, "Event/I");
  t->Branch("eventsAnalyzed", &t_eventsAnalyzed, "eventsAnalyzed/I");
  t->Branch("is_data", "std::vector<bool>", &t_is_data);
  t->Branch("is_emul", "std::vector<bool>", &t_is_emul);
  t->Branch("endcap", "std::vector<int>", &t_endcap);
  t->Branch("station", "std::vector<int>", &t_station);
  t->Branch("ring", "std::vector<int>", &t_ring);
  t->Branch("chamber", "std::vector<int>", &t_chamber);
  t->Branch("isValid", "std::vector<bool>", &t_isValid);
  t->Branch("quality", "std::vector<int>", &t_quality);
  t->Branch("keyWG", "std::vector<int>", &t_keyWG);
  t->Branch("bx", "std::vector<int>", &t_bx);
  t->Branch("trknmb", "std::vector<int>", &t_trknmb);
  t->Branch("accel", "std::vector<int>", &t_accel);
  t->Branch("collisionB", "std::vector<int>", &t_collisionB);
  return t;
}
TTree* GEMCSCTriggerPrimitivesReader::bookTTree_CLCT(){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("CLCT_tree", "GEMCSCTriggerPrimitivesReader");
  t->Branch("RUN", &t_RUN, "RUN/I");
  t->Branch("Event", &t_Event, "Event/I");
  t->Branch("eventsAnalyzed", &t_eventsAnalyzed, "eventsAnalyzed/I");
  t->Branch("is_data", "std::vector<bool>", &t_is_data);
  t->Branch("is_emul", "std::vector<bool>", &t_is_emul);
  t->Branch("endcap", "std::vector<int>", &t_endcap);
  t->Branch("station", "std::vector<int>", &t_station);
  t->Branch("ring", "std::vector<int>", &t_ring);
  t->Branch("chamber", "std::vector<int>", &t_chamber);
  t->Branch("isValid", "std::vector<bool>", &t_isValid);
  t->Branch("quality", "std::vector<int>", &t_quality);
  t->Branch("pattern", "std::vector<int>", &t_pattern);
  t->Branch("run3pattern", "std::vector<int>", &t_run3pattern);
  t->Branch("slope", "std::vector<int>", &t_slope);
  t->Branch("fractional_slope", "std::vector<int>", &t_fractional_slope);
  t->Branch("stripType", "std::vector<int>", &t_stripType);
  t->Branch("bend", "std::vector<int>", &t_bend);
  t->Branch("quartStripBit", "std::vector<int>", &t_quartStripBit);
  t->Branch("getEighthStripBit", "std::vector<int>", &t_getEighthStripBit);
  t->Branch("cfeb", "std::vector<int>", &t_cfeb);
  t->Branch("trknmb", "std::vector<int>", &t_trknmb);
  t->Branch("getKeyStrip", "std::vector<int>", &t_getKeyStrip);
  t->Branch("fractional_strip", "std::vector<int>", &t_fractional_strip);
  t->Branch("full_bx", "std::vector<int>", &t_full_bx);
  t->Branch("comp_code", "std::vector<int>", &t_comp_code);
  t->Branch("bx", "std::vector<int>", &t_bx);
  t->Branch("strip", "std::vector<int>", &t_strip);

  return t;

}
TTree* GEMCSCTriggerPrimitivesReader::bookTTree_LCT(){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("LCT_tree", "GEMCSCTriggerPrimitivesReader");
  t->Branch("RUN", &t_RUN, "RUN/I");
  t->Branch("Event", &t_Event, "Event/I");
  t->Branch("eventsAnalyzed", &t_eventsAnalyzed, "eventsAnalyzed/I");
  t->Branch("is_data", "std::vector<bool>", &t_is_data);
  t->Branch("is_emul", "std::vector<bool>", &t_is_emul);
  t->Branch("endcap", "std::vector<int>", &t_endcap);
  t->Branch("station", "std::vector<int>", &t_station);
  t->Branch("ring", "std::vector<int>", &t_ring);
  t->Branch("chamber", "std::vector<int>", &t_chamber);
  t->Branch("isValid", "std::vector<bool>", &t_isValid);
  t->Branch("quality", "std::vector<int>", &t_quality);
  t->Branch("keyWG", "std::vector<int>", &t_keyWG);
  t->Branch("strip", "std::vector<int>", &t_strip);
  t->Branch("quadStrip", "std::vector<int>", &t_quadStrip);
  t->Branch("eightStrip", "std::vector<int>", &t_eightStrip);
  t->Branch("stripType", "std::vector<int>", &t_stripType);
  t->Branch("bend", "std::vector<int>", &t_bend);
  t->Branch("slope", "std::vector<int>", &t_slope);
  t->Branch("bx", "std::vector<int>", &t_bx);
  t->Branch("pattern", "std::vector<int>", &t_pattern);
  t->Branch("run3pattern", "std::vector<int>", &t_run3pattern);
  return t;
}

DEFINE_FWK_MODULE(GEMCSCTriggerPrimitivesReader);