#ifndef D0DirectedFlowCorrelatorBase_
#define D0DirectedFlowCorrelatorBase_


#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <map>
#include <sstream>


#include <TMath.h>
#include <TH1D.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TNtuple.h>
#include <TGraph.h>
#include <TComplex.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TF1.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//////////////////////////////////////////////
// CMSSW user include files
#include "DataFormats/Common/interface/DetSetAlgorithm.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
// Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

//

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

// Vertex significance
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

// Root include files
#include "TTree.h"
//
// Track Matching and fake rate calculations     
//#include "RiceHIG/V0Analysis/interface/V0Validator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#define PI 3.1415926
#define K0sMass 0.497614
#define LambdaMass 1.115683

using namespace std;
using namespace reco;
using namespace edm;


class D0DirectedFlowCorrelator : public edm::EDAnalyzer {
   public:
      explicit D0DirectedFlowCorrelator(const edm::ParameterSet&);
      ~D0DirectedFlowCorrelator();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual TComplex q_vector(double n, double p, double w, double phi);

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
      edm::EDGetTokenT<reco::TrackCollection> trackSrc_;
      edm::EDGetTokenT<CaloTowerCollection> towerSrc_;
      edm::EDGetTokenT<reco::GenParticleCollection> genSrc_;
      
      edm::EDGetTokenT<reco::Centrality> centralityToken_;
      edm::EDGetTokenT<int> centralityBinToken_;

      edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> recoVertexCompositeCandidateCollection_Token_;
      edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token1_;
      edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token2_;

      edm::InputTag vertexName_;
      edm::InputTag trackName_;
      edm::InputTag towerName_;
      edm::InputTag genName_;

      //correction table
      TH2D* effTable[5];

      TH1D* Ntrk;
      TH1D* vtxZ;
      TH1D* trkPhi;
      TH1D* hfPhi;
      TH1D* trkPt;
      TH1D* trk_eta;
      TH1D* cbinHist;

      TH1D* c2_cb_plus;
      TH1D* c2_ac_plus;
      TH1D* c2_cb_minus;
      TH1D* c2_ac_minus;
      
      TH1D* c2_ab;

      TH1D* c2_a_real;
      TH1D* c2_b_real;
      TH1D* c2_c_plus_real;
      TH1D* c2_c_minus_real;

      TH1D* c2_a_imag;
      TH1D* c2_b_imag;
      TH1D* c2_c_plus_imag;
      TH1D* c2_c_minus_imag;

      TH1D* c2_v1[20][2][2];
      TH1D* c2_trk_accept[20][2][2];

      TH1D* c2_d0obs_v1[20][3][2];
      TH1D* c2_d0obs_trk_accept[20][3][2];

      TH1D* c2_d0bkg_v1[20][3][2];
      TH1D* c2_d0bkg_trk_accept[20][3][2];

      TH1D* D0Mass_Hist[20][3];

      TH1D* C_1_YY[20];
      TH1D* C_1_YmY[20];

      TH1D* C_2_YmY[20];
      TH1D* C_3_YmY[20];

      int Nmin_;
      int Nmax_;

      int eff_;

      double etaTracker_;
      double gapValue_;
      double etaLowHF_;
      double etaHighHF_;
      
      double vzLow_;
      double vzHigh_;
      
      double ptLow_;
      double ptHigh_;

      double offlineptErr_;
      double offlineDCA_;
      double offlineChi2_;
      double offlinenhits_;

      bool useCentrality_;
      bool doEffCorrection_;
      bool useEtaGap_;
      bool doBothSide_;
      bool doPixelReco_;
      bool doHiReco_;

      std::vector<double> etaBins_;
      std::vector<double> rapidityBins_;
      std::vector<double> ptBins_;
      std::vector<double> centBins_;

      double D0MassHigh_;
      double D0MassLow_;
      double D0EtaHigh_;
      double D0EtaLow_;
      double D0PtHigh_;
      double D0PtLow_;
      double D0YHigh_;
      double D0YLow_;
      double D0DcaHigh_;
      double D0DcaLow_;
      double D0VtxProbHigh_;
      double D0VtxProbLow_;
      double D03DAngleHigh_;
      double D03DAngleLow_;
      double D0DlosHigh_;
      double D0DlosLow_;

      double TrkPtHigh_;
      double TrkPtLow_;
      double TrkEtaHigh_;
      double TrkEtaLow_;
      double TrkChiOverNLayerHigh_;
      double TrkChiOverNLayerLow_;
      double TrkPtErrOverPtHigh_;
      double TrkPtErrOverPtLow_;
      int TrkNHitHigh_;
      int TrkNHitLow_;

};



#endif