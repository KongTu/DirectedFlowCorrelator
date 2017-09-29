// -*- C++ -*-
//
// Package:    DirectedFlowCorrelator/DirectedFlowCorrelator
// Class:      DirectedFlowCorrelator
// 
/**\class DirectedFlowCorrelator DirectedFlowCorrelator.cc DirectedFlowCorrelator/DirectedFlowCorrelator/plugins/DirectedFlowCorrelator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhoudunming Tu
//         Created:  Wed, 20 Sep 2017 22:44:00 GMT
//
//


#include "DirectedFlowCorrelator/DirectedFlowCorrelator/interface/DirectedFlowCorrelatorBase.h"


DirectedFlowCorrelator::DirectedFlowCorrelator(const edm::ParameterSet& iConfig)

{
  trackName_  =  iConfig.getParameter<edm::InputTag>("trackName");
  vertexName_ =  iConfig.getParameter<edm::InputTag>("vertexName");
  towerName_ =  iConfig.getParameter<edm::InputTag>("towerName");

  trackSrc_ = consumes<reco::TrackCollection>(trackName_);
  vertexSrc_ = consumes<reco::VertexCollection>(vertexName_);
  towerSrc_ = consumes<CaloTowerCollection>(towerName_);

  Nmin_ = iConfig.getUntrackedParameter<int>("Nmin");
  Nmax_ = iConfig.getUntrackedParameter<int>("Nmax");

  useCentrality_ = iConfig.getUntrackedParameter<bool>("useCentrality");
  doEffCorrection_ = iConfig.getUntrackedParameter<bool>("doEffCorrection");
  useEtaGap_ = iConfig.getUntrackedParameter<bool>("useEtaGap");
  doBothSide_ = iConfig.getUntrackedParameter<bool>("doBothSide");

  eff_ = iConfig.getUntrackedParameter<int>("eff");

  etaTracker_ = iConfig.getUntrackedParameter<double>("etaTracker");
  gapValue_ = iConfig.getUntrackedParameter<double>("gapValue");
  
  etaLowHF_ = iConfig.getUntrackedParameter<double>("etaLowHF");
  etaHighHF_ = iConfig.getUntrackedParameter<double>("etaHighHF");

  vzLow_ = iConfig.getUntrackedParameter<double>("vzLow");
  vzHigh_ = iConfig.getUntrackedParameter<double>("vzHigh");
  
  ptLow_ = iConfig.getUntrackedParameter<double>("ptLow");
  ptHigh_ = iConfig.getUntrackedParameter<double>("ptHigh");

  offlineptErr_ = iConfig.getUntrackedParameter<double>("offlineptErr", 0.0);
  offlineDCA_ = iConfig.getUntrackedParameter<double>("offlineDCA", 0.0);
  offlineChi2_ = iConfig.getUntrackedParameter<double>("offlineChi2", 0.0);
  offlinenhits_ = iConfig.getUntrackedParameter<double>("offlinenhits", 0.0);

  etaBins_ = iConfig.getUntrackedParameter<std::vector<double>>("etaBins");
  ptBins_ = iConfig.getUntrackedParameter<std::vector<double>>("ptBins");
  centBins_ = iConfig.getUntrackedParameter<std::vector<double>>("centBins");

}


DirectedFlowCorrelator::~DirectedFlowCorrelator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DirectedFlowCorrelator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexSrc_,vertices);
  double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
  double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
  const reco::Vertex & vtx = (*vertices)[0];
  bestvz = vtx.z(); 
  bestvx = vtx.x(); 
  bestvy = vtx.y();
  bestvzError = vtx.zError(); 
  bestvxError = vtx.xError(); 
  bestvyError = vtx.yError();

  //first selection; vertices
  if( fabs(bestvz) < vzLow_ || fabs(bestvz) > vzHigh_ ) return;

  vtxZ->Fill( bestvz );

  Handle<CaloTowerCollection> towers;
  iEvent.getByToken(towerSrc_, towers);

  Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackSrc_, tracks);

  int nTracks = 0;
  for(unsigned it = 0; it < tracks->size(); it++){

     const reco::Track & trk = (*tracks)[it];
  
     math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);

        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt() > 0.1 ) continue;
        if(fabs(dzvtx/dzerror) > 3.0) continue;
        if(fabs(dxyvtx/dxyerror) > 3.0) continue;
        if(trk.pt() < 0.4 || fabs(trk.eta()) > 2.4) continue;
        nTracks++;//count multiplicity

  }

  if( !useCentrality_ ) if( nTracks < Nmin_ || nTracks >= Nmax_ ) return;

  double etHFtowerSumPlus = 0.0;
  double etHFtowerSumMinus = 0.0;
  double etHFtowerSum = 0.0;
  
  if( useCentrality_ ){

    for( unsigned i = 0; i<towers->size(); ++ i){
       const CaloTower & tower = (*towers)[ i ];
       double eta = tower.eta();
       bool isHF = tower.ietaAbs() > 29;
          if(isHF && eta > 0){
            etHFtowerSumPlus += tower.pt();
          }
          if(isHF && eta < 0){
            etHFtowerSumMinus += tower.pt();
          }
    }
    etHFtowerSum=etHFtowerSumPlus + etHFtowerSumMinus;

    int bin = -1;
    for(int j=0; j<200; j++){
      if( etHFtowerSum >= centBins_[j] ){
         bin = j; break;
      }
    }

    int hiBin = bin;
    if( hiBin < Nmin_ || hiBin >= Nmax_ ) return;
    cbinHist->Fill( hiBin );

  }

  Ntrk->Fill( nTracks );

  const int NetaBins = etaBins_.size() - 1 ;

  TComplex  Q_n3_1_HFplus, Q_n3_1_HFminus, Q_0_1_HFplus, Q_0_1_HFminus;
  //HF towers loop to fill the towers' Q-vectors:
  TComplex Q_n3_trk_minus, Q_0_trk_minus, Q_n3_trk_plus, Q_0_trk_plus;
  //charge independent, |eta|<1.0
  TComplex Q_n1_1[NetaBins][2], Q_0_1[NetaBins][2];

  for(unsigned i = 0; i < towers->size(); ++i){

          const CaloTower & hit= (*towers)[i];

          double caloEta = hit.eta();
          double caloPhi = hit.phi();
          double w = hit.hadEt( vtx.z() ) + hit.emEt( vtx.z() );

          hfPhi->Fill(caloPhi, w);
          
          if( caloEta < etaHighHF_ && caloEta > etaLowHF_ ){
            
              Q_n3_1_HFplus += q_vector(-1, 1, w, caloPhi);
              Q_0_1_HFplus += q_vector(0, 1, w, caloPhi);

          }
          else if( caloEta < -etaLowHF_ && caloEta > -etaHighHF_ ){

              Q_n3_1_HFminus += q_vector(-1, 1, -w, caloPhi);
              Q_0_1_HFminus += q_vector(0, 1, w, caloPhi); //normalization needs to be positive in order to be not cancel out

          }
          else{continue;}
  }

  //track loop to fill charged particles Q-vectors
  for(unsigned it = 0; it < tracks->size(); it++){

    const reco::Track & trk = (*tracks)[it];

    math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

    double dzvtx = trk.dz(bestvtx);
    double dxyvtx = trk.dxy(bestvtx);
    double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
    double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError); 
    double nhits = trk.numberOfValidHits();
    double chi2n = trk.normalizedChi2();
    double nlayers = trk.hitPattern().trackerLayersWithMeasurement();
    chi2n = chi2n/nlayers;
    double nPixelLayers = trk.hitPattern().pixelLayersWithMeasurement();//only pixel layers
    double phi = trk.phi();
    double trkEta = trk.eta();

    double weight = 1.0;
    if( doEffCorrection_ ) { 
      weight = 1.0/effTable[eff_]->GetBinContent( effTable[eff_]->FindBin(trk.eta(), trk.pt()) );
    }

    if(!trk.quality(reco::TrackBase::highPurity)) continue;
    if(fabs(trk.ptError())/trk.pt() > offlineptErr_ ) continue;
    if(fabs(dzvtx/dzerror) > offlineDCA_) continue;
    if(fabs(dxyvtx/dxyerror) > offlineDCA_) continue;
    if(chi2n > offlineChi2_ ) continue;
    if(nhits < offlinenhits_ ) continue;
    if( nPixelLayers <= 0 ) continue;
    if(trk.pt() < ptLow_ || trk.pt() > ptHigh_ ) continue;
    if(fabs(trkEta) > etaTracker_ ) continue;

    trkPhi->Fill(phi, weight);
    trkPt->Fill(trk.pt(), weight);
    trk_eta->Fill(trkEta, weight);

    if( trkEta < -2.0 && trkEta > -2.4 ){
   
      Q_n3_trk_minus += q_vector(+1, 1, -weight, phi);//for scalar product in tracker
      Q_0_trk_minus += q_vector(0, 1, weight, phi);
    }
    if( trkEta < 2.4 && trkEta > 2.0 ){
   
      Q_n3_trk_plus += q_vector(+1, 1, weight, phi);//for scalar product in tracker
      Q_0_trk_plus += q_vector(0, 1, weight, phi);
    }
    
    for(int eta = 0; eta < NetaBins; eta++){
      if( trkEta > etaBins_[eta] && trkEta < etaBins_[eta+1] ){

        if( trk.charge() == +1 ){//positive charge

          //3p:
          Q_n1_1[eta][0] += q_vector(+1, 1, weight, phi);
          Q_0_1[eta][0] += q_vector(0, 1, weight, phi);

        }
        if( trk.charge() == -1 ){//negative charge

          Q_n1_1[eta][1] += q_vector(+1, 1, weight, phi);
          Q_0_1[eta][1] += q_vector(0, 1, weight, phi);
        }
      }
    }//end of eta dimension
  }

/*
event average v1
*/

//resolution factor
  TComplex N_2_trk, D_2_trk;

  N_2_trk = Q_n3_trk_plus*Q_n3_1_HFplus;
  D_2_trk = Q_0_trk_plus*Q_0_1_HFplus;

  c2_cb_plus->Fill( N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re());

  N_2_trk = Q_n3_trk_plus*Q_n3_1_HFminus;
  D_2_trk = Q_0_trk_plus*Q_0_1_HFminus;

  c2_ac_plus->Fill( N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re() );

  N_2_trk = Q_n3_trk_minus*Q_n3_1_HFplus;
  D_2_trk = Q_0_trk_minus*Q_0_1_HFplus;

  c2_cb_minus->Fill( N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re());

  N_2_trk = Q_n3_trk_minus*Q_n3_1_HFminus;
  D_2_trk = Q_0_trk_minus*Q_0_1_HFminus;

  c2_ac_minus->Fill( N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re() );

  N_2_trk = Q_n3_1_HFplus*TComplex::Conjugate(Q_n3_1_HFminus);
  D_2_trk = Q_0_1_HFplus*Q_0_1_HFminus;

  c2_ab->Fill( N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re() );

  c2_a_real->Fill( Q_n3_1_HFminus.Re()/Q_0_1_HFminus.Re(), Q_0_1_HFminus.Re() );
  c2_b_real->Fill( Q_n3_1_HFplus.Re()/Q_0_1_HFplus.Re(), Q_0_1_HFplus.Re() );
  c2_c_real->Fill( Q_n3_trk_plus.Re()/Q_0_trk_plus.Re(), Q_0_trk_plus.Re() );

  c2_a_imag->Fill( Q_n3_1_HFminus.Im()/Q_0_1_HFminus.Re(), Q_0_1_HFminus.Re() );
  c2_b_imag->Fill( Q_n3_1_HFplus.Im()/Q_0_1_HFplus.Re(), Q_0_1_HFplus.Re() );
  c2_c_imag->Fill( Q_n3_trk_plus.Im()/Q_0_trk_plus.Re(), Q_0_trk_plus.Re() );


//numerator
  for(int eta = 0; eta < NetaBins; eta++){
    for(int charge = 0; charge < 2; charge++){

      TComplex N_v1_A_SP, D_v1_A_SP, N_v1_B_SP, D_v1_B_SP;

      N_v1_A_SP = Q_n1_1[eta][charge]*Q_n3_1_HFminus;
      D_v1_A_SP = Q_0_1[eta][charge]*Q_0_1_HFminus;

      double V1_A = N_v1_A_SP.Re()/D_v1_A_SP.Re();

      N_v1_B_SP = Q_n1_1[eta][charge]*Q_n3_1_HFplus;
      D_v1_B_SP = Q_0_1[eta][charge]*Q_0_1_HFplus;

      double V1_B = N_v1_B_SP.Re()/D_v1_B_SP.Re();

      c2_v1[eta][charge][0]->Fill( V1_A, D_v1_A_SP.Re() );
      c2_v1[eta][charge][1]->Fill( V1_B, D_v1_B_SP.Re() );

      c2_trk_accept[eta][charge]->Fill(Q_n1_1[eta][charge].Re()/Q_0_1[eta][charge].Re(), Q_0_1[eta][charge].Re());

    }
  }

//
// /*
// EbyE analysis, in progress.
// */
//   TComplex N_1_SP, D_1_SP, N_2_SP, D_2_SP, N_3_SP, D_3_SP;

//   N_1_SP = Q_n3_1_HFminus*Q_n3_trk_plus;
//   D_1_SP = Q_0_1_HFminus*Q_0_trk_plus;

//   N_2_SP = Q_n3_1_HFplus*Q_n3_trk_plus;
//   D_2_SP = Q_0_1_HFplus*Q_0_trk_plus;

//   N_3_SP = Q_n3_1_HFplus*TComplex::Conjugate(Q_n3_1_HFminus);
//   D_3_SP = Q_0_1_HFplus*Q_0_1_HFminus;

//   double t1 = N_1_SP.Re()/D_1_SP.Re();
//   double t2 = N_2_SP.Re()/D_2_SP.Re();
//   double t3 = N_3_SP.Re()/D_3_SP.Re();

//   double Res_A = sqrt(t1*t3/t2);
//   double Res_B = sqrt(t2*t3/t1);

//   double v1[NetaBins][2];//EbyE v1 in different eta slices and charge

//   for(int eta = 0; eta < NetaBins; eta++){
//     for(int charge = 0; charge < 2; charge++){

//       TComplex N_v1_A_SP, D_v1_A_SP, N_v1_B_SP, D_v1_B_SP;

//       N_v1_A_SP = Q_n1_1[eta][charge]*Q_n3_1_HFminus;
//       D_v1_A_SP = Q_0_1[eta][charge]*Q_0_1_HFminus;

//       double V1_A = N_v1_A_SP.Re()/D_v1_A_SP.Re();

//       N_v1_B_SP = Q_n1_1[eta][charge]*Q_n3_1_HFplus;
//       D_v1_B_SP = Q_0_1[eta][charge]*Q_0_1_HFplus;

//       double V1_B = N_v1_B_SP.Re()/D_v1_B_SP.Re();

//       v1[eta][charge] = (V1_A/Res_A + V1_B/Res_B)/2.0;
//     }
//   }

//   //EbyE v1 is calculated for different eta and charge
//   //Now calculate the A correlator. 

//   for(int eta = 0; eta < NetaBins/2; eta++){

//     double A_1_pm_YY = v1[eta][0] - v1[eta][1];
//     double A_1_pm_YmY = v1[eta][0] - v1[NetaBins-eta-1][1];

//     double A_1_pp_YmY = v1[eta][0] - v1[NetaBins-eta-1][0];
//     double A_1_mm_YmY = v1[eta][1] - v1[NetaBins-eta-1][1];

//     double C_1 = A_1_pm_YY*A_1_pm_YY;
//     C_1_YY[eta]->Fill(C_1, Q_0_1[eta][0]);

//     double C_1_1 = A_1_pm_YmY*A_1_pm_YmY;
//     C_1_YmY[eta]->Fill(C_1_1, Q_0_1[eta][0]);

//     double C_2 = (v1[eta][0]+v1[NetaBins-eta-1][1])*(v1[eta][0]+v1[NetaBins-eta-1][1]);
//     C_2_YmY[eta]->Fill(C_2, Q_0_1[eta][0]);

//     double C_3 = A_1_pp_YmY*A_1_mm_YmY;
//     C_3_YmY[eta]->Fill(C_3, Q_0_1[eta][0]);

//   }  

}
// ------------ method called once each job just before starting event loop  ------------
void 
DirectedFlowCorrelator::beginJob()
{
  edm::Service<TFileService> fs;
    
  TH1D::SetDefaultSumw2();

  const int NetaBins = etaBins_.size() - 1 ;

  double ptBinsArray[100];
  const int Nptbins = ptBins_.size() - 1;
  for(unsigned i = 0; i < ptBins_.size(); i++){
    ptBinsArray[i] = ptBins_[i];
  }

  edm::FileInPath fip1("DirectedFlowCorrelator/DirectedFlowCorrelator/data/Hydjet_eff_mult_v1.root");
  TFile f1(fip1.fullPath().c_str(),"READ");
  for(int i = 0; i < 5; i++){
     effTable[i] = (TH2D*)f1.Get(Form("rTotalEff3D_%d",i));
  }

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",5000,0,5000);
  vtxZ = fs->make<TH1D>("vtxZ",";vz", 400,-20,20);
  cbinHist = fs->make<TH1D>("cbinHist",";cbin",200,0,200);
  trkPhi = fs->make<TH1D>("trkPhi", ";#phi", 700, -3.5, 3.5);
  hfPhi = fs->make<TH1D>("hfPhi", ";#phi", 700, -3.5, 3.5);
  trkPt = fs->make<TH1D>("trkPt", ";p_{T}(GeV)", Nptbins,ptBinsArray);
  trk_eta = fs->make<TH1D>("trk_eta", ";#eta", 50,-2.5,2.5);

  // for(int eta = 0; eta < NetaBins/2; eta++){

  //   C_1_YY[eta] = fs->make<TH1D>(Form("C_1_YY_%d",eta),";C_1_YY", 100,-1,1);
  //   C_1_YmY[eta] = fs->make<TH1D>(Form("C_1_YmY_%d",eta),";C_1_YmY", 100,-1,1);
  //   C_2_YmY[eta] = fs->make<TH1D>(Form("C_2_YmY_%d",eta),";C_2_YmY", 100,-1,1);
  //   C_3_YmY[eta] = fs->make<TH1D>(Form("C_3_YmY_%d",eta),";C_3_YmY", 100,-1,1);
  // }

  for(int eta = 0; eta < NetaBins; eta++){
    for(int charge = 0; charge < 2; charge++){

      c2_trk_accept[eta][charge] = fs->make<TH1D>(Form("c2_trk_accept_%d_%d",eta,charge), ";c1", 1,-1,1);
      for(int dir = 0; dir < 2; dir++){

        c2_v1[eta][charge][dir] = fs->make<TH1D>(Form("c2_v1_%d_%d_%d",eta,charge,dir),";c1", 1,-1,1);
      }
    }
  }


  c2_ab = fs->make<TH1D>("c2_ab",";c2_ab", 1,-1,1);
  
  c2_ac_plus = fs->make<TH1D>("c2_ac_plus",";c2_ac_plus", 1,-1,1);
  c2_cb_plus = fs->make<TH1D>("c2_cb_plus",";c2_cb_plus", 1,-1,1);

  c2_ac_minus = fs->make<TH1D>("c2_ac_minus",";c2_ac_minus", 1,-1,1);
  c2_cb_minus = fs->make<TH1D>("c2_cb_minus",";c2_cb_minus", 1,-1,1);

  c2_a_real = fs->make<TH1D>("c2_a_real",";c2_a_real", 1,-1,1);
  c2_b_real = fs->make<TH1D>("c2_b_real",";c2_b_real", 1,-1,1);
  c2_c_real = fs->make<TH1D>("c2_c_real",";c2_c_real", 1,-1,1);

  c2_a_imag = fs->make<TH1D>("c2_a_imag",";c2_a_imag", 1,-1,1);
  c2_b_imag = fs->make<TH1D>("c2_b_imag",";c2_b_imag", 1,-1,1);
  c2_c_imag = fs->make<TH1D>("c2_c_imag",";c2_c_imag", 1,-1,1);

}

TComplex 
DirectedFlowCorrelator::q_vector(double n, double p, double w, double phi) 
{
  double term1 = pow(w,p);
  TComplex e(1, n*phi, 1);
  return term1*e;
}
// ------------ method called once each job just after ending the event loop  ------------
void 
DirectedFlowCorrelator::endJob() 
{
}
void 
DirectedFlowCorrelator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DirectedFlowCorrelator::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DirectedFlowCorrelator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DirectedFlowCorrelator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DirectedFlowCorrelator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DirectedFlowCorrelator);
