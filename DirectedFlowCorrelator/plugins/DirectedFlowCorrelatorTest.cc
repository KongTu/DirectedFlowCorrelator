// -*- C++ -*-
//
// Package:    DirectedFlowCorrelatorTest/DirectedFlowCorrelatorTest
// Class:      DirectedFlowCorrelatorTest
// 
/**\class DirectedFlowCorrelatorTest DirectedFlowCorrelatorTest.cc DirectedFlowCorrelatorTest/DirectedFlowCorrelatorTest/plugins/DirectedFlowCorrelatorTest.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhoudunming Tu
//         Created:  Wed, 20 Sep 2017 22:44:00 GMT
//
//


#include "DirectedFlowCorrelator/DirectedFlowCorrelator/interface/DirectedFlowCorrelatorTestBase.h"


DirectedFlowCorrelatorTest::DirectedFlowCorrelatorTest(const edm::ParameterSet& iConfig)

{
  trackName_  =  iConfig.getParameter<edm::InputTag>("trackName");
  vertexName_ =  iConfig.getParameter<edm::InputTag>("vertexName");
  towerName_ =  iConfig.getParameter<edm::InputTag>("towerName");

  trackSrc_ = consumes<reco::TrackCollection>(trackName_);
  vertexSrc_ = consumes<reco::VertexCollection>(vertexName_);
  towerSrc_ = consumes<CaloTowerCollection>(towerName_);

  centralityToken_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
  centralityBinToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinSrc"));

  Nmin_ = iConfig.getUntrackedParameter<int>("Nmin");
  Nmax_ = iConfig.getUntrackedParameter<int>("Nmax");

  useCentrality_ = iConfig.getUntrackedParameter<bool>("useCentrality");
  doEffCorrection_ = iConfig.getUntrackedParameter<bool>("doEffCorrection");
  useEtaGap_ = iConfig.getUntrackedParameter<bool>("useEtaGap");
  doBothSide_ = iConfig.getUntrackedParameter<bool>("doBothSide");
  doPixelReco_ = iConfig.getUntrackedParameter<bool>("doPixelReco");
  
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

  phiCos_ = iConfig.getUntrackedParameter<std::vector<double>>("phiCos");
  phiSin_ = iConfig.getUntrackedParameter<std::vector<double>>("phiSin");

  etaBins_ = iConfig.getUntrackedParameter<std::vector<double>>("etaBins");
  ptBins_ = iConfig.getUntrackedParameter<std::vector<double>>("ptBins");
  centBins_ = iConfig.getUntrackedParameter<std::vector<double>>("centBins");

}


DirectedFlowCorrelatorTest::~DirectedFlowCorrelatorTest()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DirectedFlowCorrelatorTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

    if( doPixelReco_ ){

      edm::Handle<reco::Centrality> centrality;
      iEvent.getByToken(centralityToken_, centrality);
      edm::Handle<int> cbin;
      iEvent.getByToken(centralityBinToken_, cbin);
      int hiBin = *cbin;

      if( hiBin < Nmin_ || hiBin >= Nmax_ ) return;
      cbinHist->Fill( hiBin );

    }
    else{

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

  }

  Ntrk->Fill( nTracks );

  const int NetaBins = etaBins_.size() - 1 ;

  //HF towers loop to fill the towers' Q-vectors:
  TComplex  Q_n3_1_HFplus, Q_n3_1_HFminus, Q_0_1_HFplus, Q_0_1_HFminus;
  
  //HF both side
  TComplex Q_n3_1_HFcombined, Q_0_1_HFcombined;
  
  //charge independent, |eta|<0.8
  TComplex Q_n3_trk_minus, Q_0_trk_minus, Q_n3_trk_plus, Q_0_trk_plus;

  //Psi_2 event plane resolution variables:
  TComplex Q_n1_Psi2_minus, Q_0_Psi2_minus, Q_n1_Psi2_plus, Q_0_Psi2_plus;
  
  TComplex Q_n1_1[NetaBins][2], Q_0_1[NetaBins][2];


  double HF_Psi_1_cosine = 0.0;
  double HF_Psi_1_sine = 0.0;

  for(unsigned i = 0; i < towers->size(); ++i){

          const CaloTower & hit= (*towers)[i];

          double caloEta = hit.eta();
          double caloPhi = hit.phi();
          double w = hit.hadEt( vtx.z() ) + hit.emEt( vtx.z() );

          hfPhi->Fill(caloPhi, w);
  
          if( caloEta < etaHighHF_ && caloEta > etaLowHF_ ){
            
              Q_n3_1_HFplus += q_vector(-1, 1, w, caloPhi);
              Q_0_1_HFplus += q_vector(0, 1, w, caloPhi);

              Q_n3_1_HFcombined += q_vector(-1, 1, w, caloPhi);
              Q_0_1_HFcombined += q_vector(0, 1, w, caloPhi);

              HF_Psi_1_sine += -w*sin( 1*caloPhi );
              HF_Psi_1_cosine += -w*cos( 1*caloPhi );

          }
          else if( caloEta < -etaLowHF_ && caloEta > -etaHighHF_ ){

              Q_n3_1_HFminus += q_vector(-1, 1, -w, caloPhi);
              Q_0_1_HFminus += q_vector(0, 1, w, caloPhi); //normalization needs to be positive in order to be not cancel out

              Q_n3_1_HFcombined += q_vector(-1, 1, -w, caloPhi);
              Q_0_1_HFcombined += q_vector(0, 1, w, caloPhi);

              HF_Psi_1_sine += w*sin( 1*caloPhi );
              HF_Psi_1_cosine += w*cos( 1*caloPhi );

          }
          else{continue;}
  }

  //Psi_1 in HF:
  HF_Psi_1_sine = HF_Psi_1_sine - (-0.1744);
  HF_Psi_1_cosine = HF_Psi_1_cosine - 0.4749;
  double Psi_1 = TMath::ATan(HF_Psi_1_sine/HF_Psi_1_cosine)/1;
  Psi_1_cos->Fill(HF_Psi_1_cosine);
  Psi_1_sin->Fill(HF_Psi_1_sine);

  double TRK_Psi_2_sine = 0.0;
  double TRK_Psi_2_cosine = 0.0;
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
    
    if( !doPixelReco_){ if(fabs(dxyvtx/dxyerror) > offlineDCA_) continue; }
    if(chi2n > offlineChi2_ ) continue;

    if( doPixelReco_ ){ if(nhits != 3 && nhits != 4 && nhits != 5 && nhits != 6) continue;}
    else{ if(nhits < offlinenhits_ ) continue; if( nPixelLayers <= 0 ) continue;}
    
    if(trk.pt() < ptLow_ || trk.pt() > ptHigh_ ) continue;
    if(fabs(trkEta) > etaTracker_ ) continue;

    trkPhi->Fill(phi, weight);
    trkPt->Fill(trk.pt(), weight);
    trk_eta->Fill(trkEta, weight);

    //For Psi_2 plane resoution:
    if( trkEta < 0.0 && trkEta > -1.0 ){

      Q_n1_Psi2_minus += q_vector(+2, 1, weight, phi);
      Q_0_Psi2_minus += q_vector(0, 1, weight, phi); 

    }
    if( trkEta < 1.0 && trkEta > 0.0 ){

      Q_n1_Psi2_plus += q_vector(+2, 1, weight, phi);
      Q_0_Psi2_plus += q_vector(0, 1, weight, phi);

    }
    //For Psi_2 event plane angle:
    if( trkEta > -1.0 && trkEta < 1.0){

      TRK_Psi_2_sine += weight*sin( 2*phi );
      TRK_Psi_2_cosine += weight*cos( 2*phi );
    }

    if( trkEta < 0.0 && trkEta > -0.8 ){
   
      Q_n3_trk_minus += q_vector(+1, 1, -weight, phi);//for scalar product in tracker
      Q_0_trk_minus += q_vector(0, 1, weight, phi);

      Q_n3_trk_plus += q_vector(+1, 1, -weight, phi);//for scalar product in tracker
      Q_0_trk_plus += q_vector(0, 1, weight, phi);
    }
    if( trkEta < 0.8 && trkEta > 0.0 ){
    
      Q_n3_trk_minus += q_vector(+1, 1, weight, phi);//for scalar product in tracker
      Q_0_trk_minus += q_vector(0, 1, weight, phi);

      Q_n3_trk_plus += q_vector(+1, 1, weight, phi);//for scalar product in tracker
      Q_0_trk_plus += q_vector(0, 1, weight, phi);
    }
    
    for(int eta = 0; eta < NetaBins; eta++){
      if( trkEta > etaBins_[eta] && trkEta < etaBins_[eta+1] ){

          double delta_phi = phi - Psi_1;
          if( delta_phi > PI ) {
            
            delta_phi = 2*PI - delta_phi;
            //if( delta_phi > PI/2.0 ) delta_phi = PI-delta_phi;

          }
          if( delta_phi < -PI ) {
            delta_phi = delta_phi + 2*PI;

            //if(delta_phi > PI/2.0) delta_phi = PI-delta_phi;
          }
              
          if( trk.charge() == +1 ){//positive charge

            //3p:
            Q_n1_1[eta][0] += q_vector(+1, 1, weight, phi);
            Q_0_1[eta][0] += q_vector(0, 1, weight, phi);

            delta_phi_positive[eta]->Fill(delta_phi);
          }
          if( trk.charge() == -1 ){//negative charge

            Q_n1_1[eta][1] += q_vector(+1, 1, weight, phi);
            Q_0_1[eta][1] += q_vector(0, 1, weight, phi);
            
            delta_phi_negative[eta]->Fill(delta_phi);

          }
        }
      
    }//end of eta dimension
  }


  //Psi_2 event plane angle:

  TRK_Psi_2_sine = TRK_Psi_2_sine - 0.0423;
  TRK_Psi_2_cosine = TRK_Psi_2_cosine - 0.1018;
  double Psi_2 = TMath::ATan(TRK_Psi_2_sine/TRK_Psi_2_cosine)/2;
  Psi_2_cos->Fill(TRK_Psi_2_cosine);
  Psi_2_sin->Fill(TRK_Psi_2_sine);

  //Two terms are calculated 
  //<cos(phi+Psi_1-2*Psi_2)>
  //<cos(2*(Psi_1-Psi_2))
  double term_1[NetaBins];
  double term_1_weight[NetaBins];

  double term_2 = 0.0;
  double term_2_weight = 0.0;

  for(int i = 0; i < NetaBins; i++){

    term_1[i] = 0.0;
    term_1_weight[i] = 0.0;
  }

  //track loop to calculate the correlator:
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
    
    if( !doPixelReco_){ if(fabs(dxyvtx/dxyerror) > offlineDCA_) continue; }
    if(chi2n > offlineChi2_ ) continue;

    if( doPixelReco_ ){ if(nhits != 3 && nhits != 4 && nhits != 5 && nhits != 6) continue;}
    else{ if(nhits < offlinenhits_ ) continue; if( nPixelLayers <= 0 ) continue;}
    
    if(trk.pt() < ptLow_ || trk.pt() > ptHigh_ ) continue;
    if(fabs(trkEta) > etaTracker_ ) continue;

    for(int eta = 0; eta < NetaBins; eta++){
      if( trkEta > etaBins_[eta] && trkEta < etaBins_[eta+1] ){

        double cos_phi_angle = cos(phi)-phiCos_[eta];
        double sin_phi_angle = sin(phi)-phiSin_[eta];

        term_1[eta] += weight*(cos_phi_angle*cos(Psi_1-2*Psi_2) - sin_phi_angle*sin(Psi_1-2*Psi_2));
        term_1_weight[eta] += weight;

        Phi_Average_cos[eta]->Fill( cos_phi_angle );
        Phi_Average_sin[eta]->Fill( sin_phi_angle );
      
      }
    }
    
  }

  for(int eta = 0; eta < NetaBins; eta++){

    Phi_Psi_1_Psi_2[eta]->Fill(term_1[eta]/term_1_weight[eta], term_1_weight[eta]); 
  }
  
 
  term_2 = cos(2*Psi_1-2*Psi_2);
  term_2_weight = nTracks; //Use track multiplicity
  Psi_1_Psi_2->Fill(term_2, term_2_weight);

  TComplex N_2_trk, D_2_trk;

  N_2_trk = Q_n1_Psi2_plus*TComplex::Conjugate(Q_n1_Psi2_minus);
  D_2_trk = Q_0_Psi2_plus*Q_0_Psi2_minus;

  Psi_2_trk_reso->Fill(N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re());
  Psi_2_trk_accept_real[0]->Fill(Q_n1_Psi2_minus.Re()/Q_0_Psi2_minus.Re(), Q_0_Psi2_minus.Re());
  Psi_2_trk_accept_imag[0]->Fill(Q_n1_Psi2_minus.Im()/Q_0_Psi2_minus.Re(), Q_0_Psi2_minus.Re());

  Psi_2_trk_accept_real[1]->Fill(Q_n1_Psi2_plus.Re()/Q_0_Psi2_plus.Re(), Q_0_Psi2_plus.Re());
  Psi_2_trk_accept_imag[1]->Fill(Q_n1_Psi2_plus.Im()/Q_0_Psi2_plus.Re(), Q_0_Psi2_plus.Re());

/*
event average v1
*/

//resolution factor
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
  c2_c_plus_real->Fill( Q_n3_trk_plus.Re()/Q_0_trk_plus.Re(), Q_0_trk_plus.Re() );
  c2_c_minus_real->Fill( Q_n3_trk_minus.Re()/Q_0_trk_minus.Re(), Q_0_trk_minus.Re() );

  c2_a_imag->Fill( Q_n3_1_HFminus.Im()/Q_0_1_HFminus.Re(), Q_0_1_HFminus.Re() );
  c2_b_imag->Fill( Q_n3_1_HFplus.Im()/Q_0_1_HFplus.Re(), Q_0_1_HFplus.Re() );
  c2_c_plus_imag->Fill( Q_n3_trk_plus.Im()/Q_0_trk_plus.Re(), Q_0_trk_plus.Re() );
  c2_c_minus_imag->Fill( Q_n3_trk_minus.Im()/Q_0_trk_minus.Re(), Q_0_trk_minus.Re() );

  N_2_trk = Q_n3_1_HFcombined*TComplex::Conjugate(Q_n3_1_HFcombined);
  D_2_trk = Q_0_1_HFcombined*Q_0_1_HFcombined;

  c2_ab_combined->Fill( N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re() );
  c2_ab_real->Fill( Q_n3_1_HFcombined.Re()/Q_0_1_HFcombined.Re(), Q_0_1_HFcombined.Re() );
  c2_ab_imag->Fill( Q_n3_1_HFcombined.Im()/Q_0_1_HFcombined.Re(), Q_0_1_HFcombined.Re() );


//numerator
  for(int eta = 0; eta < NetaBins; eta++){
    for(int charge = 0; charge < 2; charge++){

      TComplex N_v1_A_SP, D_v1_A_SP, N_v1_B_SP, D_v1_B_SP;
      TComplex N_v1_AB_SP, D_v1_AB_SP;

      N_v1_A_SP = Q_n1_1[eta][charge]*Q_n3_1_HFminus;
      D_v1_A_SP = Q_0_1[eta][charge]*Q_0_1_HFminus;

      double V1_A = N_v1_A_SP.Re()/D_v1_A_SP.Re();

      N_v1_B_SP = Q_n1_1[eta][charge]*Q_n3_1_HFplus;
      D_v1_B_SP = Q_0_1[eta][charge]*Q_0_1_HFplus;

      double V1_B = N_v1_B_SP.Re()/D_v1_B_SP.Re();

      N_v1_AB_SP = Q_n1_1[eta][charge]*Q_n3_1_HFcombined;
      D_v1_AB_SP = Q_0_1[eta][charge]*Q_0_1_HFcombined;

      double V1_AB = N_v1_AB_SP.Re()/D_v1_AB_SP.Re();

      c2_v1[eta][charge][0]->Fill( V1_A, D_v1_A_SP.Re() );
      c2_v1[eta][charge][1]->Fill( V1_B, D_v1_B_SP.Re() );
      c2_v1[eta][charge][2]->Fill( V1_AB, D_v1_AB_SP.Re() );

      c2_trk_accept[eta][charge][0]->Fill(Q_n1_1[eta][charge].Re()/Q_0_1[eta][charge].Re(), Q_0_1[eta][charge].Re());
      c2_trk_accept[eta][charge][1]->Fill(Q_n1_1[eta][charge].Im()/Q_0_1[eta][charge].Re(), Q_0_1[eta][charge].Re());

    }
  }

}
// ------------ method called once each job just before starting event loop  ------------
void 
DirectedFlowCorrelatorTest::beginJob()
{
  edm::Service<TFileService> fs;
    
  TH1D::SetDefaultSumw2();

  const int NetaBins = etaBins_.size() - 1 ;

  double ptBinsArray[100];
  const int Nptbins = ptBins_.size() - 1;
  for(unsigned i = 0; i < ptBins_.size(); i++){
    ptBinsArray[i] = ptBins_[i];
  }

  if( !doPixelReco_ ){
    edm::FileInPath fip1("DirectedFlowCorrelator/DirectedFlowCorrelator/data/Hydjet_eff_mult_v1.root");
    TFile f1(fip1.fullPath().c_str(),"READ");
    for(int i = 0; i < 5; i++){
       effTable[i] = (TH2D*)f1.Get(Form("rTotalEff3D_%d",i));
    }
  }
  else{
    edm::FileInPath fip2("DirectedFlowCorrelator/DirectedFlowCorrelator/data/EffCorrectionsPixel_TT_pt_0_10_v2.root");
    TFile f2(fip2.fullPath().c_str(),"READ");
     
    effTable[0] = (TH2D*)f2.Get("Eff_50_100");
    effTable[1] = (TH2D*)f2.Get("Eff_30_50");
    effTable[2] = (TH2D*)f2.Get("Eff_10_30");
    effTable[3] = (TH2D*)f2.Get("Eff_5_10");
    effTable[4] = (TH2D*)f2.Get("Eff_0_5");
  }

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",5000,0,5000);
  vtxZ = fs->make<TH1D>("vtxZ",";vz", 400,-20,20);
  cbinHist = fs->make<TH1D>("cbinHist",";cbin",200,0,200);
  trkPhi = fs->make<TH1D>("trkPhi", ";#phi", 700, -3.5, 3.5);
  hfPhi = fs->make<TH1D>("hfPhi", ";#phi", 700, -3.5, 3.5);
  trkPt = fs->make<TH1D>("trkPt", ";p_{T}(GeV)", Nptbins,ptBinsArray);
  trk_eta = fs->make<TH1D>("trk_eta", ";#eta", 50,-2.5,2.5);

  for(int eta = 0; eta < NetaBins; eta++){
    for(int charge = 0; charge < 2; charge++){
      for(int dir = 0; dir < 3; dir++){

        c2_v1[eta][charge][dir] = fs->make<TH1D>(Form("c2_v1_%d_%d_%d",eta,charge,dir),";c1", 1,-1,1);
        c2_trk_accept[eta][charge][dir] = fs->make<TH1D>(Form("c2_trk_accept_%d_%d_%d",eta,charge,dir), ";c1", 1,-1,1);

      }
    }
  }

  c2_ab = fs->make<TH1D>("c2_ab",";c2_ab", 1,-1,1);
  c2_ab_combined = fs->make<TH1D>("c2_ab_combined",";c2_ab_combined", 1,-1,1);

  c2_ac_plus = fs->make<TH1D>("c2_ac_plus",";c2_ac_plus", 1,-1,1);
  c2_cb_plus = fs->make<TH1D>("c2_cb_plus",";c2_cb_plus", 1,-1,1);

  c2_ac_minus = fs->make<TH1D>("c2_ac_minus",";c2_ac_minus", 1,-1,1);
  c2_cb_minus = fs->make<TH1D>("c2_cb_minus",";c2_cb_minus", 1,-1,1);

  c2_a_real = fs->make<TH1D>("c2_a_real",";c2_a_real", 1,-1,1);
  c2_b_real = fs->make<TH1D>("c2_b_real",";c2_b_real", 1,-1,1);
  c2_c_plus_real = fs->make<TH1D>("c2_c_plus_real",";c2_c_plus_real", 1,-1,1);
  c2_c_minus_real = fs->make<TH1D>("c2_c_minus_real",";c2_c_minus_real", 1,-1,1);
  c2_ab_real = fs->make<TH1D>("c2_ab_real",";c2_ab_real", 1,-1,1);

  c2_a_imag = fs->make<TH1D>("c2_a_imag",";c2_a_imag", 1,-1,1);
  c2_b_imag = fs->make<TH1D>("c2_b_imag",";c2_b_imag", 1,-1,1);
  c2_c_plus_imag = fs->make<TH1D>("c2_c_plus_imag",";c2_c_plus_imag", 1,-1,1);
  c2_c_minus_imag = fs->make<TH1D>("c2_c_minus_imag",";c2_c_minus_imag", 1,-1,1);
  c2_ab_imag = fs->make<TH1D>("c2_ab_imag",";c2_ab_imag", 1,-1,1);

  for(int eta = 0; eta < NetaBins; eta++){

    Phi_Psi_1_Psi_2[eta] = fs->make<TH1D>(Form("Phi_Psi_1_Psi_2_%d",eta),";Phi_Psi_1_Psi_2", 1,-1,1);
    Phi_Average_cos[eta] = fs->make<TH1D>(Form("Phi_Average_cos_%d",eta),";Phi_Average_cos", 1,-1,1);
    Phi_Average_sin[eta] = fs->make<TH1D>(Form("Phi_Average_sin_%d",eta),";Phi_Average_sin", 1,-1,1);
  }

  Psi_1_Psi_2 = fs->make<TH1D>("Psi_1_Psi_2",";Psi_1_Psi_2", 1,-1,1);
  Psi_2_trk_reso = fs->make<TH1D>("Psi_2_trk_reso",";Psi_2_trk_reso", 1,-1,1);

  for(int charge = 0; charge < 2; charge++){

    Psi_2_trk_accept_real[charge] = fs->make<TH1D>(Form("Psi_2_trk_accept_real_%d",charge),";Psi_2_trk_accept_real", 1,-1,1);
    Psi_2_trk_accept_imag[charge] = fs->make<TH1D>(Form("Psi_2_trk_accept_imag_%d",charge),";Psi_2_trk_accept_imag", 1,-1,1);
  }

  Psi_2_sin = fs->make<TH1D>("Psi_2_sin",";Psi_2_sin", 1,-10000,10000);
  Psi_2_cos = fs->make<TH1D>("Psi_2_cos",";Psi_2_cos", 1,-10000,10000);
  Psi_1_sin = fs->make<TH1D>("Psi_1_sin",";Psi_1_sin", 1,-10000,10000);
  Psi_1_cos = fs->make<TH1D>("Psi_1_cos",";Psi_1_cos", 1,-10000,10000);
  
  for(int eta = 0; eta < NetaBins; eta++){

    delta_phi_positive[eta] = fs->make<TH1D>(Form("delta_phi_positive_%d",eta),";#Delta#phi", 20,-PI,PI);
    delta_phi_negative[eta] = fs->make<TH1D>(Form("delta_phi_negative_%d",eta),";#Delta#phi", 20,-PI,PI);
  }

}
TComplex 
DirectedFlowCorrelatorTest::q_vector(double n, double p, double w, double phi) 
{
  double term1 = pow(w,p);
  TComplex e(1, n*phi, 1);
  return term1*e;
}
// ------------ method called once each job just after ending the event loop  ------------
void 
DirectedFlowCorrelatorTest::endJob() 
{
}
void 
DirectedFlowCorrelatorTest::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DirectedFlowCorrelatorTest::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DirectedFlowCorrelatorTest::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DirectedFlowCorrelatorTest::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DirectedFlowCorrelatorTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DirectedFlowCorrelatorTest);
