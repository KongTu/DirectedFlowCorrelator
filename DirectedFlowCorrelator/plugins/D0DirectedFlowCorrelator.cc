// -*- C++ -*-
//
// Package:    D0DirectedFlowCorrelator/D0DirectedFlowCorrelator
// Class:      D0DirectedFlowCorrelator
// 
/**\class D0DirectedFlowCorrelator D0DirectedFlowCorrelator.cc D0DirectedFlowCorrelator/D0DirectedFlowCorrelator/plugins/D0DirectedFlowCorrelator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhoudunming Tu
//         Created:  Wed, 20 Sep 2017 22:44:00 GMT
//
//


#include "DirectedFlowCorrelator/DirectedFlowCorrelator/interface/D0DirectedFlowCorrelatorBase.h"


D0DirectedFlowCorrelator::D0DirectedFlowCorrelator(const edm::ParameterSet& iConfig)

{
  trackName_  =  iConfig.getParameter<edm::InputTag>("trackName");
  vertexName_ =  iConfig.getParameter<edm::InputTag>("vertexName");
  towerName_ =  iConfig.getParameter<edm::InputTag>("towerName");

  trackSrc_ = consumes<reco::TrackCollection>(trackName_);
  vertexSrc_ = consumes<reco::VertexCollection>(vertexName_);
  towerSrc_ = consumes<CaloTowerCollection>(towerName_);

  centralityToken_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
  centralityBinToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinSrc"));

  recoVertexCompositeCandidateCollection_Token_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("D0src"));
  Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
  Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxTruncated40"));

  Nmin_ = iConfig.getUntrackedParameter<int>("Nmin");
  Nmax_ = iConfig.getUntrackedParameter<int>("Nmax");

  useCentrality_ = iConfig.getUntrackedParameter<bool>("useCentrality");
  doEffCorrection_ = iConfig.getUntrackedParameter<bool>("doEffCorrection");
  doD0EffCorrection_ = iConfig.getUntrackedParameter<bool>("doD0EffCorrection");
  useEtaGap_ = iConfig.getUntrackedParameter<bool>("useEtaGap");
  doBothSide_ = iConfig.getUntrackedParameter<bool>("doBothSide");
  doPixelReco_ = iConfig.getUntrackedParameter<bool>("doPixelReco");
  doHiReco_ = iConfig.getUntrackedParameter<bool>("doHiReco");
  doTight_ = iConfig.getUntrackedParameter<bool>("doTight");
  doLoose_ = iConfig.getUntrackedParameter<bool>("doLoose");

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
  rapidityBins_ = iConfig.getUntrackedParameter<std::vector<double>>("rapidityBins");
  ptBins_ = iConfig.getUntrackedParameter<std::vector<double>>("ptBins");
  centBins_ = iConfig.getUntrackedParameter<std::vector<double>>("centBins");
  massBins_ = iConfig.getUntrackedParameter<std::vector<double>>("massBins");

  D0MassHigh_ = iConfig.getUntrackedParameter<double>("D0MassHigh", 1.92);
  D0MassLow_ = iConfig.getUntrackedParameter<double>("D0MassLow", 1.82);
  D0EtaHigh_ = iConfig.getUntrackedParameter<double>("D0EtaHigh", 2.4);
  D0EtaLow_ = iConfig.getUntrackedParameter<double>("D0EtaLow", -2.4);
  D0PtHigh_ = iConfig.getUntrackedParameter<double>("D0PtHigh", 40.0);
  D0PtLow_ = iConfig.getUntrackedParameter<double>("D0PtLow", 1.0);
  D0YHigh_ = iConfig.getUntrackedParameter<double>("D0YHigh", 1.8);
  D0YLow_ = iConfig.getUntrackedParameter<double>("D0YLow", -1.8);
  NSigmas_ = iConfig.getUntrackedParameter<double>("NSigmas", 3.0);
  
  D0DcaHigh_ = iConfig.getUntrackedParameter<std::vector<double>>("D0DcaHigh");
  D0VtxProbLow_ = iConfig.getUntrackedParameter<std::vector<double>>("D0VtxProbLow");
  D03DAngleHigh_ = iConfig.getUntrackedParameter<std::vector<double>>("D03DAngleHigh");
  D0DlosLow_ = iConfig.getUntrackedParameter<std::vector<double>>("D0DlosLow");
  TrkPtLow_ = iConfig.getUntrackedParameter<std::vector<double>>("TrkPtLow");
  D0MassWidth_ = iConfig.getUntrackedParameter<std::vector<double>>("D0MassWidth");

  TrkEtaHigh_ = iConfig.getUntrackedParameter<double>("TrkEtaHigh", 2.4);
  TrkEtaLow_ = iConfig.getUntrackedParameter<double>("TrkEtaLow", -2.4);
  TrkChiOverNLayerHigh_ = iConfig.getUntrackedParameter<double>("TrkChiOverNLayerHigh", 0.15);
  TrkPtErrOverPtHigh_ = iConfig.getUntrackedParameter<double>("TrkPtErrOverPtHigh", 0.1);
  TrkNHitLow_ = iConfig.getUntrackedParameter<int>("NHitLow", 11);

}


D0DirectedFlowCorrelator::~D0DirectedFlowCorrelator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
D0DirectedFlowCorrelator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

    if( doPixelReco_ || doHiReco_ ){

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

  const int NetaBins = etaBins_.size() - 1;
  const int NyBins = rapidityBins_.size() - 1;
  const int NmassBins = massBins_.size() - 1;

  TComplex Q_n3_1_HFplus, Q_n3_1_HFminus, Q_0_1_HFplus, Q_0_1_HFminus, Q_n3_1_HFplusANDminus, Q_0_1_HFplusANDminus;
  //HF towers loop to fill the towers' Q-vectors:
  TComplex Q_n3_trk_minus, Q_0_trk_minus, Q_n3_trk_plus, Q_0_trk_plus;
  //charge independent, |eta|<1.0
  TComplex Q_n1_1[NetaBins][2], Q_0_1[NetaBins][2];
  //D0 Q-vectors for both obs and bkg
  TComplex Q_D0obs_n1_1[NyBins][3], Q_D0obs_0_1[NyBins][3], Q_D0bkg_n1_1[NyBins][3], Q_D0bkg_0_1[NyBins][3];
  //D0 mass-fit method, no separation of obs and bkg
  TComplex Q_D0_n1_1[NyBins][30][3], Q_D0_0_1[NyBins][30][3];

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

              Q_n3_1_HFplusANDminus += q_vector(-1, 1, w, caloPhi);
              Q_0_1_HFplusANDminus += q_vector(0, 1, w, caloPhi);

              HF_Psi_1_sine += w*sin( 1*caloPhi );
              HF_Psi_1_cosine += w*cos( 1*caloPhi );

          }
          else if( caloEta < -etaLowHF_ && caloEta > -etaHighHF_ ){

              Q_n3_1_HFminus += q_vector(-1, 1, -w, caloPhi);
              Q_0_1_HFminus += q_vector(0, 1, w, caloPhi); //normalization needs to be positive in order to be not cancel out

              Q_n3_1_HFplusANDminus += q_vector(-1, 1, -w, caloPhi);
              Q_0_1_HFplusANDminus += q_vector(0, 1, w, caloPhi);

              HF_Psi_1_sine += -w*sin( 1*caloPhi );
              HF_Psi_1_cosine += -w*cos( 1*caloPhi );
          }
          else{continue;}
  }

  double Psi_1 = TMath::ATan(HF_Psi_1_sine/HF_Psi_1_cosine)/1;
  Psi_1_cos->Fill(HF_Psi_1_cosine);
  Psi_1_sin->Fill(HF_Psi_1_sine);

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
    else{ 
      if(nhits < offlinenhits_ ) continue; 
      if( nPixelLayers <= 0 ) continue;
    }
    
    if(trk.pt() < ptLow_ || trk.pt() > ptHigh_ ) continue;
    if(fabs(trkEta) > etaTracker_ ) continue;

    trkPhi->Fill(phi, weight);
    trkPt->Fill(trk.pt(), weight);
    trk_eta->Fill(trkEta, weight);

    if( trkEta < -1.0 && trkEta > -2.4 ){
   
      Q_n3_trk_minus += q_vector(+1, 1, -weight, phi);//for scalar product in tracker
      Q_0_trk_minus += q_vector(0, 1, weight, phi);
    }
    if( trkEta < 2.4 && trkEta > 1.0 ){
   
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
//end of track loop

/*
D0 candiates' loop
 */

  edm::Handle<reco::VertexCompositeCandidateCollection> D0candidates;
  iEvent.getByToken(recoVertexCompositeCandidateCollection_Token_,D0candidates);
  if(!D0candidates.isValid()) return;
  const reco::VertexCompositeCandidateCollection* D0 = D0candidates.product();       

    for(unsigned it=0; it<D0->size(); ++it) {
        
      const reco::VertexCompositeCandidate & trk = (*D0)[it];
      double secvz = -999.9, secvx = -999.9, secvy = -999.9;
      secvz = trk.vz();
      secvx = trk.vx();
      secvy = trk.vy();
      double px = trk.px();
      double py = trk.py();
      double pz = trk.pz();
      TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
      TVector3 secvec(px,py,pz);

      //vtxChi2
      double vtxChi2 = trk.vertexChi2();
      double ndf = trk.vertexNdof();
      double VtxProb = TMath::Prob(vtxChi2,ndf);

      //PAngle
      double agl_abs = secvec.Angle(ptosvec);

      //Decay length 3D
      typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
      typedef ROOT::Math::SVector<double, 3> SVector3;
      typedef ROOT::Math::SVector<double, 6> SVector6;

      SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
      SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);

      double dl = ROOT::Math::Mag(distanceVector);
      double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;        
      double dlos = dl/dlerror;
      double d0dca = dl*sin(agl_abs);

      double y_D0 = trk.rapidity();
      double pt = trk.pt();
      double phi = trk.phi();
      double mass = trk.mass();
      double eta = trk.eta();

      if (pt < D0PtLow_) continue;
      if (eta < D0EtaLow_ || eta > D0EtaHigh_) continue;
      
      int rap_index = 0;

      if( fabs(y_D0) > 0.0 && fabs(y_D0) < 0.6 ){rap_index = 0;}
      if( fabs(y_D0) > 0.6 && fabs(y_D0) < 1.2 ){rap_index = 1;}
      if( fabs(y_D0) > 1.2 && fabs(y_D0) < 2.0 ){rap_index = 2;}

      if (d0dca > D0DcaHigh_[rap_index]) continue;
      if (VtxProb < D0VtxProbLow_[rap_index]) continue;
      if (agl_abs > D03DAngleHigh_[rap_index]) continue;
      if (dlos < D0DlosLow_[rap_index]) continue;

      const reco::Candidate * d1 = trk.daughter(0);
      const reco::Candidate * d2 = trk.daughter(1);
              
      auto dau1 = d1->get<reco::TrackRef>();
      auto dau2 = d2->get<reco::TrackRef>();

      //trk quality       
      bool trkquality1 = dau1->quality(reco::TrackBase::highPurity);
      bool trkquality2 = dau2->quality(reco::TrackBase::highPurity);
      //track pt
      double pt1 = d1->pt();
      double pt2 = d2->pt();
      //track eta
      double eta1 = d1->eta();
      double eta2 = d2->eta();
      //track charge
      int charge1 = d1->charge();
      int charge2 = d2->charge();
      // //mass
      double mass1 = d1->mass();
      double mass2 = d2->mass();
      //track pT error
      double ptErr1 = dau1->ptError();
      double ptErr2 = dau2->ptError();
      //trkNHits
      int nhit1 = dau1->numberOfValidHits();
      int nhit2 = dau2->numberOfValidHits();

      //track nlayer
      double nlayer1 = dau1->hitPattern().trackerLayersWithMeasurement();
      double nlayer2 = dau2->hitPattern().trackerLayersWithMeasurement();

      //track Chi2
      double trkChi1 = dau1->normalizedChi2();
      double trkChi2 = dau2->normalizedChi2();

      if (trkChi1/nlayer1 > TrkChiOverNLayerHigh_) continue;
      if (trkChi2/nlayer2 > TrkChiOverNLayerHigh_) continue;

      if (!trkquality1) continue;
      if (!trkquality2) continue;

      if (pt1 < TrkPtLow_[rap_index]) continue;
      if (pt2 < TrkPtLow_[rap_index]) continue;

      if (eta1 < TrkEtaLow_ || eta1 > TrkEtaHigh_) continue;
      if (eta2 < TrkEtaLow_ || eta2 > TrkEtaHigh_) continue;

      if (ptErr1/pt1 > TrkPtErrOverPtHigh_) continue;
      if (ptErr2/pt2 > TrkPtErrOverPtHigh_) continue;

      if (nhit1 < TrkNHitLow_) continue;
      if (nhit2 < TrkNHitLow_) continue;

      double weight_D0 = 1.0;

      for(int rap = 0; rap < NyBins; rap++){
        for(int imass = 0; imass < NmassBins; imass++){
          if( y_D0 > rapidityBins_[rap] && y_D0 < rapidityBins_[rap+1] && mass > massBins_[imass] && mass < massBins_[imass+1]){

            if( doD0EffCorrection_ && !doTight_ && !doLoose_ ){

               int index = 0;
               if( rap == 0 || rap == 5 ) index = 2;
               if( rap == 1 || rap == 4 ) index = 1;
               if( rap == 2 || rap == 3 ) index = 0;
  
                if( pt > 3.0 && pt < 30.0 ){
                  weight_D0 = 1.0/(f1[index]->Eval(trk.pt()));
                }
                else if(pt < 3.0){
                  weight_D0 = 1.0/(f1[index]->Eval(3.0));
                }
                else if(pt > 30.0){
                  weight_D0 = 1.0/(f1[index]->Eval(30.0));
                }
            }
            if( doD0EffCorrection_ && doTight_ && !doLoose_ ){

               int index = 0;
               if( rap == 0 || rap == 5 ) index = 2;
               if( rap == 1 || rap == 4 ) index = 1;
               if( rap == 2 || rap == 3 ) index = 0;
  
                if( pt > 3.0 && pt < 30.0 ){
                  weight_D0 = 1.0/(f2[index]->Eval(trk.pt()));
                }
                else if(pt < 3.0){
                  weight_D0 = 1.0/(f2[index]->Eval(3.0));
                }
                else if(pt > 30.0){
                  weight_D0 = 1.0/(f2[index]->Eval(30.0));
                }
            }
            if( doD0EffCorrection_ && !doTight_ && doLoose_ ){

               int index = 0;
               if( rap == 0 || rap == 5 ) index = 2;
               if( rap == 1 || rap == 4 ) index = 1;
               if( rap == 2 || rap == 3 ) index = 0;
  
                if( pt > 3.0 && pt < 30.0 ){
                  weight_D0 = 1.0/(f33[index]->Eval(trk.pt()));
                }
                else if(pt < 3.0){
                  weight_D0 = 1.0/(f33[index]->Eval(3.0));
                }
                else if(pt > 30.0){
                  weight_D0 = 1.0/(f33[index]->Eval(30.0));
                }
            }

          //signal+bkg region altogether
          //D0
            if( charge1 == +1 && mass1 < 0.14 && mass1 > 0.13 ){
              
              Q_D0_n1_1[rap][imass][0] += q_vector(+1, 1, weight_D0, phi);
              Q_D0_0_1[rap][imass][0] += q_vector(0, 1, weight_D0, phi);
            }
          //D0bar
            if( charge2 == -1 && mass2 < 0.14 && mass2 > 0.13 ){

              Q_D0_n1_1[rap][imass][1] += q_vector(+1, 1, weight_D0, phi);
              Q_D0_0_1[rap][imass][1] += q_vector(0, 1, weight_D0, phi);
            }
          //inclusive D0
  
            Q_D0_n1_1[rap][imass][2] += q_vector(+1, 1, weight_D0, phi);
            Q_D0_0_1[rap][imass][2] += q_vector(0, 1, weight_D0, phi);
            
          }
        }
      }

      for(int rap = 0; rap < NyBins; rap++){
        if( y_D0 > rapidityBins_[rap] && y_D0 < rapidityBins_[rap+1] ){

        //Fill mass in each y
        if( charge1 == +1 && mass1 < 0.14 && mass1 > 0.13 ){

          D0Mass_Hist[rap][0]->Fill( mass );
        }
        if( charge2 == -1 && mass2 < 0.14 && mass2 > 0.13 ){

          D0Mass_Hist[rap][1]->Fill( mass );
        }

        D0Mass_Hist[rap][2]->Fill(mass);

        //signal region
          //D0
            if( charge1 == +1 && mass1 < 0.14 && mass1 > 0.13 && mass > (D0Mass-NSigmas_*D0MassWidth_[rap]) && mass < (D0Mass+NSigmas_*D0MassWidth_[rap]) ){
              
              Q_D0obs_n1_1[rap][0] += q_vector(+1, 1, weight_D0, phi);
              Q_D0obs_0_1[rap][0] += q_vector(0, 1, weight_D0, phi);
            }
          //D0bar
            if( charge2 == -1 && mass2 < 0.14 && mass2 > 0.13 && mass > (D0Mass-NSigmas_*D0MassWidth_[rap]) && mass < (D0Mass+NSigmas_*D0MassWidth_[rap]) ){

              Q_D0obs_n1_1[rap][1] += q_vector(+1, 1, weight_D0, phi);
              Q_D0obs_0_1[rap][1] += q_vector(0, 1, weight_D0, phi);
            }
          //inclusive D0
            if( mass > (D0Mass-NSigmas_*D0MassWidth_[rap]) && mass < (D0Mass+NSigmas_*D0MassWidth_[rap]) ){

              Q_D0obs_n1_1[rap][2] += q_vector(+1, 1, weight_D0, phi);
              Q_D0obs_0_1[rap][2] += q_vector(0, 1, weight_D0, phi);
            }
        //bkg region
          //D0
            if( charge1 == +1 && mass1 < 0.14 && mass1 > 0.13 && (mass < (D0Mass-NSigmas_*D0MassWidth_[rap]) || mass > (D0Mass+NSigmas_*D0MassWidth_[rap])) ){

              Q_D0bkg_n1_1[rap][0] += q_vector(+1, 1, weight_D0, phi);
              Q_D0bkg_0_1[rap][0] += q_vector(0, 1, weight_D0, phi);

            }
          //D0bar
            if( charge2 == -1 && mass2 < 0.14 && mass2 > 0.13 && (mass < (D0Mass-NSigmas_*D0MassWidth_[rap]) || mass > (D0Mass+NSigmas_*D0MassWidth_[rap])) ){

              Q_D0bkg_n1_1[rap][1] += q_vector(+1, 1, weight_D0, phi);
              Q_D0bkg_0_1[rap][1] += q_vector(0, 1, weight_D0, phi);
            }
          //inclusive D0
            if( mass < (D0Mass-NSigmas_*D0MassWidth_[rap]) || mass > (D0Mass+NSigmas_*D0MassWidth_[rap]) ){

              Q_D0bkg_n1_1[rap][2] += q_vector(+1, 1, weight_D0, phi);
              Q_D0bkg_0_1[rap][2] += q_vector(0, 1, weight_D0, phi);

            }
      }
    }

    for(int rap = 0; rap < NyBins; rap++){
        if( y_D0 > rapidityBins_[rap] && y_D0 < rapidityBins_[rap+1] ){
          
          double delta_phi = phi - Psi_1;
          if( delta_phi > PI ) delta_phi = 2*PI - delta_phi;
          if( delta_phi < -PI ) delta_phi = delta_phi + 2*PI;

          //Fill mass in each y
          if( charge1 == +1 && mass1 < 0.14 && mass1 > 0.13 ){

            D0Mass_Hist_DeltaPhi[rap][0]->Fill(delta_phi, mass );
          }
          if( charge2 == -1 && mass2 < 0.14 && mass2 > 0.13 ){

            D0Mass_Hist_DeltaPhi[rap][1]->Fill(delta_phi, mass );
          }

          D0Mass_Hist_DeltaPhi[rap][2]->Fill(delta_phi, mass);
        }
    }
 
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
  c2_c_plus_real->Fill( Q_n3_trk_plus.Re()/Q_0_trk_plus.Re(), Q_0_trk_plus.Re() );
  c2_c_minus_real->Fill( Q_n3_trk_minus.Re()/Q_0_trk_minus.Re(), Q_0_trk_minus.Re() );

  c2_a_imag->Fill( Q_n3_1_HFminus.Im()/Q_0_1_HFminus.Re(), Q_0_1_HFminus.Re() );
  c2_b_imag->Fill( Q_n3_1_HFplus.Im()/Q_0_1_HFplus.Re(), Q_0_1_HFplus.Re() );
  c2_c_plus_imag->Fill( Q_n3_trk_plus.Im()/Q_0_trk_plus.Re(), Q_0_trk_plus.Re() );
  c2_c_minus_imag->Fill( Q_n3_trk_minus.Im()/Q_0_trk_minus.Re(), Q_0_trk_minus.Re() );

//combine two HF into one plane:

  N_2_trk = Q_n3_1_HFplusANDminus*TComplex::Conjugate(Q_n3_1_HFplusANDminus);
  D_2_trk = Q_0_1_HFplusANDminus*Q_0_1_HFplusANDminus;

  c2_ab_one->Fill( N_2_trk.Re()/D_2_trk.Re(), D_2_trk.Re());

  c2_ab_one_real->Fill( Q_n3_1_HFplusANDminus.Re()/Q_0_1_HFplusANDminus.Re(), Q_0_1_HFplusANDminus.Re() );
  c2_ab_one_imag->Fill( Q_n3_1_HFplusANDminus.Im()/Q_0_1_HFplusANDminus.Re(), Q_0_1_HFplusANDminus.Re() );

//numerator charged particle
  for(int eta = 0; eta < NetaBins; eta++){
    for(int charge = 0; charge < 2; charge++){

      TComplex N_v1_A_SP, D_v1_A_SP, N_v1_B_SP, D_v1_B_SP, N_v1_AB_SP, D_v1_AB_SP;

      N_v1_A_SP = Q_n1_1[eta][charge]*Q_n3_1_HFminus;
      D_v1_A_SP = Q_0_1[eta][charge]*Q_0_1_HFminus;

      double V1_A = N_v1_A_SP.Re()/D_v1_A_SP.Re();

      N_v1_B_SP = Q_n1_1[eta][charge]*Q_n3_1_HFplus;
      D_v1_B_SP = Q_0_1[eta][charge]*Q_0_1_HFplus;

      double V1_B = N_v1_B_SP.Re()/D_v1_B_SP.Re();

      //use one plane
      N_v1_AB_SP = Q_n1_1[eta][charge]*Q_n3_1_HFplusANDminus;
      D_v1_AB_SP = Q_0_1[eta][charge]*Q_0_1_HFplusANDminus;

      double V1_AB = N_v1_AB_SP.Re()/D_v1_AB_SP.Re();

      c2_v1[eta][charge][0]->Fill( V1_A, D_v1_A_SP.Re() );
      c2_v1[eta][charge][1]->Fill( V1_B, D_v1_B_SP.Re() );
      c2_v1[eta][charge][2]->Fill( V1_AB, D_v1_AB_SP.Re() );

      c2_trk_accept[eta][charge][0]->Fill(Q_n1_1[eta][charge].Re()/Q_0_1[eta][charge].Re(), Q_0_1[eta][charge].Re());
      c2_trk_accept[eta][charge][1]->Fill(Q_n1_1[eta][charge].Im()/Q_0_1[eta][charge].Re(), Q_0_1[eta][charge].Re());

    }
  }

//numerator D0
  for(int rap = 0; rap < NyBins; rap++){
    for(int charge = 0; charge < 3; charge++){
      
      //obs:
      TComplex N_v1_A_SP, D_v1_A_SP, N_v1_B_SP, D_v1_B_SP, N_v1_AB_SP, D_v1_AB_SP;

      N_v1_A_SP = Q_D0obs_n1_1[rap][charge]*Q_n3_1_HFminus;
      D_v1_A_SP = Q_D0obs_0_1[rap][charge]*Q_0_1_HFminus;

      double V1_A = N_v1_A_SP.Re()/D_v1_A_SP.Re();

      N_v1_B_SP = Q_D0obs_n1_1[rap][charge]*Q_n3_1_HFplus;
      D_v1_B_SP = Q_D0obs_0_1[rap][charge]*Q_0_1_HFplus;

      double V1_B = N_v1_B_SP.Re()/D_v1_B_SP.Re();
     
      //use one plane
      N_v1_AB_SP = Q_D0obs_n1_1[rap][charge]*Q_n3_1_HFplusANDminus;
      D_v1_AB_SP = Q_D0obs_0_1[rap][charge]*Q_0_1_HFplusANDminus;

      double V1_AB = N_v1_AB_SP.Re()/D_v1_AB_SP.Re();

      c2_d0obs_v1[rap][charge][0]->Fill( V1_A, D_v1_A_SP.Re() );
      c2_d0obs_v1[rap][charge][1]->Fill( V1_B, D_v1_B_SP.Re() );
      c2_d0obs_v1[rap][charge][2]->Fill( V1_AB, D_v1_AB_SP.Re() );

      c2_d0obs_trk_accept[rap][charge][0]->Fill(Q_D0obs_n1_1[rap][charge].Re()/Q_D0obs_0_1[rap][charge].Re(), Q_D0obs_0_1[rap][charge].Re());
      c2_d0obs_trk_accept[rap][charge][1]->Fill(Q_D0obs_n1_1[rap][charge].Im()/Q_D0obs_0_1[rap][charge].Re(), Q_D0obs_0_1[rap][charge].Re());

      //bkg:
      N_v1_A_SP = Q_D0bkg_n1_1[rap][charge]*Q_n3_1_HFminus;
      D_v1_A_SP = Q_D0bkg_0_1[rap][charge]*Q_0_1_HFminus;

      V1_A = N_v1_A_SP.Re()/D_v1_A_SP.Re();

      N_v1_B_SP = Q_D0bkg_n1_1[rap][charge]*Q_n3_1_HFplus;
      D_v1_B_SP = Q_D0bkg_0_1[rap][charge]*Q_0_1_HFplus;

      V1_B = N_v1_B_SP.Re()/D_v1_B_SP.Re();

      //use one plane
      N_v1_AB_SP = Q_D0bkg_n1_1[rap][charge]*Q_n3_1_HFplusANDminus;
      D_v1_AB_SP = Q_D0bkg_0_1[rap][charge]*Q_0_1_HFplusANDminus;

      V1_AB = N_v1_AB_SP.Re()/D_v1_AB_SP.Re();

      c2_d0bkg_v1[rap][charge][0]->Fill( V1_A, D_v1_A_SP.Re() );
      c2_d0bkg_v1[rap][charge][1]->Fill( V1_B, D_v1_B_SP.Re() );
      c2_d0bkg_v1[rap][charge][2]->Fill( V1_AB, D_v1_AB_SP.Re() );

      c2_d0bkg_trk_accept[rap][charge][0]->Fill(Q_D0bkg_n1_1[rap][charge].Re()/Q_D0bkg_0_1[rap][charge].Re(), Q_D0bkg_0_1[rap][charge].Re());
      c2_d0bkg_trk_accept[rap][charge][1]->Fill(Q_D0bkg_n1_1[rap][charge].Im()/Q_D0bkg_0_1[rap][charge].Re(), Q_D0bkg_0_1[rap][charge].Re());

      
    }
  }

  //numerator D0 mass-fit method
  for(int rap = 0; rap < NyBins; rap++){
    for(int imass = 0; imass < NmassBins; imass++){
      for(int charge = 0; charge < 3; charge++){
        
        TComplex N_v1_A_SP, D_v1_A_SP, N_v1_B_SP, D_v1_B_SP, N_v1_AB_SP, D_v1_AB_SP;

        N_v1_A_SP = Q_D0_n1_1[rap][imass][charge]*Q_n3_1_HFminus;
        D_v1_A_SP = Q_D0_0_1[rap][imass][charge]*Q_0_1_HFminus;

        double V1_A = N_v1_A_SP.Re()/D_v1_A_SP.Re();

        N_v1_B_SP = Q_D0_n1_1[rap][imass][charge]*Q_n3_1_HFplus;
        D_v1_B_SP = Q_D0_0_1[rap][imass][charge]*Q_0_1_HFplus;

        double V1_B = N_v1_B_SP.Re()/D_v1_B_SP.Re();
       
        //use one plane
        N_v1_AB_SP = Q_D0_n1_1[rap][imass][charge]*Q_n3_1_HFplusANDminus;
        D_v1_AB_SP = Q_D0_0_1[rap][imass][charge]*Q_0_1_HFplusANDminus;

        double V1_AB = N_v1_AB_SP.Re()/D_v1_AB_SP.Re();

        c2_d0_v1[rap][imass][charge][0]->Fill( V1_A, D_v1_A_SP.Re() );
        c2_d0_v1[rap][imass][charge][1]->Fill( V1_B, D_v1_B_SP.Re() );
        c2_d0_v1[rap][imass][charge][2]->Fill( V1_AB, D_v1_AB_SP.Re() );

        c2_d0_trk_accept[rap][imass][charge][0]->Fill(Q_D0_n1_1[rap][imass][charge].Re()/Q_D0_0_1[rap][imass][charge].Re(), Q_D0_0_1[rap][imass][charge].Re());
        c2_d0_trk_accept[rap][imass][charge][1]->Fill(Q_D0_n1_1[rap][imass][charge].Im()/Q_D0_0_1[rap][imass][charge].Re(), Q_D0_0_1[rap][imass][charge].Re());
        
      }
    }
  }


}
// ------------ method called once each job just before starting event loop  ------------
void 
D0DirectedFlowCorrelator::beginJob()
{
  edm::Service<TFileService> fs;
    
  TH1D::SetDefaultSumw2();

  const int NetaBins = etaBins_.size() - 1 ;
  const int NyBins = rapidityBins_.size() - 1;
  const int NmassBins = massBins_.size() - 1;

  double ptBinsArray[100];
  const int Nptbins = ptBins_.size() - 1;
  for(unsigned i = 0; i < ptBins_.size(); i++){
    ptBinsArray[i] = ptBins_[i];
  }

  if( !doPixelReco_ && !doHiReco_ ){
    edm::FileInPath fip1("DirectedFlowCorrelator/DirectedFlowCorrelator/data/Hydjet_eff_mult_v1.root");
    TFile f1(fip1.fullPath().c_str(),"READ");
    for(int i = 0; i < 5; i++){
       effTable[i] = (TH2D*)f1.Get(Form("rTotalEff3D_%d",i));
    }
  }
  if( doHiReco_ && !doPixelReco_ ){
    edm::FileInPath fip3("DirectedFlowCorrelator/DirectedFlowCorrelator/data/PbPb_MB_TT_5TeV_v2.root");
    TFile f3(fip3.fullPath().c_str(),"READ");

    effTable[0] = (TH2D*)f3.Get("rTotalEff3D_50_100");
    effTable[1] = (TH2D*)f3.Get("rTotalEff3D_30_50");
    effTable[2] = (TH2D*)f3.Get("rTotalEff3D_10_30");
    effTable[3] = (TH2D*)f3.Get("rTotalEff3D_5_10");
    effTable[4] = (TH2D*)f3.Get("rTotalEff3D_0_5");

    edm::FileInPath fip4("DirectedFlowCorrelator/DirectedFlowCorrelator/data/D0_eff_table.root");
    TFile f4(fip4.fullPath().c_str(),"READ");

    d0EffTable[0] = (TH1D*)f4.Get("hist_D0_eff_1");
    d0EffTable[1] = (TH1D*)f4.Get("hist_D0_eff_2");
    d0EffTable[2] = (TH1D*)f4.Get("hist_D0_eff_3");

    f1[0] = (TF1*)f4.Get("f1");
    f1[1] = (TF1*)f4.Get("f2");
    f1[2] = (TF1*)f4.Get("f3");

    edm::FileInPath fip5("DirectedFlowCorrelator/DirectedFlowCorrelator/data/D0_eff_table_tight.root");
    TFile f5(fip5.fullPath().c_str(),"READ");

    d0EffTable_tight[0] = (TH1D*)f5.Get("hist_D0_eff_1");
    d0EffTable_tight[1] = (TH1D*)f5.Get("hist_D0_eff_2");
    d0EffTable_tight[2] = (TH1D*)f5.Get("hist_D0_eff_3");

    f2[0] = (TF1*)f5.Get("f1");
    f2[1] = (TF1*)f5.Get("f2");
    f2[2] = (TF1*)f5.Get("f3");

    edm::FileInPath fip7("DirectedFlowCorrelator/DirectedFlowCorrelator/data/D0_eff_table_loose.root");
    TFile f7(fip7.fullPath().c_str(),"READ");

    d0EffTable_loose[0] = (TH1D*)f7.Get("hist_D0_eff_1");
    d0EffTable_loose[1] = (TH1D*)f7.Get("hist_D0_eff_2");
    d0EffTable_loose[2] = (TH1D*)f7.Get("hist_D0_eff_3");

    f33[0] = (TF1*)f7.Get("f1");
    f33[1] = (TF1*)f7.Get("f2");
    f33[2] = (TF1*)f7.Get("f3");

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

  for(int rap = 0; rap < NyBins; rap++){
    for(int charge = 0; charge < 3; charge++){
      for(int dir = 0; dir < 3; dir++){

        c2_d0obs_v1[rap][charge][dir] = fs->make<TH1D>(Form("c2_d0obs_v1_%d_%d_%d",rap,charge,dir),";c1", 1,-1,1);
        c2_d0obs_trk_accept[rap][charge][dir] = fs->make<TH1D>(Form("c2_d0obs_trk_accept_%d_%d_%d",rap,charge,dir), ";c1", 1,-1,1);
        
        c2_d0bkg_v1[rap][charge][dir] = fs->make<TH1D>(Form("c2_d0bkg_v1_%d_%d_%d",rap,charge,dir),";c1", 1,-1,1);
        c2_d0bkg_trk_accept[rap][charge][dir] = fs->make<TH1D>(Form("c2_d0bkg_trk_accept_%d_%d_%d",rap,charge,dir), ";c1", 1,-1,1);

      }
    }
  }

  for(int rap = 0; rap < NyBins; rap++){
    for(int imass = 0; imass < NmassBins; imass++){
      for(int charge = 0; charge < 3; charge++){
        for(int dir = 0; dir < 3; dir++){

          c2_d0_v1[rap][imass][charge][dir] = fs->make<TH1D>(Form("c2_d0_v1_%d_%d_%d_%d",rap,imass,charge,dir),";c1", 1,-1,1);
          c2_d0_trk_accept[rap][imass][charge][dir] = fs->make<TH1D>(Form("c2_d0_trk_accept_%d_%d_%d_%d",rap,imass,charge,dir), ";c1", 1,-1,1);
          
        }
      }
    }
  }

  for(int rap = 0; rap < NyBins; rap++){
    for(int charge = 0; charge < 3; charge++){

      D0Mass_Hist[rap][charge] = fs->make<TH1D>(Form("D0Mass_Hist_%d_%d",rap,charge),";mass",300,1.7,2.0);
      D0Mass_Hist_DeltaPhi[rap][charge] = fs->make<TH2D>(Form("D0Mass_Hist_DeltaPhi_%d_%d",rap,charge),";#Deltaphi;mass",20,-PI,PI,300,1.7,2.0);
    }
  }


  c2_ab = fs->make<TH1D>("c2_ab",";c2_ab", 1,-1,1);
  
  c2_ac_plus = fs->make<TH1D>("c2_ac_plus",";c2_ac_plus", 1,-1,1);
  c2_cb_plus = fs->make<TH1D>("c2_cb_plus",";c2_cb_plus", 1,-1,1);

  c2_ac_minus = fs->make<TH1D>("c2_ac_minus",";c2_ac_minus", 1,-1,1);
  c2_cb_minus = fs->make<TH1D>("c2_cb_minus",";c2_cb_minus", 1,-1,1);

  c2_a_real = fs->make<TH1D>("c2_a_real",";c2_a_real", 1,-1,1);
  c2_b_real = fs->make<TH1D>("c2_b_real",";c2_b_real", 1,-1,1);
  c2_c_plus_real = fs->make<TH1D>("c2_c_plus_real",";c2_c_real", 1,-1,1);
  c2_c_minus_real = fs->make<TH1D>("c2_c_minus_real",";c2_c_real", 1,-1,1);

  c2_a_imag = fs->make<TH1D>("c2_a_imag",";c2_a_imag", 1,-1,1);
  c2_b_imag = fs->make<TH1D>("c2_b_imag",";c2_b_imag", 1,-1,1);
  c2_c_plus_imag = fs->make<TH1D>("c2_c_plus_imag",";c2_c_imag", 1,-1,1);
  c2_c_minus_imag = fs->make<TH1D>("c2_c_minus_imag",";c2_c_imag", 1,-1,1);

  c2_ab_one = fs->make<TH1D>("c2_ab_one",";c2_ab_one", 1,-1,1);

  c2_ab_one_real = fs->make<TH1D>("c2_ab_one_real",";c2_ab_one_real", 1,-1,1);
  c2_ab_one_imag = fs->make<TH1D>("c2_ab_one_imag",";c2_ab_one_imag", 1,-1,1);

  Psi_1_sin = fs->make<TH1D>("Psi_1_sin",";Psi_1_sin", 1,-10000,10000);
  Psi_1_cos = fs->make<TH1D>("Psi_1_cos",";Psi_1_cos", 1,-10000,10000);

}

TComplex 
D0DirectedFlowCorrelator::q_vector(double n, double p, double w, double phi) 
{
  double term1 = pow(w,p);
  TComplex e(1, n*phi, 1);
  return term1*e;
}
// ------------ method called once each job just after ending the event loop  ------------
void 
D0DirectedFlowCorrelator::endJob() 
{
}
void 
D0DirectedFlowCorrelator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
D0DirectedFlowCorrelator::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
D0DirectedFlowCorrelator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
D0DirectedFlowCorrelator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
D0DirectedFlowCorrelator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}

//define this as a plug-in
DEFINE_FWK_MODULE(D0DirectedFlowCorrelator);
