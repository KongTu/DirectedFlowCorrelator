#include <iostream>
#include <vector>
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "TTree.h"
#include "TChain.h"
#include "TF1.h"

using namespace std;

void GetPromptNonPrompt_D0()
{
    
    TChain* t= new TChain("demo/D0para");
    t->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_MC/PromptD0_pt_1_Centrality/*.root");
    //t->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/PromptD0_PbPb_PtHat_MC/pt0/*.root");
    // t->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/PromptD0_PbPb_PtHat_MC/pt2/*.root");
    // t->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/PromptD0_PbPb_PtHat_MC/pt4/*.root");
    // t->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/PromptD0_PbPb_PtHat_MC/pt10/*.root");
    // t->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/PromptD0_PbPb_PtHat_MC/pt20/*.root");
    // t->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/PromptD0_PbPb_PtHat_MC/pt30/*.root");
    //t->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_MC/EPOS_pPb_embed_D0_prompt_matched_deltaR0p5_pt1p2_eta2p4_OFFICIAL/*.root");

    TChain* t1= new TChain("demo/D0para");
    t1->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_MC/NonPromptD0_pt_1_Centrality/*.root");
    // t1->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/PromptD0_PbPb_PtHat_MC/pt0/*.root");
    // t1->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/PromptD0_PbPb_PtHat_MC/pt2/*.root");
    // t1->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/PromptD0_PbPb_PtHat_MC/pt4/*.root");
    // t1->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/PromptD0_PbPb_PtHat_MC/pt10/*.root");
    // t1->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/PromptD0_PbPb_PtHat_MC/pt20/*.root");
    // t1->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/PromptD0_PbPb_PtHat_MC/pt30/*.root");

    //double pt_bins[10] = {2.0,2.5,3.0,4.0,5.0,8.0,12.0,15.0,20.0,30.0};
    double pt_bins[2] = {3.,30.};
    //double dca_bins[14] = {0.0,0.002,0.004,0.006,0.008,0.012,0.022,0.028,0.0281,0.0367,0.0476,0.07,0.0758,0.0817};
    double dca_bins[15] = {0.0,0.001,0.0023,0.0039,0.0059,0.008,0.0118,0.016,0.0214,0.0281,0.0367,0.0476,0.07,0.0758,0.0817};

    TH1D* hist_prompt[10];
    TH1D* hist_nonprompt[10];
    for(int iy = 0; iy < 1; iy++){
        hist_prompt[iy] = new TH1D(Form("hist_prompt_%d",iy),Form("hist_prompt_%d",iy), 14, dca_bins);
        hist_nonprompt[iy] = new TH1D(Form("hist_nonprompt_%d",iy),Form("hist_nonprompt_%d",iy), 14, dca_bins);
    }

    TFile* fileopen = new TFile("./vtxRatio.root");
    TH1D* weight_vtx = (TH1D*) fileopen->Get("ratio");

    TF1* funCentralityWeight = new TF1("funCentralityWeight", "exp((1.637)+(-0.0212332)*x+(-5.1822e-05)*x*x)", 0, 200);

    double D0_pt = 2.0;
    double ptdau = 1.5;
    
    // double final_PAngle[] = {0.12,0.08,0.12};
    // double final_VtxProb[] = {0.04,0.0,0.08};
    // double final_DLoS[] = {3.5,5.0,5.0};
    
    double final_PAngle[] = {0.12};
    double final_VtxProb[] = {0.08};
    double final_DLoS[] = {5.0};

    float mass;
    float pt;
    float ETT;
    int hiBin;
    float VtxProb;
    float PAngle;
    float DLoS;
    float pTD1;
    float pTD2;
    float vtxZ;
    float y;
    float D0DCA;
    int NHitD1;
    int NHitD2;
    float pTerrD1;
    float pTerrD2;
    float EtaD1;
    float EtaD2;
    bool isSwap;
    bool isPrompt;
    bool matchGEN;
    
    t->SetBranchAddress("mass",&mass);
    t->SetBranchAddress("pT",&pt);
    t->SetBranchAddress("hiBin",&hiBin);
    t->SetBranchAddress("VtxProb",&VtxProb);
    t->SetBranchAddress("3DPointingAngle",&PAngle);
    t->SetBranchAddress("3DDecayLengthSignificance",&DLoS);
    t->SetBranchAddress("D0DCA", &D0DCA);
    t->SetBranchAddress("pTD1",&pTD1);
    t->SetBranchAddress("pTD2",&pTD2);
    t->SetBranchAddress("vtxZ",&vtxZ);
    t->SetBranchAddress("y",&y);
    t->SetBranchAddress("NHitD1",&NHitD1);
    t->SetBranchAddress("NHitD2",&NHitD2);
    t->SetBranchAddress("pTerrD1",&pTerrD1);
    t->SetBranchAddress("pTerrD2",&pTerrD2);
    t->SetBranchAddress("EtaD1",&EtaD1);
    t->SetBranchAddress("EtaD2",&EtaD2);
    t->SetBranchAddress("isSwap",&isSwap);
    t->SetBranchAddress("isPrompt",&isPrompt);
    t->SetBranchAddress("matchGEN",&matchGEN);

    Int_t nentries = t->GetEntries();
    cout << "prompt entry: " << nentries << endl;
    
    for (Int_t i=0;i<nentries;i++)
    {
        if(i%500000==0) cout<<"process "<<i<<endl;
        t->GetEntry(i);
        
        if( pt < D0_pt ) continue;
        if(pTD1<ptdau || pTD2<ptdau) continue;
        if(fabs(vtxZ)>15) continue;
        if(NHitD1<11 || NHitD2<11) continue;
        if(fabs(pTerrD1/pTD1)>0.1 || fabs(pTerrD2/pTD2)>0.1) continue;
        if(fabs(EtaD1)>2.2 || fabs(EtaD2)>2.2) continue;
        if( !matchGEN ) continue;
        if( !isPrompt ) continue;
        if( fabs(y) > 1.2 && fabs(y) < 2.0 ) continue;

        double weight = weight_vtx->FindBin(vtxZ);
        double centrality_weight = funCentralityWeight->Eval(hiBin);
        double total_weight = weight*centrality_weight;

        for(int ipt = 0; ipt < 1; ipt++){
            if( pt > pt_bins[ipt] && pt < pt_bins[ipt+1] ){
                
                if(PAngle > final_PAngle[0]) continue;
                if(VtxProb < final_VtxProb[0]) continue;
                if(DLoS < final_DLoS[0]) continue;

                hist_prompt[ipt]->Fill(D0DCA);
            }
        }
    }
    
    float mass_nonprompt;
    float pt_nonprompt;
    int hiBin_nonprompt;
    float ETT_nonprompt;
    float VtxProb_nonprompt;
    float PAngle_nonprompt;
    float DLoS_nonprompt;
    float pTD1_nonprompt;
    float pTD2_nonprompt;
    float vtxZ_nonprompt;
    float y_nonprompt;
    float D0DCA_nonprompt;
    int NHitD1_nonprompt;
    int NHitD2_nonprompt;
    float pTerrD1_nonprompt;
    float pTerrD2_nonprompt;
    float EtaD1_nonprompt;
    float EtaD2_nonprompt;
    bool isSwap_nonprompt;
    bool isPrompt_nonprompt;
    bool matchGEN_nonprompt;
    
    t1->SetBranchAddress("mass",&mass_nonprompt);
    t1->SetBranchAddress("pT",&pt_nonprompt);
    t1->SetBranchAddress("ETT",&ETT_nonprompt);
    t1->SetBranchAddress("hiBin",&hiBin_nonprompt);
    t1->SetBranchAddress("VtxProb",&VtxProb_nonprompt);
    t1->SetBranchAddress("3DPointingAngle",&PAngle_nonprompt);
    t1->SetBranchAddress("3DDecayLengthSignificance",&DLoS_nonprompt);
    t1->SetBranchAddress("D0DCA", &D0DCA_nonprompt);
    t1->SetBranchAddress("pTD1",&pTD1_nonprompt);
    t1->SetBranchAddress("pTD2",&pTD2_nonprompt);
    t1->SetBranchAddress("vtxZ",&vtxZ_nonprompt);
    t1->SetBranchAddress("y",&y_nonprompt);
    t1->SetBranchAddress("NHitD1",&NHitD1_nonprompt);
    t1->SetBranchAddress("NHitD2",&NHitD2_nonprompt);
    t1->SetBranchAddress("pTerrD1",&pTerrD1_nonprompt);
    t1->SetBranchAddress("pTerrD2",&pTerrD2_nonprompt);
    t1->SetBranchAddress("EtaD1",&EtaD1_nonprompt);
    t1->SetBranchAddress("EtaD2",&EtaD2_nonprompt);
    t1->SetBranchAddress("isSwap",&isSwap_nonprompt);
    t1->SetBranchAddress("isPrompt",&isPrompt_nonprompt);
    t1->SetBranchAddress("matchGEN",&matchGEN_nonprompt);

    Int_t nentries_nonprompt = t1->GetEntries();
   
    cout << "Nonprompt entry: " << nentries_nonprompt << endl;

    for (Int_t i=0;i<nentries_nonprompt;i++)
    {
        if(i%50000==0) cout<<"nonprompt process "<<i<<endl;
        t1->GetEntry(i);
        
        if( pt_nonprompt < D0_pt ) continue;
        if(pTD1_nonprompt<ptdau || pTD2_nonprompt<ptdau) continue;
        if(fabs(vtxZ_nonprompt)>15) continue;
        if(NHitD1_nonprompt<11 || NHitD2_nonprompt<11) continue;
        if(fabs(pTerrD1_nonprompt/pTD1_nonprompt)>0.1 || fabs(pTerrD2_nonprompt/pTD2_nonprompt)>0.1) continue;
        if(fabs(EtaD1_nonprompt)>2.2 || fabs(EtaD2_nonprompt)>2.2) continue;
        if(!matchGEN_nonprompt) continue;
        if(isPrompt_nonprompt) continue;
        if( fabs(y_nonprompt) > 1.2 && fabs(y_nonprompt) < 2.0 ) continue;

        double weight = weight_vtx->FindBin(vtxZ_nonprompt);
        double centrality_weight = funCentralityWeight->Eval(hiBin_nonprompt);
        double total_weight = weight*centrality_weight;

        for(int ipt = 0; ipt < 1; ipt++){
            if( pt_nonprompt> pt_bins[ipt] && pt_nonprompt< pt_bins[ipt+1] ){
                
                if(PAngle_nonprompt > final_PAngle[0]) continue;
                if(VtxProb_nonprompt < final_VtxProb[0]) continue;
                if(DLoS_nonprompt < final_DLoS[0]) continue;

                hist_nonprompt[ipt]->Fill(D0DCA_nonprompt);
            }
        }
    }

    TFile ofile("./D0_promptNonPrompt_MC_rap2.root","RECREATE");

    for(int iy = 0; iy < 1; iy++){
        hist_prompt[iy]->Write();
        hist_nonprompt[iy]->Write();
    }


}
