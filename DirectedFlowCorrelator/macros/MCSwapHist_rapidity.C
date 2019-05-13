#include <iostream>

#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "TString.h"
#include "TChain.h"

#include <vector>

void MCSwapHist_rapidity()
{
    
    TChain* D0tree= new TChain("demo/D0para");
    D0tree->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_MC/PromptD0_pt_1/*.root");

    TFile ofile("../systematics/ReconstructionCuts/hMass_MCSwap_reweightZvtx_rapidity_pt1to40_PbPbMC_ext_eta2p4_dau1p5_loose.root","RECREATE");
    
    double D0_pt = 2.0;
    double dau_pt = 1.5;
    // double final_PAngle[] = {0.12,0.08,0.12,0.12,0.08,0.12};
    // double final_VtxProb[] = {0.08,0.0,0.04,0.04,0.0,0.08};
    // double final_DLoS[] = {5.0,5.0,3.5,3.5,5.0,5.0};
    // double final_d0dca[] = {0.008,0.12,0.02,0.02,0.012,0.008};

    //tight cuts
    // double final_PAngle[] = {0.08,0.06,0.08,0.08,0.06,0.08};
    // double final_VtxProb[] = {0.12,0.04,0.08,0.08,0.04,0.12};
    // double final_DLoS[] = {6.0,6.0,5.0,5.0,6.0,6.0};
    // double final_d0dca[] = {0.008,0.12,0.02,0.02,0.012,0.008};

    //loose cuts
    double final_PAngle[] = {0.16,0.12,0.16,0.16,0.12,0.16};
    double final_VtxProb[] = {0.04,0.0,0.0,0.0,0.0,0.04};
    double final_DLoS[] = {3.5,3.5,3.0,3.0,3.5,3.5};
    double final_d0dca[] = {0.008,0.12,0.02,0.02,0.012,0.008};

    float mass;
    float pt;
    float eta;
    float y;
    bool isSwap;
    bool matchGEN;
    float VtxProb;
    float PointingAngle;
    float DecayLengthSignificance;
    int nhit1;
    int nhit2;
    float pt1;
    float pt2;
    float pterr1;
    float pterr2;
    float eta1;
    float eta2;
    float d0dca;
    float Zvtx;
    
    D0tree->SetBranchAddress("mass",&mass);
    D0tree->SetBranchAddress("pT",&pt);
    D0tree->SetBranchAddress("eta",&eta);
    D0tree->SetBranchAddress("y",&y);
    D0tree->SetBranchAddress("isSwap",&isSwap);
    D0tree->SetBranchAddress("matchGEN",&matchGEN);
    D0tree->SetBranchAddress("VtxProb",&VtxProb);
    D0tree->SetBranchAddress("3DPointingAngle",&PointingAngle);
    D0tree->SetBranchAddress("3DDecayLengthSignificance",&DecayLengthSignificance);
    D0tree->SetBranchAddress("NHitD1",&nhit1);
    D0tree->SetBranchAddress("NHitD2",&nhit2);
    D0tree->SetBranchAddress("pTD1",&pt1);
    D0tree->SetBranchAddress("pTD2",&pt2);
    D0tree->SetBranchAddress("pTerrD1",&pterr1);
    D0tree->SetBranchAddress("pTerrD2",&pterr2);
    D0tree->SetBranchAddress("EtaD1",&eta1);
    D0tree->SetBranchAddress("EtaD2",&eta2);
    D0tree->SetBranchAddress("D0DCA",&d0dca);
    D0tree->SetBranchAddress("vtxZ",&Zvtx);
    
    unsigned long long int nentries = D0tree->GetEntries();

    TH1D* hmass[9];
    TH1D* hmass_swap[9];
    TH1D* hmass_all[9];
    double ptbin[18] = {-2.0,-1.2,-0.6,0,0.6,1.2,2.0};
    //double ptbin[10] = {1.4,2,2.6,3.3,4,5,6,7,8,10};
    
    for(int i=0;i<6;i++)
    {
        hmass[i] = new TH1D(Form("mass_rapidity%d",i),Form("mass_rapidity%d",i),60,1.7,2);
        hmass_all[i] = new TH1D(Form("mass_all_rapidity%d",i),Form("mass_all_rapidity%d",i),60,1.7,2);
        hmass_swap[i] = new TH1D(Form("mass_swap_rapidity%d",i),Form("mass_swap_rapidity%d",i),60,1.7,2);
    }
    
    for(unsigned long long int i=0;i<nentries;i++)
    {
        D0tree->GetEntry(i);
        if(i%1000000==0) cout<<i<<" / "<<nentries<<endl;
        if(!matchGEN) continue;
        if(fabs(Zvtx)>15) continue;
        if(pt<D0_pt) continue;
        if(fabs(y)>2.0) continue;
        if(nhit1<11 || nhit2<11) continue;
        if(fabs(pterr1/pt1)>0.1 || fabs(pterr2/pt2)>0.1) continue;
        if(fabs(eta1)>2.4 || fabs(eta2)>2.4) continue;
        if(pt1<dau_pt || pt2<dau_pt) continue;

        for(int j=0;j<6;j++)
        {
            if(y<ptbin[j] || y>ptbin[j+1]) continue;
            if(VtxProb < final_VtxProb[j]) continue;
            if(PointingAngle > final_PAngle[j]) continue;
            if(d0dca > final_d0dca[j]) continue;
            if(DecayLengthSignificance < final_DLoS[j]) continue;

            double weight = 1.0;
            
            if(!isSwap) hmass[j]->Fill(mass,weight);
            if(isSwap) hmass_swap[j]->Fill(mass,weight);
            hmass_all[j]->Fill(mass,weight);
        }
    }
    
    for(int i=0;i<6;i++)
    {
        hmass[i]->Write();
        hmass_swap[i]->Write();
        hmass_all[i]->Write();
    }
    
}
