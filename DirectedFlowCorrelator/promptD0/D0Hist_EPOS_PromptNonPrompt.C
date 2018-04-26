#include <iostream>
#include <vector>
#include "TCanvas.h"
#include "TH1.h"
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

using namespace std;

void D0Hist_EPOS_PromptNonPrompt()
{
    
    TChain* t= new TChain("demo/D0para");
    t->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_MC/PromptD0_pt_1_update/*.root");
    TFile ofile("./D0masshist_EPOS_PromptNonPrompt_rap2.root","RECREATE");
    
    double rap_bins[7] = {-2.0,-1.2,-0.6,0.0,0.6,1.2,2.0};
    //double pt_bins[10] = {2.0,2.5,3.0,4.0,5.0,8.0,12.0,15.0,20.0,30.0}; 
    double pt_bins[2] = {3.,30.};
    //double dca_bins[15] = {0.0,0.002,0.004,0.006,0.008,0.012,0.022,0.028,0.0281,0.0367,0.0476,0.07,0.0758,0.0817, 0.09};
    double dca_bins[15] = {0.0,0.001,0.0023,0.0039,0.0059,0.008,0.0118,0.016,0.0214,0.0281,0.0367,0.0476,0.07,0.0758,0.0817};

    double D0_pt = 2.0;
    double dau_pt = 1.5;
    // double final_PAngle[] = {0.12,0.08,0.12};
    // double final_VtxProb[] = {0.04,0.0,0.08};
    // double final_DLoS[] = {3.5,5.0,5.0};

    double final_PAngle[] = {0.12};
    double final_VtxProb[] = {0.08};
    double final_DLoS[] = {5.0};

    TH1D* masshist[9][18];
    TH1D* masshist_all[9][18];

    for(int iAgl=0;iAgl<1;iAgl++)
    {
        for(int iVP=0;iVP<14;iVP++)
        {
            masshist[iAgl][iVP] = new TH1D(Form("mass_%d_%d_",iAgl,iVP),Form("mass_%d_%d",iAgl,iVP),60,1.7,2.0);
            masshist_all[iAgl][iVP] = new TH1D(Form("mass_all_%d_%d_",iAgl,iVP),Form("mass_all_%d_%d",iAgl,iVP),60,1.7,2.0);

        }
    }
    
    float mass;
    float pt;
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
    bool matchGEN;
    
    t->SetBranchAddress("mass",&mass);
    t->SetBranchAddress("pT",&pt);
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
    t->SetBranchAddress("matchGEN",&matchGEN);
    
    Int_t nentries = t->GetEntries();
    
    for (Int_t i=0;i<nentries;i++)
    {
        if(i%50000==0) cout<<"process "<<i<<endl;
        t->GetEntry(i);
        
        if(pt<D0_pt) continue;
        if(pTD1<dau_pt || pTD2<dau_pt) continue;
        if(fabs(vtxZ)>15) continue;
        if(NHitD1<11 || NHitD2<11) continue;
        if(fabs(pTerrD1/pTD1)>0.1 || fabs(pTerrD2/pTD2)>0.1) continue;
        if(fabs(EtaD1)>2.2 || fabs(EtaD2)>2.2) continue;
        if(!matchGEN) continue;
        if( fabs(y) < 1.2 || fabs(y) > 2.0  ) continue;        
        
        for(int ipt = 0; ipt < 1; ipt++){
            for(int idca = 0; idca < 14; idca++){
                if( D0DCA > dca_bins[idca] && D0DCA < dca_bins[idca+1] ){      
                    if( pt > pt_bins[ipt] && pt < pt_bins[ipt+1] ){            

                        if(PAngle > final_PAngle[0]) continue;
                        if(VtxProb < final_VtxProb[0]) continue;
                        if(DLoS < final_DLoS[0]) continue;

                        if( !isSwap ) masshist[ipt][idca]->Fill(mass);
                        masshist_all[ipt][idca]->Fill(mass);

                    }
                }
            }
        }
    }
    
    for(int iAgl=0;iAgl<1;iAgl++)
    {
        for(int iVP=0;iVP<14;iVP++)
        {
            masshist[iAgl][iVP]->Write();
            masshist_all[iAgl][iVP]->Write();
        }
    }




}
