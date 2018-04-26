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

void D0Hist_PromptNonPrompt_Data()
{
    // TFile* f1 = TFile::Open("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/D0_Ntuple_HI_MB6_skim.root");
    TChain* t= new TChain("D0para");
    t->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/D0_Ntuple_HI_MB6_skim.root");
    t->Add("/Volumes/NormalDrive/D_Ntuple/PbPb_5TeV_data/D0_Ntuple_HI_MB7_skim.root");

    TFile ofile("./D0masshist_MB67_data_PromptNonPrompt_rap2.root","RECREATE");
    
    //TTree* t;
    //f1->GetObject("D0para",t);

    double rap_bins[7] = {-2.0,-1.2,-0.6,0.0,0.6,1.2,2.0};
    //double pt_bins[10] = {2.0,2.5,3.0,4.0,5.0,8.0,12.0,15.0,20.0,30.0}; 
    double pt_bins[2] = {3.0,30.0};
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

    TH1D* masshist[9][14];
    
    for(int iAgl=0;iAgl<1;iAgl++)
    {
        for(int iVP=0;iVP<14;iVP++)
        {
            masshist[iAgl][iVP] = new TH1D(Form("mass_%d_%d",iAgl,iVP),Form("mass_%d_%d",iAgl,iVP),60,1.7,2.0);
 
        }
    }
    
    float mass;
    float pt;
    float eta;
    float rap;
    float VtxProb;
    float PAngle;
    float DLoS;
    float D0DCA;
    float pTD1;
    float pTD2;
    float EtaD1;
    float EtaD2;

    t->SetBranchAddress("mass",&mass);
    t->SetBranchAddress("pT",&pt);
    t->SetBranchAddress("eta", &eta);
    t->SetBranchAddress("y", &rap);
    t->SetBranchAddress("VtxProb",&VtxProb);
    t->SetBranchAddress("3DPointingAngle",&PAngle);
    t->SetBranchAddress("3DDecayLengthSignificance",&DLoS);
    t->SetBranchAddress("d0dca", &D0DCA);
    t->SetBranchAddress("EtaD1",&EtaD1);
    t->SetBranchAddress("EtaD2",&EtaD2);
    t->SetBranchAddress("pTD1",&pTD1);
    t->SetBranchAddress("pTD2",&pTD2);

    unsigned long long int nentries = t->GetEntries();
    cout << "total number of entries: " << nentries << endl;
    unsigned long long int first_nentries = 1000000000;
    
    for (unsigned long long int i=0;i<nentries;i++)
    {
        if(i%500000==0) cout<<"process "<<i<<endl;
        t->GetEntry(i);
        
        if(pt<D0_pt) continue;
        if(pTD1<dau_pt || pTD2<dau_pt) continue;
        if(fabs(EtaD1)>2.2 || fabs(EtaD2)>2.2) continue;
        if( fabs(rap) < 1.2 || fabs(rap) > 2.0  ) continue;
        
        for(int ipt = 0; ipt < 1; ipt++){
            for(int idca = 0; idca < 14; idca++){
                if( D0DCA > dca_bins[idca] && D0DCA < dca_bins[idca+1] ){
                    if( pt > pt_bins[ipt] && pt < pt_bins[ipt+1] ){

                        if(PAngle > final_PAngle[0]) continue;
                        if(VtxProb < final_VtxProb[0]) continue;
                        if(DLoS < final_DLoS[0]) continue;

                        masshist[ipt][idca]->Fill(mass);
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
        }
    }

}
