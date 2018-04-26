//#include "makeMultiPanelCanvas.C"
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
#include "TNtuple.h"
#include "TFitter.h"
#include "TFitResult.h"

#include <vector>

TH1D* hD0DcaMCPSignal;
TH1D* hD0DcaMCNPSignal;
TH1D* hD0DcaData;

Double_t funMix(Double_t* x_, Double_t* para)
{
    float x = x_[0];
    float APrompt = para[0];
    float ANonPrompt = para[1];
    float promptYield = 0;
    float nonPromptYield = 0;
    
    promptYield = hD0DcaMCPSignal->GetBinContent(hD0DcaMCPSignal->GetXaxis()->FindBin(x));
    nonPromptYield = hD0DcaMCNPSignal->GetBinContent(hD0DcaMCNPSignal->GetXaxis()->FindBin(x));
    
    return APrompt*promptYield+ANonPrompt*nonPromptYield;
}

Double_t funNonPrompt(Double_t* x_, Double_t* para)
{
    float x = x_[0];
    float APrompt = para[0];
    float ANonPrompt = para[1];
    float nonPromptYield = 0;
    nonPromptYield = hD0DcaMCNPSignal->GetBinContent(hD0DcaMCNPSignal->GetXaxis()->FindBin(x));
    return ANonPrompt*nonPromptYield;
}

Double_t funPrompt(Double_t* x_, Double_t* para)
{
    float x = x_[0];
    float APrompt = para[0];
    float ANonPrompt = para[1];
    float nonPromptYield = 0;
    float promptYield = hD0DcaMCPSignal->GetBinContent(hD0DcaMCPSignal->GetXaxis()->FindBin(x));
    //nonPromptYield = hD0DcaMCNPSignal->GetBinContent(hD0DcaMCNPSignal->GetXaxis()->FindBin(x));
    return APrompt*promptYield;
}

double dca_bins[8] = {0.001,0.0023,0.0039,0.0059,0.008,0.0118,0.016,0.0214};

void fitDCA_largeDCA()
{
    double promptCut = 1.0;
    
    TFile* file0 = TFile::Open("bFeedDownPbPbMBMC.hist.root");
    //TFile* file1 = TFile::Open("Yield_DCA_pPb_pT_EPOS_nonprompt_droplowNbin_reweightZvtx.root");
    
    TFile* file2[8];
    file2[0] = TFile::Open("D0_DCA_data_rap0.root");
    file2[1] = TFile::Open("D0_DCA_data_rap1.root");
    file2[2] = TFile::Open("D0_DCA_data_rap2.root");
    
    TCanvas* c[3];
    TF1* fMix[3];
    //TF1* fNP[3];
    
    double dca_bins_hist[9] = {0.001,0.0023,0.0039,0.0059,0.008,0.0118,0.016,0.0214,0.03};
    double ptbin[] = {-2.0,-1.2,-0.6,0.0,0.6,1.2,2.0};
    double PromptFrac[8];
    double PromptFracErr[8];
    
    double PromptFrac_cut[8];
    double PromptFracErr_cut[8];
    
    double cutfraction[8];
    TPad* pad[10][2];
    TLine* Line = new TLine(0,1,0.036,1);
    Line->SetLineStyle(7);
    
    TGraphErrors* ratio[10];
    
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);

    TH1D* hist = new TH1D("hist","",100,0,0.0367);
    hist->SetLineWidth(0);
    hist->SetStats(kFALSE);
    hist->GetXaxis()->SetTitle("D^{0} DCA (cm)");
    hist->GetYaxis()->SetTitle("data/fit");
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleOffset(1);
    hist->GetYaxis()->SetTitleOffset(1);
    hist->GetXaxis()->SetLabelOffset(0.007);
    hist->GetYaxis()->SetLabelOffset(0.007);
    hist->GetXaxis()->SetTitleSize(0.075);
    hist->GetYaxis()->SetTitleSize(0.075);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->SetMinimum(0);
    hist->SetMaximum(2);
    hist->GetYaxis()->SetRangeUser(0,2.5);
    hist->GetXaxis()->SetRangeUser(0,0.0367);
    
    for(int i=0;i<1;i++)
    {
        c[i] = new TCanvas(Form("c_%d",i),Form("c_%d",i),400,500);
        c[i]->Range(0,0,1,1);
        
        pad[i][0] = new TPad(Form("pad0_%d",i),Form("pad0_%d",i),0.0,0.346,1,1);
        pad[i][1] = new TPad(Form("pad1_%d",i),Form("pad1_%d",i),0.0,0.0,1,0.346);

        pad[i][0]->SetLeftMargin(0.18);
        pad[i][0]->SetRightMargin(0.043);
        pad[i][0]->SetTopMargin(0.06);
        pad[i][0]->SetBottomMargin(0);
        
        pad[i][1]->SetLeftMargin(0.18);
        pad[i][1]->SetRightMargin(0.043);
        pad[i][1]->SetTopMargin(0);
        pad[i][1]->SetBottomMargin(0.145);
        
        pad[i][0]->Draw();
        pad[i][1]->Draw();
        
        pad[i][0]->SetLogy();
        pad[i][0]->cd();
        

        TH3D* hMCPSignal = (TH3D*)file0->Get("hMCPSignal");
        TH3D* hMCNPSignal = (TH3D*)file0->Get("hMCNPSignal");

        hD0DcaMCPSignal = (TH1D*)hMCPSignal->ProjectionY("hD0DcaMCPSignal",1,1000,1,1000);
        hD0DcaMCNPSignal = (TH1D*)hMCNPSignal->ProjectionY("hD0DcaMCNPSignal",1,1000,1,1000);
        
        //TH2D* data_2d = (TH2D*)file2[i]->Get("hist_Data");
        //hD0DcaData = (TH1D*)data_2d->ProjectionY(Form("DCA_data_%d",i),i+1,i+1);
        hD0DcaData = (TH1D*) file2[i]->Get("hist_Data");


        hD0DcaData->SetTitle("");
        hD0DcaData->GetXaxis()->SetRangeUser(0,0.03);
        hD0DcaData->GetYaxis()->SetTitle("dN / d(D^{0} DCA) (cm^{-1})");
        hD0DcaData->GetXaxis()->SetTitle("D^{0} DCA (cm)");
        hD0DcaData->GetXaxis()->CenterTitle();
        hD0DcaData->GetYaxis()->CenterTitle();
        hD0DcaData->SetMinimum(2e1);
        hD0DcaData->SetMaximum(4e5);
        
        hD0DcaData->SetMarkerStyle(20);
        hD0DcaData->Draw("PE");


        //hD0DcaMCPSignal->Draw("same");
        //hD0DcaMCNPSignal->Draw("");

        double integralTotalYield = hD0DcaData->Integral(1,hD0DcaData->GetXaxis()->GetNbins());
        cout<<"Data Total yield: "<<integralTotalYield<<endl;
        fMix[i] = new TF1(Form("fMix_%d",i),&funMix, -0.5, 0.5, 2);
        fMix[i]->SetParameters(0.5*integralTotalYield,0.5*integralTotalYield);
        fMix[i]->SetParLimits(0,0,0.5*integralTotalYield);
        fMix[i]->SetParLimits(1,0,0.5*integralTotalYield);
  
        fMix[i]->SetLineColor(2);
        fMix[i]->SetFillColorAlpha(kRed,0.36);
        fMix[i]->SetFillStyle(1001);
 
        //Set fit range
        float fitRangeL = 0.000;
        float fitRangeH = 0.025;//0.036

        int fitStatus = 1;
        TFitResultPtr fitResult;
        double fitPrecision = 1.e-10;
        while(fitStatus)
        {
            TFitter::SetPrecision(fitPrecision);
            fMix[i]->SetParameters(0.5*integralTotalYield,0.5*integralTotalYield);
            fMix[i]->SetParError(0,0.5*integralTotalYield);
            fMix[i]->SetParError(1,0.5*integralTotalYield);
            fitResult = hD0DcaData->Fit(Form("fMix_%d",i),"E SNQ0", "", fitRangeL, fitRangeH);
            fitStatus = fitResult->Status();
            cout<<"fit precision: "<<TFitter::GetPrecision()<<"   status: "<<fitStatus<<endl;
            if(fitStatus)
                fitPrecision *= 5;
        }
        
        cout<<"============== do main fit ============"<<endl;
        fMix[i]->SetParameters(0.5*integralTotalYield,0.5*integralTotalYield);
        fMix[i]->SetParError(0,0.5*integralTotalYield);
        fMix[i]->SetParError(1,0.5*integralTotalYield);
        fMix[i]->SetNpx(10000);//10000
        fitResult = hD0DcaData->Fit(Form("fMix_%d",i),"E S0", "", fitRangeL, fitRangeH);
        fMix[i]->SetRange(0,fitRangeH);
        fMix[i]->Draw("same");
        //hD0DcaData->GetFunction(Form("fMix_%d",i))->Draw("flsame");
        fitStatus = fitResult->Status();
        cout<<"fit precision: "<<TFitter::GetPrecision()<<"   status: "<<fitStatus<<endl;
        
        TF1* fNP = new TF1("fNP",&funNonPrompt, 0., 0.5, 2);
        fNP->SetParameters(fMix[i]->GetParameter(0),fMix[i]->GetParameter(1));
        fNP->SetRange(0,0.35);
        fNP->SetLineColor(4);
        fNP->SetFillStyle(1001);
        fNP->SetFillColorAlpha(kBlue-6,0.8);
        fNP->SetNpx(10000);
        fNP->Draw("same");
        
        TF1* fP = new TF1("fP",&funPrompt, 0., 0.5, 2);
        fP->SetParameters(fMix[i]->GetParameter(0),fMix[i]->GetParameter(1));
        fP->SetRange(0,0.35);
        fP->SetLineColor(2);
        fP->SetFillColorAlpha(kRed,0.36);
        fP->SetFillStyle(1001);
        
        hD0DcaData->Draw("PEsame");

        double promptDYield = fMix[i]->GetParameter(0);
        double promptDYieldErr = fMix[i]->GetParError(0);
        double nonpromptDYield = fMix[i]->GetParameter(1);
        double nonpromptDYieldErr = fMix[i]->GetParError(1);
        
        PromptFrac[i] = promptDYield/(promptDYield+nonpromptDYield);
        PromptFracErr[i] = promptDYieldErr/(promptDYield+nonpromptDYield);
            
        for(int j = 0; j < 8; j++){

            promptCut = dca_bins[j];
            double promptDYield_cut = fMix[i]->GetParameter(0)*(fP->Integral(0,promptCut)/fP->Integral(0,fitRangeH));
            double promptDYieldErr_cut = fMix[i]->GetParError(0)*(fP->Integral(0,promptCut)/fP->Integral(0,fitRangeH));
            double nonpromptDYield_cut = fMix[i]->GetParameter(1)*(fNP->Integral(0,promptCut)/fNP->Integral(0,fitRangeH));
            double nonpromptDYieldErr_cut = fMix[i]->GetParError(1)*(fNP->Integral(0,promptCut)/fNP->Integral(0,fitRangeH));

            cutfraction[j] = fMix[i]->Integral(0,promptCut)/fMix[i]->Integral(0,fitRangeH);
            
            PromptFrac_cut[j] = promptDYield_cut/(promptDYield_cut+nonpromptDYield_cut);
            PromptFracErr_cut[j] = promptDYieldErr_cut/(promptDYield_cut+nonpromptDYield_cut);
        }

        

        
        TLegend* leg4 = new TLegend(0.56,0.56,0.84,0.79);
        leg4->SetBorderSize(0);
        leg4->SetTextSize(0.045);
        leg4->SetTextFont(42);
        leg4->SetFillStyle(0);
        leg4->AddEntry(hD0DcaData,"Data","pl");
        leg4->AddEntry(fMix[i],"Prompt D^{0}","f");
        leg4->AddEntry(fNP,"B to D^{0}","f");
        leg4->Draw("same");
        
        tex->DrawLatex(0.4,0.88,Form("0 < |y| < 0.6, 2.0 < p_{T} < 30.0 GeV/c",ptbin[i],ptbin[i+1]));
        tex->DrawLatex(0.4,0.8,Form("Prompt frac. = %.2f #pm %.2f",PromptFrac[i],PromptFracErr[i]));
        
        double max = hD0DcaData->GetMaximum();
        TLine* l = new TLine(promptCut,0,promptCut,max);
        l->SetLineStyle(7);
        //l->Draw("LSAME");
        
        pad[i][1]->cd();
        double DCAx[10];
        double ratioFit[10];
        double ratioFitErr[10];
        float Chi2 =0;
        double pullDCA[10];
        int ndf = 9-1;

        for(int k=0;k<10;k++)
        {
            DCAx[k] = hD0DcaData->GetBinCenter(k+1);
            ratioFit[k] = hD0DcaData->GetBinContent(k+1)/fMix[i]->Eval(DCAx[k]);
            ratioFitErr[k] = hD0DcaData->GetBinError(k+1)/fMix[i]->Eval(DCAx[k]);
            
            pullDCA[k] = (hD0DcaData->GetBinContent(k+1) - fMix[i]->Eval(DCAx[k]))/hD0DcaData->GetBinError(k+1);
            Chi2 = Chi2+pullDCA[k]*pullDCA[k];
        }
        
        ratio[i] = new TGraphErrors(10,DCAx,ratioFit,0,ratioFitErr);
        
        hist->Draw();
         ratio[i]->SetMarkerStyle(20);
        ratio[i]->Draw("PESAME");
        Line->Draw("LSAME");
        
        pad[i][0]->cd();
        //tex->DrawLatex(0.4,0.5,Form("Chi2/ndf=%.0f/%d",Chi2,ndf));
        
        c[i]->Print(Form("plots/largeDCA_promptfit_%d_ratio.pdf",i));
        
        TFile ofile("prompt_frac_test.root","RECREATE");
        
        TH1D* hPromptFrac = new TH1D("PromptFrac","PromptFrac",8,dca_bins_hist);
        TH1D* hPromptFrac_cut = new TH1D("PromptFrac_cut","PromptFrac_cut",8,dca_bins_hist);
        
        for(int j=0;j<1;j++)
        {
            hPromptFrac->SetBinContent(j+1,PromptFrac[j]);
            hPromptFrac->SetBinError(j+1,PromptFracErr[j]);
        }
        for(int j=0;j<8;j++)
        {
            hPromptFrac_cut->SetBinContent(j+1,PromptFrac_cut[j]);
            hPromptFrac_cut->SetBinError(j+1,PromptFracErr_cut[j]);
        }
        
        hPromptFrac->Write();
        hPromptFrac_cut->Write();

    }
    
    for(int i=0;i<1;i++)
    {
        cout<<"pt"<<i<<", left fraction "<<cutfraction[i]<<endl;
    }
    //c->Print("plots/EPOSDCA/DCA_promptfit.pdf");
}

