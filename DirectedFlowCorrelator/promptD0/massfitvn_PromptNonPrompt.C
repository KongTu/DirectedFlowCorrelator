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
#include "TGaxis.h"
#include "TString.h"

#include <vector>

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"

void massfitvn_PromptNonPrompt()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    
    double fit_range_low = 1.7;
    double fit_range_high = 2.0;
    double D0_mass = 1.8648;
    double histmassbinsize = (2.0-1.7)/60.;
    TFile* file0 = TFile::Open("./D0masshist_EPOS_PromptNonPrompt_rap2.root");
    TFile* file1 = TFile::Open("./D0masshist_MB67_data_PromptNonPrompt_rap2.root");
    
    TFile ofile("D0_DCA_data_rap2.root","RECREATE");
    
    double rap_bins[4] = {0.0, 0.6, 1.2, 2.0};
    //double pt_bins[10] = {2.0,2.5,3.0,4.0,5.0,8.0,12.0,15.0,20.0,30.0}; 
    double pt_bins[2] = {3.,30.0}; 
    double dca_bins[15] = {0.0,0.001,0.0023,0.0039,0.0059,0.008,0.0118,0.016,0.0214,0.0281,0.0367,0.0476,0.07,0.0758,0.0817};
    //double dca_bins[15] = {0.0,0.002,0.004,0.006,0.008,0.012,0.022,0.028,0.0281,0.0367,0.0476,0.07,0.0758,0.0817, 0.09};
    
    double dca_binwidths[14];
    for(int i = 0; i < 14; i++){
        dca_binwidths[i] = dca_bins[i+1]-dca_bins[i];
        cout << "width " << dca_binwidths[i] << endl;
    }

    // TH2D* hist_Data = new TH2D("hist_Data","hist_Data",3, rap_bins, 14, dca_bins);
    //TH2D* hist_Data = new TH2D("hist_Data","hist_Data",9, pt_bins, 14, dca_bins);

    TH1D* hist_Data = new TH1D("hist_Data","hist_Data", 14, dca_bins);


    TCanvas* c = new TCanvas("c","c",400,400);
        c->cd(1)->SetTopMargin(0.06);
        c->cd(1)->SetLeftMargin(0.18);
        c->cd(1)->SetRightMargin(0.043);
        c->cd(1)->SetBottomMargin(0.145);
        c->cd(2)->SetTopMargin(0.06);
        c->cd(2)->SetLeftMargin(0.18);
        c->cd(2)->SetRightMargin(0.043);
        c->cd(2)->SetBottomMargin(0.145);
    
    TLatex* tex = new TLatex;
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.045);
    tex->SetLineWidth(2);
 
    TLatex* texCMS = new TLatex;
    texCMS->SetNDC();
    texCMS->SetTextFont(42);
    texCMS->SetTextSize(0.05);
    texCMS->SetTextAlign(12);
    
    double sigsig[9][18];
    double sigyield[9][18];
    double totyield[9][18];
    double sigwidth[9][18];
    double error[9][18];

    TF1* f[9][18];

    for(int iAgl = 0; iAgl<1; iAgl++){
        for(int iVP = 0; iVP < 14; iVP++){
            
            f[iAgl][iVP] = new TF1(Form("f_%d_%d",iAgl,iVP),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + [9] + [10]*x + [11]*x*x + [12]*x*x*x", fit_range_low, fit_range_high);

        }
    }

    for(int iAgl=0;iAgl<1;iAgl++)
    {
        for(int iVP=0;iVP<10;iVP++)
        {
           
            TH1D* h_mc_match_signal = (TH1D*)file0->Get(Form("mass_%d_%d_",iAgl,iVP));
            TH1D* h_mc_match_all = (TH1D*)file0->Get(Form("mass_all_%d_%d_",iAgl,iVP));
            
            TH1D* h_data = (TH1D*)file1->Get(Form("mass_%d_%d",iAgl,iVP));
            h_data->SetMinimum(0);
            h_data->SetMarkerSize(0.8);
            h_data->SetMarkerStyle(20);
            h_data->SetLineWidth(1);
            h_data->SetOption("e");
            h_data->GetXaxis()->SetRangeUser(1.7,2);
            h_data->GetXaxis()->SetTitle("m_{#piK} (GeV/c^{2})");
            h_data->GetYaxis()->SetTitle("Entries / 5 MeV");
            h_data->GetXaxis()->CenterTitle();
            h_data->GetYaxis()->CenterTitle();
            h_data->GetXaxis()->SetTitleOffset(1.3);
            h_data->GetYaxis()->SetTitleOffset(2);
            h_data->GetXaxis()->SetLabelOffset(0.007);
            h_data->GetYaxis()->SetLabelOffset(0.007);
            h_data->GetXaxis()->SetTitleSize(0.045);
            h_data->GetYaxis()->SetTitleSize(0.045);
            h_data->GetXaxis()->SetTitleFont(42);
            h_data->GetYaxis()->SetTitleFont(42);
            h_data->GetXaxis()->SetLabelFont(42);
            h_data->GetYaxis()->SetLabelFont(42);
            h_data->GetXaxis()->SetLabelSize(0.04);
            h_data->GetYaxis()->SetLabelSize(0.04);
            
            h_data->GetXaxis()->SetNoExponent(true);
            ((TGaxis*)h_data->GetXaxis())->SetMaxDigits(7);
            
            h_data->SetMaximum(h_data->GetMaximum()*2.0);
            
            c->cd(1);
            /*The full fitting function is constructed as follow
             [0] is signal + swap yield;
             [1] is common mean of double gaussian;
             [2] is signal gaussian 1 sigma;
             [3] is signal gaussian 2 sigma;
             [4] is fractional signal gaussian 1 yield; 1-[4] is fractional signal gaussian 2 yield;
             [5] is fractional double gaussian signal yield, 1-[5] is fractional swap yield;
             [6] is a factor to let width of the gaussians to vary in data;
             [7] is swap gaussian sigma;
             [8] is swap gaussian mean;
             [9-12] is 3rd order poly parameters
             */
            
            // f[iAgl][iVP] = new TF1(Form("f_%d_%d",iAgl,iVP),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + [9] + [10]*x + [11]*x*x + [12]*x*x*x", fit_range_low, fit_range_high);
            
            f[iAgl][iVP]->SetLineColor(2);
            f[iAgl][iVP]->SetLineWidth(1);
            
            //first fit MC signal, swap and poly bkg set to 0
            
            f[iAgl][iVP]->SetParameter(0,1.);
            f[iAgl][iVP]->SetParameter(1,D0_mass);
            f[iAgl][iVP]->SetParameter(2,0.03);//0.03
            f[iAgl][iVP]->SetParameter(3,0.005);//0.005
            f[iAgl][iVP]->SetParameter(4,0.1);//0.1
            f[iAgl][iVP]->FixParameter(5,1);//1
            f[iAgl][iVP]->FixParameter(6,0.5); //always 0 in MC
            if(iVP < 6) f[iAgl][iVP]->FixParameter(7,0.0); //does not really mater here as yield is fix to 0
            else f[iAgl][iVP]->FixParameter(7,0.1);
            f[iAgl][iVP]->FixParameter(8,D0_mass); //does not really mater here as yield is fix to 0
            f[iAgl][iVP]->FixParameter(9,0.0);//0
            f[iAgl][iVP]->FixParameter(10,0.0);//0
            f[iAgl][iVP]->FixParameter(11,0.0);//0
            f[iAgl][iVP]->FixParameter(12,0.0);//0
            
            f[iAgl][iVP]->SetParLimits(2,0.01,0.1);
            f[iAgl][iVP]->SetParLimits(3,0.001,0.05);
            f[iAgl][iVP]->SetParLimits(4,0,1);
            f[iAgl][iVP]->SetParLimits(5,0,1);
            
            f[iAgl][iVP]->FixParameter(1,1.8648); //for first few attempt fix mean of gaussian to get reasonable estimation of other pars; later open it up
            h_mc_match_signal->Fit(Form("f_%d_%d",iAgl,iVP),"q","",fit_range_low,fit_range_high);
            h_mc_match_signal->Fit(Form("f_%d_%d",iAgl,iVP),"q","",fit_range_low,fit_range_high);
            f[iAgl][iVP]->ReleaseParameter(1); //now let gaussian mean float
            h_mc_match_signal->Fit(Form("f_%d_%d",iAgl,iVP),"L q","",fit_range_low,fit_range_high);
            h_mc_match_signal->Fit(Form("f_%d_%d",iAgl,iVP),"L q","",fit_range_low,fit_range_high);
            h_mc_match_signal->Fit(Form("f_%d_%d",iAgl,iVP),"L m","",fit_range_low,fit_range_high);
            
            //now fix signal double gaussian mean, sigma and gaus1,gaus2 yield ratio
            f[iAgl][iVP]->FixParameter(1,f[iAgl][iVP]->GetParameter(1));
            f[iAgl][iVP]->FixParameter(2,f[iAgl][iVP]->GetParameter(2));
            f[iAgl][iVP]->FixParameter(3,f[iAgl][iVP]->GetParameter(3));
            f[iAgl][iVP]->FixParameter(4,f[iAgl][iVP]->GetParameter(4));
            
            //now release swap bkg parameters to fit signal+swap MC
            f[iAgl][iVP]->ReleaseParameter(5);
            f[iAgl][iVP]->ReleaseParameter(7);
            f[iAgl][iVP]->ReleaseParameter(8);
            
            f[iAgl][iVP]->SetParameter(7,0.1);
            f[iAgl][iVP]->SetParameter(8,D0_mass);
            
            //fit signal+swap MC
            h_mc_match_all->Fit(Form("f_%d_%d",iAgl,iVP),"L q","",fit_range_low,fit_range_high);
            h_mc_match_all->Fit(Form("f_%d_%d",iAgl,iVP),"L q","",fit_range_low,fit_range_high);
            h_mc_match_all->Fit(Form("f_%d_%d",iAgl,iVP),"L q","",fit_range_low,fit_range_high);
            h_mc_match_all->Fit(Form("f_%d_%d",iAgl,iVP),"L q","",fit_range_low,fit_range_high);
            h_mc_match_all->Fit(Form("f_%d_%d",iAgl,iVP),"L m","",fit_range_low,fit_range_high);
            
            //now fix swap bkg parameters to fit data
            f[iAgl][iVP]->FixParameter(5,f[iAgl][iVP]->GetParameter(5));
            f[iAgl][iVP]->FixParameter(7,f[iAgl][iVP]->GetParameter(7));
            f[iAgl][iVP]->FixParameter(8,f[iAgl][iVP]->GetParameter(8));
            
            //now release poly bkg pars
            f[iAgl][iVP]->ReleaseParameter(9);
            f[iAgl][iVP]->ReleaseParameter(10);
            f[iAgl][iVP]->ReleaseParameter(11);
            f[iAgl][iVP]->ReleaseParameter(12);
            
            //now fit data
            h_data->Fit(Form("f_%d_%d",iAgl,iVP),"q","",fit_range_low,fit_range_high);
            h_data->Fit(Form("f_%d_%d",iAgl,iVP),"q","",fit_range_low,fit_range_high);
            f[iAgl][iVP]->ReleaseParameter(1); //allow data to have different mass peak mean than MC
            f[iAgl][iVP]->ReleaseParameter(6); //allow data to have different peak width than MC
            f[iAgl][iVP]->SetParameter(6,0);
            f[iAgl][iVP]->SetParLimits(6,-1,1);
            h_data->Fit(Form("f_%d_%d",iAgl,iVP),"L q","",fit_range_low,fit_range_high);
            h_data->Fit(Form("f_%d_%d",iAgl,iVP),"L q","",fit_range_low,fit_range_high);
            h_data->Fit(Form("f_%d_%d",iAgl,iVP),"L q","",fit_range_low,fit_range_high);
            h_data->Fit(Form("f_%d_%d",iAgl,iVP),"L m","",fit_range_low,fit_range_high);
            
            //draw D0 signal separately
            TF1* f1 = new TF1(Form("f_sig_%d_%d",iAgl,iVP),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
            f1->SetLineColor(kOrange-3);
            f1->SetLineWidth(1);
            f1->SetLineStyle(2);
            f1->SetFillColorAlpha(kOrange-3,0.3);
            f1->SetFillStyle(1001);
            f1->FixParameter(0,f[iAgl][iVP]->GetParameter(0));
            f1->FixParameter(1,f[iAgl][iVP]->GetParameter(1));
            f1->FixParameter(2,f[iAgl][iVP]->GetParameter(2));
            f1->FixParameter(3,f[iAgl][iVP]->GetParameter(3));
            f1->FixParameter(4,f[iAgl][iVP]->GetParameter(4));
            f1->FixParameter(5,f[iAgl][iVP]->GetParameter(5));
            f1->FixParameter(6,f[iAgl][iVP]->GetParameter(6));
            f1->Draw("LSAME");
            
            //draw swap bkg separately
            TF1* f2 = new TF1(Form("f_swap_%d_%d",iAgl,iVP),"[0]*((1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])))", fit_range_low, fit_range_high);
            f2->SetLineColor(kGreen+4);
            f2->SetLineWidth(1);
            f2->SetLineStyle(1);
            f2->SetFillColorAlpha(kGreen+4,0.3);
            f2->SetFillStyle(1001);
            f2->FixParameter(0,f[iAgl][iVP]->GetParameter(0));
            f2->FixParameter(5,f[iAgl][iVP]->GetParameter(5));
            f2->FixParameter(6,f[iAgl][iVP]->GetParameter(6));
            f2->FixParameter(7,f[iAgl][iVP]->GetParameter(7));
            f2->FixParameter(8,f[iAgl][iVP]->GetParameter(8));
            f2->Draw("LSAME");
            
            //draw poly bkg separately
            TF1* f3 = new TF1(Form("f_bkg_%d_%d",iAgl,iVP),"[9] + [10]*x + [11]*x*x + [12]*x*x*x", fit_range_low, fit_range_high);
            f3->SetLineColor(4);
            f3->SetLineWidth(1);
            f3->SetLineStyle(2);
            f3->FixParameter(9,f[iAgl][iVP]->GetParameter(9));
            f3->FixParameter(10,f[iAgl][iVP]->GetParameter(10));
            f3->FixParameter(11,f[iAgl][iVP]->GetParameter(11));
            f3->FixParameter(12,f[iAgl][iVP]->GetParameter(12));
            f3->Draw("LSAME");
            
            tex->DrawLatex(0.22,0.86,"Cent.30-80%");
            tex->DrawLatex(0.22,0.80,Form("%.1f < |y| < %.1f",pt_bins[iAgl],pt_bins[iAgl+1]));
            tex->DrawLatex(0.22,0.74,"p_{T,D^{0}} > 1.0");
            
            texCMS->DrawLatex(.18,.97,"#font[61]{CMS} #it{Preliminary}");
            texCMS->DrawLatex(0.62,0.97, "#scale[0.8]{PbPb #sqrt{s_{NN}} = 5.02 TeV}");
            
            TLegend* leg = new TLegend(0.65,0.58,0.81,0.9,NULL,"brNDC");
            leg->SetBorderSize(0);
            leg->SetTextSize(0.045);
            leg->SetTextFont(42);
            leg->SetFillStyle(0);
            leg->AddEntry(h_data,"data","p");
            leg->AddEntry(f[iAgl][iVP],"Fit","L");
            leg->AddEntry(f1,"D^{0}+#bar{D^{#lower[0.2]{0}}} Signal","f");
            leg->AddEntry(f2,"K-#pi swap","f");
            leg->AddEntry(f3,"Combinatorial","l");
            leg->Draw("SAME");
            
            double sigma1 = f[iAgl][iVP]->GetParameter(2);
            double sigma2 = f[iAgl][iVP]->GetParameter(3);
            double sig1frac = f[iAgl][iVP]->GetParameter(4);
            double sigmatot = sqrt(sig1frac*sigma1*sigma1 + (1-sig1frac)*sigma2*sigma2);
            double mean = f[iAgl][iVP]->GetParameter(1);
            
            //sigma[iAgl][iVP][iDL][idau][ipt] = sqrt(sig1frac*sigma1*sigma1 + (1-sig1frac)*sigma2*sigma2);
            
            sigwidth[iAgl][iVP] = sigmatot;
            sigyield[iAgl][iVP] = (f1->Integral(mean-3*sigmatot,mean+3*sigmatot)+f2->Integral(mean-3*sigmatot,mean+3*sigmatot))/histmassbinsize;
            totyield[iAgl][iVP] = (f1->Integral(mean-3*sigmatot,mean+3*sigmatot)+f2->Integral(mean-3*sigmatot,mean+3*sigmatot)+f3->Integral(mean-3*sigmatot,mean+3*sigmatot))/histmassbinsize;
            
            error[iAgl][iVP] = f[iAgl][iVP]->GetParError(0);
            //if( error[iAgl][iVP] > 100 ) error[iAgl][iVP] = 0.1*sigyield[iAgl][iVP];

            cout << "error " << iAgl << " " << iVP << " " << error[iAgl][iVP]/dca_binwidths[iVP] << endl;            
            //cout << "yield " << iAgl << " " << iVP << " " << sigyield[iAgl][iVP]/dca_binwidths[iVP] << endl;
            cout << "yield " << iAgl << " " << iVP << " " << (f[iAgl][iVP]->GetParameter(0))/dca_binwidths[iVP] << endl;

            sigsig[iAgl][iVP] = sigyield[iAgl][iVP]/sqrt(totyield[iAgl][iVP]);
            
            // if(iVP == 4 ){
            //     hist_Data->SetBinContent(iAgl+1, iVP+1, sigyield[iAgl][iVP]/dca_binwidths[iVP]);
            //     hist_Data->SetBinError(iAgl+1, iVP+1, 0.1*sigyield[iAgl][iVP]/dca_binwidths[iVP]);
            // }
            // if(iVP == 7 ){
            //     hist_Data->SetBinContent(iAgl+1, iVP+1, 117509.);
            //     hist_Data->SetBinError(iAgl+1, iVP+1, 2000);
            // }
            // if(iVP == 9 ){
            //     hist_Data->SetBinContent(iAgl+1, iVP+1, 56042.5);
            //     hist_Data->SetBinError(iAgl+1, iVP+1, 500);
            // }
            //else if(iVP != 4){
            hist_Data->SetBinContent(iVP+1, (sigyield[iAgl][iVP]*histmassbinsize)/dca_binwidths[iVP]);
            hist_Data->SetBinError(iVP+1, sqrt(sigyield[iAgl][iVP]*histmassbinsize)/dca_binwidths[iVP] ) ;
            //}
            double sigsigtmp = sigsig[iAgl][iVP];
            
            tex->DrawLatex(0.22,0.68,Form("Sigsig: %.2f",sigsigtmp));
            
            c->Print(Form("../promptD0/mass_fit/massfit_pt%d_dca%d.pdf",iAgl,iVP));

            /*The full fitting function is constructed as follow
             [0] is signal + swap yield;
             [1] is common mean of double gaussian;
             [2] is signal gaussian 1 sigma;
             [3] is signal gaussian 2 sigma;
             [4] is fractional signal gaussian 1 yield; 1-[4] is fractional signal gaussian 2 yield;
             [5] is fractional double gaussian signal yield, 1-[5] is fractional swap yield;
             [6] is a factor to let width of the gaussians to vary in data;
             [7] is swap gaussian sigma;
             [8] is swap gaussian mean;
             [9-12] is 3rd order poly parameters
             */
             
        }
    }

    hist_Data->Write();

    for(int i = 0; i < 1; i++){
        for(int j = 0; j < 10; j++){
            cout << "error " << i << " " << j << " " << error[i][j]/sigyield[i][j]/dca_binwidths[j] << endl;
        }
    }
    
    
}
