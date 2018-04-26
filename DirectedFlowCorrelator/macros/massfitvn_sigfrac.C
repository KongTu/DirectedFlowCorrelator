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

void massfitvn_sigfrac()
{
    // gStyle->SetOptStat(0);
    // gStyle->SetOptTitle(0);
    
    double fit_range_low = 1.7;
    double fit_range_high = 2.0;
    double D0_mass = 1.865;
    double histmassbinsize = (2.0-1.7)/60.;
    TFile* file0 = TFile::Open("hMass_MCSwap_reweightZvtx_rapidity_pt1to40_PbPbMC_ext_eta2p4_dau1p5_loose.root");
    TFile* file1 = TFile::Open("../rootfiles/D0_CDDF_v15.root");
    
    //TFile ofile("FitScan_D0para.root","RECREATE");
    
    double ybin[] = {-2.0,-1.2,-0.6,0,0.6,1.2,2.0};
    double ptbin[14] = {1.5,2.4,3,3.5,4.2,5,6,7,8,10};
    double d0mean[14] = {1.865,1.865,1.865,1.865,1.865,1.865,1.865,1.865,1.865,1.865,1.865,1.865,1.865};
    double d0sigma[14] = {0.0123,0.0123,0.0123,0.0123,0.0123,0.0128,0.0122,0.0120,0.0116,0.0116,0.0120,0.0118,0.0118};
    
    double sigsig[14];
    double sigyield[14];
    double totyield[14];
    double sigfrac_SB[14];
    
    TCanvas* c[10];
    // for(int i=0;i<10;i++)
    // {
    //     c[i] = new TCanvas(Form("c_%d",i),Form("c_%d",i),800,400);
    //     c[i]->Divide(2,1);
    // }
    
    // for(int i=0;i<9;i++)
    // {
    //     c[i]->cd(1)->SetTopMargin(0.06);
    //     c[i]->cd(1)->SetLeftMargin(0.18);
    //     c[i]->cd(1)->SetRightMargin(0.043);
    //     c[i]->cd(1)->SetBottomMargin(0.145);
    //     c[i]->cd(2)->SetTopMargin(0.06);
    //     c[i]->cd(2)->SetLeftMargin(0.18);
    //     c[i]->cd(2)->SetRightMargin(0.043);
    //     c[i]->cd(2)->SetBottomMargin(0.145);
        
    // }
    
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

    double f_obs_sig[6];
    double f_sb_sig[6];
    double f_obs_swap[6];
    double f_sb_swap[6];

    double f_obs[6];
    double f_sb[6];
    for(int ipt=0;ipt<6;ipt++)
    {

        c[ipt] = new TCanvas(Form("c_%d",ipt),Form("c_%d",ipt),600,600);
                        TH1D* h_mc_match_signal = (TH1D*)file0->Get(Form("mass_rapidity%d",ipt));
                        TH1D* h_mc_match_all = (TH1D*)file0->Get(Form("mass_all_rapidity%d",ipt));
                        
                        TH1D* h_data = (TH1D*)file1->Get(Form("ana_loose/D0Mass_Hist_%d_2",ipt));
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

                        h_data->SetMaximum(h_data->GetMaximum()*1.5);

                        c[ipt]->cd(1);
                        gPad->SetLeftMargin(0.1);
                        gPad->SetBottomMargin(0.1);
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
                        
                        TF1* f = new TF1(Form("f_%d",ipt),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])))+(1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + [9] + [10]*x + [11]*x*x + [12]*x*x*x", fit_range_low, fit_range_high);
                        f->SetLineColor(2);
                        f->SetLineWidth(1);
                        
                        //first fit MC signal, swap and poly bkg set to 0
                        
                        f->SetParameter(0,100.);
                        f->SetParameter(1,D0_mass);
                        f->SetParameter(2,0.03);
                        f->SetParameter(3,0.005);
                        f->SetParameter(4,0.1);
                        
                        f->FixParameter(5,1);
                        f->FixParameter(6,0); //always 0 in MC
                        f->FixParameter(7,0.1); //does not really mater here as yield is fix to 0
                        f->FixParameter(8,D0_mass); //does not really mater here as yield is fix to 0
                        f->FixParameter(9,0);
                        f->FixParameter(10,0);
                        f->FixParameter(11,0);
                        f->FixParameter(12,0);
                        
                        f->SetParLimits(2,0.01,0.1);
                        f->SetParLimits(3,0.001,0.05);
                        f->SetParLimits(4,0,1);
                        f->SetParLimits(5,0,1);
                        
                        f->FixParameter(1,1.8648); //for first few attempt fix mean of gaussian to get reasonable estimation of other pars; later open it up
                        h_mc_match_signal->Fit(Form("f_%d",ipt),"q","",fit_range_low,fit_range_high);
                        h_mc_match_signal->Fit(Form("f_%d",ipt),"q","",fit_range_low,fit_range_high);
                        f->ReleaseParameter(1); //now let gaussian mean float
                        h_mc_match_signal->Fit(Form("f_%d",ipt),"L q","",fit_range_low,fit_range_high);
                        h_mc_match_signal->Fit(Form("f_%d",ipt),"L q","",fit_range_low,fit_range_high);
                        h_mc_match_signal->Fit(Form("f_%d",ipt),"L m","",fit_range_low,fit_range_high);
                        
                        //now fix signal double gaussian mean, sigma and gaus1,gaus2 yield ratio
                        f->FixParameter(1,f->GetParameter(1));
                        f->FixParameter(2,f->GetParameter(2));
                        f->FixParameter(3,f->GetParameter(3));
                        f->FixParameter(4,f->GetParameter(4));
                        
                        //now release swap bkg parameters to fit signal+swap MC
                        f->ReleaseParameter(5);
                        f->ReleaseParameter(7);
                        f->ReleaseParameter(8);
                        
                        f->SetParameter(7,0.1);
                        f->SetParameter(8,D0_mass);
                        
                        //fit signal+swap MC
                        h_mc_match_all->Fit(Form("f_%d",ipt),"L q","",fit_range_low,fit_range_high);
                        h_mc_match_all->Fit(Form("f_%d",ipt),"L q","",fit_range_low,fit_range_high);
                        h_mc_match_all->Fit(Form("f_%d",ipt),"L q","",fit_range_low,fit_range_high);
                        h_mc_match_all->Fit(Form("f_%d",ipt),"L q","",fit_range_low,fit_range_high);
                        h_mc_match_all->Fit(Form("f_%d",ipt),"L m","",fit_range_low,fit_range_high);
                        
                        //now fix swap bkg parameters to fit data
                        f->FixParameter(5,f->GetParameter(5));
                        f->FixParameter(7,f->GetParameter(7));
                        f->FixParameter(8,f->GetParameter(8));
                        
                        //now release poly bkg pars
                        f->ReleaseParameter(9);
                        f->ReleaseParameter(10);
                        f->ReleaseParameter(11);
                        f->ReleaseParameter(12);
                        
                        //now fit data
                        h_data->Fit(Form("f_%d",ipt),"q","",fit_range_low,fit_range_high);
                        h_data->Fit(Form("f_%d",ipt),"q","",fit_range_low,fit_range_high);
                        f->ReleaseParameter(1); //allow data to have different mass peak mean than MC
                        f->ReleaseParameter(6); //allow data to have different peak width than MC
                        f->SetParameter(6,0);
                        f->SetParLimits(6,-1,1);
                        h_data->Fit(Form("f_%d",ipt),"L q","",fit_range_low,fit_range_high);
                        h_data->Fit(Form("f_%d",ipt),"L q","",fit_range_low,fit_range_high);
                        h_data->Fit(Form("f_%d",ipt),"L q","",fit_range_low,fit_range_high);
                        h_data->Fit(Form("f_%d",ipt),"L m","",fit_range_low,fit_range_high);
                        
                        //draw D0 signal separately
                        TF1* f1 = new TF1(Form("f_sig_%d",ipt),"[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6]))+(1-[4])*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
                        f1->SetLineColor(kOrange-3);
                        f1->SetLineWidth(1);
                        f1->SetLineStyle(2);
                        f1->SetFillColorAlpha(kOrange-3,0.3);
                        f1->SetFillStyle(1001);
                        f1->FixParameter(0,f->GetParameter(0));
                        f1->FixParameter(1,f->GetParameter(1));
                        f1->FixParameter(2,f->GetParameter(2));
                        f1->FixParameter(3,f->GetParameter(3));
                        f1->FixParameter(4,f->GetParameter(4));
                        f1->FixParameter(5,f->GetParameter(5));
                        f1->FixParameter(6,f->GetParameter(6));
                        f1->Draw("LSAME");
                        
                        //draw swap bkg separately
                        TF1* f2 = new TF1(Form("f_swap_%d",ipt),"[0]*((1-[5])*TMath::Gaus(x,[8],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6])))", fit_range_low, fit_range_high);
                        f2->SetLineColor(kGreen+4);
                        f2->SetLineWidth(1);
                        f2->SetLineStyle(1);
                        f2->SetFillColorAlpha(kGreen+4,0.3);
                        f2->SetFillStyle(1001);
                        f2->FixParameter(0,f->GetParameter(0));
                        f2->FixParameter(5,f->GetParameter(5));
                        f2->FixParameter(6,f->GetParameter(6));
                        f2->FixParameter(7,f->GetParameter(7));
                        f2->FixParameter(8,f->GetParameter(8));
                        f2->Draw("LSAME");
                        
                        //draw poly bkg separately
                        TF1* f3 = new TF1(Form("f_bkg_%d",ipt),"[9] + [10]*x + [11]*x*x + [12]*x*x*x", fit_range_low, fit_range_high);
                        f3->SetLineColor(4);
                        f3->SetLineWidth(1);
                        f3->SetLineStyle(2);
                        f3->FixParameter(9,f->GetParameter(9));
                        f3->FixParameter(10,f->GetParameter(10));
                        f3->FixParameter(11,f->GetParameter(11));
                        f3->FixParameter(12,f->GetParameter(12));
                        f3->Draw("LSAME");
                        
                        tex->DrawLatex(0.22,0.86,"Cent. 30-80%");
                        tex->DrawLatex(0.22,0.80,Form("%.1f < y < %.1f",ybin[ipt],ybin[ipt+1]));
                        tex->DrawLatex(0.22,0.74,"p_{T} > 2");
                        
                        texCMS->DrawLatex(.18,.97,"#font[61]{CMS} #it{Preliminary}");
                        texCMS->DrawLatex(0.62,0.97, "#scale[0.8]{PbPb #sqrt{s_{NN}} = 5.02 TeV}");
                        
                        TLegend* leg = new TLegend(0.65,0.58,0.81,0.9,NULL,"brNDC");
                        leg->SetBorderSize(0);
                        leg->SetTextSize(0.045);
                        leg->SetTextFont(42);
                        leg->SetFillStyle(0);
                        leg->AddEntry(h_data,"data","p");
                        leg->AddEntry(f,"Fit","L");
                        leg->AddEntry(f1,"D^{0} Signal","f");
                        leg->AddEntry(f2,"K-#pi swap","f");
                        leg->AddEntry(f3,"Combinatorial","l");
                        leg->Draw("SAME");

                        double sigma1 = f->GetParameter(2);
                        double sigma2 = f->GetParameter(3);
                        double sig1frac = f->GetParameter(4);
                        double sigmatot = d0sigma[ipt];
                        double mean = d0mean[ipt];
                        
                        //sigma[iAgl][iVP][iDL][idau][ipt] = sqrt(sig1frac*sigma1*sigma1 + (1-sig1frac)*sigma2*sigma2);
                        
                        //sigyield[ipt] = (f1->Integral(mean-2*sigmatot,mean+2*sigmatot)+f2->Integral(mean-2*sigmatot,mean+2*sigmatot))/histmassbinsize;
                        
                        //2 and 3 sig mass window
                        sigyield[ipt] = (f1->Integral(mean-3*sigmatot,mean+3*sigmatot)+f2->Integral(mean-3*sigmatot,mean+3*sigmatot))/histmassbinsize;
                        totyield[ipt] = (f1->Integral(mean-3*sigmatot,mean+3*sigmatot)+f2->Integral(mean-3*sigmatot,mean+3*sigmatot)+f3->Integral(mean-3*sigmatot,mean+3*sigmatot))/histmassbinsize;
                        
                        sigsig[ipt] = sigyield[ipt]/sqrt(totyield[ipt]);
        
                        double sigyield_SB = (f1->Integral(1.7,mean-3*sigmatot) + f1->Integral(mean+3*sigmatot,2.0)+f2->Integral(1.7,mean-3*sigmatot) + f2->Integral(mean+3*sigmatot,2.0))/histmassbinsize;
                        double totyield_SB = (f1->Integral(1.7,mean-3*sigmatot) + f1->Integral(mean+3*sigmatot,2.0) + f2->Integral(1.7,mean-3*sigmatot) + f2->Integral(mean+3*sigmatot,2.0) + f3->Integral(1.7,mean-3*sigmatot) + f3->Integral(mean+3*sigmatot,2.0))/histmassbinsize;
                        sigfrac_SB[ipt] = sigyield_SB/totyield_SB;
        
                        double sigsigtmp = sigsig[ipt];

                        
                        // hard cut
                        // sigyield[ipt] = (f1->Integral(1.82,1.92)+f2->Integral(1.82,1.92))/histmassbinsize;
                        // totyield[ipt] = (f1->Integral(1.82,1.92)+f2->Integral(1.82,1.92)+f3->Integral(1.82,1.92))/histmassbinsize;
                        
                        // sigsig[ipt] = sigyield[ipt]/totyield[ipt];
        
                        // double sigyield_SB = (f1->Integral(1.7,1.82) + f1->Integral(1.92,2.0) + f2->Integral(1.7,1.82) + f2->Integral(1.92,2.0))/histmassbinsize;
                        // double totyield_SB = (f1->Integral(1.7,1.82) + f1->Integral(1.92,2.0) + f2->Integral(1.7,1.82) + f2->Integral(1.92,2.0) + f3->Integral(1.7,1.82) + f3->Integral(1.92,2.0))/histmassbinsize;
                        // sigfrac_SB[ipt] = sigyield_SB/totyield_SB;
        
                        //double sigsigtmp = sigsig[ipt];
                        
                        tex->DrawLatex(0.22,0.68,Form("Sigfrac: %.6f",sigsigtmp));
                       
                        f_obs[ipt] = sigsigtmp;

                        tex->DrawLatex(0.22,0.6,Form("Sigfrac_SB: %.6f",sigfrac_SB[ipt]));
                        
                        f_sb[ipt] = sigfrac_SB[ipt];



                        //split sig and swap
                        // sigyield[ipt] = (f1->Integral(mean-3*sigmatot,mean+3*sigmatot))/histmassbinsize;
                        // double swap_yield = (f2->Integral(mean-3*sigmatot,mean+3*sigmatot))/histmassbinsize;
                        // totyield[ipt] = (f1->Integral(mean-3*sigmatot,mean+3*sigmatot)+f2->Integral(mean-3*sigmatot,mean+3*sigmatot)+f3->Integral(mean-3*sigmatot,mean+3*sigmatot))/histmassbinsize;
                        
                        // sigsig[ipt] = sigyield[ipt]/totyield[ipt];
                        // double swap_swap = swap_yield/totyield[ipt];

                        // double sigyield_SB = (f1->Integral(1.7,mean-3*sigmatot) + f1->Integral(mean+3*sigmatot,2.0))/histmassbinsize;
                        // double swapyield_SB = (f2->Integral(1.7,mean-3*sigmatot) + f2->Integral(mean+3*sigmatot,2.0))/histmassbinsize;
                        // double totyield_SB = (f1->Integral(1.7,mean-3*sigmatot) + f1->Integral(mean+3*sigmatot,2.0) + f2->Integral(1.7,mean-3*sigmatot) + f2->Integral(mean+3*sigmatot,2.0) + f3->Integral(1.7,mean-3*sigmatot) + f3->Integral(mean+3*sigmatot,2.0))/histmassbinsize;
                        
                        // sigfrac_SB[ipt] = sigyield_SB/totyield_SB;
                        // double swap_swap_SB = swapyield_SB/totyield_SB;


                        // double sigsigtmp = sigsig[ipt];
                        
                        // tex->DrawLatex(0.22,0.68,Form("Sigfrac: %.6f",sigsigtmp));
                       
                        // f_obs_sig[ipt] = sigsigtmp;
                        // f_obs_swap[ipt] = swap_swap;

                        // tex->DrawLatex(0.22,0.6,Form("Sigfrac_SB: %.6f",sigfrac_SB[ipt]));
                        
                        // f_sb_sig[ipt] = sigfrac_SB[ipt];
                        // f_sb_swap[ipt] = swap_swap_SB;





                        //c->Print(Form("plots/massfit_pt%d_Agl%d_VP%d_DL%d_daupt%d.pdf",ipt,iAgl,iVP,iDL,idau));
                        
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

    cout << "f_obs[]={";
    for(int i = 0; i < 6; i++){
        cout << f_obs[i] << ",";
    }
    cout << "}" << endl;

    cout << "f_sb[]={";
    for(int i = 0; i < 6; i++){
        cout << f_sb[i] << ",";
    }
    cout << "}" << endl;

    // cout << "f_obs_sig[]={";
    // for(int i = 0; i < 6; i++){
    //     cout << f_obs_sig[i] << ",";
    // }
    // cout << "}" << endl;

    // cout << "f_sb_sig[]={";
    // for(int i = 0; i < 6; i++){
    //     cout << f_sb_sig[i] << ",";
    // }
    // cout << "}" << endl;

    // cout << "f_obs_swap[]={";
    // for(int i = 0; i < 6; i++){
    //     cout << f_obs_swap[i] << ",";
    // }
    // cout << "}" << endl;

    // cout << "f_sb_swap[]={";
    // for(int i = 0; i < 6; i++){
    //     cout << f_sb_swap[i] << ",";
    // }
    // cout << "}" << endl;


}
