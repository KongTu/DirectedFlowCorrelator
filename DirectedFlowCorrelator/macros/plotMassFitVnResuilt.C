#include "RiceStyle.h"

using namespace std;

void plotMassFitVnResuilt(){

	TFile* file1 = new TFile("./v1vsy_combined_v9.root" );
	TFile* file2 = new TFile("./v1vsy_test_v9.root" );

	TGraphErrors* gr1_d0d0bar = (TGraphErrors*) file1->Get("v1vsy");
	TGraphErrors* gr1_bkg = (TGraphErrors*) file1->Get("v1bkgvsy");

	TGraphErrors* gr2_d0 = (TGraphErrors*) file2->Get("v1vsy");
	TGraphErrors* gr2_dobar = (TGraphErrors*) file2->Get("v1vsy_anti");
	TGraphErrors* gr2_diff = (TGraphErrors*) file2->Get("v1vsy_diff");
	TGraphErrors* gr2_bkg = (TGraphErrors*) file2->Get("v1bkgvsy");

	TCanvas* c1 = new TCanvas("c1","c1",600,600);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	//gStyle->SetPadBorderMode(0.1);
	//gStyle->SetOptTitle(0);

	TH1D* base2 = makeHist("base2", "", "y", "v^{odd}_{1}(y)", 100,-3.0,3.0,kBlack);
	base2->GetYaxis()->SetRangeUser(-0.2, 0.2);
	base2->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base2,1.1,1.25);

	base2->GetYaxis()->SetTitleOffset(1.3);
	base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.6);
	base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.6);
	base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.6);
	base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.6);
	base2->GetXaxis()->SetNdivisions(4,6,0);
	base2->GetYaxis()->SetNdivisions(4,6,0);

	base2->Draw();

	gr2_d0->SetMarkerStyle(20);
	gr2_d0->SetLineColor(kRed);
	gr2_d0->SetMarkerColor(kRed);
	gr2_d0->SetMarkerSize(1.5);

	gr2_dobar->SetMarkerStyle(20);
	gr2_dobar->SetLineColor(kBlue);
	gr2_dobar->SetMarkerColor(kBlue);
	gr2_dobar->SetMarkerSize(1.5);

	// gr2_d0->Draw("Psame");
	// gr2_dobar->Draw("Psame");

	gr1_d0d0bar->SetMarkerStyle(20);
	gr1_d0d0bar->SetMarkerColor(kRed);
	gr1_d0d0bar->SetLineColor(kBlue);
	gr1_d0d0bar->SetMarkerSize(1.4);

	gr1_bkg->SetMarkerStyle(20);
	gr1_bkg->SetMarkerColor(kBlue);
	gr1_bkg->SetMarkerSize(1.4);

	gr1_d0d0bar->Draw("Psame");
	//gr1_bkg->Draw("Psame");

	TLatex* r42 = new TLatex(0.2, 0.85, "PbPb 5.02 TeV");
    r42->SetNDC();
    r42->SetTextSize(23);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);

    TLatex* r43 = new TLatex(0.65,0.91, "CMS");
    r43->SetNDC();
    r43->SetTextSize(0.04);
    
    TLatex* r44 = new TLatex(0.75,0.91, "Preliminary");
    r44->SetNDC();
    r44->SetTextSize(22);
    r44->SetTextFont(53);

    TLatex* r45 = new TLatex(0.2, 0.8, "Cent.0-30%");
    r45->SetNDC();
    r45->SetTextSize(23);
    r45->SetTextFont(43);
    r45->SetTextColor(kBlack);

    TLatex* r47 = new TLatex(0.2, 0.75, "p_{T,D^{0}} > 2.0 GeV");
    r47->SetNDC();
    r47->SetTextSize(23);
    r47->SetTextFont(43);
    r47->SetTextColor(kBlack);

    r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");
    r45->Draw("same");
    r47->Draw("same");

	TLegend *w4 = new TLegend(0.23,0.20,0.70,0.26);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(22);
    w4->SetTextFont(45);

    w4->AddEntry(gr2_d0, "D^{0}  ", "P");
    w4->AddEntry(gr2_dobar, "#bar{D^{0}}","P");
    //w4->AddEntry(gr1_d0d0bar, "D^{0}+#bar{D^{0}}","P");
    w4->Draw("same");

	TCanvas* c2 = new TCanvas("c2","c2",600,600);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	TH1D* base3 = makeHist("base3", "Pb-going", "y", "v^{D^{0}}_{1}#minusv^{#bar{D^{0}}}_{1}", 100,-3.0,3.0,kBlack);
	base3->GetYaxis()->SetRangeUser(-0.3, 0.3.);
	base3->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base3,1.1,1.25);

	base3->GetYaxis()->SetTitleOffset(1.3);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.6);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.6);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.6);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.6);
	base3->GetXaxis()->SetNdivisions(4,6,0);
	base3->GetYaxis()->SetNdivisions(4,6,0);

	
	
	gr2_diff->Fit("pol1");
	gr2_diff->SetMarkerStyle(20);
	gr2_diff->SetMarkerSize(1.6);
	gr2_diff->SetMarkerColor(kBlack);
	gr2_diff->SetLineColor(kBlack);
	//gr2_diff->Draw("Psame");
	base3->Draw();
	gr2_diff->Draw("Psame");

	TF1* func1 = (TF1*) gr2_diff->GetFunction("pol1");
	func1->SetLineColor(kGreen-3);
	func1->SetMarkerColor(kGreen-3);
	double slope = func1->GetParameter(1);
	double slope_error = func1->GetParError(1);
	TH1D *hint1 = new TH1D("hint1","Fitted gaussian with .95 conf.band", 1000, -3, 3.0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint1, 0.68);
     //Now the "hint1" histogram has the fitted function values as the
     //bin contents and the confidence intervals as bin errors
    hint1->SetStats(kFALSE);
    hint1->SetFillColor(kGreen-3);
    hint1->SetMarkerColor(kGreen-3);
    hint1->SetFillStyle(1001);
    hint1->SetFillColorAlpha(kGreen-3,0.4);
    hint1->Draw("e3 same");

	r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");
    r45->Draw("same");
    r47->Draw("same");

    TLatex* r48 = new TLatex(0.2, 0.7, Form("Slope: %.3f +/- %.3f",slope, slope_error ));
    r48->SetNDC();
    r48->SetTextSize(23);
    r48->SetTextFont(43);
    r48->SetTextColor(kBlack);
    r48->Draw("same");

    c1->Print("../plots/v1_combineD0D0bar_0_30.pdf");
    //c2->Print("../plots/v1_D0D0bar_diffdiff_2p0.pdf");


}