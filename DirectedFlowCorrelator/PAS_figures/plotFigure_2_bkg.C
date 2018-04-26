#include "../macros/RiceStyle.h"

using namespace std;


void plotFigure_2_bkg(){

    double total_inclusive_D0D0BAR_sys = sqrt(0.005*0.005 + 0.02*0.02 + 0.008*0.008 + 0.01*0.01 ); 
    double total_inclusive_D0_sys = sqrt(0.005*0.005 + 0.025*0.025 + 0.034*0.034 + 0.01*0.01 ); 
    double total_inclusive_D0BAR_sys = sqrt(0.005*0.005 + 0.035*0.035 + 0.015*0.015 + 0.01*0.01 ); 

	gStyle->SetErrorX(0);

	TFile* file1 = new TFile("../macros/V1_odd_chargedParticles.root");
	TH1D* cV1_plus = (TH1D*) file1->Get("hist1_2");
	TH1D* cV1_minus = (TH1D*) file1->Get("hist1_3");
	TH1D* cV1_comb = (TH1D*) file1->Get("charge_inde");
	TH1D* cV1_diff = (TH1D*) file1->Get("ratio");

    cV1_plus->SetBinContent(1,100); cV1_plus->SetBinContent(8,100);
    cV1_minus->SetBinContent(1,100); cV1_minus->SetBinContent(8,100);
    cV1_comb->SetBinContent(1,100); cV1_comb->SetBinContent(8,100);

	TFile* file2 = new TFile("../macros/v1vsy_combined_v13.root" );
	TFile* file3 = new TFile("../macros/v1vsy_test_v13_bkgseparate.root" );

	TGraphErrors* gr1_d0d0bar = (TGraphErrors*) file2->Get("v1vsy");
	TGraphErrors* gr1_bkg = (TGraphErrors*) file2->Get("v1bkgvsy");

	TGraphErrors* gr2_d0 = (TGraphErrors*) file3->Get("v1vsy");
	TGraphErrors* gr2_dobar = (TGraphErrors*) file3->Get("v1vsy_anti");
	TGraphErrors* gr2_diff = (TGraphErrors*) file3->Get("v1vsy_diff");
    TGraphErrors* gr2_bkg = (TGraphErrors*) file3->Get("v1bkgvsy");
    TGraphErrors* gr2_antibkg = (TGraphErrors*) file3->Get("v2antibkg");


	TCanvas* c1 = new TCanvas("c1","c1",600,600);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	//gStyle->SetPadBorderMode(0.1);
	//gStyle->SetOptTitle(0);

	TH1D* base2 = makeHist("base2", "", "y^{D^{0}},#eta^{ch}", "v^{odd}_{1}", 100,-2.2,2.2,kBlack);
	base2->GetYaxis()->SetRangeUser(-0.17, 0.17);
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

    //drawBoxTGraph(gr1_d0d0bar, 6, total_inclusive_D0D0BAR_sys, false, true);

    //drawBox(cV1_comb, 0.06, true);

	cV1_comb->SetMarkerStyle(25);
	//cV1_comb->Draw("Psame");

	gr2_bkg->SetMarkerStyle(20);
	gr2_bkg->SetMarkerColor(kRed);
	gr2_bkg->SetLineColor(kRed);
	gr2_bkg->SetMarkerSize(1.4);

	gr2_antibkg->SetMarkerStyle(20);
	gr2_antibkg->SetMarkerColor(kBlue);
	gr2_antibkg->SetMarkerSize(1.4);

	gr2_bkg->Draw("Psame");
	gr2_antibkg->Draw("Psame");

	TLatex* r42 = new TLatex(0.2, 0.85, "PbPb 5.02 TeV");
    r42->SetNDC();
    r42->SetTextSize(23);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);

    TLatex* r43 = new TLatex(0.62,0.91, "CMS");
    r43->SetNDC();
    r43->SetTextSize(0.04);
    
    TLatex* r44 = new TLatex(0.71,0.91, "Preliminary");
    r44->SetNDC();
    r44->SetTextSize(22);
    r44->SetTextFont(53);

    TLatex* r45 = new TLatex(0.2, 0.8, "Cent.30-80%");
    r45->SetNDC();
    r45->SetTextSize(23);
    r45->SetTextFont(43);
    r45->SetTextColor(kBlack);

    TLatex* r47 = new TLatex(0.2, 0.75, "p_{T} > 2.0 GeV");
    r47->SetNDC();
    r47->SetTextSize(23);
    r47->SetTextFont(43);
    r47->SetTextColor(kBlack);

    r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");
    r45->Draw("same");
    r47->Draw("same");

	TLegend *w4 = new TLegend(0.23,0.20,0.6,0.32);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(22);
    w4->SetTextFont(45);

    w4->AddEntry(gr2_bkg, "D^{0} bkg  ", "P");
    w4->AddEntry(gr2_antibkg, "#bar{D} bkg","P");
    //w4->AddEntry(gr1_d0d0bar, "D^{0}+#bar{D^{0}}","P");
    w4->Draw("same");

    c1->Print("Figure_Sup_1.pdf");




}