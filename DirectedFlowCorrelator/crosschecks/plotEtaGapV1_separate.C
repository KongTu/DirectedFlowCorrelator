#include "RiceStyle.h"

using namespace std;


void plotEtaGapV1_separate(){

	gStyle->SetErrorX(0);

	TFile* file[4];
	TH1D* v1_odd_minus[4];
	TH1D* v1_odd_plus[4];
	for(int i = 0; i < 4; i++){

		file[i] = new TFile(Form("V1_separate_chargedParticles_ana%d.root",i+1));
		v1_odd_minus[i] = (TH1D*) file[i]->Get("hist1_2");
		v1_odd_plus[i] = (TH1D*) file[i]->Get("hist1_3");
	}

	TCanvas* c1 = new TCanvas("c1","c1",700,700);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	TH1D* base2 = makeHist("base2", "Pb-going", "#eta", "v^{odd}_{1}(#eta)", 100,-3.0,3.0,kBlack);
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

	v1_odd_minus[0]->SetMarkerStyle(20);
	v1_odd_minus[0]->SetMarkerColor(kRed);
	v1_odd_minus[0]->SetMarkerSize(1.4);

	v1_odd_minus[1]->SetMarkerStyle(20);
	v1_odd_minus[1]->SetMarkerColor(kBlue);
	v1_odd_minus[1]->SetMarkerSize(1.4);

	v1_odd_minus[2]->SetMarkerStyle(20);
	v1_odd_minus[2]->SetMarkerColor(kBlack);
	v1_odd_minus[2]->SetMarkerSize(1.4);

	v1_odd_minus[3]->SetMarkerStyle(20);
	v1_odd_minus[3]->SetMarkerColor(kGreen-2);
	v1_odd_minus[3]->SetMarkerSize(1.4);

	v1_odd_minus[0]->Draw("Psame");
	v1_odd_minus[1]->Draw("Psame");
	v1_odd_minus[2]->Draw("Psame");
    v1_odd_minus[3]->Draw("Psame");

    v1_odd_plus[0]->SetMarkerStyle(21);
	v1_odd_plus[0]->SetMarkerColor(kRed);
	v1_odd_plus[0]->SetMarkerSize(1.4);

	v1_odd_plus[1]->SetMarkerStyle(21);
	v1_odd_plus[1]->SetMarkerColor(kBlue);
	v1_odd_plus[1]->SetMarkerSize(1.4);

	v1_odd_plus[2]->SetMarkerStyle(21);
	v1_odd_plus[2]->SetMarkerColor(kBlack);
	v1_odd_plus[2]->SetMarkerSize(1.4);

	v1_odd_plus[3]->SetMarkerStyle(21);
	v1_odd_plus[3]->SetMarkerColor(kGreen-2);
	v1_odd_plus[3]->SetMarkerSize(1.4);

	v1_odd_plus[0]->Draw("Psame");
	v1_odd_plus[1]->Draw("Psame");
	v1_odd_plus[2]->Draw("Psame");
    v1_odd_plus[3]->Draw("Psame");

    TLegend *w4 = new TLegend(0.23,0.2,0.50,0.4);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(24);
    w4->SetTextFont(45);
    w4->AddEntry(v1_odd_plus[0], "HF 3.0-3.5", "P");
    w4->AddEntry(v1_odd_plus[1], "HF 3.5-4.0 ","P");
    w4->AddEntry(v1_odd_plus[2], "HF 4.0-4.5", "P");
    w4->AddEntry(v1_odd_plus[3], "HF 4.5-5.0","P");
    w4->Draw("same");


}