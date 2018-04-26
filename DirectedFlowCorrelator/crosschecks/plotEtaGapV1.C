#include "RiceStyle.h"

using namespace std;


void plotEtaGapV1(){

	gStyle->SetErrorX(0);

	TFile* file[4];
	TH1D* v1_odd[4];
	for(int i = 0; i < 4; i++){

		file[i] = new TFile(Form("V1_odd_chargedParticles_ana%d.root",i+1));
		v1_odd[i] = (TH1D*) file[i]->Get("hist1_2");
	}

	TCanvas* c1 = new TCanvas("c1","c1",600,600);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	TH1D* base2 = makeHist("base2", "Pb-going", "#eta", "v^{odd}_{1}(#eta)", 100,-3.0,3.0,kBlack);
	base2->GetYaxis()->SetRangeUser(-0.1, 0.1);
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

	v1_odd[0]->SetMarkerStyle(24);
	v1_odd[0]->SetMarkerColor(kRed);
	v1_odd[0]->SetMarkerSize(1.4);

	v1_odd[1]->SetMarkerStyle(24);
	v1_odd[1]->SetMarkerColor(kBlue);
	v1_odd[1]->SetMarkerSize(1.4);

	v1_odd[2]->SetMarkerStyle(20);
	v1_odd[2]->SetMarkerColor(kRed);
	v1_odd[2]->SetMarkerSize(1.4);

	v1_odd[3]->SetMarkerStyle(20);
	v1_odd[3]->SetMarkerColor(kBlue);
	v1_odd[3]->SetMarkerSize(1.4);

	v1_odd[0]->Draw("Psame");
	v1_odd[1]->Draw("Psame");
	v1_odd[2]->Draw("Psame");
    v1_odd[3]->Draw("Psame");

    TLegend *w4 = new TLegend(0.23,0.2,0.50,0.4);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(24);
    w4->SetTextFont(45);
    w4->AddEntry(v1_odd[0], "HF 3.0-3.5", "P");
    w4->AddEntry(v1_odd[1], "HF 3.5-4.0 ","P");
    w4->AddEntry(v1_odd[2], "HF 4.0-4.5", "P");
    w4->AddEntry(v1_odd[3], "HF 4.5-5.0","P");
    w4->Draw("same");

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

    c1->Print("v1_gap_study.pdf");

}