#include "RiceStyle.h"

using namespace std;

void plotPromptFraction_Cut(){
	
	TFile* file1 = new TFile("prompt_frac_test.root");
	TH1D* hist1 = (TH1D*) file1->Get("PromptFrac_cut");

	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);

	hist1->GetXaxis()->SetTitle("DCA");
	hist1->GetYaxis()->SetTitle("Prompt Fraction");
	hist1->GetYaxis()->SetRangeUser(0.3,1.1);

	hist1->SetMarkerStyle(20);
	hist1->SetMarkerColor(kBlack);

	hist1->SetTitle("");
	hist1->Draw("P");

	TFile* file2 = new TFile("prompt_frac_test_1.root");
	TH1D* hist2 = (TH1D*) file2->Get("PromptFrac_cut");

	hist2->GetXaxis()->SetTitle("DCA");
	hist2->GetYaxis()->SetTitle("Prompt Fraction");

	hist2->SetMarkerStyle(20);
	hist2->SetMarkerColor(kRed);

	hist2->Draw("Psame");

	TLegend *w4 = new TLegend(0.23,0.20,0.6,0.32);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(22);
    w4->SetTextFont(45);

    w4->AddEntry(hist1, "0 < |y| < 0.6  ", "P");
    w4->AddEntry(hist2, "0 < |y| < 1.2","P");
    w4->Draw("same");
}