#include "RiceStyle.h"

using namespace std;

void plotDCA(){

	TFile* file1 = new TFile("./D0_promptNonPrompt_MC_rap0_pt_vtxZreweight.root");
	TFile* file2 = new TFile("./D0_promptNonPrompt_MC_rap0_pt_vtxZreweight_0p3Matching.root");

	TH1D* mass_default = (TH1D*) file1->Get("hist_prompt_0");
	TH1D* mass_tight = (TH1D*) file2->Get("hist_prompt_0");

	mass_default->SetMarkerStyle(24);
	mass_tight->SetMarkerStyle(20);
	mass_tight->SetMarkerColor(kRed);

	
	mass_default->DrawNormalized("P");
	mass_tight->DrawNormalized("Psame");

	TLegend *w4 = new TLegend(0.53,0.50,0.80,0.66);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(22);
    w4->SetTextFont(45);

    w4->AddEntry(mass_default, "#Deltap_{T} < 0.5 #DeltaR < 0.5","P");
    w4->AddEntry(mass_tight, "#Deltap_{T} < 0.3 #DeltaR < 0.3","P");
    w4->Draw("same");
	

}