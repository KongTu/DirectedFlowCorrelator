#include "RiceStyle.h"

using namespace std;

void plotMass(){

	TFile* file1 = new TFile("./mass_default.root");
	TFile* file2 = new TFile("./mass_tight.root");

	TH1D* mass_default = (TH1D*) file1->Get("htemp");
	TH1D* mass_tight = (TH1D*) file2->Get("htemp");

	mass_default->SetMarkerStyle(24);
	mass_tight->SetMarkerStyle(20);
	mass_tight->SetMarkerColor(kRed);

	mass_default->Draw("P");
	mass_tight->Draw("Psame");


	

}