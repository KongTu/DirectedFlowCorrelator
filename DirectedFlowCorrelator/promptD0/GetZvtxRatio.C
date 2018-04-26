#include "RiceStyle.h"

using namespace std;

void GetZvtxRatio(){

	TFile* file_data = new TFile("../rootfiles/D0_CDDF_v5plusv7_HI.root");
	TH1D* data = (TH1D*) file_data->Get("ana/vtxZ");

	TFile* file_mc = new TFile("./vtxZ_mc.root");
	TH1D* mc = (TH1D*) file_mc->Get("hist");

	cout << "data bins: " << data->GetNbinsX() << endl;
	cout << "mc bins: " << mc->GetNbinsX() << endl;

	TH1D* ratio = (TH1D*) data->Clone("ratio");

	ratio->Divide( mc );

	ratio->Draw();

	TFile f("vtxRatio.root","RECREATE");

	ratio->Write();
	

}