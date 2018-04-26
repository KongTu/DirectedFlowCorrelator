#include "RiceStyle.h"

using namespace std;

void plotMixedHarmonics(){

	TFile* file = new TFile("../rootfiles/CDDF_v23.root");

	TH1D* threePart[8];
	TH1D* phi_cos[8];
	TH1D* phi_sin[8];

	TH1D* Res_Psi_1_Psi_2;
	TH1D* Res_2;
	TH1D* Res_2_real_0;
	TH1D* Res_2_imag_0;
	TH1D* Res_2_real_1;
	TH1D* Res_2_imag_1;

	for(int eta = 0; eta < 8; eta++){

		threePart[eta] = (TH1D*) file->Get(Form("ana/Phi_Psi_1_Psi_2_%d",eta));
		
	}
	for(int eta = 0; eta < 8; eta++){

		phi_cos[eta] = (TH1D*) file->Get(Form("ana/Phi_Average_cos_%d",eta));
		phi_sin[eta] = (TH1D*) file->Get(Form("ana/Phi_Average_sin_%d",eta));
	}

	Res_Psi_1_Psi_2 = (TH1D*) file->Get("ana/Psi_1_Psi_2");

	Res_2 = (TH1D*) file->Get("ana/Psi_2_trk_reso");
	Res_2_real_0 = (TH1D*) file->Get("ana/Psi_2_trk_accept_real_0");
	Res_2_imag_0 = (TH1D*) file->Get("ana/Psi_2_trk_accept_imag_0");
	Res_2_real_1 = (TH1D*) file->Get("ana/Psi_2_trk_accept_real_1");
	Res_2_imag_1 = (TH1D*) file->Get("ana/Psi_2_trk_accept_imag_1");

	cout << "cos: {";
	for(int eta = 0; eta < 8; eta++){

		cout << phi_cos[eta]->GetMean() << ", ";
	}
	cout << endl;

	cout << "sin: {";
	for(int eta = 0; eta < 8; eta++){

		cout << phi_sin[eta]->GetMean() << ", ";
	}
	cout << endl;

/*
Res Psi_1 and Psi_2
*/

	// double Res_Psi_12_left = Res_Psi_1_Psi_2[0]->GetMean();
	// double Res_Psi_12_right = Res_Psi_1_Psi_2[1]->GetMean();
	
	//double Res_2_Tracker = Res_2->GetMean() - (Res_2_real_0->GetMean()*Res_2_real_1->GetMean()) - (Res_2_imag_0->GetMean()*Res_2_imag_1->GetMean());
	//Res_2_Tracker = sqrt(Res_2_Tracker);
	// double resolution = sqrt( Res_Psi_12 );
	// cout << "resolution " << resolution << endl;

	double etabins[] = {-2.4,-2.0,-1.2,-0.6,0.0,0.6,1.2,2.0,2.4};
	TH1D* hist1[5];
	for(int sign = 0; sign < 2; sign++){		
		hist1[sign] = new TH1D(Form("hist1_%d",sign), "", 8, etabins);
	}

	for(int eta = 0; eta < 8; eta++){

		double value1 = threePart[eta]->GetMean();
		double error1 = threePart[eta]->GetMeanError();

		cout << "eta " << eta+1 << " " << value1 << endl;

		// double total = value1/resolution;
		// double total_error = error1/resolution;

		hist1[0]->SetBinContent(eta+1, value1);
		hist1[0]->SetBinError(eta+1, error1);

		
	}

	TCanvas* c1 = new TCanvas("c1","c1",700,700);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	TH1D* base2 = makeHist("base2", "Pb-going", "#eta", "v^{odd}_{1}(#eta)", 100,-3.0,3.0,kBlack);
	base2->GetYaxis()->SetRangeUser(-0.05, 0.05);
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

	hist1[0]->SetMarkerStyle(25);
	hist1[0]->SetMarkerColor(kRed);
	hist1[0]->SetMarkerSize(1.4);
	hist1[0]->Draw("Psame");

	// hist1[1]->SetMarkerStyle(24);
	// hist1[1]->SetMarkerColor(kBlue);
	// hist1[1]->SetMarkerSize(1.4);
	// hist1[1]->Draw("Psame");


}