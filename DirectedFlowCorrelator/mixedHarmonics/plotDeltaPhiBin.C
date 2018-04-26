#include "RiceStyle.h"

using namespace std;

void plotDeltaPhiBin(){
	
	TFile* file = new TFile("../rootfiles/CDDF_v21.root");

	TH1D* c2_ab = (TH1D*) file->Get("ana/c2_ab");

	TH1D* c2_a_real = (TH1D*) file->Get("ana/c2_a_real");
	TH1D* c2_b_real = (TH1D*) file->Get("ana/c2_b_real");
	TH1D* c2_a_imag = (TH1D*) file->Get("ana/c2_a_imag");
	TH1D* c2_b_imag = (TH1D*) file->Get("ana/c2_b_imag");

	TH1D* delta_phi_positive[8][2];
	TH1D* delta_phi_negative[8][2];

	for(int eta = 0; eta < 8; eta++){
		for(int side = 0; side < 2; side++){
			
			delta_phi_positive[eta][side] = (TH1D*) file->Get(Form("ana/delta_phi_positive_%d_%d",eta,side));
			delta_phi_negative[eta][side] = (TH1D*) file->Get(Form("ana/delta_phi_negative_%d_%d",eta,side));
		}		
	}

	double Qa_real = c2_a_real->GetMean();
	double Qb_real = c2_b_real->GetMean();
	double Qa_imag = c2_a_imag->GetMean();
	double Qb_imag = c2_b_imag->GetMean();
	double Qab = c2_ab->GetMean();
	
	double QaQb_accept = (Qa_real * Qb_real + Qa_imag * Qb_imag);
	double resolution =  sqrt( Qab - QaQb_accept );
	
	double etabins[] = {-2.4,-2.0,-1.2,-0.6,0.0,0.6,1.2,2.0,2.4};
	TH1D* hist1[5];
	for(int sign = 0; sign < 5; sign++){		
		hist1[sign] = new TH1D(Form("hist1_%d",sign), "", 8, etabins);
	}

	TF1 * f1[8][2];
	for(int eta = 0; eta < 8; eta++){
		for(int side = 0; side < 2; side++){

			f1[eta][side] = new TF1(Form("f1_%d_%d",eta,side), "[0]+2*[1]*TMath::Cos(1*x[0]) + 2*[2]*TMath::Cos(2*x[0])");
			f1[eta][side]->SetRange(0,3.1415926);
			delta_phi_positive[eta][side]->Scale(1/delta_phi_positive[eta][side]->Integral());
			delta_phi_positive[eta][side]->SetStats(kFALSE);
			delta_phi_positive[eta][side]->Fit(Form("f1_%d_%d",eta,side), "ER");

			double v1 = f1[eta][side]->GetParameter(1);
			//v1 = v1/resolution;
			hist1[side]->SetBinContent(eta+1, v1);	
		}
	}

	TCanvas* c1 = new TCanvas("c1","c1",700,700);
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

	hist1[0]->SetMarkerStyle(25);
	hist1[0]->SetMarkerColor(kRed);
	hist1[0]->SetMarkerSize(1.4);

	hist1[1]->SetMarkerStyle(24);
	hist1[1]->SetMarkerColor(kBlue);
	hist1[1]->SetMarkerSize(1.4);

	hist1[0]->Draw("Psame");
	hist1[1]->Draw("Psame");


}