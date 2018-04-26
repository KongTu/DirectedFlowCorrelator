#include "RiceStyle.h"

using namespace std;

void plotCombinedHFQvectors(){

	gStyle->SetErrorX(0);
	
	TFile* file = new TFile("../rootfiles/CDDF_v23.root");

	TH1D* c2_cb_plus = (TH1D*) file->Get("ana/c2_cb_plus");
	TH1D* c2_ac_plus = (TH1D*) file->Get("ana/c2_ac_plus");
	TH1D* c2_cb_minus = (TH1D*) file->Get("ana/c2_cb_minus");
	TH1D* c2_ac_minus = (TH1D*) file->Get("ana/c2_ac_minus");

	TH1D* c2_ab = (TH1D*) file->Get("ana/c2_ab");
	TH1D* c2_ab_combined = (TH1D*) file->Get("ana/c2_ab_combined");

	TH1D* c2_a_real = (TH1D*) file->Get("ana/c2_a_real");
	TH1D* c2_b_real = (TH1D*) file->Get("ana/c2_b_real");
	TH1D* c2_c_plus_real = (TH1D*) file->Get("ana/c2_c_plus_real");
	TH1D* c2_c_minus_real = (TH1D*) file->Get("ana/c2_c_minus_real");
	TH1D* c2_ab_real = (TH1D*) file->Get("ana/c2_ab_real");

	TH1D* c2_a_imag = (TH1D*) file->Get("ana/c2_a_imag");
	TH1D* c2_b_imag = (TH1D*) file->Get("ana/c2_b_imag");
	TH1D* c2_c_plus_imag = (TH1D*) file->Get("ana/c2_c_plus_imag");
	TH1D* c2_c_minus_imag = (TH1D*) file->Get("ana/c2_c_minus_imag");
	TH1D* c2_ab_imag = (TH1D*) file->Get("ana/c2_ab_imag");

	TH1D* c2_v1[20][2][3];
	TH1D* c2_trk_accept[20][2][3];

	for(int eta = 0; eta < 8; eta++){
		for(int charge = 0; charge < 2; charge++){
			for(int dir = 0; dir < 3; dir++){

				c2_v1[eta][charge][dir] = (TH1D*) file->Get(Form("ana/c2_v1_%d_%d_%d",eta,charge,dir)); 
				c2_trk_accept[eta][charge][dir] = (TH1D*) file->Get(Form("ana/c2_trk_accept_%d_%d_%d",eta,charge,dir));

			}
		}
	}

	double etabins[] = {-2.4,-2.0,-1.2,-0.6,0.0,0.6,1.2,2.0,2.4};
	//double etabins[] = {-2.4,-2.0,-1.4,-1.0,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,1.0,1.4,2.0,2.4};
	TH1D* hist1[5];
	for(int sign = 0; sign < 5; sign++){		
		hist1[sign] = new TH1D(Form("hist1_%d",sign), "", 8, etabins);
	}

	double Qab = c2_ab_combined->GetMean();
	double Qab_real = c2_ab_real->GetMean();
	double Qab_imag = c2_ab_imag->GetMean();

	double QaQb_accept = (Qab_real * Qab_real + Qab_imag * Qab_imag);
	double resolution =  sqrt(2)*sqrt( Qab - QaQb_accept );
	cout << "acceptance: " << QaQb_accept << endl;
	cout << "resolution: " << resolution << endl;

	double Qab = c2_ab->GetMean();
	double Qa_real = c2_a_real->GetMean();
	double Qb_real = c2_b_real->GetMean();
	double Qa_imag = c2_a_imag->GetMean();
	double Qb_imag = c2_b_imag->GetMean();

	double QaQb_accept = (Qa_real * Qb_real + Qa_imag * Qb_imag);
	double resolution_separate =  sqrt( Qab - QaQb_accept );
	cout << "acceptance: " << QaQb_accept << endl;
	cout << "resolution: " << resolution_separate << endl;

/*
two sub-events
 */
	for(int eta = 0; eta < 8; eta++){

		double t1 = c2_v1[eta][0][2]->GetMean();
		double t1_error = c2_v1[eta][0][2]->GetMeanError();

		double Q1_real = c2_trk_accept[eta][0][0]->GetMean();
		double Q1_imag = c2_trk_accept[eta][0][1]->GetMean();

		t1 = t1 - (Q1_real*Qab_real-Q1_imag*Qab_imag);

		double total = t1/resolution_separate;
		double total_error = t1_error/resolution_separate;

		hist1[2]->SetBinContent(eta+1, total);
		hist1[2]->SetBinError(eta+1, total_error);
	
	}

/*
two sub-events
 */
	for(int eta = 0; eta < 8; eta++){

		double t1 = c2_v1[eta][1][0]->GetMean();
		double t1_error = c2_v1[eta][1][0]->GetMeanError();

		double Q1_real = c2_trk_accept[eta][1][0]->GetMean();
		double Q1_imag = c2_trk_accept[eta][1][1]->GetMean();

		t1 = t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);

		double t2 = c2_v1[eta][1][1]->GetMean();
		double t2_error = c2_v1[eta][1][1]->GetMeanError();

		t2 = t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

		double total = (t1/resolution_separate + t2/resolution_separate)/2.0;
		double total_error = (1/(2*resolution_separate))*sqrt(t1_error*t1_error + t2_error*t2_error);

		hist1[3]->SetBinContent(eta+1, total);
		hist1[3]->SetBinError(eta+1, total_error);
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

	hist1[2]->SetMarkerStyle(20);
	hist1[2]->SetMarkerColor(kRed);
	hist1[2]->SetMarkerSize(1.4);

	hist1[2]->Draw("Psame");

	hist1[3]->SetMarkerStyle(20);
	hist1[3]->SetMarkerColor(kBlue);
	hist1[3]->SetMarkerSize(1.4);

	hist1[3]->Draw("Psame");


}