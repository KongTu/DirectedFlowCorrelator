#include "RiceStyle.h"

using namespace std;

void plotMixedHarmonics_Qvectors(){
	
	gStyle->SetErrorX(0);

	TFile* file = new TFile("../rootfiles/CDDF_v25.root");

	TH1D* threePart[8][2];
	TH1D* phi_cos[8];
	TH1D* phi_sin[8];

	TH1D* Res_Psi_1_Psi_2;
	TH1D* Res_2;
	TH1D* Res_2_real_0;
	TH1D* Res_2_imag_0;
	TH1D* Res_2_real_1;
	TH1D* Res_2_imag_1;

	for(int eta = 0; eta < 8; eta++){
		for(int charge = 0; charge < 2; charge++){
			
			threePart[eta][charge] = (TH1D*) file->Get(Form("ana/c2_v1_mixed_%d_%d",eta, charge));

		}
	}
	// for(int eta = 0; eta < 8; eta++){

	// 	phi_cos[eta] = (TH1D*) file->Get(Form("ana/Phi_Average_cos_%d",eta));
	// 	phi_sin[eta] = (TH1D*) file->Get(Form("ana/Phi_Average_sin_%d",eta));
	// }

	TH1D* c2_ab = (TH1D*) file->Get("ana/c2_ab");
	TH1D* c2_a_real = (TH1D*) file->Get("ana/c2_a_real");
	TH1D* c2_b_real = (TH1D*) file->Get("ana/c2_b_real");
	TH1D* c2_a_imag = (TH1D*) file->Get("ana/c2_a_imag");
	TH1D* c2_b_imag = (TH1D*) file->Get("ana/c2_b_imag");

	double Qab = c2_ab->GetMean();
	double Qa_real = c2_a_real->GetMean();
	double Qb_real = c2_b_real->GetMean();
	double Qa_imag = c2_a_imag->GetMean();
	double Qb_imag = c2_b_imag->GetMean();

	double QaQb_accept = (Qa_real * Qb_real + Qa_imag * Qb_imag);
	double resolution =  sqrt( Qab - QaQb_accept );
	cout << "acceptance: " << QaQb_accept << endl;
	cout << "resolution: " << resolution << endl;

	Res_2 = (TH1D*) file->Get("ana/Psi_2_trk_reso");
	Res_2_real_0 = (TH1D*) file->Get("ana/Psi_2_trk_accept_real_0");
	Res_2_imag_0 = (TH1D*) file->Get("ana/Psi_2_trk_accept_imag_0");
	Res_2_real_1 = (TH1D*) file->Get("ana/Psi_2_trk_accept_real_1");
	Res_2_imag_1 = (TH1D*) file->Get("ana/Psi_2_trk_accept_imag_1");

	double Q_Res_2 = Res_2->GetMean();
	double Q_Res_2_real_0 = Res_2_real_0->GetMean();
	double Q_Res_2_real_1 = Res_2_real_1->GetMean();
	double Q_Res_2_imag_0 = Res_2_imag_0->GetMean();
	double Q_Res_2_imag_1 = Res_2_imag_1->GetMean();

	double Q_Res_2_correction = (Q_Res_2_real_0 * Q_Res_2_real_1 + Q_Res_2_imag_0 * Q_Res_2_imag_1);
	double resolution_2 = sqrt( Q_Res_2 - Q_Res_2_correction);
	cout << "resolution 2: " << resolution_2 << endl;
	cout << "total resolution: " << resolution*resolution_2 << endl;

	double etabins[] = {-2.4,-2.0,-1.2,-0.8,0.0,0.8,1.2,2.0,2.4};
	TH1D* hist1[5];
	for(int sign = 0; sign < 2; sign++){		
		hist1[sign] = new TH1D(Form("hist1_%d",sign), "", 8, etabins);
	}

	for(int eta = 0; eta < 8; eta++){

		double value1 = threePart[eta][0]->GetMean();
		double error1 = threePart[eta][0]->GetMeanError();

		cout << "eta " << eta+1 << " " << value1 << endl;

		double total = value1/(resolution*resolution_2);
		double total_error = error1/(resolution*resolution_2);

		cout << "eta " << eta+1 << " " << total << endl;

		hist1[0]->SetBinContent(eta+1, total);
		hist1[0]->SetBinError(eta+1, total_error);
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
	hist1[0]->Draw("Psame");

	TFile* file1[4];
	TH1D* v1_odd[4];
	for(int i = 0; i < 4; i++){

		file1[i] = new TFile(Form("../crosschecks/V1_odd_chargedParticles_ana%d.root",i+1));
		v1_odd[i] = (TH1D*) file1[i]->Get("hist1_2");
	}

	v1_odd[3]->SetMarkerStyle(20);
	v1_odd[3]->SetMarkerColor(kBlue);
	v1_odd[3]->SetMarkerSize(1.4);
	v1_odd[3]->Draw("Psame");

	TLegend *w4 = new TLegend(0.23,0.2,0.50,0.3);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(24);
    w4->SetTextFont(45);
   	w4->AddEntry(v1_odd[3], "Symmetric sub-events","P");
    w4->AddEntry(hist1[0], "Mixed harmonics","P");
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

    c1->Print("mixedHarmonics.pdf");

    TFile outfile("mixedHarmonics_tracker_psi2.root","RECREATE");
    hist1[0]->Write();

}