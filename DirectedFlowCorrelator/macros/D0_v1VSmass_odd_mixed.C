#include "RiceStyle.h"
#include "LinearSolver.h"

using namespace std;

void D0_v1VSmass_odd_mixed(){

	gStyle->SetErrorX(0);
	
/*
Loading all the histograms
 */

	TFile* file = new TFile("../rootfiles/D0_CDDF_v17.root");

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

	Res_2 = (TH1D*) file->Get("ana/Psi_2_HF_reso");
	Res_2_real_0 = (TH1D*) file->Get("ana/Psi_2_HF_accept_real_0");
	Res_2_imag_0 = (TH1D*) file->Get("ana/Psi_2_HF_accept_imag_0");
	Res_2_real_1 = (TH1D*) file->Get("ana/Psi_2_HF_accept_real_1");
	Res_2_imag_1 = (TH1D*) file->Get("ana/Psi_2_HF_accept_imag_1");

	double Q_Res_2 = Res_2->GetMean();
	double Q_Res_2_real_0 = Res_2_real_0->GetMean();
	double Q_Res_2_real_1 = Res_2_real_1->GetMean();
	double Q_Res_2_imag_0 = Res_2_imag_0->GetMean();
	double Q_Res_2_imag_1 = Res_2_imag_1->GetMean();

	double Q_Res_2_correction = (Q_Res_2_real_0 * Q_Res_2_real_1 + Q_Res_2_imag_0 * Q_Res_2_imag_1);
	double resolution_2 = sqrt( Q_Res_2 - Q_Res_2_correction);
	cout << "resolution 2: " << resolution_2 << endl;
	cout << "total resolution: " << resolution*resolution_2 << endl;


	TH1D* c2_d0_v1[20][30][3][3];
	TH1D* c2_d0_trk_accept[20][30][3][2];
	TH1D* D0Mass_Hist[20][3];

	for(int eta = 0; eta < 6; eta++){
		for(int charge = 0; charge < 3; charge++){
			
			D0Mass_Hist[eta][charge] = (TH1D*) file->Get(Form("ana/D0Mass_Hist_%d_%d",eta,charge));
			
			for(int imass = 0; imass < 12; imass++){	
				for(int dir = 0; dir < 3; dir++){

					c2_d0_v1[eta][imass][charge][dir] = (TH1D*) file->Get(Form("ana/c2_d0_v1_%d_%d_%d_%d",eta,imass,charge,dir));
			
				}
			}
		}
	}

/*
end of loading
 */
	double etabins[] = {-2.4,-1.4,-0.7,0.0,0.7,1.4,2.4};
	double ybins[] = {-2.0,-1.2,-0.6,0.0,0.6,1.2,2.0};
	double massbins[] = {1.7,1.73,1.76,1.79,1.81,1.83,1.85,1.87,1.89,1.91,1.94,1.97,2.0};
	TH1D* hist1[4];
	for(int sign = 0; sign < 4; sign++){		
		hist1[sign] = new TH1D(Form("hist1_%d",sign), "", 6, etabins);
	}

/*
D0 v1 
*/
	TH1D* hist_d0_1[6];
	TH1D* hist_d0bar_1[6];
	TH1D* hist_d0d0bar_1[6];
	for(int rap = 0; rap < 6; rap++){		
		hist_d0_1[rap] = new TH1D(Form("hist_d0_1_%d",rap), "", 12, massbins);
		hist_d0bar_1[rap] = new TH1D(Form("hist_d0bar_1_%d",rap), "", 12, massbins);
		hist_d0d0bar_1[rap] = new TH1D(Form("hist_d0d0bar_1_%d",rap), "", 12, massbins);
	}

	for(int rap = 0; rap < 6; rap++){
		for(int imass = 0; imass < 12; imass++){

		//D0 v1 vs mass

			double d0_t1 = c2_d0_v1[rap][imass][0][2]->GetMean();
			double d0_t1_error = c2_d0_v1[rap][imass][0][2]->GetMeanError();

			double total_d0 = (d0_t1/(resolution*resolution_2));
			double total_d0_error = (d0_t1_error/(resolution*resolution_2));

			hist_d0_1[rap]->SetBinContent(imass+1, total_d0);
			hist_d0_1[rap]->SetBinError(imass+1, total_d0_error);

		//D0bar v1 vs mass

			d0_t1 = c2_d0_v1[rap][imass][1][2]->GetMean();
			d0_t1_error = c2_d0_v1[rap][imass][1][2]->GetMeanError();

			double total_d0bar = (d0_t1/(resolution*resolution_2));
			double total_d0bar_error = (d0_t1_error/(resolution*resolution_2));

			hist_d0bar_1[rap]->SetBinContent(imass+1, total_d0bar);
			hist_d0bar_1[rap]->SetBinError(imass+1, total_d0bar_error);

		//D0 inclusive v1 vs mass
			d0_t1 = c2_d0_v1[rap][imass][2][2]->GetMean();
			d0_t1_error = c2_d0_v1[rap][imass][2][2]->GetMeanError();

			double total = (d0_t1/(resolution*resolution_2));
			double total_error = (d0_t1_error/(resolution*resolution_2));

			hist_d0d0bar_1[rap]->SetBinContent(imass+1, total);
			hist_d0d0bar_1[rap]->SetBinError(imass+1, total_error);

		}
	}


	TCanvas* c1 = new TCanvas("c1","c1",600,600);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	//gStyle->SetPadBorderMode(0.1);
	//gStyle->SetOptTitle(0);

	TH1D* base2 = makeHist("base2", "", "mass", "v^{odd}_{1}(y)", 100,1.7,2.0,kBlack);
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

	hist_d0_1[5]->SetMarkerStyle(20);
	hist_d0_1[5]->SetLineColor(kRed);
	hist_d0_1[5]->SetMarkerColor(kRed);
	hist_d0_1[5]->SetMarkerSize(1.5);

	hist_d0_1[5]->Draw("Psame");

	hist_d0bar_1[5]->SetMarkerStyle(20);
	hist_d0bar_1[5]->SetLineColor(kBlue);
	hist_d0bar_1[5]->SetMarkerColor(kBlue);
	hist_d0bar_1[5]->SetMarkerSize(1.5);

	hist_d0bar_1[5]->Draw("Psame");

	TLatex* r42 = new TLatex(0.2, 0.85, "PbPb 5.02 TeV");
    r42->SetNDC();
    r42->SetTextSize(23);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);

    TLatex* r43 = new TLatex(0.65,0.91, "CMS");
    r43->SetNDC();
    r43->SetTextSize(0.04);
    
    TLatex* r44 = new TLatex(0.75,0.91, "Preliminary");
    r44->SetNDC();
    r44->SetTextSize(22);
    r44->SetTextFont(53);

    TLatex* r45 = new TLatex(0.2, 0.8, "Cent.30-80%");
    r45->SetNDC();
    r45->SetTextSize(23);
    r45->SetTextFont(43);
    r45->SetTextColor(kBlack);

    TLatex* r47 = new TLatex(0.2, 0.75, "p_{T,D^{0}} > 2.0 GeV");
    r47->SetNDC();
    r47->SetTextSize(23);
    r47->SetTextFont(43);
    r47->SetTextColor(kBlack);

    r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");
    r45->Draw("same");
    r47->Draw("same");

	TLegend *w4 = new TLegend(0.23,0.20,0.70,0.26);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(22);
    w4->SetTextFont(45);

    TFile f1("output_v1VSmass_v17.root","RECREATE");

    for(int rap = 0; rap < 6; rap++){

    	hist_d0_1[rap]->Write();
    	hist_d0bar_1[rap]->Write();
    	hist_d0d0bar_1[rap]->Write();

    }

    for(int eta = 0; eta < 6; eta++){
		for(int charge = 0; charge < 3; charge++){

			D0Mass_Hist[eta][charge]->Write();
		}
	}



}