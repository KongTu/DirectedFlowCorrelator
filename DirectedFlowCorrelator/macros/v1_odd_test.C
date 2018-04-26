#include "RiceStyle.h"

using namespace std;

void v1_odd_test(){

	gStyle->SetErrorX(0);
	
	TFile* file = new TFile("../rootfiles/CDDF_v22.root");

	TH1D* c2_cb_plus = (TH1D*) file->Get("ana4/c2_cb_plus");
	TH1D* c2_ac_plus = (TH1D*) file->Get("ana4/c2_ac_plus");
	TH1D* c2_cb_minus = (TH1D*) file->Get("ana4/c2_cb_minus");
	TH1D* c2_ac_minus = (TH1D*) file->Get("ana4/c2_ac_minus");

	TH1D* c2_ab = (TH1D*) file->Get("ana4/c2_ab");

	TH1D* c2_a_real = (TH1D*) file->Get("ana4/c2_a_real");
	TH1D* c2_b_real = (TH1D*) file->Get("ana4/c2_b_real");
	TH1D* c2_c_plus_real = (TH1D*) file->Get("ana4/c2_c_plus_real");
	TH1D* c2_c_minus_real = (TH1D*) file->Get("ana4/c2_c_minus_real");

	TH1D* c2_a_imag = (TH1D*) file->Get("ana4/c2_a_imag");
	TH1D* c2_b_imag = (TH1D*) file->Get("ana4/c2_b_imag");
	TH1D* c2_c_plus_imag = (TH1D*) file->Get("ana4/c2_c_plus_imag");
	TH1D* c2_c_minus_imag = (TH1D*) file->Get("ana4/c2_c_minus_imag");

	TH1D* c2_v1[20][2][2];
	TH1D* c2_trk_accept[20][2][2];

	for(int eta = 0; eta < 8; eta++){
		for(int charge = 0; charge < 2; charge++){
			for(int dir = 0; dir < 2; dir++){

				c2_v1[eta][charge][dir] = (TH1D*) file->Get(Form("ana4/c2_v1_%d_%d_%d",eta,charge,dir)); 
				c2_trk_accept[eta][charge][dir] = (TH1D*) file->Get(Form("ana4/c2_trk_accept_%d_%d_%d",eta,charge,dir));

			}
		}
	}

	double etabins[] = {-2.4,-2.0,-1.2,-0.6,0.0,0.6,1.2,2.0,2.4};
	//double etabins[] = {-2.4,-2.0,-1.4,-1.0,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,1.0,1.4,2.0,2.4};
	TH1D* hist1[5];
	for(int sign = 0; sign < 5; sign++){		
		hist1[sign] = new TH1D(Form("hist1_%d",sign), "", 8, etabins);
	}

/*
3-sub events
*/
	double Qa_real = c2_a_real->GetMean();
	double Qb_real = c2_b_real->GetMean();
	double Qc_plus_real = c2_c_plus_real->GetMean();
	double Qc_minus_real = c2_c_minus_real->GetMean();

	double Qa_imag = c2_a_imag->GetMean();
	double Qb_imag = c2_b_imag->GetMean();
	double Qc_plus_imag = c2_c_plus_imag->GetMean();
	double Qc_minus_imag = c2_c_minus_imag->GetMean();

	double QaQb_Real = Qa_real*Qb_real;
	double QaQb_Imag = Qa_imag*Qb_imag;

//plus
	double QaQc_plus_Real = Qa_real*Qc_plus_real;
	double QaQc_plus_Imag = Qa_imag*Qc_plus_imag;

	double QbQc_plus_Real = Qb_real*Qc_plus_real;
	double QbQc_plus_Imag = Qb_imag*Qc_plus_imag;

//minus
	double QaQc_minus_Real = Qa_real*Qc_minus_real;
	double QaQc_minus_Imag = Qa_imag*Qc_minus_imag;

	double QbQc_minus_Real = Qb_real*Qc_minus_real;
	double QbQc_minus_Imag = Qb_imag*Qc_minus_imag;
	

	double Qab = c2_ab->GetMean();
	double Qac_plus = c2_ac_plus->GetMean();
	double Qbc_plus = c2_cb_plus->GetMean();

	double Qac_minus = c2_ac_minus->GetMean();
	double Qbc_minus = c2_cb_minus->GetMean();

	cout << "Qab " << (Qab-QaQb_Real-QaQb_Imag) << endl;
	cout << "Qac_plus " << (Qac_plus-QaQc_plus_Real+QaQc_plus_Imag) << endl;
	cout << "Qbc_plus " << (Qbc_plus-QbQc_plus_Real+QbQc_plus_Imag) << endl;
	cout << "Qac_minus " << (Qac_minus-QaQc_minus_Real+QaQc_minus_Imag) << endl;
	cout << "Qbc_minus " << (Qbc_minus-QbQc_minus_Real+QbQc_minus_Imag) << endl;


	double Res_HFminus_TrackerPlus = sqrt(  fabs((Qac_plus-QaQc_plus_Real+QaQc_plus_Imag)*(Qab-QaQb_Real-QaQb_Imag)/(Qbc_plus-QbQc_plus_Real+QbQc_plus_Imag )) );
	double Res_HFplus_TrackerPlus = sqrt( fabs((Qbc_plus-QbQc_plus_Real+QbQc_plus_Imag)*(Qab-QaQb_Real-QaQb_Imag)/(Qac_plus-QaQc_plus_Real+QaQc_plus_Imag)) );

	// double Res_HFminus_TrackerPlus = sqrt( (Qac_plus)*(Qab)/(Qbc_plus)  );
	// double Res_HFplus_TrackerPlus = sqrt(  (Qbc_plus)*(Qab)/(Qac_plus) );

	// double Res_HFplus_TrackerMinus = sqrt( Qbc_minus*Qab/Qac_minus );
	// double Res_HFminus_TrackerMinus = sqrt( Qac_minus*Qab/Qbc_minus );

	double Res_HFplus_TrackerMinus = sqrt( fabs((Qbc_minus-QbQc_minus_Real+QbQc_minus_Imag)*(Qab-QaQb_Real-QaQb_Imag)/(Qac_minus-QaQc_minus_Real+QaQc_minus_Imag)) );
	double Res_HFminus_TrackerMinus = sqrt( fabs((Qac_minus-QaQc_minus_Real+QaQc_minus_Imag)*(Qab-QaQb_Real-QaQb_Imag)/(Qbc_minus-QbQc_minus_Real+QbQc_minus_Imag)) );


	cout << "Res_HFminus_TrackerPlus " << Res_HFminus_TrackerPlus << endl;
	cout << "Res_HFminus_TrackerMinus " << Res_HFminus_TrackerMinus << endl;
	cout << "Res_HFplus_TrackerPlus " << Res_HFplus_TrackerPlus << endl;
	cout << "Res_HFplus_TrackerMinus " << Res_HFplus_TrackerMinus << endl;
//*******************
/*
Three sub-events
 */
	for(int eta = 0; eta < 8; eta++){

		if(eta < 4){

			double t1 = c2_v1[eta][0][0]->GetMean();
			double t1_error = c2_v1[eta][0][0]->GetMeanError();
			
			double Q1_real = c2_trk_accept[eta][0][0]->GetMean();
			double Q1_imag = c2_trk_accept[eta][0][1]->GetMean();

			t1 = t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);

			double t2 = c2_v1[eta][0][1]->GetMean();
			double t2_error = c2_v1[eta][0][1]->GetMeanError();

			t2 = t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

			double total = (t1/Res_HFminus_TrackerPlus + t2/Res_HFplus_TrackerPlus)/2.0;
			double total_error = 0.5*sqrt(t1_error*t1_error + t2_error*t2_error);

			hist1[0]->SetBinContent(eta+1, total);
			hist1[0]->SetBinError(eta+1, total_error);

			double t1 = c2_v1[eta][1][0]->GetMean();
			double t1_error = c2_v1[eta][1][0]->GetMeanError();

			double Q1_real = c2_trk_accept[eta][1][0]->GetMean();
			double Q1_imag = c2_trk_accept[eta][1][1]->GetMean();

			t1 = t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);

			double t2 = c2_v1[eta][1][1]->GetMean();
			double t2_error = c2_v1[eta][1][1]->GetMeanError();

			t2 = t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

			double total = (t1/Res_HFminus_TrackerPlus + t2/Res_HFplus_TrackerPlus)/2.0;
			double total_error = 0.5*sqrt(t1_error*t1_error + t2_error*t2_error);

			hist1[1]->SetBinContent(eta+1, total);
			hist1[1]->SetBinError(eta+1, total_error);


		}
		else{

			double t1 = c2_v1[eta][0][0]->GetMean();
			double t1_error = c2_v1[eta][0][0]->GetMeanError();

			double Q1_real = c2_trk_accept[eta][0][0]->GetMean();
			double Q1_imag = c2_trk_accept[eta][0][1]->GetMean();

			t1 = t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);
			
			double t2 = c2_v1[eta][0][1]->GetMean();
			double t2_error = c2_v1[eta][0][1]->GetMeanError();

			t2 = t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

			double total = (t1/Res_HFminus_TrackerMinus + t2/Res_HFplus_TrackerMinus)/2.0;
			double total_error = 0.5*sqrt(t1_error*t1_error + t2_error*t2_error);

			hist1[0]->SetBinContent(eta+1, total);
			hist1[0]->SetBinError(eta+1, total_error);

			double t1 = c2_v1[eta][1][0]->GetMean();
			double t1_error = c2_v1[eta][1][0]->GetMeanError();

			double Q1_real = c2_trk_accept[eta][1][0]->GetMean();
			double Q1_imag = c2_trk_accept[eta][1][1]->GetMean();

			t1 = t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);

			double t2 = c2_v1[eta][1][1]->GetMean();
			double t2_error = c2_v1[eta][1][1]->GetMeanError();

			t2 = t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

			double total = (t1/Res_HFminus_TrackerMinus + t2/Res_HFplus_TrackerMinus)/2.0;
			double total_error = 0.5*sqrt(t1_error*t1_error + t2_error*t2_error);

			hist1[1]->SetBinContent(eta+1, total);
			hist1[1]->SetBinError(eta+1, total_error);
		}
	}

	double QaQb_accept = (Qa_real * Qb_real + Qa_imag * Qb_imag);
	double resolution =  sqrt( Qab - QaQb_accept );
	cout << "acceptance: " << QaQb_accept << endl;
	cout << "resolution: " << resolution << endl;

/*
two sub-events
 */
	for(int eta = 0; eta < 8; eta++){

			double t1 = c2_v1[eta][0][0]->GetMean();
			double t1_error = c2_v1[eta][0][0]->GetMeanError();

			double Q1_real = c2_trk_accept[eta][0][0]->GetMean();
			double Q1_imag = c2_trk_accept[eta][0][1]->GetMean();

			t1 = t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);

			double t2 = c2_v1[eta][0][1]->GetMean();
			double t2_error = c2_v1[eta][0][1]->GetMeanError();
			
			t2 = t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

			double total = (t1/resolution + t2/resolution)/2.0;
			double total_error = (1/(2*resolution))*sqrt(t1_error*t1_error + t2_error*t2_error);

			hist1[2]->SetBinContent(eta+1, total);
			hist1[2]->SetBinError(eta+1, total_error);

			double t1 = c2_v1[eta][1][0]->GetMean();
			double t1_error = c2_v1[eta][1][0]->GetMeanError();

			double Q1_real = c2_trk_accept[eta][1][0]->GetMean();
			double Q1_imag = c2_trk_accept[eta][1][1]->GetMean();

			t1 = t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);

			double t2 = c2_v1[eta][1][1]->GetMean();
			double t2_error = c2_v1[eta][1][1]->GetMeanError();

			t2 = t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

			double total = (t1/resolution + t2/resolution)/2.0;
			double total_error = (1/(2*resolution))*sqrt(t1_error*t1_error + t2_error*t2_error);

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

	hist1[0]->SetMarkerStyle(25);
	hist1[0]->SetMarkerColor(kRed);
	hist1[0]->SetMarkerSize(1.4);

	hist1[1]->SetMarkerStyle(25);
	hist1[1]->SetMarkerColor(kBlue);
	hist1[1]->SetMarkerSize(1.4);

	hist1[2]->SetMarkerStyle(20);
	hist1[2]->SetMarkerColor(kRed);
	hist1[2]->SetMarkerSize(1.4);

	hist1[3]->SetMarkerStyle(20);
	hist1[3]->SetMarkerColor(kBlue);
	hist1[3]->SetMarkerSize(1.4);

	TH1D* charge_inde = (TH1D*) hist1[2]->Clone("charge_inde");
	charge_inde->Add(hist1[3], +1);
	charge_inde->Scale(0.5);
	//charge_inde->SetMarkerSize(1.6);
	charge_inde->SetMarkerColor(kBlue);
	charge_inde->SetLineColor(kBlue);

	//hist1[0]->Draw("Psame");
	//hist1[1]->Draw("Psame");

	hist1[2]->Draw("Psame");
    hist1[3]->Draw("Psame");

	//charge_inde->Draw("Psame");
	//hist1[3]->Draw("Psame");

	TFile* file1 = new TFile("v1Kong.root");
	TH1D* v1odd = (TH1D*) file1->Get("30_50/v1oddSUB2");

	v1odd->SetMarkerColor(kRed);
    v1odd->SetMarkerStyle(25);
	//v1odd->Draw("same");

	double bins[] = {-2.2,-1.8,-1.4,-1.0,-0.6,-0.2,0.2,0.6,1.0,1.4,1.8,2.2};
	double values_1[] = {0.00387749, 0.00896082, 0.00968393, 0.00792812, 0.00502962, 0.00176095, -0.00212748, -0.0054055, -0.00815283, -0.0102848, -0.00953422, -0.00452408};
	double errors_1[] = {0.000206731, 0.000200724, 0.000179953, 0.000162642, 0.000150315, 0.000146505, 0.000146353, 0.000149666, 0.000162489, 0.000180968, 0.000200114, 0.000210354};

	double values_2[] = {0.00720323, 0.012818, 0.0134219, 0.0105773, 0.00698389, 0.00252322, -0.0025126, -0.00726318, -0.011089, -0.0135221, -0.0130803, -0.00804368};

	TGraphErrors* gr = new TGraphErrors(12);
	for(int i = 0; i < 12; i++){
		gr->SetPoint(i, bins[i], (values_1[i]+values_2[i])/2.0);
		gr->SetPointError(i, 0.0, errors_1[i]);
	}

	gr->SetMarkerStyle(24);
	gr->SetMarkerColor(kBlack);
	gr->SetMarkerSize(1.4);

	//gr->Draw("Psame");

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

    TLatex* r45 = new TLatex(0.2, 0.8, "Cent.30-50%");
    r45->SetNDC();
    r45->SetTextSize(23);
    r45->SetTextFont(43);
    r45->SetTextColor(kBlack);

    TLatex* r47 = new TLatex(0.2, 0.75, "0.3 < p_{T} < 3.0 GeV");
    r47->SetNDC();
    r47->SetTextSize(23);
    r47->SetTextFont(43);
    r47->SetTextColor(kBlack);

    r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");
    r45->Draw("same");
    r47->Draw("same");

	TLegend *w4 = new TLegend(0.23,0.2,0.50,0.4);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(24);
    w4->SetTextFont(45);
    //w4->AddEntry(hist1[0], "Kong Pixel 3-sub +", "P");
    //w4->AddEntry(hist1[1], "Charge inclusive 3-sub ","P");
    //w4->AddEntry(hist1[2], "h^{+}", "P");
    //w4->AddEntry(hist1[3], "h^{-}","P");
    w4->AddEntry(charge_inde, "Charge inclusive 2-sub","P");
    //w4->AddEntry(v1odd, "Will's 2-sub Charge-indep", "P");
    w4->Draw("same");

	TCanvas* c2 = new TCanvas("c2","c2",700,700);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	TH1D* base3 = makeHist("base3", "Pb-going", "#eta", "v^{odd}_{1,+}#minusv^{odd}_{1,#minus}", 100,-3.0,3.0,kBlack);
	base3->GetYaxis()->SetRangeUser(-0.001, 0.001);
	base3->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base3,1.1,1.25);

	base3->GetYaxis()->SetTitleOffset(1.3);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.6);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.6);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.6);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.6);
	base3->GetXaxis()->SetNdivisions(4,6,0);
	base3->GetYaxis()->SetNdivisions(4,6,0);


	TH1D* ratio = (TH1D*) hist1[2]->Clone("ratio");

	ratio->Add( hist1[3], -1 );
	ratio->SetMarkerStyle(20);
	ratio->SetLineColor(kBlack);
	ratio->SetMarkerColor(kBlack);
	//ratio->Fit("pol1");

	base3->Draw();

	ratio->Draw("Psame");

	// TF1* func1 = (TF1*) ratio->GetFunction("pol1");
	// func1->SetLineStyle(2);
	// func1->SetLineColor(kGreen);
	// double slope = func1->GetParameter(1);
	// double slope_error = func1->GetParError(1);
	// TH1D *hint1 = new TH1D("hint1","Fitted gaussian with .95 conf.band", 1000, -3, 3.0);
 //    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint1, 0.68);
 //     //Now the "hint1" histogram has the fitted function values as the
 //     //bin contents and the confidence intervals as bin errors
 //    hint1->SetStats(kFALSE);
 //    hint1->SetFillColor(kGreen-3);
 //    hint1->SetMarkerColor(kGreen-3);
 //    hint1->SetFillStyle(1001);
 //    hint1->SetFillColorAlpha(kGreen-3,0.4);
 //    hint1->Draw("e3 same");
	// TLatex* r48 = new TLatex(0.2, 0.7, Form("Slope: %.6f +/- %.6f",slope, slope_error ));
 //    r48->SetNDC();
 //    r48->SetTextSize(23);
 //    r48->SetTextFont(43);
 //    r48->SetTextColor(kBlack);
 //    r48->Draw("same");

	r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");
    r45->Draw("same");
    r47->Draw("same");

    // TFile f("../crosschecks/V1_odd_chargedParticles_ana4.root","RECREATE");
    // hist1[2]->Write();
    // hist1[3]->Write();
    // charge_inde->Write();
    // ratio->Write();

    // c1->Print("../plots/v1_first.pdf");
    // c2->Print("../plots/v1_diff.pdf");


}