#include "RiceStyle.h"
#include "LinearSolver.h"

using namespace std;

// double f_obs[] = {0.0941031,0.186697,0.414384,0.40759,0.181193,0.0926985};
// double f_sb[] = {0.0165897,0.0354349,0.0881079,0.0862269,0.0345133,0.0164944};

// double f_obs_sig[]={0.0715718,0.140314,0.294066,0.291345,0.13919,0.0740226};
// double f_sb_sig[]={0.00111912,0.000417825,3.55151e-05,6.23979e-05,0.000439751,0.00127018};
// double f_obs_swap[]={0.018402,0.0448815,0.115502,0.11206,0.0430628,0.0192499};
// double f_sb_swap[]={0.0148678,0.034836,0.0859537,0.085172,0.0342337,0.0154975};


//v4
// double f_obs[]={0.516179,0.609793,0.723054,0.71735,0.623363,0.506904};
// double f_sb[]={0.13998,0.15446,0.191713,0.188131,0.162887,0.131097};

// double f_obs_sig[]={0.391226,0.457678,0.531815,0.524143,0.465996,0.372659};
// double f_sb_sig[]={0.0153115,0.00476865,0.000929062,0.000925871,0.00551677,0.0118153};
// double f_obs_swap[]={0.124953,0.152115,0.19124,0.193207,0.157367,0.134245};
// double f_sb_swap[]={0.124669,0.149692,0.190784,0.187206,0.15737,0.119282};

//v5
// double f_obs[]={0.528726,0.613236,0.422516,0.406017,0.626313,0.527612};
// double f_sb[]={0.145269,0.179533,0.083349,0.0737874,0.175505,0.152898};

// double f_obs_sig[]={0.38369,0.466746,0.309781,0.302316,0.488343,0.394301};
// double f_sb_sig[]={0.0119856,0.00428895,0.000148885,0.000317198,0.00733542,0.016279};
// double f_obs_swap[]={0.145036,0.14649,0.112735,0.103701,0.13797,0.133312};
// double f_sb_swap[]={0.133283,0.175245,0.0832001,0.0734702,0.16817,0.136619};

//v5+v7
double f_obs[]={0.533844,0.609318,0.412744,0.404786,0.616582,0.525866};
double f_sb[]={0.144471,0.177617,0.0807023,0.0736628,0.171656,0.154892};

double f_obs_sig[]={0.386281,0.465018,0.301731,0.301974,0.480418,0.394693};
double f_sb_sig[]={0.011275,0.00448473,0.000126761,0.000340304,0.00703724,0.0179306};
double f_obs_swap[]={0.147563,0.1443,0.111013,0.102812,0.136164,0.131173};
double f_sb_swap[]={0.133196,0.173132,0.0805756,0.0733225,0.164619,0.136962};

//v6
// double f_obs[]={0.163448,0.412486,0.483638,0.472148,0.406099,0.165195};
// double f_sb[]={0.0300882,0.097723,0.111672,0.105621,0.0898693,0.0347131};

// double f_obs_sig[]={0.141058,0.326773,0.352637,0.3514,0.324065,0.131222};
// double f_sb_sig[]={0.00385635,0.00293768,2.42017e-05,0.000166429,0.00282148,0.00403884};
// double f_obs_swap[]={0.0223908,0.0857123,0.131001,0.120749,0.0820336,0.0339731};
// double f_sb_swap[]={0.0262319,0.0947853,0.111647,0.105455,0.0870478,0.0306743};

/*
order of the fraction_mc
0: f_obs_sig_d0
1: f_obs_swap_d0
2: f_obs_sig_d0bar
3: f_obs_swap_d0bar
4: f_sb_sig_d0
5: f_sb_swap_d0
6: f_sb_sig_d0bar
7: f_sb_swap_d0bar

order of data v1
0: v1_obs_d0
1: v1_obs_d0bar
2: v1_side_d0
3: v1_side_d0bar
*/


void D0_v1_odd(){

	gStyle->SetErrorX(0);
	
/*
Loading all the histograms
 */

	TFile* file = new TFile("../rootfiles/D0_CDDF_v5plusv7_HI.root");

	TH1D* c2_cb_plus = (TH1D*) file->Get("ana/c2_cb_plus");
	TH1D* c2_ac_plus = (TH1D*) file->Get("ana/c2_ac_plus");;
	TH1D* c2_cb_minus = (TH1D*) file->Get("ana/c2_cb_minus");;
	TH1D* c2_ac_minus = (TH1D*) file->Get("ana/c2_ac_minus");;

	TH1D* c2_ab = (TH1D*) file->Get("ana/c2_ab");;

	TH1D* c2_a_real = (TH1D*) file->Get("ana/c2_a_real");;
	TH1D* c2_b_real = (TH1D*) file->Get("ana/c2_b_real");;
	TH1D* c2_c_plus_real = (TH1D*) file->Get("ana/c2_c_plus_real");;
	TH1D* c2_c_minus_real = (TH1D*) file->Get("ana/c2_c_minus_real");;

	TH1D* c2_a_imag = (TH1D*) file->Get("ana/c2_a_imag");;
	TH1D* c2_b_imag = (TH1D*) file->Get("ana/c2_b_imag");;
	TH1D* c2_c_plus_imag = (TH1D*) file->Get("ana/c2_c_plus_imag");;
	TH1D* c2_c_minus_imag = (TH1D*) file->Get("ana/c2_c_minus_imag");;

	TH1D* c2_v1[20][2][2];
	TH1D* c2_trk_accept[20][2][2];

	for(int eta = 0; eta < 6; eta++){
		for(int charge = 0; charge < 2; charge++){
			for(int dir = 0; dir < 2; dir++){

				c2_v1[eta][charge][dir] = (TH1D*) file->Get(Form("ana/c2_v1_%d_%d_%d",eta,charge,dir)); 
				c2_trk_accept[eta][charge][dir] = (TH1D*) file->Get(Form("ana/c2_trk_accept_%d_%d_%d",eta,charge,dir));

			}
		}
	}

	TH1D* c2_d0obs_v1[20][3][2];
	TH1D* c2_d0obs_trk_accept[20][3][2];

	TH1D* c2_d0bkg_v1[20][3][2];
	TH1D* c2_d0bkg_trk_accept[20][3][2];

	TH1D* D0Mass_Hist[20][3];

	for(int eta = 0; eta < 6; eta++){
		for(int charge = 0; charge < 3; charge++){
			
			D0Mass_Hist[eta][charge] = (TH1D*) file->Get(Form("ana/D0Mass_Hist_%d_%d",eta,charge));
			for(int dir = 0; dir < 2; dir++){

				c2_d0obs_v1[eta][charge][dir] = (TH1D*) file->Get(Form("ana/c2_d0obs_v1_%d_%d_%d",eta,charge,dir));
				c2_d0obs_trk_accept[eta][charge][dir] = (TH1D*) file->Get(Form("ana/c2_d0obs_trk_accept_%d_%d_%d",eta,charge,dir));
				c2_d0bkg_v1[eta][charge][dir] = (TH1D*) file->Get(Form("ana/c2_d0bkg_v1_%d_%d_%d",eta,charge,dir));
				c2_d0bkg_trk_accept[eta][charge][dir] = (TH1D*) file->Get(Form("ana/c2_d0bkg_trk_accept_%d_%d_%d",eta,charge,dir));

			}
		}
	}

/*
end of loading
 */
	double etabins[] = {-2.4,-1.4,-0.7,0.0,0.7,1.4,2.4};
	double ybins[] = {-2.0,-1.2,-0.6,0.0,0.6,1.2,2.0};
	TH1D* hist1[4];
	for(int sign = 0; sign < 4; sign++){		
		hist1[sign] = new TH1D(Form("hist1_%d",sign), "", 6, etabins);
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

	double Res_HFminus_TrackerPlus = sqrt(  fabs((Qac_plus-QaQc_plus_Real+QaQc_plus_Imag)*(Qab-QaQb_Real-QaQb_Imag)/(Qbc_plus-QbQc_plus_Real+QbQc_plus_Imag )) );
	double Res_HFplus_TrackerPlus = sqrt( fabs((Qbc_plus-QbQc_plus_Real+QbQc_plus_Imag)*(Qab-QaQb_Real-QaQb_Imag)/(Qac_plus-QaQc_plus_Real+QaQc_plus_Imag)) );

	double Res_HFplus_TrackerMinus = sqrt( fabs((Qbc_minus-QbQc_minus_Real+QbQc_minus_Imag)*(Qab-QaQb_Real-QaQb_Imag)/(Qac_minus-QaQc_minus_Real+QaQc_minus_Imag)) );
	double Res_HFminus_TrackerMinus = sqrt( fabs((Qac_minus-QaQc_minus_Real+QaQc_minus_Imag)*(Qab-QaQb_Real-QaQb_Imag)/(Qbc_minus-QbQc_minus_Real+QbQc_minus_Imag)) );

	cout << "Res_HFminus_TrackerPlus " << Res_HFminus_TrackerPlus << endl;
	cout << "Res_HFminus_TrackerMinus " << Res_HFminus_TrackerMinus << endl;
	cout << "Res_HFplus_TrackerPlus " << Res_HFplus_TrackerPlus << endl;
	cout << "Res_HFplus_TrackerMinus " << Res_HFplus_TrackerMinus << endl;
//*******************

	double QaQb_accept = (Qa_real * Qb_real + Qa_imag * Qb_imag);
	double Qab_error = c2_ab->GetMeanError();
	double resolution =  sqrt( Qab - QaQb_accept );
	cout << "acceptance: " << QaQb_accept << endl;
	cout << "resolution: " << resolution << endl;
	cout << "error: " << sqrt(Qab_error) << endl;

/*
two sub-events
 */
	// for(int eta = 0; eta < 6; eta++){

	// 		double t1 = c2_v1[eta][0][0]->GetMean();
	// 		double t1_error = c2_v1[eta][0][0]->GetMeanError();

	// 		double Q1_real = c2_trk_accept[eta][0][0]->GetMean();
	// 		double Q1_imag = c2_trk_accept[eta][0][1]->GetMean();

	// 		t1 = t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);

	// 		double t2 = c2_v1[eta][0][1]->GetMean();
	// 		double t2_error = c2_v1[eta][0][1]->GetMeanError();
			
	// 		t2 = t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

	// 		double total = (t1/resolution + t2/resolution)/2.0;
	// 		double total_error = (1/(2*resolution))*sqrt(t1_error*t1_error + t2_error*t2_error);

	// 		hist1[2]->SetBinContent(eta+1, total);
	// 		hist1[2]->SetBinError(eta+1, total_error);

	// 		double t1 = c2_v1[eta][1][0]->GetMean();
	// 		double t1_error = c2_v1[eta][1][0]->GetMeanError();

	// 		double Q1_real = c2_trk_accept[eta][1][0]->GetMean();
	// 		double Q1_imag = c2_trk_accept[eta][1][1]->GetMean();

	// 		t1 = t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);

	// 		double t2 = c2_v1[eta][1][1]->GetMean();
	// 		double t2_error = c2_v1[eta][1][1]->GetMeanError();

	// 		t2 = t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

	// 		double total = (t1/resolution + t2/resolution)/2.0;
	// 		double total_error = (1/(2*resolution))*sqrt(t1_error*t1_error + t2_error*t2_error);

	// 		hist1[3]->SetBinContent(eta+1, total);
	// 		hist1[3]->SetBinError(eta+1, total_error);
		
	// }

/*
D0 v1 obs and bkg, charge independent
*/
	TH1D* hist_d0_1[6];//obs_d0, obs_d0bar, side_d0, side_d0bar, obs_d0_inclusive, side_d0_inclusive
	TH1D* hist_d0d0bar_1[6];//sig_d0, sig_d0bar, bkg_d0, bkg_d0bar
	for(int sign = 0; sign < 6; sign++){		
		hist_d0_1[sign] = new TH1D(Form("hist_d0_1_%d",sign), "", 6, ybins);
		hist_d0d0bar_1[sign] = new TH1D(Form("hist_d0d0bar_1_%d",sign), "", 6, ybins);

	}

	for(int rap = 0; rap < 6; rap++){

		vector<double> frac_mc;
		frac_mc.push_back( f_obs_sig[rap] );
		frac_mc.push_back( f_obs_swap[rap] );
		frac_mc.push_back( f_obs_sig[rap] );//use d0 as d0bar
		frac_mc.push_back( f_obs_swap[rap] );//use d0 as d0bar

		frac_mc.push_back( f_sb_sig[rap] );
		frac_mc.push_back( f_sb_swap[rap] );
		frac_mc.push_back( f_sb_sig[rap] );//use d0 as d0bar
		frac_mc.push_back( f_sb_swap[rap] );//use d0 as d0bar
	//obs D0

		double d0obs_t1 = c2_d0obs_v1[rap][0][0]->GetMean();
		double d0obs_t1_error = c2_d0obs_v1[rap][0][0]->GetMeanError();

		double Q1_real = c2_d0obs_trk_accept[rap][0][0]->GetMean();
		double Q1_imag = c2_d0obs_trk_accept[rap][0][1]->GetMean();

		d0obs_t1 = d0obs_t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);

		double d0obs_t2 = c2_d0obs_v1[rap][0][1]->GetMean();
		double d0obs_t2_error = c2_d0obs_v1[rap][0][1]->GetMeanError();

		d0obs_t2 = d0obs_t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

		double total_obs_d0 = (d0obs_t1/resolution + d0obs_t2/resolution)/2.0;
		double total_obs_d0_error = (1/(2*resolution))*sqrt(d0obs_t1_error*d0obs_t1_error + d0obs_t2_error*d0obs_t2_error);

		hist_d0_1[0]->SetBinContent(rap+1, total_obs_d0);
		hist_d0_1[0]->SetBinError(rap+1, total_obs_d0_error);

	//obs D0bar

		d0obs_t1 = c2_d0obs_v1[rap][1][0]->GetMean();
		d0obs_t1_error = c2_d0obs_v1[rap][1][0]->GetMeanError();

		Q1_real = c2_d0obs_trk_accept[rap][1][0]->GetMean();
		Q1_imag = c2_d0obs_trk_accept[rap][1][1]->GetMean();

		d0obs_t1 = d0obs_t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);

		d0obs_t2 = c2_d0obs_v1[rap][1][1]->GetMean();
		d0obs_t2_error = c2_d0obs_v1[rap][1][1]->GetMeanError();

		d0obs_t2 = d0obs_t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

		double total_obs_d0bar = (d0obs_t1/resolution + d0obs_t2/resolution)/2.0;
		double total_obs_d0bar_error = (1/(2*resolution))*sqrt(d0obs_t1_error*d0obs_t1_error + d0obs_t2_error*d0obs_t2_error);

		hist_d0_1[1]->SetBinContent(rap+1, total_obs_d0bar);
		hist_d0_1[1]->SetBinError(rap+1, total_obs_d0bar_error);

	//bkg D0:

		double d0bkg_t1 = c2_d0bkg_v1[rap][0][0]->GetMean();
		double d0bkg_t1_error = c2_d0bkg_v1[rap][0][0]->GetMeanError();

		Q1_real = c2_d0bkg_trk_accept[rap][0][0]->GetMean();
		Q1_imag = c2_d0bkg_trk_accept[rap][0][1]->GetMean();

		d0bkg_t1 = d0bkg_t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);

		double d0bkg_t2 = c2_d0bkg_v1[rap][0][1]->GetMean();
		double d0bkg_t2_error = c2_d0bkg_v1[rap][0][1]->GetMeanError();

		d0bkg_t2 = d0bkg_t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

		double total_bkg_d0 = (d0bkg_t1/resolution + d0bkg_t2/resolution)/2.0;
		double total_bkg_d0_error = (1/(2*resolution))*sqrt(d0bkg_t1_error*d0bkg_t1_error + d0bkg_t2_error*d0bkg_t2_error);

		hist_d0_1[2]->SetBinContent(rap+1, total_bkg_d0);
		hist_d0_1[2]->SetBinError(rap+1, total_bkg_d0_error);

	//bkg D0bar:

		d0bkg_t1 = c2_d0bkg_v1[rap][1][0]->GetMean();
		d0bkg_t1_error = c2_d0bkg_v1[rap][1][0]->GetMeanError();

		Q1_real = c2_d0bkg_trk_accept[rap][1][0]->GetMean();
		Q1_imag = c2_d0bkg_trk_accept[rap][1][1]->GetMean();

		d0bkg_t1 = d0bkg_t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);

		d0bkg_t2 = c2_d0bkg_v1[rap][1][1]->GetMean();
		d0bkg_t2_error = c2_d0bkg_v1[rap][1][1]->GetMeanError();

		d0bkg_t2 = d0bkg_t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

		double total_bkg_d0bar = (d0bkg_t1/resolution + d0bkg_t2/resolution)/2.0;
		double total_bkg_d0bar_error = (1/(2*resolution))*sqrt(d0bkg_t1_error*d0bkg_t1_error + d0bkg_t2_error*d0bkg_t2_error);
		
		hist_d0_1[3]->SetBinContent(rap+1, total_bkg_d0bar);
		hist_d0_1[3]->SetBinError(rap+1, total_bkg_d0bar_error);

		vector<double> v1_data;
		v1_data.push_back(  total_obs_d0  );
		v1_data.push_back(  total_obs_d0bar  );
		v1_data.push_back(  total_bkg_d0  );
		v1_data.push_back(  total_bkg_d0bar  );

		vector<double> v1_data_error;
		v1_data_error.push_back(  total_obs_d0_error  );
		v1_data_error.push_back(  total_obs_d0bar_error  );
		v1_data_error.push_back(  total_bkg_d0_error  );
		v1_data_error.push_back(  total_bkg_d0bar_error  );

	//obs D0 inclusive
		d0obs_t1 = c2_d0obs_v1[rap][2][0]->GetMean();
		d0obs_t1_error = c2_d0obs_v1[rap][2][0]->GetMeanError();

		Q1_real = c2_d0obs_trk_accept[rap][2][0]->GetMean();
		Q1_imag = c2_d0obs_trk_accept[rap][2][1]->GetMean();

		d0obs_t1 = d0obs_t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);

		d0obs_t2 = c2_d0obs_v1[rap][2][1]->GetMean();
		d0obs_t2_error = c2_d0obs_v1[rap][2][1]->GetMeanError();

		d0obs_t2 = d0obs_t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

		double total_obs = (d0obs_t1/resolution + d0obs_t2/resolution)/2.0;
		double total_obs_error = (1/(2*resolution))*sqrt(d0obs_t1_error*d0obs_t1_error + d0obs_t2_error*d0obs_t2_error);

		hist_d0_1[4]->SetBinContent(rap+1, total_obs);
		hist_d0_1[4]->SetBinError(rap+1, total_obs_error);


	//bkg D0 inclusive:
		d0bkg_t1 = c2_d0bkg_v1[rap][2][0]->GetMean();
		d0bkg_t1_error = c2_d0bkg_v1[rap][2][0]->GetMeanError();

		Q1_real = c2_d0bkg_trk_accept[rap][2][0]->GetMean();
		Q1_imag = c2_d0bkg_trk_accept[rap][2][1]->GetMean();

		d0bkg_t1 = d0bkg_t1 - (Q1_real*Qa_real-Q1_imag*Qa_imag);

		d0bkg_t2 = c2_d0bkg_v1[rap][2][1]->GetMean();
		d0bkg_t2_error = c2_d0bkg_v1[rap][2][1]->GetMeanError();

		d0bkg_t2 = d0bkg_t2 - (Q1_real*Qb_real-Q1_imag*Qb_imag);

		double total_bkg = (d0bkg_t1/resolution + d0bkg_t2/resolution)/2.0;
		double total_bkg_error = (1/(2*resolution))*sqrt(d0bkg_t1_error*d0bkg_t1_error + d0bkg_t2_error*d0bkg_t2_error);

		hist_d0_1[5]->SetBinContent(rap+1, total_bkg);
		hist_d0_1[5]->SetBinError(rap+1, total_bkg_error);
	
	//v1 sig and bkg inclusive
		double v1b = (total_bkg*f_obs[rap] - total_obs*f_sb[rap])/(f_obs[rap] - f_sb[rap]);
		double v1be = sqrt((total_bkg_error*f_obs[rap])**2 + (total_obs_error*f_sb[rap])**2)/(f_obs[rap] - f_sb[rap]);

		double v1t = (total_obs - v1b*(1-f_obs[rap]))/f_obs[rap];
		double v1te = sqrt((v1be*(1-f_obs[rap]))**2 + total_obs_error**2)/f_obs[rap];

		hist_d0d0bar_1[4]->SetBinContent(rap+1, v1t);
		hist_d0d0bar_1[4]->SetBinError(rap+1, v1te);

		hist_d0d0bar_1[5]->SetBinContent(rap+1, v1b);
		hist_d0d0bar_1[5]->SetBinError(rap+1, v1be);

    //v1 d0 and d0bar separately

        vector<double> final_v1 = GetDmesonSigBkg(frac_mc, v1_data, v1_data_error);

        hist_d0d0bar_1[0]->SetBinContent(rap+1, final_v1[0]);
        hist_d0d0bar_1[0]->SetBinError(rap+1, final_v1[4]);

        cout << "d0 error: " << final_v1[4] << endl;


        hist_d0d0bar_1[1]->SetBinContent(rap+1, final_v1[1]);
      	hist_d0d0bar_1[1]->SetBinError(rap+1, final_v1[5]);

        cout << "d0bar error: " << final_v1[5] << endl;


        hist_d0d0bar_1[2]->SetBinContent(rap+1, final_v1[2]);
        hist_d0d0bar_1[2]->SetBinError(rap+1, final_v1[6]);

        // cout << "error D0: "

        hist_d0d0bar_1[3]->SetBinContent(rap+1, final_v1[3]);
        hist_d0d0bar_1[3]->SetBinError(rap+1, final_v1[7]);
	}


	TCanvas* c1 = new TCanvas("c1","c1",600,600);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	//gStyle->SetPadBorderMode(0.1);
	//gStyle->SetOptTitle(0);

	TH1D* base2 = makeHist("base2", "", "y", "v^{odd}_{1}(y)", 100,-3.0,3.0,kBlack);
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

	hist_d0d0bar_1[0]->SetMarkerStyle(20);
	hist_d0d0bar_1[0]->SetLineColor(kRed);
	hist_d0d0bar_1[0]->SetMarkerColor(kRed);
	hist_d0d0bar_1[0]->SetMarkerSize(1.5);

	hist_d0d0bar_1[1]->SetMarkerStyle(20);
	hist_d0d0bar_1[1]->SetLineColor(kBlue);
	hist_d0d0bar_1[1]->SetMarkerColor(kBlue);
	hist_d0d0bar_1[1]->SetMarkerSize(1.5);

	hist_d0d0bar_1[0]->Draw("Psame");
	hist_d0d0bar_1[1]->Draw("Psame");

	hist_d0d0bar_1[4]->SetMarkerStyle(20);
	hist_d0d0bar_1[4]->SetMarkerColor(kRed);
	hist_d0d0bar_1[4]->SetLineColor(kBlue);
	hist_d0d0bar_1[4]->SetMarkerSize(1.4);

	hist_d0d0bar_1[5]->SetMarkerStyle(20);
	hist_d0d0bar_1[5]->SetMarkerColor(kBlue);
	hist_d0d0bar_1[5]->SetMarkerSize(1.4);

	//hist_d0d0bar_1[4]->Draw("Psame");
	//hist_d0d0bar_1[5]->Draw("Psame");

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

    w4->AddEntry(hist_d0d0bar_1[0], "D^{0}  ", "P");
    w4->AddEntry(hist_d0d0bar_1[1], "#bar{D^{0}}","P");
    //w4->AddEntry(hist_d0d0bar_1[4], "D^{0}+#bar{D^{0}}","P");
    w4->Draw("same");

    TCanvas* c2 = new TCanvas("c2","c2",600,600);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	TH1D* base3 = makeHist("base3", "Pb-going", "y", "v^{D^{0}}_{1}#minusv^{#bar{D^{0}}}_{1}", 100,-3.0,3.0,kBlack);
	base3->GetYaxis()->SetRangeUser(-0.3, 0.3.);
	base3->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base3,1.1,1.25);

	base3->GetYaxis()->SetTitleOffset(1.3);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.6);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.6);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.6);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.6);
	base3->GetXaxis()->SetNdivisions(4,6,0);
	base3->GetYaxis()->SetNdivisions(4,6,0);

	TH1D* ratio = (TH1D*) hist_d0d0bar_1[0]->Clone("ratio");

	ratio->Add( hist_d0d0bar_1[1], -1 );
	
	ratio->Fit("pol1");
	ratio->SetMarkerStyle(20);
	ratio->SetMarkerColor(kBlack);
	ratio->SetLineColor(kBlack);
	//ratio->Draw("Psame");
	base3->Draw();
	ratio->Draw("Psame");

	TF1* func1 = (TF1*) ratio->GetFunction("pol1");
	func1->SetLineColor(kGreen-3);
	func1->SetMarkerColor(kGreen-3);
	double slope = func1->GetParameter(1);
	double slope_error = func1->GetParError(1);
	TH1D *hint1 = new TH1D("hint1","Fitted gaussian with .95 conf.band", 1000, -3, 3.0);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint1, 0.68);
     //Now the "hint1" histogram has the fitted function values as the
     //bin contents and the confidence intervals as bin errors
    hint1->SetStats(kFALSE);
    hint1->SetFillColor(kGreen-3);
    hint1->SetMarkerColor(kGreen-3);
    hint1->SetFillStyle(1001);
    hint1->SetFillColorAlpha(kGreen-3,0.4);
    hint1->Draw("e3 same");

	r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");
    r45->Draw("same");
    r47->Draw("same");

    TLatex* r48 = new TLatex(0.2, 0.7, Form("Slope: %.3f +/- %.3f",slope, slope_error ));
    r48->SetNDC();
    r48->SetTextSize(23);
    r48->SetTextFont(43);
    r48->SetTextColor(kBlack);
    r48->Draw("same");

    c1->Print("../plots/v1_D0_30_80_update1.pdf");
    c2->Print("../plots/v1_D0_diff_30_80_update1.pdf");



}