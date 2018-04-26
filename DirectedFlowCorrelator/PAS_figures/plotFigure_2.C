#include "../macros/RiceStyle.h"

using namespace std;

double ybins_center[] = {-1.6,-0.8,-0.3,0.3,0.8,1.6};

void plotFigure_2(){

    // double total_inclusive_D0D0BAR_sys = sqrt(0.005*0.005 + 0.02*0.02 + 0.008*0.008 + 0.01*0.01 ); 
    // double total_inclusive_D0_sys = sqrt(0.005*0.005 + 0.025*0.025 + 0.034*0.034 + 0.01*0.01 ); 
    // double total_inclusive_D0BAR_sys = sqrt(0.005*0.005 + 0.035*0.035 + 0.015*0.015 + 0.01*0.01 ); 

    double total_D0_sys_slope = 0.048;
    double total_D0bar_sys_slope = 0.053;
    double total_D0D0bar_sys_slope = 0.018;
    double total_diff_sys_slope = 0.091;

	gStyle->SetErrorX(0);

	TFile* file1 = new TFile("../macros/V1_odd_chargedParticles.root");
	TH1D* cV1_plus = (TH1D*) file1->Get("hist1_2");
	TH1D* cV1_minus = (TH1D*) file1->Get("hist1_3");
	TH1D* cV1_comb = (TH1D*) file1->Get("charge_inde");
	TH1D* cV1_diff = (TH1D*) file1->Get("ratio");

    cV1_plus->SetBinContent(1,100); cV1_plus->SetBinContent(8,100);
    cV1_minus->SetBinContent(1,100); cV1_minus->SetBinContent(8,100);
    cV1_comb->SetBinContent(1,100); cV1_comb->SetBinContent(8,100);

	TFile* file2 = new TFile("../macros/v1vsy_combined_v17.root" );
	TFile* file3 = new TFile("../macros/v1vsy_test_v17.root" );

	TGraphErrors* gr1_d0d0bar = (TGraphErrors*) file2->Get("v1vsy");
	TGraphErrors* gr1_bkg = (TGraphErrors*) file2->Get("v1bkgvsy");

	TGraphErrors* gr2_d0 = (TGraphErrors*) file3->Get("v1vsy");
	TGraphErrors* gr2_dobar = (TGraphErrors*) file3->Get("v1vsy_anti");
	TGraphErrors* gr2_diff = (TGraphErrors*) file3->Get("v1vsy_diff");
	TGraphErrors* gr2_bkg = (TGraphErrors*) file3->Get("v1bkgvsy");


	TCanvas* c1 = new TCanvas("c1","c1",600,600);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	//gStyle->SetPadBorderMode(0.1);
	//gStyle->SetOptTitle(0);

	TH1D* base2 = makeHist("base2", "", "y^{D^{0}},#eta^{ch}", "v^{odd}_{1}", 100,-2.2,2.2,kBlack);
	base2->GetYaxis()->SetRangeUser(-0.17, 0.17);
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

    //drawBoxTGraph(gr1_d0d0bar, 6, total_inclusive_D0D0BAR_sys, false, true);
    TBox* box1[6];
    for(int eta = 0; eta < 6; eta++){

        double xe,ye;
        double width = 0.05;
        double sys = total_D0D0bar_sys_slope*ybins_center[eta];
        gr1_d0d0bar->GetPoint(eta,xe,ye);

        box1[eta] = new TBox(xe-width,ye-sys,xe+width,ye+sys);
        box1[eta]->SetFillColor(kRed);
        box1[eta]->SetFillColorAlpha(kGray+2,0.4);
        box1[eta]->SetFillStyle(1001);
        box1[eta]->SetLineWidth(0);
        box1[eta]->SetLineColor(kRed);
        box1[eta]->Draw("SAME");
    }
    //

    drawBox(cV1_comb, 0.06, true);

	cV1_comb->SetMarkerStyle(25);
	cV1_comb->Draw("Psame");

	//gr1_d0d0bar->Fit("pol0","","",0,2.0);
	gr1_d0d0bar->SetMarkerStyle(20);
	gr1_d0d0bar->SetMarkerColor(kRed);
	gr1_d0d0bar->SetLineColor(kRed);
	gr1_d0d0bar->SetMarkerSize(1.4);

	gr1_bkg->SetMarkerStyle(20);
	gr1_bkg->SetMarkerColor(kBlue);
	gr1_bkg->SetMarkerSize(1.4);

	gr1_d0d0bar->Draw("Psame");
	//gr1_bkg->Draw("Psame");

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

	TLegend *w4 = new TLegend(0.23,0.20,0.6,0.32);
    w4->SetLineColor(kWhite);
    w4->SetFillColor(0);
    w4->SetTextSize(22);
    w4->SetTextFont(45);

    w4->AddEntry(gr1_d0d0bar, "D^{0}+#bar{D^{0}}  ", "P");
    w4->AddEntry(cV1_comb, "h^{+}+h^{#minus}","P");
    //w4->AddEntry(gr1_d0d0bar, "D^{0}+#bar{D^{0}}","P");
    w4->Draw("same");


    TCanvas* c2 = new TCanvas("c2","c2",600,600);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);

	TH1D* base1 = (TH1D*)base2->Clone("base1");
	base1->GetYaxis()->SetRangeUser(-0.29,0.29);
	base1->Draw();

	cV1_plus->SetMarkerStyle(25);
	cV1_plus->SetLineColor(kRed);
	cV1_plus->SetMarkerColor(kRed);
	cV1_plus->SetMarkerSize(1.5);

	cV1_minus->SetMarkerStyle(25);
	cV1_minus->SetLineColor(kBlue);
	cV1_minus->SetMarkerColor(kBlue);
	cV1_minus->SetMarkerSize(1.5);

	gr2_d0->SetMarkerStyle(20);
	gr2_d0->SetLineColor(kRed);
	gr2_d0->SetMarkerColor(kRed);
	gr2_d0->SetMarkerSize(1.5);

	gr2_dobar->SetMarkerStyle(20);
	gr2_dobar->SetLineColor(kBlue);
	gr2_dobar->SetMarkerColor(kBlue);
	gr2_dobar->SetMarkerSize(1.5);

    // drawBoxTGraph(gr2_d0, 6, total_inclusive_D0_sys, false, true);
    // drawBoxTGraph(gr2_dobar, 6, total_inclusive_D0BAR_sys, false, true);
   
    TBox* box1[6];
    TBox* box2[6];
    for(int eta = 0; eta < 6; eta++){

        double xe,ye;
        double width = 0.05;
        double sys = total_D0_sys_slope*ybins_center[eta];
        gr2_d0->GetPoint(eta,xe,ye);

        box1[eta] = new TBox(xe-width,ye-sys,xe+width,ye+sys);
        box1[eta]->SetFillColor(kRed);
        box1[eta]->SetFillColorAlpha(kGray+2,0.4);
        box1[eta]->SetFillStyle(1001);
        box1[eta]->SetLineWidth(0);
        box1[eta]->SetLineColor(kRed);
        box1[eta]->Draw("SAME");

        gr2_dobar->GetPoint(eta,xe,ye);
        double sys = total_D0bar_sys_slope*ybins_center[eta];

        box2[eta] = new TBox(xe-width,ye-sys,xe+width,ye+sys);
        box2[eta]->SetFillColor(kRed);
        box2[eta]->SetFillColorAlpha(kGray+2,0.4);
        box2[eta]->SetFillStyle(1001);
        box2[eta]->SetLineWidth(0);
        box2[eta]->SetLineColor(kRed);
        box2[eta]->Draw("SAME");
    }

    //
    drawBox(cV1_plus, 0.06, true);
    drawBox(cV1_minus, 0.06, true);

	gr2_d0->Draw("Psame");
	gr2_dobar->Draw("Psame");
	cV1_plus->Draw("Psame");
	cV1_minus->Draw("Psame");

 	r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");
    r45->Draw("same");
    r47->Draw("same");

    TLegend *w5 = new TLegend(0.23,0.18,0.70,0.35);
    w5->SetLineColor(kWhite);
    w5->SetFillColor(0);
    w5->SetTextSize(22);
    w5->SetTextFont(45);

    w5->AddEntry(gr2_d0, "D^{0}", "P");
    w5->AddEntry(gr2_dobar, "#bar{D^{0}}","P");
    w5->AddEntry(cV1_plus, "h^{+}", "P");
    w5->AddEntry(cV1_minus, "h^{#minus}","P");
    w5->Draw("same");

	TCanvas* c3 = new TCanvas("c3","c3",600,600);
    gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.13);
	gStyle->SetPadBorderMode(0.1);
	gStyle->SetOptTitle(0);

	TH1D* base3 = makeHist("base3", "Pb-going", "y^{D^{0}},#eta^{ch}", "v^{D^{0}/h^{+}}_{1}#minusv^{#bar{D^{0}}/h^{#minus}}_{1}", 100,-2.2,2.2,kBlack);
	base3->GetYaxis()->SetRangeUser(-0.37, 0.37);
	base3->GetXaxis()->SetTitleColor(kBlack);
	
	fixedFontHist1D(base3,1.1,1.25);

	base3->GetYaxis()->SetTitleOffset(1.3);
	base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.6);
	base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.6);
	base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.6);
	base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.6);
	base3->GetXaxis()->SetNdivisions(4,6,0);
	base3->GetYaxis()->SetNdivisions(4,6,0);


	TF1* f1 = new TF1("f1","[0]*x[0]",-2,2);
    f1->SetLineStyle(2);
    f1->SetLineColor(kRed);
	gr2_diff->Fit("f1");
	gr2_diff->SetMarkerStyle(20);
	gr2_diff->SetMarkerSize(1.6);
	gr2_diff->SetMarkerColor(kRed);
	gr2_diff->SetLineColor(kRed);
	
	TF1* func1 = (TF1*) gr2_diff->GetFunction("f1");
	double slope1 = func1->GetParameter(0);
	double slope_error1 = func1->GetParError(0);
    // double intercept1 = func1->GetParameter(1);
    // double intercept_error1 = func1->GetParError(1);
	// TH1D *hint1 = new TH1D("hint1","Fitted gaussian with .95 conf.band", 100, -2, 2);
 //    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint1, 0.68);
 //     //Now the "hint1" histogram has the fitted function values as the
 //     //bin contents and the confidence intervals as bin errors
 //    hint1->SetStats(kFALSE);
 //    hint1->SetFillColor(kGreen);
 //    hint1->SetMarkerColor(kGreen);
 // 	hint1->SetMarkerStyle(0);
 //    hint1->SetFillStyle(1001);
 //    hint1->SetFillColorAlpha(kGreen,0.4);

 //    TH1D *hint1_1 = new TH1D("hint1_1","Fitted gaussian with .95 conf.band", 100, -2, 2);
 //    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint1_1, 0.95);
 //     //Now the "hint1_1" histogram has the fitted function values as the
 //     //bin contents and the confidence intervals as bin errors
 //    hint1_1->SetStats(kFALSE);
 //    hint1_1->SetFillColor(kGray);
 //    hint1_1->SetMarkerColor(kGray);
 // 	hint1_1->SetMarkerStyle(0);
 //    hint1_1->SetFillStyle(1001);
 //    hint1_1->SetFillColorAlpha(kGray,0.4);
	
	TF1* f2 = new TF1("f2","[0]*x[0]",-2,2);
    f2->SetLineStyle(2);
    f2->SetLineColor(kBlue);
	cV1_diff->Fit("f2");
    cV1_diff->SetMarkerStyle(24);
	cV1_diff->SetMarkerColor(kBlue);

	TF1* func2 = (TF1*) cV1_diff->GetFunction("f2");
	double slope2 = func2->GetParameter(0);
	double slope_error2 = func2->GetParError(0);
    // double intercept2 = func2->GetParameter(1);
    // double intercept_error2 = func2->GetParError(1);
	// TH1D *hint2 = new TH1D("hint2","Fitted gaussian with .95 conf.band", 100, -2.4, 2.4);
 //    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint2, 0.68);
 //     //Now the "hint2" histogram has the fitted function values as the
 //     //bin contents and the confidence intervals as bin errors
 //    hint2->SetStats(kFALSE);
 //    hint2->SetFillColor(kBlue);
 //    hint2->SetMarkerColor(kBlue);
 //    hint2->SetMarkerStyle(0);
 //    hint2->SetFillStyle(1001);
 //    hint2->SetFillColorAlpha(kBlue,0.4);

    cV1_diff->SetBinContent(1,100); cV1_diff->SetBinContent(8,100);

    base3->Draw();

    TBox*box1[6];
    TBox*box2[6];

    double v1_default_value[6];
    double v1_default_error[6];

    double v1_value[6][2];
    double v1_value_error[6][2];

    // for(int i = 0; i < 6; i++){

    //     double xe = 0.05;
    //     double x1,value1;
    //     gr2_diff->GetPoint(i,x1,value1);
    //     double ye = sqrt(total_inclusive_D0_sys*total_inclusive_D0_sys - total_inclusive_D0BAR_sys*total_inclusive_D0BAR_sys);

    //     v1_value[i][0] = value1+ye;
    //     v1_value[i][1] = value1-ye;
    //     v1_value_error[i][0] = gr2_diff->GetErrorY(i);
    //     v1_value_error[i][1] = gr2_diff->GetErrorY(i);

    //     v1_default_value[i] = value1;
    //     v1_default_error[i] = gr2_diff->GetErrorY(i);

    //     box1[i] = new TBox(x1-xe,value1-ye,x1+xe,value1+ye);
    //     box1[i]->SetFillColorAlpha(kGray+2,0.4);
    //     box1[i]->SetFillStyle(1001);
    //     box1[i]->SetLineWidth(0);
    //     box1[i]->SetLineColor(kRed);
    //     box1[i]->Draw("SAME");

    //     double xe = 0.05;
    //     double x1,value1;
    //     x1 = cV1_diff->GetBinCenter(i+2);
    //     value1 = cV1_diff->GetBinContent(i+2);
    //     double ye = sqrt(2)*value1*0.06;

    //     box2[i] = new TBox(x1-xe,value1-ye,x1+xe,value1+ye);
    //     box2[i]->SetFillColorAlpha(kGray+2,0.4);
    //     box2[i]->SetFillStyle(1001);
    //     box2[i]->SetLineWidth(0);
    //     box2[i]->SetLineColor(kRed);
    //     box2[i]->Draw("SAME");
    // }

    TBox* box1[6];
    for(int eta = 0; eta < 6; eta++){

        double xe,ye;
        double width = 0.05;
        double sys = total_diff_sys_slope*ybins_center[eta];
        gr2_diff->GetPoint(eta,xe,ye);

        v1_value[eta][0] = ye+sys;
        v1_value[eta][1] = ye-sys;
        v1_value_error[eta][0] = gr2_diff->GetErrorY(eta);
        v1_value_error[eta][1] = gr2_diff->GetErrorY(eta);

        box1[eta] = new TBox(xe-width,ye-sys,xe+width,ye+sys);
        box1[eta]->SetFillColor(kRed);
        box1[eta]->SetFillColorAlpha(kGray+2,0.4);
        box1[eta]->SetFillStyle(1001);
        box1[eta]->SetLineWidth(0);
        box1[eta]->SetLineColor(kRed);
        box1[eta]->Draw("SAME");
    }

//test the fit in chi2
/*
    double chi2_set[10000];
    double slope_value[10000];
    double rap_bins[6] = {-1.6,-0.9,-0.3,0.3,0.9,1.6};
    double A = 0.0001;
    double y_expected[6];

    for(int count = 0; count < 10000; count++){
        double chi2 = 0.;
        for(int i = 0; i < 6; i++){

            y_expected[i] = A*rap_bins[i];
            chi2 = chi2 + ((v1_default_value[i] - y_expected[i])*(v1_default_value[i] - y_expected[i]))/(v1_default_error[i]*v1_default_error[i]);
        }
        chi2_set[count] = chi2;
        slope_value[count] = A;
        A = A + 0.0001;

    }

    double min =10000.0;
    double min_slope = 0.0;
    for(int count = 0; count < 10000; count++){

        if( min > chi2_set[count]) 
        {
            min = chi2_set[count];
            min_slope = slope_value[count];
        }
    }

    cout << "min chi2" << min << endl;
    cout << "correspoinding slope " << min_slope << endl;

    double chi2_set_error[20000];
    double min_error_set[20000];
    double derror = -1.0;
    for(int count = 0; count < 20000; count++){
        double chi2 = 0.;
        for(int i = 0; i < 6; i++){

            y_expected[i] = (0.0646+derror)*rap_bins[i];
            chi2 = chi2 + ((v1_default_value[i] - y_expected[i])*(v1_default_value[i] - y_expected[i]))/(v1_default_error[i]*v1_default_error[i]);
        }
        chi2_set_error[count] = chi2;
        min_error_set[count] = derror;
        derror += 0.0001;
    }

    double min_dchi2 =20000.0;
    double min_error = 0.0;

    TH1D* hist_error = new TH1D("hist1","hist1",20000,-1.0,1.0);
    for(int count = 0; count < 20000; count++){

        double diff = chi2_set_error[count] - min;
        if(  diff > 3.99 && diff < 4.01 ) 
        {
            min_dchi2 = chi2_set_error[count];
            min_error = min_error_set[count];
        }
        hist_error->SetBinContent(count+1, diff);
    }

    cout << "min chi2 error" << min_dchi2 << endl;
    cout << "correspoinding error " << min_error << endl;

*/


//systematics on slope

    // TGraphErrors* gr_fit[64];

    // double x_bin[6] = {-1.6,-0.9,-0.3,0.3,0.9,1.6};

    // for(int index = 0; index < 64; index++){
    //     gr_fit[index] = new TGraphErrors(6);
    // }
    // int index = 0;    
    // for(int i = 0; i < 2; i++){
    //     for(int j = 0; j < 2; j++){
    //         for(int k = 0; k < 2; k++){
    //             for(int l = 0; l < 2; l++){
    //                 for(int m = 0; m < 2; m++){
    //                     for(int p = 0; p < 2; p++){
    //                         gr_fit[index]->SetPoint(0, x_bin[0], v1_value[0][i]);
    //                         gr_fit[index]->SetPoint(1, x_bin[1], v1_value[1][j]);
    //                         gr_fit[index]->SetPoint(2, x_bin[2], v1_value[2][k]);
    //                         gr_fit[index]->SetPoint(3, x_bin[3], v1_value[3][l]);
    //                         gr_fit[index]->SetPoint(4, x_bin[4], v1_value[4][m]);
    //                         gr_fit[index]->SetPoint(5, x_bin[5], v1_value[5][p]);

    //                         index++;
    //                     }
    //                 }
    //             }   
    //         }
    //     }
    // }

    // double slope_para[64];
    // TF1* f3[64];
    // for(int index = 0; index < 64; index++){

    //     f3[index] = new TF1(Form("f3_%d",index),"[0]*x[0]",-2,2);
    //     gr_fit[index]->Fit(Form("f3_%d",index));
    //     TF1* myFunc = gr_fit[index]->GetFunction(Form("f3_%d",index));
    //     slope_para[index] = fabs(myFunc->GetParameter(0) - 0.065);
    // }

    // TH1D* hist_slope = new TH1D("hist_slope","hist_slope",100,0,0.5);
    // double max = 0.0;
    // for(int index = 0; index < 64; index++){

    //     hist_slope->Fill(slope_para[index]);

    //     //if( max < slope_para[index] ) max = slope_para[index];
    // }

    // cout << "max " << max << endl;


    // hint1->Draw("e4 same");
    // hint1_1->Draw("e4 same");
	gr2_diff->Draw("Psame");
	cV1_diff->Draw("Psame");

	r42->Draw("same");
    r43->Draw("same");
    r44->Draw("same");
    r45->Draw("same");
    r47->Draw("same");

    TLatex* r48 = new TLatex(0.35, 0.23, Form(" slope: %.3f #pm %.3f #pm %.3f ",slope1, slope_error1,0.053 ));
    r48->SetNDC();
    r48->SetTextSize(20);
    r48->SetTextFont(43);
    r48->SetTextColor(kBlack);
    r48->Draw("same");
    // TLatex* r48 = new TLatex(0.35, 0.23, Form(" intercept: %.3f #pm %.3f #pm %.3f ",intercept1, intercept_error1,0.027 ));
    // r48->SetNDC();
    // r48->SetTextSize(20);
    // r48->SetTextFont(43);
    // r48->SetTextColor(kBlack);
    // r48->Draw("same");
    TLatex* r49 = new TLatex(0.35, 0.17, Form(" slope: %.5f #pm %.5f #pm %.5f ",slope2, slope_error2,0.0 ));
    r49->SetNDC();
    r49->SetTextSize(20);
    r49->SetTextFont(43);
    r49->SetTextColor(kBlack);
    r49->Draw("same");
    // TLatex* r49 = new TLatex(0.35, 0.13, Form(" intercept: %.5f #pm %.5f #pm %.5f ",intercept2, intercept_error2,0.0 ));
    // r49->SetNDC();
    // r49->SetTextSize(20);
    // r49->SetTextFont(43);
    // r49->SetTextColor(kBlack);
    // r49->Draw("same");

    TLegend *w6 = new TLegend(0.17,0.16,0.34,0.26);
    w6->SetLineColor(kWhite);
    w6->SetFillColor(0);
    w6->SetTextSize(22);
    w6->SetTextFont(45);

    w6->AddEntry(gr2_diff, "D^{0}#minus#bar{D^{0}}  ", "P");
    w6->AddEntry(cV1_diff, "h^{+}#minush^{#minus}","P");
    w6->Draw("same");

    c1->Print("Figure_2_a_mixedHarmonics.pdf");
    c2->Print("Figure_2_b_mixedHarmonics.pdf");
    c3->Print("Figure_3_mixedHarmonics.pdf");

    // TCanvas* c4 = new TCanvas();
    // hist_error->Draw();


}