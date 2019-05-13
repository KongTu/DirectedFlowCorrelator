void perisubrefalltestD0true_swap(){
TH1::SetDefaultSumw2();

    TFile *_file0 = TFile::Open("D0v2_vspt_nosub_NassFit_refpt033_obs.root");
    TFile *_file1 = TFile::Open("D0v2_vspt_nosub_NassFit_refpt033_bkg.root");
    
//    double sigfrks[13] = {0.051,0.1,0.211,0.304,0.376,0.424,0.480,0.526};
//    double sigfrks_SB[13] = {0.000113,0.000145,0.000721,0.000440,0.000770,0.001028,0.00081,0.00094};
//    double sigfrks[13] = {0.0399,0.077,0.177,0.254015,0.320527,0.363829,0.419494,0.4699};

    double sigfrks_SB[13] = {0.000182,0.000336,0.000547,0.000895,0.000864,0.000892,0.0014,0.001165};
    double sigfrks[13] = {0.055729,0.100339,0.180785,0.271435,0.321383,0.365555,0.419614,0.466810};
    
    //read vn values
    TGraphErrors* v2obs_ks = _file0->Get("D0v2_obs_GplusPP");
    TGraphErrors* v2obs_ks_sub = _file0->Get("D0v2sub_obs_GplusPP");
    
    TGraphErrors* v2bkg_ks = _file1->Get("D0v2_bkg_GplusPP");
    TGraphErrors* v2bkg_ks_sub = _file1->Get("D0v2sub_bkg_GplusPP");

    TGraphErrors* v2obs_ks_KET = _file0->Get("D0v2_obs_GplusPP_KET");
    TGraphErrors* v2obs_ks_sub_KET = _file0->Get("D0v2sub_obs_GplusPP_KET");
    
    TGraphErrors* v2bkg_ks_KET = _file1->Get("D0v2_bkg_GplusPP_KET");
    TGraphErrors* v2bkg_ks_sub_KET = _file1->Get("D0v2sub_bkg_GplusPP_KET");
    
    double* v2ks_OB = v2obs_ks->GetY();
    double* v2eks_OB = v2obs_ks->GetEY();
    double* v2ks_SB = v2bkg_ks->GetY();
    double* v2eks_SB = v2bkg_ks->GetEY();
    
    double* ptks = v2obs_ks->GetX();

    double* KETks = v2obs_ks_KET->GetX();

    //bkg subtraction
    double v2t[13];
    double v2te[13];
    double v2b[13];
    double v2be[13];
    double v2t_ncq[13];
    double v2te_ncq[13];
    double KET_ncq[13];
    
    for(int i=0;i<8;i++)
    {
        v2b[i] = (v2ks_SB[i]*sigfrks[i] - v2ks_OB[i]*sigfrks_SB[i])/(sigfrks[i] - sigfrks_SB[i]);
        v2be[i] = sqrt((v2eks_SB[i]*sigfrks[i])**2 + (v2eks_OB[i]*sigfrks_SB[i])**2)/(sigfrks[i] - sigfrks_SB[i]);
        
        v2t[i] = (v2ks_OB[i] - v2b[i]*(1-sigfrks[i]))/sigfrks[i];
        v2te[i] = sqrt((v2be[i]*(1-sigfrks[i]))**2 + v2eks_OB[i]**2)/sigfrks[i];
        
        v2t_ncq[i] = v2t[i]/2.0;
        v2te_ncq[i] = v2te[i]/2.0;

        KET_ncq[i] = KETks[i]/2.0;
    }
    
    TGraphErrors *ksv2tg = new TGraphErrors(8,ptks,v2t,0,v2te);
    TGraphErrors *ksv2bg = new TGraphErrors(8,ptks,v2b,0,v2be);
    TGraphErrors *ksv2tg_KET = new TGraphErrors(8,KET_ncq,v2t_ncq,0,v2te_ncq);
    
    ksv2tg_KET->SetName("D0v2true_KET_ncq");
  
    ksv2bg->SetName("D0v2bkg");

    ksv2tg->SetName("D0v2true");

    TFile ofile("D0v2_vspt_nosub_NassFit_highpt_pt033_true_swap.root","RECREATE");
    
    ksv2tg.Write();
    ksv2bg.Write();
    
    ksv2tg_KET.Write();
}
