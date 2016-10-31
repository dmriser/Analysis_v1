{
  
  //  gStyle->SetPalette(53);
  gStyle->SetOptStat(1000000001);
  gStyle->SetOptFit(0001);

  // getting file or dying
  TFile * _file = TFile::Open("gpp.root");
  if (_file) cout << " file opened " << endl;
  else {cout << " file error" << endl; exit(0);}
  
  TH1F * h_W[2][7];
  TH1F * h_dt[2][7];
  TH2F * h_dt_p[2][7];

  // setting up to get histos 
  string type[2] = {"data","GSIM"};
  string histos[7] = {"all","s1","s2","s3","s4","s5","s6"};
  Color_t colors[2] = {kCyan, kGreen -4};
  string fitParameters[6] = {"A-Gauss","Mu-Gauss","Sigma-Gauss","a-pol2","b-pol2","c-pol2"};

  // getting histos and setting titles 
  for (int itype=0; itype<2; itype++)
    for (int ihist = 0; ihist < 7; ihist++)
      {
	h_W[itype][ihist] = (TH1F*) _file.Get(Form("h_W_%s_%s",type[itype].c_str(),histos[ihist].c_str()));
	h_W[itype][ihist]->SetFillColor(colors[itype]);
	h_W[itype][ihist]->SetTitle(Form("W for %s of %s",histos[ihist].c_str(),type[itype].c_str()));
	h_W[itype][ihist]->GetXaxis()->SetTitle("W [GeV/c]");

	h_dt[itype][ihist] = (TH1F*) _file.Get(Form("h_dt_%s_%s",type[itype].c_str(),histos[ihist].c_str()));
	h_dt[itype][ihist]->SetFillColor(colors[itype]);
	h_dt[itype][ihist]->SetTitle(Form("dt for %s of %s",histos[ihist].c_str(),type[itype].c_str()));
	h_dt[itype][ihist]->GetXaxis()->SetTitle("dt [ns]");

	h_dt_p[itype][ihist] = (TH2F*) _file.Get(Form("h_dt_p_%s_%s",type[itype].c_str(),histos[ihist].c_str()));
	h_dt_p[itype][ihist]->SetTitle(Form("dt vs p for %s of %s",histos[ihist].c_str(),type[itype].c_str()));
	h_dt_p[itype][ihist]->GetXaxis()->SetTitle("p [GeV/c]");
	h_dt_p[itype][ihist]->GetYaxis()->SetTitle("dt [ns]");

      }

  // ----------------- defining fits ------------------------                                                                                 
  TF1 * f_W[2][7];
  TF1 * f_dt[2][7];

  // setting up limits of W fit by sector
  //                             all    s1   s2     s3    s4    s5     s6
  double limit_W_left[2][7]  = {{0.84, 0.84, 0.83, 0.84, 0.84, 0.84, 0.86},  // data
				{0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80}}; // GSIM
  double limit_W_right[2][7] = {{1.01, 1.01, 1.02, 1.01, 1.02, 1.01, 1.02},
				{1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04}};
  /*  double limit_W_left[2][7]  = {{0.84, 0.84, 0.83, 0.84, 0.84, 0.84, 0.86},  // data
				{0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80}}; // GSIM
  double limit_W_right[2][7] = {{1.02, 1.02, 1.03, 1.02, 1.02, 1.02, 1.02},
				{1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04}};
  */
  for (int itype = 0; itype < 2; itype++)
    for (int ihist = 0; ihist < 7; ihist++)
      {
        f_W[itype][ihist] = new TF1(Form("f_W_%s_%s",type[itype].c_str(),histos[ihist].c_str()),"[0]*exp(-0.5*((x-[1])/[2])^2)",limit_W_left[itype][ihist],limit_W_right[itype][ihist]);
	f_dt[itype][ihist] = new TF1(Form("f_dt_%s_%s",type[itype].c_str(),histos[ihist].c_str()),"[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*x^2 + [4]*x + [5]",-1.5,1.5);

	f_W[itype][ihist]->SetParameter(0,h_W[itype][ihist]->GetMaximum());
	f_W[itype][ihist]->SetParameter(1,h_W[itype][ihist]->GetMean());
	f_W[itype][ihist]->SetParameter(2,h_W[itype][ihist]->GetRMS());
      
	f_dt[itype][ihist]->SetParameter(0,h_dt[itype][ihist]->GetMaximum());
	f_dt[itype][ihist]->SetParameter(1,h_dt[itype][ihist]->GetMean());
	f_dt[itype][ihist]->SetParameter(2,h_dt[itype][ihist]->GetRMS());
	f_dt[itype][ihist]->SetParameter(3,-0.01);
	f_dt[itype][ihist]->SetParameter(4,0.01);
	f_dt[itype][ihist]->SetParameter(5,0.01);
	
	f_W[itype][ihist] ->SetParNames("A","mu","sigma");
	f_dt[itype][ihist]->SetParNames("A","mu","sigma","a","b","c");

      }

  for(int itype = 0; itype < 2; itype++)
    for (int ihist = 0; ihist < 7; ihist++)
      {
	h_W[itype][ihist]->Fit(f_W[itype][ihist],"rq");
	h_dt[itype][ihist]->Fit(f_dt[itype][ihist],"rq");
      }

  // ------------------ end fits --------------------------                                                                                 


  // ----------- print fit params -----------------------
  cout << "\n > Values for W Fits " << endl;

  cout.width(12);
  cout << "type";

  cout.width(12);
  cout << "parameter";

  for(int ihist = 0; ihist < 7; ihist++)
    { 
      cout.width(12);
      cout << histos[ihist];
    }  

  cout << endl;

  for (int itype = 0; itype < 2; itype++)
    {
      for (int iparam = 0; iparam < 3; iparam++)
	{     
	  cout.width(12);
	  cout << type[itype];
	  cout.width(12);
	  cout << fitParameters[iparam];
	  for(int ihist = 0; ihist < 7; ihist++)
	    {
	      cout.width(12);
	      cout << f_W[itype][ihist]->GetParameter(iparam);
	    }
	  cout << endl;
	}
    }

  cout << "\n > Values for dt Fits " << endl;

  cout.width(12);
  cout << "type";

  cout.width(12);
  cout << "parameter";

  for(int ihist = 0; ihist < 7; ihist++)
    { 
      cout.width(12);
      cout << histos[ihist];
    }  

  cout << endl;

  for (int itype = 0; itype < 2; itype++)
    {
      for (int iparam = 0; iparam < 6; iparam++)
	{     
	  cout.width(12);
	  cout << type[itype];
	  cout.width(12);
	  cout << fitParameters[iparam];

	  for(int ihist = 0; ihist < 7; ihist++)
	    {
	      cout.width(12);
	      cout << f_dt[itype][ihist]->GetParameter(iparam);
	    }
	  cout << endl;
	}
    }


  // --------- end print fit params ---------------------

  // --------- start drawing ---------------------------

  TCanvas * c1 = new TCanvas("c1","",800,800);
  c1->Divide(2,2);

  c1->cd(1);
  h_W[0][0]->Draw();

  c1->cd(2);
  h_W[1][0]->Draw();

  c1->cd(3);
  h_dt[0][0]->Draw();

  c1->cd(4);
  h_dt[1][0]->Draw();

  c1->Print("gpp_summary.png");

  TCanvas * c2 = new TCanvas("c2","",1200,800);
  c2->Divide(3,2);
  
  for(int i=1; i<7; i++)
    {
      c2->cd(i);
      h_W[0][i]->Draw();
    }

  c2->Print("W_data.png");

  TCanvas * c3 = new TCanvas("c3","",1200,800);
  c3->Divide(3,2);
  
  for(int i=1; i<7; i++)
    {
      c3->cd(i);
      h_W[1][i]->Draw();
    }

  c3->Print("W_GSIM.png");

  TCanvas * c4 = new TCanvas("c4","",1200,800);
  c4->Divide(3,2);
  
  for(int i=1; i<7; i++)
    {
      c4->cd(i);
      h_dt[0][i]->Draw();
    }

  c4->Print("dt_data.png");

  TCanvas * c5 = new TCanvas("c5","",1200,800);
  c5->Divide(3,2);
  
  for(int i=1; i<7; i++)
    {
      c5->cd(i);
      h_dt[1][i]->Draw();
    }

  c5->Print("dt_GSIM.png");

  TCanvas * c6 = new TCanvas("c6","",800,400);
  c6->Divide(2,1);
  c6->cd(1);
  h_dt_p[0][0]->Draw("colz");
  c6->cd(2);
  h_dt_p[1][0]->Draw("colz");

  // ----------------- end drawing ------------------
}
