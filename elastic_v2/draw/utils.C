 void change_sector()
{
  SECT++;
  if (SECT == 7) SECT = 0;
  cout << " Sector changed to " << SECT << endl;
}

void change_type()
{

  string sType[2] = {"GSIM","Data."};

  TYPE++;
  if (TYPE == 2) TYPE = 0;
  cout << " Type changed to " << sType[TYPE] << endl;
}

void change_view()
{
  VIEW = !VIEW;
  string v[2] = {"one at a time.","all."};
  cout << " View changed to " << v[VIEW] << endl; 
}


void do_xsection()
{

  /* 
     xs = xs_measure * (N_GEN_RAD/N_REC_RAD) * (RAD CORR)
     
     RAD CORR = (NORAD GEN)/(RAD GEN)

   */ 

  cout << "\n Doing cross section... ";
  cout << " Acceptance calculated with radiative effects... ";

  for (int t=0; t<2; t++)
    for (int s=0; s<7; s++)
      {
	histos.h2_acc[t][s] = (TH2F*) histos.h2_rec[t][s]->Clone();
	histos.h2_acc[t][s]->Divide(histos.h2_gen[t][s]);       
      }

  /*
  //! Setting Model Errors To 0 
  for (int b=0; b<thetaBins.number(); b++)
    {
      histos.h1_xs_model[0]->SetBinError(b+1,0.0);
      histos.h1_xs_model[1]->SetBinError(b+1,0.0);
    }
  */

  for (int s=0; s<7; s++)
    {

      histos.h1_rc = (TH1F*) histos.h1_xs_model[1]->Clone();
      histos.h1_rc->Divide( histos.h1_xs_model[0] );
      //      if (s == 0) for (int t=0; t<thetaBins.number(); t++){ histos.h1_rc->SetBinError(t+1, 0.00); }

      // clone data hits 
      histos.h1_xs[s]     = new TH1D(Form("h1_xs_S%d",s),Form("Cross Section for S%d",s),thetaBins.number(),thetaBins.min(),thetaBins.max());
      //      histos.h1_xs_raw[s] = new TH1D(Form("h1_xs_raw_S%d",s),Form("Raw Cross Section for S%d",s),thetaBins.number(),thetaBins.min(),thetaBins.max());
      histos.h2_xs[s]     = (TH2F*) histos.h2_rec[0][s]->Clone();
      histos.h2_xs[s]             ->SetName(Form("h2_xs_S%d",s));
      histos.h2_xs[s]             ->Divide(histos.h2_scale);
      histos.h2_xs[s]             ->Scale(1/(thetaBins.width()*relPhiBins.width()));         // per unit angle in degrees 
      if (s > 0) histos.h2_xs[s]  ->Scale(3.14159/3);                                        // Integrated over phi means 60 deg per sector. 
      if (s == 0) histos.h2_xs[s] ->Scale(1/6.284);                                          // Divide out factor of 6 sectors going into 1 bin.
      histos.h2_xs[s]             ->Divide(histos.h2_acc[1][s]);                             // Doing cross section with radiative 
      histos.h2_xs[s]             ->ProjectionY(Form("h1_xs_raw_S%d",s),1,relPhiBins.number());  // somehow throwing off the error bars 
      histos.h2_xs[s]             ->ProjectionY(Form("h1_xs_S%d",s),1,relPhiBins.number());  // somehow throwing off the error bars 
      histos.h1_xs[s]             ->Scale((double) 1/relPhiBins.number());                   // because we summed over the phi bins
      //      histos.h1_xs_raw[s]         ->Scale((double) 1/relPhiBins.number());                   // because we summed over the phi bins
      //      histos.h1_xs[s]             ->Divide(histos.h1_rc);
      
      /*
      // Manually setting errors for xs
      for (int t=0; t<thetaBins.number(); t++)
	{
	  double err2 = 0.0;
	  for (int p=0; p<relPhiBins.number(); p++)
	    {
	      err2 += (histos.h2_xs[s]->GetBinError(p+1,t+1))*(histos.h2_xs[s]->GetBinError(p+1,t+1));	      
	    }
	  histos.h1_xs[s]->SetBinError(t+1,sqrt(err2));
	}
      */

      // Ratio plots
      histos.h1_xs_ratio[0][s] = (TH1F*) histos.h1_xs[s]     ->Clone();
      histos.h1_xs_ratio[1][s] = (TH1F*) histos.h1_xs[s] ->Clone();
      histos.h1_xs_ratio[0][s]->SetName(Form(" Cross Section Ratio for S%d (NO RAD)",s));
      histos.h1_xs_ratio[1][s]->SetName(Form(" Cross Section Ratio for S%d ( RAD )",s));
      histos.h1_xs_ratio[0][s]->Divide(histos.h1_xs_model[0]);
      histos.h1_xs_ratio[1][s]->Divide(histos.h1_xs_model[1]);

      double RAT_MIN = 0.1;
      double RAT_MAX = 5.0;
      /*
      // Manually doing errors for xs ratio 
      for (int t=0; t<thetaBins.number(); t++)
	{
	  double r0_err = histos.h1_xs[s]->GetBinError(t+1)/histos.h1_xs_model[0]->GetBinContent(t+1);
	  double r1_err = histos.h1_xs[s]->GetBinError(t+1)/histos.h1_xs_model[1]->GetBinContent(t+1);
	  histos.h1_xs_ratio[0][s]->SetBinError(t+1,r0_err);
	  histos.h1_xs_ratio[1][s]->SetBinError(t+1,r1_err);

	  if (histos.h1_xs_ratio[0][s]->GetBinContent(t+1)>RAT_MAX || 
	      histos.h1_xs_ratio[0][s]->GetBinContent(t+1)<RAT_MIN )
	    {
	      histos.h1_xs_ratio[0][s]->SetBinContent(t+1, 0.0);
	    }
	
	  if (histos.h1_xs_ratio[1][s]->GetBinContent(t+1)>RAT_MAX || 
	      histos.h1_xs_ratio[1][s]->GetBinContent(t+1)<RAT_MIN )
	    {
	      histos.h1_xs_ratio[1][s]->SetBinContent(t+1, 0.0);
	    }	
	}
      */
    }
  
  cout << " Done! " << endl; 
}

void print_all()
{
  bool pVIEW = VIEW;

  PRINT = "png";  // change if you want other extension. 
  VIEW  = true;    // printing supported for combined view only.
  SECT  = 0;  
  TYPE  = 0;

  for (int t=0; t<3; t++)
    {
      show_xs();
      show_acc();
      show_rec();
      show_gen();
      show_rec_theta();
      show_gen_theta();
      show_xs_ratio();
      show_rc();
      show_fc();
      show_W_comparison();
      show_xs_raw();

      // goto next type
      change_type();
  }

  // set values back 
  PRINT = "";
  VIEW = pVIEW;

}

void print_error_report()
{
    cout << " Error summary " << endl;
    cout.width(12); cout << " bin ";
    cout.width(12); cout << " h2_scale ";
    cout.width(12); cout << " h2_acc ";
    cout.width(12); cout << " h2_rec ";
    cout.width(12); cout << " h2_gen ";
    cout.width(12); cout << " h2_xs ";
    cout.width(12); cout << " data " << endl;

  for (int t=0; t<thetaBins.number(); t++)
    for (int p=0; p<relPhiBins.number(); p++)
      {
      cout.width(12); cout << "(" << p << "," << t << ")";
      cout.width(12); cout << histos.h2_scale     ->GetBinError(p,t);
      cout.width(12); cout << histos.h2_acc[1][1] ->GetBinError(p,t);
      cout.width(12); cout << histos.h2_rec[1][1] ->GetBinError(p,t);
      cout.width(12); cout << histos.h2_gen[1][1] ->GetBinError(p,t);
      cout.width(12); cout << histos.h2_xs[1]     ->GetBinError(p,t);
      cout.width(12); cout << histos.h2_gen[2][1] ->GetBinError(p,t) << endl;
      }
}

void do_slices()
{
  TH1D * h1_xs_slice[7][NPHIBINS];
  
  for(int sl=0; sl<NPHIBINS; sl++)
    for (int isect=0; isect<7; isect++)
    {
      h1_xs_slice[isect][sl] = new TH1D(Form("h1_xs_slice_S%d_%d",isect,sl),Form("h1_xs_slice_S%d_%d",isect,sl),thetaBins.number(),thetaBins.min(),thetaBins.max());
      histos.h2_xs[isect]->ProjectionY(Form("h1_xs_slice_S%d_%d",isect,sl),sl,sl+1);
    }

  TCanvas * c1 = new TCanvas("c1","",800,800);
  int n = ceil(sqrt((double)NPHIBINS));
  cout << " canvas width " << n << endl;
  c1->Divide(n,n);

  // settings errors from mother histo 
  for (int is=0; is<7; is++)
    { 
      for(int isl=0; isl<NPHIBINS; isl++)
	{
	  
	  //	  for (int tbin=0; tbin<thetaBins.number(); tbin++) h1_xs_slice[is][isl]->SetBinError(tbin,sqrt(histos.h2_xs[is]->GetBinError(isl+1,tbin+1)*histos.h2_xs[is]->GetBinError(isl+1,tbin+1)));
	  c1->cd(isl+1);
	  h1_xs_slice[is][isl]->SetMarkerStyle(4);
	  h1_xs_slice[is][isl]->Divide(histos.h1_rc);
	  h1_xs_slice[is][isl]->Divide(histos.h1_xs_model[0]);
	  //	  h1_xs_slice[is][isl]->Draw("P");
	  h1_xs_slice[is][isl]->Draw("P");
	  //	histos.h1_xs_model[0]->SetMarkerStyle(3);
	  histos.h1_xs_model[0]->SetLineColor(kBlue);
	  //	  histos.h1_xs_model[0]->Draw("Lsame");
	}
      c1->Print(Form("../img/xs_slices/h1_xs_slice_S%d.pdf",is));
    }
}

void do_acc_slices()
{
  TH1D * h1_acc_slice[7][NPHIBINS];

  TLatex lab;
  lab.SetNDC();
  lab.SetTextFont(32);
  lab.SetTextSize(0.50);
  
  for(int sl=0; sl<NPHIBINS; sl++)
    for (int isect=0; isect<7; isect++)
    {
      h1_acc_slice[isect][sl] = new TH1D(Form("h1_acc_slice_S%d_%d",isect,sl),Form("h1_acc_slice_S%d_%d",isect,sl),thetaBins.number(),thetaBins.min(),thetaBins.max());
      histos.h2_acc[1][isect]->ProjectionY(Form("h1_acc_slice_S%d_%d",isect,sl),sl,sl+1);
    }

  TCanvas * c1 = new TCanvas("c1","",800,800);
  int n = ceil(sqrt((double)NPHIBINS));
  cout << " canvas width " << n << endl;
  c1->Divide(n,n);

  // settings errors from mother histogram 
  for (int is=0; is<7; is++)
    { 
      for(int isl=0; isl<NPHIBINS; isl++)
	{
	  
	  //	  for (int tbin=0; tbin<thetaBins.number(); tbin++) h1_xs_slice[is][isl]->SetBinError(tbin,sqrt(histos.h2_xs[is]->GetBinError(isl+1,tbin+1)*histos.h2_xs[is]->GetBinError(isl+1,tbin+1)));
	  c1->cd(isl+1);
	  h1_acc_slice[is][isl]->SetMarkerStyle(2);
	  h1_acc_slice[is][isl]->Draw("P");
	  lab.DrawLatex(0.75, 0.85, Form(" #phi %f - %f deg. ",histos.h2_acc[1][is]->GetBin(isl),histos.h2_acc[1][is]->GetBin(isl+1)));
	}
      c1->Print(Form("../img/acc_slices/h1_acc_slice_S%d.pdf",is));
    }
}

void do_w_debug_plots()
{

  cout << " Doing W Debug plots... ";

  TCanvas * c1 = new TCanvas("c1","",1200,800);
  c1->Divide(3,2);
  
  for (int s=1; s<7; s++){
    c1->cd(s);
    histos.h1_W[0][s]->Scale( 1/histos.h1_W[0][s]->GetMaximum() );
    histos.h1_W[1][s]->Scale( 1/histos.h1_W[1][s]->GetMaximum() );
    histos.h1_W[2][s]->Scale( 1/histos.h1_W[2][s]->GetMaximum() );
    histos.h1_W[0][s]->SetLineColor(kGreen);
    histos.h1_W[1][s]->SetLineColor(kBlue);
    histos.h1_W[2][s]->SetLineColor(kRed);
    histos.h1_W[0][s]->Draw();
    histos.h1_W[1][s]->Draw("same");
    histos.h1_W[2][s]->Draw("same");    
  }
  c1->Print("../img/w_debug/h1_w_debug.pdf");
  cout << " Done! " << endl;  

}
