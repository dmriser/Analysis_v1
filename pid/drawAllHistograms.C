{

  TFile * file = TFile::Open("MCKeppelRad.root");
  TCanvas * c1 = new TCanvas("c1", "", 800, 800);


  // Electron Identification Histograms [2][7]                                                                                                                                                                      // [0] - without cuts,  [1] - with cuts                                                                                                                                                                           // [0] - all sectors,   [1-6] - sectors 1-6                                                                                                                                                                  
  string type[11] = {"allNegatives", "cuts", "Z_VERTEX","CC_FID","CC_PHI","CC_THETA","DC_R1_FID","DC_R3_FID","EC_FID","EC_IN_OUT","EC_SAMPLING"};
  string sect[7] = {"all", "s1", "s2", "s3", "s4", "s5", "s6"};

  TLatex lab; 
  lab.SetNDC();
  //  lab.SetTextFont(42);
  lab.SetTextSize(0.04);

  // 1-D                                                                                                                                                                                                       
  TH1D * h1_nphe[11][7];
  TH1F * h1_ec_edep_inner[11][7];
  TH1F * h1_ec_edep_outer[11][7];
  TH1F * h1_p[11][7];
  TH1F * h1_z_vertex[11][7];

  // 2-D                                                                                                                                                                                                          
  TH2F * h2_cc_theta[11][7];
  TH2F * h2_etot_p[11][7];
  TH2F * h2_ang_fid[11][7];
  TH2F * h2_ec_edep[11][7];
  TH2F * h2_dcr1_fid[11][7];
  TH2F * h2_dcr3_fid[11][7];
  TH2F * h2_ec_fid[11][7];

  // initialize                                                                                                                                                                           
  string printer = "";
                        
  for (int itype = 0; itype < 11; itype++)
    for(int isect = 0; isect < 7; isect++){
	// 1d                                                                                                                                                                                                   
	h1_nphe[itype][isect]          = (TH1D*) file.Get(Form("h1_nphe_%s_%s",type[itype].c_str(),sect[isect].c_str()));
	h1_nphe[itype][isect]->SetFillColorAlpha(85, 0.4);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h1_nphe[itype][isect]->Draw();
	printer = Form(" #rightarrow pass: %f ", (double) h1_nphe[itype][isect]->GetEntries()/h1_nphe[0][isect]->GetEntries() );
	lab.DrawLatex(0.5, 0.9, printer.c_str());
	c1->Print(Form("img/h1_nphe_%s_%s.png",type[itype].c_str(),sect[isect].c_str()));

	h1_ec_edep_inner[itype][isect] = (TH1F*) file.Get(Form("h1_ec_edep_inner_%s_%s",type[itype].c_str(),sect[isect].c_str()));
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h1_ec_edep_inner[itype][isect]->Draw();
	printer = Form(" #rightarrow pass: %f ", (double) h1_ec_edep_inner[itype][isect]->GetEntries()/h1_nphe[0][isect]->GetEntries() );
	lab.DrawLatex(0.5, 0.9, printer.c_str());
	c1->Print(Form("img/h1_ec_edep_inner_%s_%s.png",type[itype].c_str(),sect[isect].c_str())); 

	h1_ec_edep_outer[itype][isect] = (TH1F*) file.Get(Form("h1_ec_edep_outer_%s_%s",type[itype].c_str(),sect[isect].c_str()));
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h1_ec_edep_outer[itype][isect]->Draw();
	printer = Form(" #rightarrow pass: %f ", (double) h1_ec_edep_outer[itype][isect]->GetEntries()/h1_nphe[0][isect]->GetEntries() );
	lab.DrawLatex(0.5, 0.9, printer.c_str());
	c1->Print(Form("img/h1_ec_edep_outer_%s_%s.png",type[itype].c_str(),sect[isect].c_str()));

	h1_p[itype][isect]             = (TH1F*) file.Get(Form("h1_p_%s_%s",type[itype].c_str(),sect[isect].c_str()));
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h1_p[itype][isect]->SetFillColorAlpha(62, 0.4);
	h1_p[itype][isect]->Draw();
	printer = Form(" #rightarrow pass: %f ", (double) h1_p[itype][isect]->GetEntries()/h1_nphe[0][isect]->GetEntries() );
	lab.DrawLatex(0.5, 0.9, printer.c_str());
	c1->Print(Form("img/h1_p_%s_%s.png",type[itype].c_str(),sect[isect].c_str()));

	h1_z_vertex[itype][isect]      = (TH1F*) file.Get(Form("h1_z_vertex_%s_%s",type[itype].c_str(),sect[isect].c_str()));
	gPad->SetGridx(1);
	gPad->SetGridy(1);
	h1_z_vertex[itype][isect]->SetFillColorAlpha(55, 0.4);
	h1_z_vertex[itype][isect]->Draw();
	printer = Form(" #rightarrow pass: %f ", (double) h1_z_vertex[itype][isect]->GetEntries()/h1_nphe[0][isect]->GetEntries() );
	lab.DrawLatex(0.5, 0.9, printer.c_str());
	c1->Print(Form("img/h1_z_vertex_%s_%s.png",type[itype].c_str(),sect[isect].c_str()));
	
	// 2d                                                                                                                                                                                                   
	h2_cc_theta[itype][isect] = (TH2F*) file.Get(Form("h_cc_theta_%s_%s",type[itype].c_str(),sect[isect].c_str()));
	h2_cc_theta[itype][isect]->Draw("colz");
	printer = Form(" #rightarrow pass: %f ", (double) h2_cc_theta[itype][isect]->GetEntries()/h1_nphe[0][isect]->GetEntries() );
	lab.DrawLatex(0.5, 0.9, printer.c_str());
	c1->Print(Form("img/h_cc_theta_%s_%s.png",type[itype].c_str(),sect[isect].c_str()));

	h2_etot_p[itype][isect]   = (TH2F*) file.Get(Form("h_etot_p_%s_%s",type[itype].c_str(),sect[isect].c_str()));
	h2_etot_p[itype][isect]->Draw("colz");
	printer = Form(" #rightarrow pass: %f ", (double) h2_etot_p[itype][isect]->GetEntries()/h1_nphe[0][isect]->GetEntries() );
	lab.DrawLatex(0.5, 0.9, printer.c_str());
	c1->Print(Form("img/h_etot_p_%s_%s.png",type[itype].c_str(),sect[isect].c_str()));

	h2_ang_fid[itype][isect]  = (TH2F*) file.Get(Form("h_ang_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()));
	h2_ang_fid[itype][isect]->Draw("colz");
	printer = Form(" #rightarrow pass: %f ", (double) h2_ang_fid[itype][isect]->GetEntries()/h1_nphe[0][isect]->GetEntries() );
	lab.DrawLatex(0.5, 0.9, printer.c_str());
	c1->Print(Form("img/h_ang_fid_%s_%s.png",type[itype].c_str(),sect[isect].c_str()));

	h2_ec_edep[itype][isect]  = (TH2F*) file.Get(Form("h_ec_edep_%s_%s",type[itype].c_str(),sect[isect].c_str()));
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogz(); 
	h2_ec_edep[itype][isect]->Draw("colz"); 
	printer = Form(" #rightarrow pass: %f ", (double) h2_ec_edep[itype][isect]->GetEntries()/h1_nphe[0][isect]->GetEntries() );
	lab.DrawLatex(0.5, 0.9, printer.c_str());
	c1->Print(Form("img/h_ec_edep_%s_%s.png",type[itype].c_str(),sect[isect].c_str()));

	h2_dcr1_fid[itype][isect] = (TH2F*) file.Get(Form("h_dcr1_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()));
	h2_dcr1_fid[itype][isect]->Draw("colz");
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogz();
	printer = Form(" #rightarrow pass: %f ", (double) h2_dcr1_fid[itype][isect]->GetEntries()/h1_nphe[0][isect]->GetEntries() );
	lab.DrawLatex(0.5, 0.9, printer.c_str());
	c1->Print(Form("img/h_dcr1_fid_%s_%s.png",type[itype].c_str(),sect[isect].c_str()));

	h2_dcr3_fid[itype][isect] = (TH2F*) file.Get(Form("h_dcr3_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()));
	h2_dcr3_fid[itype][isect]->Draw("colz");
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogz();
	printer = Form(" #rightarrow pass: %f ", (double) h2_dcr3_fid[itype][isect]->GetEntries()/h1_nphe[0][isect]->GetEntries() );
	lab.DrawLatex(0.5, 0.9, printer.c_str());
	c1->Print(Form("img/h_dcr3_fid_%s_%s.png",type[itype].c_str(),sect[isect].c_str()));

	h2_ec_fid[itype][isect]   = (TH2F*) file.Get(Form("h_ec_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()));
	h2_ec_fid[itype][isect]->Draw("colz");
	gPad->SetGridx();
	gPad->SetGridy();
	gPad->SetLogz();
	printer = Form(" #rightarrow pass: %f ", (double) h2_ec_fid[itype][isect]->GetEntries()/h1_nphe[0][isect]->GetEntries() );
	lab.DrawLatex(0.5, 0.9, printer.c_str());
	c1->Print(Form("img/h_ec_fid_%s_%s.png",type[itype].c_str(),sect[isect].c_str()));
      }

  // ------------ end histograms ---------------  

}
