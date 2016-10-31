
{  

#include "../src/histos.h"
#include "../src/histos.cc"
  //#include "../analysisClasses/Bins.h"


  gROOT->LoadMacro("utils.C");
  gROOT->LoadMacro("show.C");
  gStyle->SetOptTitle(1);

  //    string filename = "../out/full_backup.root";
  string filename = "../out/elastic.root";
  
  elast_histos histos;
  histos.load(filename);

  // control variables 
  int TYPE  = 0;
  int SECT  = 0;
  bool VIEW = false;
  string PRINT = "png";

  string S[7] = {"ALL","S1","S2","S3","S4","S5","S6"};
  string T[3] = {"GSIM_NORAD","GSIM_RAD","DATA"};

  // retrieve binning scheme 
  BinStructure relPhiBins(histos.h2_rec[0][0]->GetNbinsX(),
			  histos.h2_rec[0][0]->GetXaxis()->GetXmin(),
			  histos.h2_rec[0][0]->GetXaxis()->GetXmax());

  BinStructure thetaBins(histos.h1_rec_theta[0][0]->GetNbinsX(),
			 histos.h1_rec_theta[0][0]->GetXaxis()->GetXmin(),
			 histos.h1_rec_theta[0][0]->GetXaxis()->GetXmax());


  const int NPHIBINS = relPhiBins.number();
  
  do_xsection();
  do_norm_theta();

  cout << endl;

  // Display to look at dist.
  bar = new TControlBar("vertical", " David Riser - Elastic ");
  bar->AddButton("","");
  bar->AddButton("Change View (Combined/Single)","change_view()");
  bar->AddButton("Change Type ",                "change_type()");
  bar->AddButton("Change Sector ",              "change_sector()");
  bar->AddButton("","");
  bar->AddButton(" Show Cross Section  ","show_xs()");
  bar->AddButton(" Show XS/Model Ratio ","show_xs_ratio()");
  bar->AddButton(" Show Acceptance  ","show_acc()");
  bar->AddButton(" Show Radiative Correction  ","show_rc()");
  bar->AddButton(" Show Reconstructed  ","show_rec()");
  bar->AddButton(" Show Generated ","show_gen()");
  bar->AddButton("","");
  bar->AddButton(" Show Reconstructed Theta ","show_rec_theta()");
  bar->AddButton(" Show Generated Theta ","show_gen_theta()");
  bar->AddButton(" Show Normalized Theta ","show_norm_theta()");
  bar->AddButton("","");
  bar->AddButton(" Show Bjorken x ","show_x()");
  bar->AddButton(" Show Mass Squared ","show_MM2()");
  bar->AddButton(" Show W ","show_W()");
  bar->AddButton(" Show #Q^2 ","show_QQ()");
  bar->AddButton(" Show EP Angle ","show_EP()");
  bar->AddButton("","");
  bar->AddButton(" Show Faraday Cup Monitor ","show_fc()");
  bar->AddButton(" Print Error Report ","print_error_report()");
  bar->AddButton(" Print All to .pdf ","print_all()");
  bar->AddButton("","");
  bar->AddButton(" Exit ",".q");
  bar->Show();

  
}
