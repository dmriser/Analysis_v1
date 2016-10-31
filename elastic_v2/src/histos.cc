#ifndef elast_histos_cxx
#define elast_histos_cxx

// c++ includes
#include <iostream>
using namespace std;

// my includes 
#include "../../analysisClasses/h22Event.h"
#include "../../analysisClasses/Bins.h"
#include "histos.h"

// ROOT includes
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"

elast_histos::elast_histos()
{
  // birth 
  w_cut = 1.10; // default 
}

elast_histos::~elast_histos()
{
  // total destruction 
}

void elast_histos::set_bins(BinStructure *tBins, BinStructure *pBins, BinStructure * mBins)
{
  thetaBins = new BinStructure(tBins->number(), tBins->min(), tBins->max());
  phiBins   = new BinStructure(pBins->number(), pBins->min(), pBins->max());
  momBins   = new BinStructure(mBins->number(), mBins->min(), mBins->max());
}

void elast_histos::init()
{

  string TYPE[2]     = {"DATA","GSIM"};
  string SECT[7]     = {"ALL","S1","S2","S3","S4","S5","S6"};

  for (int t=0; t<2; t++)
    {
       for (int s=0; s<7; s++)
	{
	  // 1-D
	  h1_W[t][s]         = new TH1F(Form("h1_W_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),Form("h1_W_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),200,0.6,1.3);
	  h1_QQ[t][s]        = new TH1F(Form("h1_QQ_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),Form("h1_QQ_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),100,0.5,3.5);

	}      
    }

  h1_rc     = new TH1F("h1_rc"," Radiative Correction",thetaBins->number(),thetaBins->min(),thetaBins->max());
  h1_fc_mon = new TH1F("h1_fc_mon","h1_fc_mon",100,0,50);

  // with and without radiation 
  string MODELTYPE[2] = {"NORAD","RAD"};
  for (int t=0; t<2; t++) h1_xs_model[t] = new TH1F(Form("h1_xs_model_%s",MODELTYPE[t].c_str()),Form("h1_xs_model_%s",MODELTYPE[t].c_str()),thetaBins->number(),thetaBins->min(),thetaBins->max());

  // Binned by Outside Parameters 
  for (int s=0; s<7; s++){
    vector<TH1D*> v; 
    for (int p=0; p<=momBins->number(); p++)
      {
	string h1_name = Form("h1_w_by_mom_bin%d_%s",p,SECT[s].c_str());
	TH1D * h1_temp = new TH1D(h1_name.c_str(), h1_name.c_str(), 100, 0.75, 1.25);
	v.push_back(h1_temp);
      }
    h1_w_by_mom.push_back(v);
  }

  // Binned by Outside Parameters 
  for (int s=0; s<7; s++){
    vector<TH1D*> v; 
    for (int p=0; p<=phiBins->number(); p++)
      {
	string h1_name = Form("h1_acc_by_phi_bin%d_%s",p,SECT[s].c_str());
	TH1D * h1_temp = new TH1D(h1_name.c_str(), h1_name.c_str(), thetaBins->number(), thetaBins->min(), thetaBins->max());
	v.push_back(h1_temp);
      }
    h1_acc_by_phi.push_back(v);
  }

  // Binned by Outside Parameters 
  for (int s=0; s<7; s++){
    vector<TH1D*> v; 
    for (int p=0; p<=phiBins->number(); p++)
      {
	string h1_name = Form("h1_xs_by_phi_bin%d_%s",p,SECT[s].c_str());
	TH1D * h1_temp = new TH1D(h1_name.c_str(), h1_name.c_str(), thetaBins->number(), thetaBins->min(), thetaBins->max());
	v.push_back(h1_temp);
      }
    h1_xs_by_phi.push_back(v);
  }

  // Binned by Outside Parameters 
  for (int s=0; s<7; s++){
    vector<TH1D*> v; 
    for (int p=0; p<=phiBins->number(); p++)
      {
	string h1_name = Form("h1_xs_ratio_by_phi_bin%d_%s",p,SECT[s].c_str());
	TH1D * h1_temp = new TH1D(h1_name.c_str(), h1_name.c_str(), thetaBins->number(), thetaBins->min(), thetaBins->max());
	v.push_back(h1_temp);
      }
    h1_xs_ratio_by_phi.push_back(v);
  }
 

  // Binned by Outside Parameters 
  for (int s=0; s<7; s++){
    vector<TH1D*> v; 
    for (int p=0; p<=phiBins->number(); p++)
      {
	string h1_name = Form("h1_rxs_by_phi_bin%d_%s",p,SECT[s].c_str());
	TH1D * h1_temp = new TH1D(h1_name.c_str(), h1_name.c_str(), thetaBins->number(), thetaBins->min(), thetaBins->max());
	v.push_back(h1_temp);
      }
    h1_rxs_by_phi.push_back(v);
  }

  // Binned by Outside Parameters 
  for (int s=0; s<7; s++){
    vector<TH1D*> v; 
    for (int p=0; p<=phiBins->number(); p++)
      {
	string h1_name = Form("h1_rxs_ratio_by_phi_bin%d_%s",p,SECT[s].c_str());
	TH1D * h1_temp = new TH1D(h1_name.c_str(), h1_name.c_str(), thetaBins->number(), thetaBins->min(), thetaBins->max());
	v.push_back(h1_temp);
      }
    h1_rxs_ratio_by_phi.push_back(v);
  }
 

  // Binned by Outside Parameters 
  for (int s=0; s<7; s++){
    vector<TH1D*> v; 
    for (int p=0; p<=phiBins->number(); p++)
      {
	string h1_name = Form("h1_gen_by_phi_bin%d_%s",p,SECT[s].c_str());
	TH1D * h1_temp = new TH1D(h1_name.c_str(), h1_name.c_str(), thetaBins->number(), thetaBins->min(), thetaBins->max());
	v.push_back(h1_temp);
      }
    h1_gen_by_phi.push_back(v);
  }
 

  // Binned by Outside Parameters 
  for (int s=0; s<7; s++){
    vector<TH1D*> v; 
    for (int p=0; p<=phiBins->number(); p++)
      {
	string h1_name = Form("h1_rec_by_phi_bin%d_%s",p,SECT[s].c_str());
	TH1D * h1_temp = new TH1D(h1_name.c_str(), h1_name.c_str(), thetaBins->number(), thetaBins->min(), thetaBins->max());
	v.push_back(h1_temp);
      }
    h1_rec_by_phi.push_back(v);
  }
 
  // Binned by Outside Parameters 
  for (int s=0; s<7; s++){
    vector<TH1D*> v; 
    for (int p=0; p<=phiBins->number(); p++)
      {
	string h1_name = Form("h1_hits_by_phi_bin%d_%s",p,SECT[s].c_str());
	TH1D * h1_temp = new TH1D(h1_name.c_str(), h1_name.c_str(), thetaBins->number(), thetaBins->min(), thetaBins->max());
	v.push_back(h1_temp);
      }
    h1_hits_by_phi.push_back(v);
  }
 
}

void elast_histos::draw_w(string name)
{
  // Doing Plots for W by P. 
  
  TCanvas * w_canvas = new TCanvas("w_canvas","w_canvas",1200,800);
  w_canvas->Print( Form("../img/%s_w_means.pdf[",name.c_str()) );
  w_canvas->Clear();
  w_canvas->Divide(3,2);
  
  for (int s=1; s<7; s++)
    {
      w_canvas->cd(s);
      g_mp_from_w[s]->SetMarkerStyle(8);
      g_mp_from_w[s]->Draw();
    }

  w_canvas->Print( Form("../img/%s_w_means.pdf",name.c_str()) );
  w_canvas->Clear();
  w_canvas->Divide(3,2);
  
  for (int s=1; s<7; s++)
    {
      w_canvas->cd(s);
      g_dmp_from_w[s]->SetMarkerStyle(4);
      g_dmp_from_w[s]->Draw();
    }

  w_canvas->Print( Form("../img/%s_w_means.pdf",name.c_str()) );
  w_canvas->Print( Form("../img/%s_w_means.pdf]",name.c_str()) );
}

void elast_histos::fill(int e_index, int index, h22Event event, ElasticEvent elasticEvent)
{
  // Data
  if (index == 0 && elasticEvent.getW() < w_cut)
    {
      // Filling Histograms 
      int s = event.dc_sect[e_index];

      // Fill All
      h1_W[index][0]        ->Fill(elasticEvent.getW());
      h1_QQ[index][0]       ->Fill(elasticEvent.getQQ());

      h1_W[index][s]        ->Fill(elasticEvent.getW());
      h1_QQ[index][s]       ->Fill(elasticEvent.getQQ());

      int b = phiBins->getBin( event.rphi(e_index) )+1;
      if (b > 0) { h1_hits_by_phi[0][b]->Fill(event.theta(e_index)); h1_hits_by_phi[s][b]->Fill(event.theta(e_index)); h1_hits_by_phi[s][0]->Fill(event.theta(e_index)); }

      b = momBins->getBin( event.p[e_index] )+1;
      if (b > 0) { h1_w_by_mom[0][b]->Fill(elasticEvent.getW()); h1_w_by_mom[s][b]->Fill(elasticEvent.getW()); h1_w_by_mom[s][0]->Fill(elasticEvent.getW()); }
    }

  // GSIM
  if (index == 1)
    {
      // Filling Histograms 
      int s   = event.dc_sect[e_index];
      int mcs = event.mcsect(e_index);

      if (elasticEvent.getW() < w_cut) {
	// Fill All
	h1_W[index][0]        ->Fill(elasticEvent.getW());
	h1_QQ[index][0]       ->Fill(elasticEvent.getQQ());
	
	h1_W[index][s]        ->Fill(elasticEvent.getW());
	h1_QQ[index][s]       ->Fill(elasticEvent.getQQ());

	int b = phiBins->getBin( event.rphi(e_index))+1;
	if (b > 0) { h1_rec_by_phi[s][b]->Fill(event.theta(e_index)); h1_rec_by_phi[0][b]->Fill(event.theta(e_index)); h1_rec_by_phi[s][0]->Fill(event.theta(e_index)); }

      }

    }
  
}

void elast_histos::fill_gen(int e_index, int index, h22Event event, ElasticEvent elasticEvent)
{

  // GSIM
  if (index == 1)
    {
      // Filling Histograms 
      int s   = event.dc_sect[e_index];
      int mcs = event.mcsect(e_index);
      int b = phiBins->getBin( event.mcrphi(0) )+1;
      if (b > 0 && mcs > 0 && mcs < 7) { h1_gen_by_phi[mcs][b]->Fill(event.mctheta[0]); h1_gen_by_phi[0][b]->Fill(event.mctheta[0]); h1_gen_by_phi[mcs][0]->Fill(event.mctheta[0]); }
    }
  
}

void elast_histos::write_and_close(string _file)
{
  TFile f(_file.c_str(),"recreate");

  h1_xs_model[0]->Write();
  h1_xs_model[1]->Write();
  h1_fc_mon->Write();
  
  for (int t=0; t<2; t++)
    for (int s=0; s<7; s++)
      {
	// 1D
	h1_W[t][s]->Write();
	h1_QQ[t][s]->Write();
      }
  
  // Writing histograms binned by outside parameters 
  for (int s=0; s<7; s++)
    for (int b=0; b<h1_acc_by_phi[s].size(); b++) // should be same as phiBins->number()
      {
	h1_gen_by_phi[s][b]      ->Write();
	h1_rec_by_phi[s][b]      ->Write();
	h1_hits_by_phi[s][b]     ->Write();
	h1_acc_by_phi[s][b]      ->Write();
	h1_xs_by_phi[s][b]       ->Write();
	h1_xs_ratio_by_phi[s][b] ->Write();
	h1_rxs_by_phi[s][b]      ->Write();
	h1_rxs_ratio_by_phi[s][b]->Write();
      }

  for (int s=0; s<7; s++)
    for (int b=0; b<h1_w_by_mom[s].size(); b++)
      h1_w_by_mom[s][b]->Write();	

  f.Write();
  f.Close();
}

void elast_histos::draw(string name, bool gif_status)
{

  // Standard 
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);

  // Latex Writer 
  TLatex lab;
  lab.SetNDC();
  lab.SetTextFont(22);

  // Formatting
  h1_xs_model[0]->SetLineColor(kBlue);
  h1_xs_model[1]->SetLineColor(kRed);

  // Setting up line for drawing straight over 1 in ratio plots 
  TLine ratio_line(thetaBins->min(),1,thetaBins->max(),1); 

  // Drawing Histograms 
  TCanvas * c1 = new TCanvas("c1","",1200,800);

  c1->Print(Form("../img/%s.pdf[",name.c_str()));
  c1->Clear();
  
  int can_size = floor(sqrt(h1_acc_by_phi[0].size()));

  for (int s=0; s<7; s++)
    {
      c1->Divide(can_size+1, can_size);
      for (int b=0; b<h1_acc_by_phi[s].size()-1; b++)
	{
	  c1->cd(b+1);

	  double phi_start = phiBins->min() + (b-1)*phiBins->width();
	  double phi_end   = phiBins->min() + b*phiBins->width();
	  if (b == 0) { phi_start = phiBins->min(); phi_end = phiBins->max(); }

	  h1_acc_by_phi[s][b]->SetMarkerStyle(8);
	  h1_acc_by_phi[s][b]->Draw("PE");
	  lab.DrawLatex(0.3, 0.02, Form(" Acceptance: Sector %d Bin %d ",s,b));
	  lab.DrawLatex(0.55, 0.8, Form(" %.2f #rightarrow %.2f #circ ",phi_start, phi_end));
	}

      c1->Print( Form("../img/%s.pdf",name.c_str()) );
      c1->Clear();
    }

  // Doing 1-D Cross Section 
  for (int s=0; s<7; s++)
    {
      c1->Divide(can_size+1, can_size);

      for (int b=0; b<h1_xs_by_phi[s].size()-1; b++)
	{
	  c1->cd(b+1);

	  double phi_start = phiBins->min() + (b-1)*phiBins->width();
	  double phi_end   = phiBins->min() + b*phiBins->width();
	  if (b == 0) { phi_start = phiBins->min(); phi_end = phiBins->max(); }

	  h1_xs_by_phi[s][b]->SetMarkerStyle(8);
	  h1_xs_by_phi[s][b]->Draw();
	  h1_xs_model[0]->Draw("Lsame");
	  lab.DrawLatex(0.3, 0.02, Form(" #sigma: Sector %d Bin %d ",s,b));
	  lab.DrawLatex(0.55, 0.8, Form(" %.2f #rightarrow %.2f #circ ",phi_start, phi_end));
	}

      c1->Print( Form("../img/%s.pdf",name.c_str()) );
      c1->Clear();
    }

  // Doing 1-D Cross Section 
  for (int s=0; s<7; s++)
    {
      c1->Divide(can_size+1, can_size);

      for (int b=0; b<h1_xs_ratio_by_phi[s].size()-1; b++)
	{
	  c1->cd(b+1);

	  double phi_start = phiBins->min() + (b-1)*phiBins->width();
	  double phi_end   = phiBins->min() + b*phiBins->width();
	  if (b == 0) { phi_start = phiBins->min(); phi_end = phiBins->max(); }

	  h1_xs_ratio_by_phi[s][b]->SetMarkerStyle(8);
	  h1_xs_ratio_by_phi[s][b]->SetMinimum(0.0);
	  h1_xs_ratio_by_phi[s][b]->SetMaximum(3.0);
	  h1_xs_ratio_by_phi[s][b]->Draw();
	  lab.DrawLatex(0.3, 0.02, Form(" #sigma/#sigma_born: Sector %d Bin %d ",s,b));
	  lab.DrawLatex(0.55, 0.8, Form(" %.2f #rightarrow %.2f #circ ",phi_start, phi_end));
	  ratio_line.Draw("same");
	}

      c1->Print( Form("../img/%s.pdf",name.c_str()) );
      c1->Clear();
    }


  // Doing 1-D Cross Section 
  for (int s=0; s<7; s++)
    {
      c1->Divide(can_size+1, can_size);

      for (int b=0; b<h1_rxs_by_phi[s].size()-1; b++)
	{
	  c1->cd(b+1);

	  double phi_start = phiBins->min() + (b-1)*phiBins->width();
	  double phi_end   = phiBins->min() + b*phiBins->width();
	  if (b == 0) { phi_start = phiBins->min(); phi_end = phiBins->max(); }

	  h1_rxs_by_phi[s][b]->SetMarkerStyle(8);
	  h1_rxs_by_phi[s][b]->Draw();
	  h1_xs_model[1]->Draw("Lsame");
	  lab.DrawLatex(0.3, 0.02, Form(" Raw #sigma: Sector %d Bin %d ",s,b));
	  lab.DrawLatex(0.55, 0.8, Form(" %.2f #rightarrow %.2f #circ ",phi_start, phi_end));
	}

      c1->Print( Form("../img/%s.pdf",name.c_str()) );
      c1->Clear();
    }


  // Doing 1-D Cross Section 
  for (int s=0; s<7; s++)
    {
      c1->Divide(can_size+1, can_size);

      for (int b=0; b<h1_rxs_ratio_by_phi[s].size()-1; b++)
	{
	  c1->cd(b+1);

	  double phi_start = phiBins->min() + (b-1)*phiBins->width();
	  double phi_end   = phiBins->min() + b*phiBins->width();
	  if (b == 0) { phi_start = phiBins->min(); phi_end = phiBins->max(); }

	  h1_rxs_ratio_by_phi[s][b]->SetMarkerStyle(8);
	  h1_rxs_ratio_by_phi[s][b]->SetMinimum(0.0);
	  h1_rxs_ratio_by_phi[s][b]->SetMaximum(3.0);
	  h1_rxs_ratio_by_phi[s][b]->Draw();
	  ratio_line.Draw("same");
	  lab.DrawLatex(0.3, 0.02, Form(" Raw #sigma/#sigma_rad: Sector %d Bin %d ",s,b));
	  lab.DrawLatex(0.55, 0.8, Form(" %.2f #rightarrow %.2f #circ ",phi_start, phi_end));
	}

      c1->Print( Form("../img/%s.pdf",name.c_str()) );
      c1->Clear();
    }

  // Doing 1-D Cross Section 
  for (int s=0; s<7; s++)
    {
      c1->Divide(can_size+1, can_size);

      for (int b=0; b<h1_gen_by_phi[s].size()-1; b++)
	{
	  c1->cd(b+1);

	  double phi_start = phiBins->min() + (b-1)*phiBins->width();
	  double phi_end   = phiBins->min() + b*phiBins->width();
	  if (b == 0) { phi_start = phiBins->min(); phi_end = phiBins->max(); }

	  h1_gen_by_phi[s][b]->SetMarkerStyle(8);
	  h1_gen_by_phi[s][b]->Draw();
	  lab.DrawLatex(0.3, 0.02, Form(" Generated: Sector %d Bin %d ",s,b));
	  lab.DrawLatex(0.55, 0.8, Form(" %.2f #rightarrow %.2f #circ ",phi_start, phi_end));
	}

      c1->Print( Form("../img/%s.pdf",name.c_str()) );
      c1->Clear();
    }

  // Doing 1-D Cross Section 
  for (int s=0; s<7; s++)
    {
      c1->Divide(can_size+1, can_size);

      for (int b=0; b<h1_rec_by_phi[s].size()-1; b++)
	{
	  c1->cd(b+1);

	  double phi_start = phiBins->min() + (b-1)*phiBins->width();
	  double phi_end   = phiBins->min() + b*phiBins->width();
	  if (b == 0) { phi_start = phiBins->min(); phi_end = phiBins->max(); }

	  h1_rec_by_phi[s][b]->SetMarkerStyle(8);
	  h1_rec_by_phi[s][b]->Draw();
	  lab.DrawLatex(0.3, 0.02, Form(" Reconstructed: Sector %d Bin %d ",s,b));
	  lab.DrawLatex(0.55, 0.8, Form(" %.2f #rightarrow %.2f #circ ",phi_start, phi_end));
	}

      c1->Print( Form("../img/%s.pdf",name.c_str()) );
      c1->Clear();
    }

  // Doing 1-D Cross Section 
  for (int s=0; s<7; s++)
    {
      c1->Divide(can_size+1, can_size);

      for (int b=0; b<h1_hits_by_phi[s].size()-1; b++)
	{
	  c1->cd(b+1);

	  double phi_start = phiBins->min() + (b-1)*phiBins->width();
	  double phi_end   = phiBins->min() + b*phiBins->width();
	  if (b == 0) { phi_start = phiBins->min(); phi_end = phiBins->max(); }

	  h1_hits_by_phi[s][b]->SetMarkerStyle(8);
	  h1_hits_by_phi[s][b]->Draw();
	  lab.DrawLatex(0.3, 0.02, Form(" Data hits: Sector %d Bin %d ",s,b));
	  lab.DrawLatex(0.55, 0.8, Form(" %.2f #rightarrow %.2f #circ ",phi_start, phi_end));
	}

      c1->Print( Form("../img/%s.pdf",name.c_str()) );
      c1->Clear();
    }


  // Doing W cut by momentum bins 
  can_size = floor(sqrt(h1_w_by_mom[0].size()));
  for (int s=0; s<7; s++)
    {
      c1->Divide(can_size+1, can_size);
      c1->SetMargin(0.15, 0.1, 0.15, 0.15);

      for (int b=0; b<h1_w_by_mom[s].size()-1; b++)
	{
	  c1->cd(b+1);

	  double p_start = momBins->min() + (b-1)*momBins->width();
	  double p_end   = momBins->min() + b*momBins->width();

	  if (b == 0) { p_start = momBins->min(); p_end = momBins->max(); }

	  h1_w_by_mom[s][b]->SetFillStyle(3004);
	  h1_w_by_mom[s][b]->SetFillColorAlpha(kBlack, 1);
	  h1_w_by_mom[s][b]->SetLineColor(kBlack);
	  h1_w_by_mom[s][b]->Draw();
	  lab.DrawLatex(0.4, 0.01, Form(" Sector %d Bin %d ",s,b));
	  lab.DrawLatex(0.55, 0.8, Form(" %.2f #rightarrow %.2f GeV/c ",p_start, p_end));
	}

      c1->Print( Form("../img/%s.pdf",name.c_str()) );
      c1->Clear();
    }

  // Doing some plots of data with mc overlayed 
  c1->Divide(3,2);

  // Normalizing so that the sum over all 6 sectors 
  // gives 6. 
  int mc_entries   = 0.0;
  int data_entries = 0.0;

  for (int s=1; s<7; s++) { mc_entries += h1_rec_by_phi[s][0]->Integral();  data_entries += h1_hits_by_phi[s][0]->Integral(); }

  cout << " Found normalization for Monte Carlo reconstructed events by sum: " << mc_entries << " direct method: " << h1_rec_by_phi[0][0]->GetEntries() << endl; 
  cout << " Found normalization for Data reconstructed events by sum: " << data_entries << " direct method: " << h1_hits_by_phi[0][0]->GetEntries() << endl; 

  for (int s=1; s<7; s++)
    {
      c1->cd(s);
      c1->SetMargin(0.1, 0.1, 0.2, 0.1);
      h1_rec_by_phi[s][0] ->Scale( (double) 6/mc_entries);
      h1_hits_by_phi[s][0]->Scale( (double) 6/data_entries);
      h1_rec_by_phi[s][0]->SetFillColorAlpha(kMagenta, 1.0);
      h1_rec_by_phi[s][0]->SetLineColor(kMagenta);
      h1_rec_by_phi[s][0]->SetFillStyle(3004);
      h1_hits_by_phi[s][0]->Draw("PE");
      h1_rec_by_phi[s][0] ->Draw("HISTsame");

      lab.DrawLatex(0.47, 0.02, Form(" Sector %d ",s) );
      lab.DrawLatex(0.65, 0.8, " #rightarrow Data ");
      lab.SetTextColor(kMagenta);
      lab.DrawLatex(0.65, 0.75, " #rightarrow Rec. ");
      lab.SetTextColor(kBlack);
    }

  c1->Print(Form("../img/%s.pdf",name.c_str()));
  c1->Print(Form("../img/%s.pdf]",name.c_str()));

  // This takes a suprisingly long time 
  if (gif_status) 
    {
      cout << "\n Drawing .gif images, please wait.. " << endl; 
      
      // Draw .gif files for events 
      TCanvas * gifCanvas = new TCanvas("gifCanvas","",1200,800);
      int gifStep = 20;
      
      gifCanvas->Divide(3,2);
      for (int b=1; b<h1_hits_by_phi[0].size(); b++)
	{
	  for (int s=1; s<7; s++)
	    {
	      gifCanvas->cd(s);
	      h1_hits_by_phi[s][b]->Draw();
	    }
	  string gifName = Form("../gif/%s_h1_hits_by_phi.gif+%d",name.c_str(),gifStep);
	  gifCanvas->Print( gifName.c_str() );
	}
      gifCanvas->Clear();
      
      // Reconstructed Events 
      gifCanvas->Divide(3,2);
      for (int b=1; b<h1_hits_by_phi[0].size(); b++)
	{
	  for (int s=1; s<7; s++)
	    {
	      gifCanvas->cd(s);
	      h1_rec_by_phi[s][b]->Draw();
	    }
	  string gifName = Form("../gif/%s_h1_rec_by_phi.gif+%d",name.c_str(),gifStep);
	  gifCanvas->Print( gifName.c_str() );
	}
      gifCanvas->Clear();
      
      // Generated Events 
      gifCanvas->Divide(3,2);
      for (int b=1; b<h1_hits_by_phi[0].size(); b++)
	{
	  for (int s=1; s<7; s++)
	    {
	      gifCanvas->cd(s);
	      h1_gen_by_phi[s][b]->Draw();
	    }
	  string gifName = Form("../gif/%s_h1_gen_by_phi.gif+%d",name.c_str(),gifStep);
	  gifCanvas->Print( gifName.c_str() );
	}
      gifCanvas->Clear();
      
      // Acceptance 
      gifCanvas->Divide(3,2);
      for (int b=1; b<h1_hits_by_phi[0].size(); b++)
	{
	  for (int s=1; s<7; s++)
	    {
	      gifCanvas->cd(s);
	      h1_acc_by_phi[s][b]->Draw();
	}
	  string gifName = Form("../gif/%s_h1_acc_by_phi.gif+%d",name.c_str(),gifStep);
	  gifCanvas->Print( gifName.c_str() );
	}
      gifCanvas->Clear();
      
      // Cross Section 
      gifCanvas->Divide(3,2);
      for (int b=1; b<h1_hits_by_phi[0].size(); b++)
	{
	  for (int s=1; s<7; s++)
	    {
	      gifCanvas->cd(s);
	      h1_xs_by_phi[s][b]->Draw();
	      h1_xs_model[0]->Draw("Lsame");
	    }
	  string gifName = Form("../gif/%s_h1_xs_by_phi.gif+%d",name.c_str(),gifStep);
	  gifCanvas->Print( gifName.c_str() );
	}
      gifCanvas->Clear();
      
      // Cross Section Ratio 
      gifCanvas->Divide(3,2);
      for (int b=1; b<h1_hits_by_phi[0].size(); b++)
	{
	  for (int s=1; s<7; s++)
	    {
	      gifCanvas->cd(s);
	      h1_xs_ratio_by_phi[s][b]->Draw();
	    }
	  string gifName = Form("../gif/%s_h1_xs_ratio_by_phi.gif+%d",name.c_str(),gifStep);
	  gifCanvas->Print( gifName.c_str() );
	}
      gifCanvas->Clear();
      
      // Raw Cross Section 
      gifCanvas->Divide(3,2);
      for (int b=1; b<h1_hits_by_phi[0].size(); b++)
	{
	  for (int s=1; s<7; s++)
	    {
	      gifCanvas->cd(s);
	      h1_rxs_by_phi[s][b]->Draw();
	      h1_xs_model[1]->Draw("Lsame");
	    }
	  string gifName = Form("../gif/%s_h1_rxs_by_phi.gif+%d",name.c_str(),gifStep);
	  gifCanvas->Print( gifName.c_str() );
	}
      gifCanvas->Clear();
      
      // Raw Cross Section Ratio 
      gifCanvas->Divide(3,2);
      for (int b=1; b<h1_hits_by_phi[0].size(); b++)
	{
	  for (int s=1; s<7; s++)
	    {
	      gifCanvas->cd(s);
	      h1_rxs_ratio_by_phi[s][b]->Draw();
	    }
	  string gifName = Form("../gif/%s_h1_rxs_ratio_by_phi.gif+%d",name.c_str(),gifStep);
	  gifCanvas->Print( gifName.c_str() );
	}
      gifCanvas->Clear();
    }
  
  
}

#endif
