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
#include "TH1.h"
#include "TH2.h"

elast_histos::elast_histos()
{
  // birth 
}

elast_histos::~elast_histos()
{
  // total destruction 
}

void elast_histos::init(BinStructure thetaBins, BinStructure relPhiBins)
{

  string TYPE[3]     = {"GSIM_NORAD","GSIM_RAD","DATA"};
  string ALT_TYPE[3] = {"GEN","REC","DATA"};
  string SECT[7]     = {"ALL","S1","S2","S3","S4","S5","S6"};

  for (int t=0; t<3; t++)
    {
      for (int s=0; s<7; s++)
	{
	  // 1-D
	  h1_W[t][s]         = new TH1F(Form("h1_W_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),Form("h1_W_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),200,0.8,4);
	  h1_QQ[t][s]        = new TH1F(Form("h1_QQ_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),Form("h1_QQ_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),100,0.5,3.5);
	  h1_x[t][s]         = new TH1F(Form("h1_x_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),Form("h1_x_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),100,0.0,1.1);
	  h1_EP_angle[t][s]  = new TH1F(Form("h1_EP_angle_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),Form("h1_EP_angle_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),100,0,180);
	  h1_MM2[t][s]       = new TH1F(Form("h1_MM2_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),Form("h1_MM2_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),100,-0.1,3.5);
	  h1_ME[t][s]       = new TH1F(Form("h1_ME_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),Form("h1_ME_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),100,-0.1,3.5);
	  h1_rec_theta[t][s] = new TH1F(Form("h1_rec_theta_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),Form("h1_rec_theta_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),thetaBins.number(),thetaBins.min(),thetaBins.max());
	  h1_gen_theta[t][s] = new TH1F(Form("h1_gen_theta_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),Form("h1_gen_theta_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),thetaBins.number(),thetaBins.min(),thetaBins.max());

	  // 2-D 
	  h2_rec[t][s]  = new TH2F(Form("h2_rec_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),Form("h2_rec_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),relPhiBins.number(),relPhiBins.min(),relPhiBins.max(),thetaBins.number(),thetaBins.min(),thetaBins.max());
	  h2_gen[t][s]  = new TH2F(Form("h2_gen_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),Form("h2_gen_%s_%s",TYPE[t].c_str(),SECT[s].c_str()),relPhiBins.number(),relPhiBins.min(),relPhiBins.max(),thetaBins.number(),thetaBins.min(),thetaBins.max());
	  h2_W_theta[t][s] = new TH2F(Form("h2_W_theta_%s_%s",ALT_TYPE[t].c_str(),SECT[s].c_str()),Form("h2_W_theta_%s_%s",ALT_TYPE[t].c_str(),SECT[s].c_str()),100,0,5,thetaBins.number(),thetaBins.min(),thetaBins.max());
	}      
    }

  h2_scale = new TH2F("h2_scale","h2_scale",relPhiBins.number(),relPhiBins.min(),relPhiBins.max(),thetaBins.number(),thetaBins.min(),thetaBins.max());
  h1_fc_mon = new TH1F("h1_fc_mon","h1_fc_mon",100,0,40);

  // with and without radiation 
  string MODELTYPE[2] = {"NORAD","RAD"};
  for (int t=0; t<2; t++) h1_xs_model[t] = new TH1F(Form("h1_xs_model_%s",MODELTYPE[t].c_str()),Form("h1_xs_model_%s",MODELTYPE[t].c_str()),thetaBins.number(),thetaBins.min(),thetaBins.max());

  h1_e_E       = new TH1F("h1_e_E"," Electron Energy ", 400, 0, 6);
  h1_gamma_E       = new TH1F("h1_gamma_E"," BH Photon Energy ", 400, 0, 6);
  h1_gamma_theta   = new TH1F("h1_gamma_theta"," BH Photon Theta ", 400, 0, 90);
  h1_gamma_phistar = new TH1F("h1_gamma_phistar"," BH Photon Phi* ", 400, 0, 200);

  h1_proton_E       = new TH1F("h1_proton_E"," Proton Energy ", 400, 0, 6);
  h1_proton_theta   = new TH1F("h1_proton_theta"," Proton Theta ", 400, 0, 90);
  h1_proton_phistar = new TH1F("h1_proton_phistar"," Proton Phi* ", 400, 0, 200);
 
}

void elast_histos::load(string _file)
{
  
  TFile *f = TFile::Open(_file.c_str());

  string ALT_TYPE[3] = {"GEN","REC","DATA"};
  string TYPE[3]     = {"GSIM_NORAD","GSIM_RAD","DATA"};
  string SECT[7]     = {"ALL","S1","S2","S3","S4","S5","S6"};

  for (int t=0; t<3; t++)
    {
      for (int s=0; s<7; s++)
	{
	  // 1-D
	  h1_W[t][s]         = (TH1F*) f->Get(Form("h1_W_%s_%s",TYPE[t].c_str(),SECT[s].c_str()));
	  h1_QQ[t][s]        = (TH1F*) f->Get(Form("h1_QQ_%s_%s",TYPE[t].c_str(),SECT[s].c_str()));
	  h1_x[t][s]         = (TH1F*) f->Get(Form("h1_x_%s_%s",TYPE[t].c_str(),SECT[s].c_str()));
	  h1_EP_angle[t][s]  = (TH1F*) f->Get(Form("h1_EP_angle_%s_%s",TYPE[t].c_str(),SECT[s].c_str()));
	  h1_MM2[t][s]       = (TH1F*) f->Get(Form("h1_MM2_%s_%s",TYPE[t].c_str(),SECT[s].c_str()));
	  h1_ME[t][s]       = (TH1F*) f->Get(Form("h1_ME_%s_%s",TYPE[t].c_str(),SECT[s].c_str()));
	  h1_rec_theta[t][s] = (TH1F*) f->Get(Form("h1_rec_theta_%s_%s",TYPE[t].c_str(),SECT[s].c_str()));
	  h1_gen_theta[t][s] = (TH1F*) f->Get(Form("h1_gen_theta_%s_%s",TYPE[t].c_str(),SECT[s].c_str()));
 
	  // 2-D 
	  h2_rec[t][s]     = (TH2F*) f->Get(Form("h2_rec_%s_%s",TYPE[t].c_str(),SECT[s].c_str()));
	  h2_gen[t][s]     = (TH2F*) f->Get(Form("h2_gen_%s_%s",TYPE[t].c_str(),SECT[s].c_str()));
	  h2_W_theta[t][s] = (TH2F*) f->Get(Form("h2_W_theta_%s_%s",ALT_TYPE[t].c_str(),SECT[s].c_str()));
	}      
    }

  h2_scale = (TH2F*) f->Get("h2_scale");
  h1_fc_mon = (TH1F*) f->Get("h1_fc_mon");

  // with and without radiation 
  string MODELTYPE[2] = {"NORAD","RAD"};
  for (int t=0; t<2; t++) h1_xs_model[t] = (TH1F*) f->Get(Form("h1_xs_model_%s",MODELTYPE[t].c_str()));

  h1_e_E           = (TH1F*) f->Get("h1_e_E");
  h1_gamma_E       = (TH1F*) f->Get("h1_gamma_E");
  h1_gamma_theta   = (TH1F*) f->Get("h1_gamma_theta");
  h1_gamma_phistar = (TH1F*) f->Get("h1_gamma_phistar");

  h1_e_E           ->SetDirectory(0);
  h1_gamma_E       ->SetDirectory(0);
  h1_gamma_theta   ->SetDirectory(0);
  h1_gamma_phistar ->SetDirectory(0);

  h1_proton_E       = (TH1F*) f->Get("h1_proton_E");
  h1_proton_theta   = (TH1F*) f->Get("h1_proton_theta");
  h1_proton_phistar = (TH1F*) f->Get("h1_proton_phistar");

  h1_proton_E       ->SetDirectory(0);
  h1_proton_theta   ->SetDirectory(0);
  h1_proton_phistar ->SetDirectory(0);

  h2_scale       ->SetDirectory(0);
  h1_xs_model[0] ->SetDirectory(0);
  h1_xs_model[1] ->SetDirectory(0);

  for (int t=0; t<3; t++)
    for (int s=0; s<7; s++)
      {
	// 1D
	h1_W[t][s]        ->SetDirectory(0);
	h1_QQ[t][s]       ->SetDirectory(0);
	h1_x[t][s]        ->SetDirectory(0);
	h1_EP_angle[t][s] ->SetDirectory(0);
	h1_MM2[t][s]      ->SetDirectory(0);
	h1_ME[t][s]      ->SetDirectory(0);
	h1_rec_theta[t][s]->SetDirectory(0);
	h1_gen_theta[t][s]->SetDirectory(0);

	// 2D
	h2_rec[t][s]     ->SetDirectory(0);
	h2_gen[t][s]     ->SetDirectory(0);
	h2_W_theta[t][s] ->SetDirectory(0);
      }



}

void elast_histos::write_and_close(string _file)
{
  TFile f(_file.c_str(),"recreate");

  h2_scale->Write();
  h1_xs_model[0]->Write();
  h1_xs_model[1]->Write();
  h1_fc_mon->Write();

  h1_e_E->Write();
  h1_gamma_E->Write();
  h1_gamma_theta->Write();
  h1_gamma_phistar->Write();
  h1_proton_E->Write();
  h1_proton_theta->Write();
  h1_proton_phistar->Write();

  for (int t=0; t<3; t++)
    for (int s=0; s<7; s++)
      {
	// 1D
	h1_W[t][s]->Write();
	h1_QQ[t][s]->Write();
	h1_x[t][s]->Write();
	h1_EP_angle[t][s]->Write();
	h1_MM2[t][s]->Write();
	h1_ME[t][s]->Write();
	h1_rec_theta[t][s]->Write();
	h1_gen_theta[t][s]->Write();

	// 2D
	h2_rec[t][s]     ->Write();
	h2_gen[t][s]     ->Write();
	h2_W_theta[t][s] ->Write();
      }

  f.Write();
  f.Close();
}

#endif
