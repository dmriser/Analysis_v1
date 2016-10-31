/* 

      allHistograms.C 

        David Riser, 
   University of Connecticut 
        Mar 11, 2016


*/

// C++ Libraries
#include <iostream>
#include <cstdlib>
using namespace std;

// My headers
#include "../../analysisClasses/Bins.h"
#include "../../analysisClasses/ElasticEvent.h"
#include "../../analysisClasses/FaradayReader.h"
#include "../../analysisClasses/h22Event.h"
#include "../../analysisClasses/h22Reader.h"
#include "../../analysisClasses/ParticleFilter.h"
#include "../../analysisClasses/HTMLPrinter.h"

// CERN Root Libraries
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"

// Nathan's headers 
#include "../../analysisClasses/programFiles/sctimeCorr.C" // for time corrections in the SC 

// Marco
#include "../../analysisClasses/momCorr/MomCorr.C"

#define speed_of_light 29.979

int main (int argc, char * argv[])
{
  int GSIM[3] = {0,1,1};                          // data (0 - false), gsim (1 - true)
  
  // get number of files from command line 
  if (argc < 2)
    {
      cout << " expected number of files as option " << endl;
      exit(0);
    }

  TFile myFile("allHistograms.root","recreate");

  // our data/GSIM setup 
  int nFiles[3] = {atoi(argv[1]),atoi(argv[1])*20,atoi(argv[1])*20};
  string files[3] = {"../fcup/goodruns.txt", "mc_norad.txt","mc_rad.txt"};
  string type[3] = {"data","GSIM_norad","GSIM_rad"};

  // Momentum Corrections 
  MomCorr_e1f * momCorr = new MomCorr_e1f();

  // setup ParticleFilter for PID
  ParticleFilter filter;


  // ------------------------- histograms ----------------------------
  string histNames[7] = {"all","s1","s2","s3","s4","s5","s6"};

  TH1F * h_W[3][7];
  TH1F * h_x[3][7];
  TH1F * h_QQ[3][7];
  TH1F * h_dt[3][7];

  TH2F * h_dt_p[3][7];
  TH2F * h_x_W[3][7];
  TH2F * h_x_QQ[3][7];
  TH2F * h_QQ_W[3][7];
  TH2F * h_ang_fid[3][7];

  for (int itype = 0; itype < 3; itype++)
    for (int ihist=0; ihist<7; ihist++)
      {
	h_W[itype][ihist]      = new TH1F(Form("h_W_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),Form("h_W_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),100,0.7,3);
	h_x[itype][ihist]      = new TH1F(Form("h_x_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),Form("h_x_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),100,0,1);
	h_QQ[itype][ihist]     = new TH1F(Form("h_QQ_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),Form("h_QQ_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),100,0,5);
	h_dt[itype][ihist]     = new TH1F(Form("h_dt_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),Form("h_dt_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),100,-2,2);
	h_dt_p[itype][ihist]   = new TH2F(Form("h_dt_p_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),Form("h_dt_p_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),250,0,5,250,-4,4);      
	h_x_W[itype][ihist]    = new TH2F(Form("h_x_W_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),Form("h_x_W_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),250,0,1,250,0,5);
	h_x_QQ[itype][ihist]   = new TH2F(Form("h_x_QQ_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),Form("h_x_QQ_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),250,0,1,250,0,5);
	h_QQ_W[itype][ihist]   = new TH2F(Form("h_QQ_W_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),Form("h_QQ_W_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),250,0,5,250,0,5);
	h_ang_fid[itype][ihist] = new TH2F(Form("h_ang_fid_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),Form("h_ang_fid_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),250,-30,30,250,0,60);
      }
  // -------------------------- end histograms -------------------------

  // setup file reader and add files
  h22Reader * reader[2];

  // loop over data, GSIM
  for (int itype = 0; itype < 3; itype++)
    {  
      reader[itype] = new h22Reader(GSIM[itype]);

      reader[itype]->AddList(files[itype],nFiles[itype]);
      reader[itype]->Init(); // set branch addresses
  
  int nEvents = reader[itype]->GetEntries();
  
  // loop over events
  for (int iEvent = 0; iEvent < nEvents; iEvent++)
    {
      reader[itype]->GetEntry(iEvent);

      // store each event locally for the duration of the loop 
      h22Event event = reader[itype]->GetEvent();

      // setting up PID 
      // for mc it doesn't matter what run number is
      string run    = reader[itype]->GetFilename();
      string srunno = run.substr(57,5);
      int runno     = atoi(srunno.c_str());

      filter.loadEvent(event,GSIM[itype],runno);

      int e_index = filter.getIndexByPID(11);
      if (e_index > -123) 
	{
	  int sector = event.dc_sect[e_index]%1000;

	  TLorentzVector electron(event.cx[e_index]*event.p[e_index],event.cy[e_index]*event.p[e_index],event.cz[e_index]*event.p[e_index],event.p[e_index]);

	  // do mom corrections on data
	  if (GSIM[itype] == 0) electron = momCorr->PcorN(electron, -1, 11);
	  ElasticEvent elastic(electron);

	  // fill "all" spots
	  h_W[itype][0]     ->Fill(elastic.getW());
	  h_x[itype][0]     ->Fill(elastic.getx());
	  h_QQ[itype][0]    ->Fill(elastic.getQQ());
	  h_x_W[itype][0]   ->Fill(elastic.getx(),elastic.getW());
	  h_x_QQ[itype][0]  ->Fill(elastic.getx(),elastic.getQQ());
	  h_QQ_W[itype][0]  ->Fill(elastic.getQQ(),elastic.getW());
	  h_ang_fid[itype][0]->Fill(event.getRelativePhi(e_index),event.getTheta(e_index));

	  // fill by sector
	  if (sector > 0)
	    { 
	      h_W[itype][sector]     ->Fill(elastic.getW());
	      h_x[itype][sector]     ->Fill(elastic.getx());
	      h_QQ[itype][sector]    ->Fill(elastic.getQQ());
	      h_x_W[itype][sector]   ->Fill(elastic.getx(),elastic.getW());
	      h_x_QQ[itype][sector]  ->Fill(elastic.getx(),elastic.getQQ());
	      h_QQ_W[itype][sector]  ->Fill(elastic.getQQ(),elastic.getW());
	      h_ang_fid[itype][sector]->Fill(event.getRelativePhi(e_index),event.getTheta(e_index));
	    }

	  // jumping through some hoops to get timing info 
	  double startTime = e_sctimeCorr(!GSIM[itype],event.sc_t[e_index],sector,event.sc_pd[e_index],runno) - event.sc_r[e_index]/speed_of_light;

	  for(int ipart=0; ipart<event.gpart; ipart++)
	    if (event.q[ipart] > 0)
	      {
		int sector     = event.dc_sect[ipart]%1000;
		Float_t cbeta  = speed_of_light*event.p[ipart]/sqrt(event.p[ipart]*event.p[ipart] + 0.938*0.938);
		Float_t t_calc = event.sc_r[ipart]/cbeta;
		Float_t t_exp  = h_sctimeCorr(!GSIM[itype], event.sc_t[ipart], sector, event.sc_pd[ipart], runno) - startTime;

		h_dt[itype][0]       ->Fill(t_calc-t_exp);
		h_dt[itype][sector]  ->Fill(t_calc-t_exp);
		h_dt_p[itype][0]     ->Fill(event.p[ipart],t_calc-t_exp);
		h_dt_p[itype][sector]->Fill(event.p[ipart],t_calc-t_exp);
	      }
	}
      cout << "\r > " << iEvent << " of " << nEvents;
    }  // end loop over events
    }

  myFile.Write();
  myFile.Close();

  return 0;
}
