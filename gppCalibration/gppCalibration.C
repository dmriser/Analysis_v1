/* 

      gppCalibration.cc 

   Code used to check gpp smearing parameters...
     ===>  a = b = c , dc smearing , check momentum    -> check W dist.
     ===>  f       , ftof smearing , check timing info -> check dt dist.

   Just run this once on all the files you have
   then use drawGPP.C to do fits/ pretty plots 
          <user> root -l drawGPP.C

        David Riser, 
   University of Connecticut 
        Feb 29, 2016


*/

// C++ Libraries
#include <iostream>
#include <cstdlib>
using namespace std;

// My headers
#include "../analysisClasses/Bins.h"
#include "../analysisClasses/ElasticEvent.h"
#include "../analysisClasses/FaradayReader.h"
#include "../analysisClasses/h22Event.h"
#include "../analysisClasses/h22Reader.h"
#include "../analysisClasses/ParticleFilter.h"

// CERN Root Libraries
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"

// Nathan's headers 
#include "../analysisClasses/programFiles/sctimeCorr.C" // for time corrections in the SC 

// Marco
#include "../analysisClasses/momCorr/MomCorr.C"

#define speed_of_light 29.979

int main (int argc, char * argv[])
{
  int GSIM[2] = {0,1};                          // data (0 - false), gsim (1 - true)
  
  // get number of files from command line 
  if (argc < 2)
    {
      cout << " expected number of files as option " << endl;
      exit(0);
    }

  TFile myFile("gpp.root","recreate");

  int nFiles = atoi(argv[1]);

  //string files[2] = {"/volatile/clas/clas12/dmriser/nathanFiles/AnalysisCode/datafiles.txt",
  //		     "mcfiles.txt"};

  string files[2] = {"firstfiles.txt", "mcfiles.txt"};

  string type[2] = {"data","GSIM"};

  // Momentum Corrections 
  MomCorr_e1f * momCorr = new MomCorr_e1f();

  // setup ParticleFilter for PID
  ParticleFilter filter;


  // ------------------------- histograms ----------------------------
  string histNames[7] = {"all","s1","s2","s3","s4","s5","s6"};

  TH1F * h_W[2][7];
  TH1F * h_dt[2][7];
  TH2F * h_dt_p[2][7];

  for (int itype = 0; itype < 2; itype++)
    for (int ihist=0; ihist<7; ihist++)
      {
	h_W[itype][ihist]    = new TH1F(Form("h_W_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),Form("h_W_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),100,0.7,1.3);
	h_dt[itype][ihist]   = new TH1F(Form("h_dt_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),Form("h_dt_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),100,-2,2);
	h_dt_p[itype][ihist] = new TH2F(Form("h_dt_p_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),Form("h_dt_p_%s_%s",type[itype].c_str(),histNames[ihist].c_str()),250,0,5,250,-4,4);
      }
  // -------------------------- end histograms -------------------------

  // setup file reader and add files
  h22Reader * reader[2];

  // loop over data, GSIM
  for (int itype = 0; itype < 2; itype++)
    {  
      reader[itype] = new h22Reader(GSIM[itype]);

      reader[itype]->AddList(files[itype],nFiles);
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
      string srunno = run.substr(68,6);
      int runno     = atoi(srunno.c_str());

      filter.loadEvent(event,GSIM[itype],runno);

      int e_index = filter.getIndexByPID(11);
      if (e_index > -123) 
	{
	  int sector = event.dc_sect[e_index]%1000;

	  TLorentzVector electron(event.cx[e_index]*event.p[e_index],event.cy[e_index]*event.p[e_index],event.cz[e_index]*event.p[e_index],event.p[e_index]);

	  // do mom corrections
	  if (GSIM[itype] == 0) electron = momCorr->PcorN(electron, -1, 11);
	  ElasticEvent elastic(electron);

	  h_W[itype][0]->Fill(elastic.getW());
	  if (sector > 0) h_W[itype][sector]->Fill(elastic.getW());

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
