/* 

      gppTOFStudy.C 

      Code used to generate curve dt vs. f param
      time of flight smearing parameter.

        David Riser, 
   University of Connecticut 
        March 11, 2016


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
  int GSIM  = 1;   // GSIM -> TRUE
  int nFiles = atoi(argv[1]); 

  // get number of files from command line 
  if (argc < 2)
    {
      cout << " expected number of files as option " << endl;
      exit(0);
    }

  TFile myFile("gppTOFStudy.root","recreate");

  string files[11]  = {"tof.0.75.txt", "tof.0.80.txt", "tof.0.85.txt",
		     "tof.0.90.txt", "tof.0.95.txt", "tof.1.00.txt", 
		     "tof.1.05.txt", "tof.1.10.txt", "tof.1.15.txt",
		     "tof.1.20.txt", "tof.1.25.txt"};

  double type[11] = {0.75,0.80,0.85,0.90,0.95,1.00,1.05,1.10,1.15,1.20,1.25};

  // setup ParticleFilter for PID
  ParticleFilter filter;


  // ------------------------- histograms ----------------------------
  string histNames[7] = {"all","s1","s2","s3","s4","s5","s6"};

  TH1F * h_dt[11][7];

  for (int itype = 0; itype < 11; itype++)
    for (int ihist=0; ihist<7; ihist++)
      {
	h_dt[itype][ihist]   = new TH1F(Form("h_dt_%f_%s",type[itype],histNames[ihist].c_str()),Form("h_dt_%f_%s",type[itype],histNames[ihist].c_str()),100,-2,2);
      }
  // -------------------------- end histograms -------------------------

  // setup file reader and add files
  h22Reader * reader[11];

  // loop over data, GSIM
  for (int itype = 0; itype < 11; itype++)
    {  
      reader[itype] = new h22Reader(GSIM);

      reader[itype]->AddList(files[itype],nFiles);
      reader[itype]->Init(); // set branch addresses
  

  int nEvents = reader[itype]->GetEntries();

  
  // loop over events
  for (int iEvent = 0; iEvent < nEvents; iEvent++)
    {
      reader[itype]->GetEntry(iEvent);

      // store each event locally for the duration of the loop 
      h22Event event = reader[itype]->GetEvent();

      // doesnt matter for MC 
      //      int runno = atoi(reader[itype]->GetFilenameChunk(68,6).c_str());
      int runno = 000000;
      filter.loadEvent(event,GSIM,runno);

      int e_index = filter.getIndexByPID(11);
      if (e_index > -123) 
	{
	  int sector = event.dc_sect[e_index]%1000;

	  // jumping through some hoops to get timing info, not that my bool value is opposite of NH convention, we pass !GSIM which means 
	  // it is GSIM not data.  
	  double startTime = e_sctimeCorr(!GSIM,event.sc_t[e_index],sector,event.sc_pd[e_index],runno) - event.sc_r[e_index]/speed_of_light;

	  for(int ipart=0; ipart<event.gpart; ipart++)
	    if (event.q[ipart] > 0)
	      {
		int sector     = event.dc_sect[ipart]%1000;
		Float_t cbeta  = speed_of_light*event.p[ipart]/sqrt(event.p[ipart]*event.p[ipart] + 0.938*0.938);
		Float_t t_calc = event.sc_r[ipart]/cbeta;
		Float_t t_exp  = h_sctimeCorr(!GSIM, event.sc_t[ipart], sector, event.sc_pd[ipart], runno) - startTime;

		h_dt[itype][0]       ->Fill(t_calc-t_exp);
		h_dt[itype][sector]  ->Fill(t_calc-t_exp);
	      }
	}
      cout << "\r > " << iEvent << " of " << nEvents;
    }  // end loop over events
    }
  
  // ------------- fitting histogams ----------------
  TF1 * f_dt[11][7];

  for (int itype = 0; itype < 11; itype++)
    for (int ihist = 0; ihist < 7; ihist++)
      {
	f_dt = new TF1(Form("f_dt_%f_%s",type[itype],histNames[ihist].c_str()),"[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*x^2 + [4]*x + [5]",-1.5,1.5);
	f_dt[itype][ihist]->SetParameter(0,h_dt[itype][ihist]->GetMaximum());
	f_dt[itype][ihist]->SetParameter(1,h_dt[itype][ihist]->GetMean());
	f_dt[itype][ihist]->SetParameter(2,h_dt[itype][ihist]->GetRMS());
	f_dt[itype][ihist]->SetParameter(3,-0.01);
	f_dt[itype][ihist]->SetParameter(4,0.01);
	f_dt[itype][ihist]->SetParameter(5,0.01);
	h_dt[type][ihist]->Fit(f_dt[itype][ihist]);
      }

  myFile.Write();
  myFile.Close();

  return 0;
}
