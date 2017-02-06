/* 

   rates.C compare fcup charge to particle rates

       David Riser, 
   University of Connecticut 
       Mar 1, 2016

 */

// C++ Libraries
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;


// My Libraries
#include "../analysisClasses/h22Event.h"
#include "../analysisClasses/h22Reader.h"
#include "../analysisClasses/Bins.h"
#include "../analysisClasses/FaradayReader.h"
#include "../analysisClasses/ParticleFilter.h"

// CERN Root Libraries
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"

int main (int argc, char * argv[])
{

  // data (0 - false), gsim (1 - true)
  int GSIM = 0; 
  
  // ---------------- file IO ---------------
  if (argc < 3)
    {
      cout << " usage: rates nfiles startfile" << endl;
      exit(0);
    }

  int nFiles    = atoi(argv[1]);
  int startFile = atoi(argv[2]);
  int ifile     = 0;
  int listIndex = 0;

  string files = "../combinedFiles.txt";
  TFile outfile("ratesFBF.root","recreate");

  ifstream list;
  list.open(files.c_str());
  // ------------- end file IO --------------

  // particle id filter
  ParticleFilter filter;

  // ----------- Graphs, Histograms, Vectors ----------------

  vector<double> fileNumber;
  vector<int> electronsPerFile;
  vector<double> dQPerFile; 
  vector<double> v_e_dQ_ratio;

  string histname[7] = {"all","s1","s2","s3","s4","s5","s6"};
  TH1F * h_nElectrons[7];
  TH1F * h_electron_dQ_ratio[7];

  for (int ihist=0; ihist<7; ihist++)
    {
      h_nElectrons[ihist]        = new TH1F(Form("h_nElectrons_%s",histname[ihist].c_str()), Form("h_nElectrons_%s",histname[ihist].c_str()), 100, 0, 1e6);
      h_electron_dQ_ratio[ihist] = new TH1F(Form("h_electron_dQ_ratio_%s",histname[ihist].c_str()), Form("h_electron_dQ_ratio_%s",histname[ihist].c_str()), 50, 0, 20); 
    }

  // --------------------------------------------------------


  // loop over files 
  while ( !list.eof() && ifile < nFiles )
    {
      // skip to starting file
      while (listIndex < startFile) 
	{
	  string tempString;
	  list >> tempString;
	  cout << " passing " << tempString << endl;
	  listIndex++;
	}
      // setup file reader and add files
      h22Reader reader(GSIM);

      string currentFullFilename;
      list >> currentFullFilename;

      reader.AddFile(currentFullFilename);
      reader.Init();

      // if we have events in file
      if (reader.GetEntries() > 0)
	{
	  int runno = atoi(reader.GetFilenameChunk(57,5).c_str());      
	  
	  // doing file-by-file mode
	  FaradayCup fcupReader(runno);
	  double dQ = fcupReader.getFcupCharge();
	  
	  // how many entries each file has 
	  int nEvents = reader.GetEntries();

	  int nElectrons[7] = {};
	  
	  // loop over events in each file 
	  for (int iEvent = 0; iEvent < nEvents; iEvent++)
	    {
	      reader.GetEntry(iEvent);
	      
	      h22Event event = reader.GetEvent();
	      filter.loadEvent(event,GSIM,runno);
	      
	      // asking for electron, proceeding if found
	      int e_index = filter.getIndexByPID(11);
	      if (e_index > -123)
		{
		  int sector = event.dc_sect[e_index];
		  
		  nElectrons[0]++;
		  nElectrons[sector]++;
		}
	      
	    }  // end loop over events

	  // ------------- status update ------------
	  cout.width(12);
	  cout << ifile;
	  cout.width(12);
	  cout << dQ;
	  cout.width(12);
	  cout << nElectrons[0];
	  cout.width(12);
	  cout << (double) nElectrons[0]/dQ << endl;
	  // -------------end status update---------
	  
	  // update the overall information before proceed to next file
	  fileNumber      .push_back((double)runno);
	  dQPerFile       .push_back(dQ);
	  electronsPerFile.push_back(nElectrons[0]);
	  if ((double)nElectrons[0]/dQ > 0) v_e_dQ_ratio    .push_back((double)nElectrons[0]/dQ);
	  
	  for (int i=0; i<7; i++) 
	    {
	      h_nElectrons[i]        ->Fill(nElectrons[i]);
	      h_electron_dQ_ratio[i] ->Fill(nElectrons[i]/dQ);
	    }     
	  
	  
	  ifile++;
	  
	} // end conditional on has file
    } // end loop on files

  cout << endl;  

  // ---------------- graph results over all files --------------
  TGraph g_electron_dQ_ratio(fileNumber.size(), &(fileNumber[0]), &v_e_dQ_ratio[0]);
  g_electron_dQ_ratio.SetMarkerStyle(7);
  g_electron_dQ_ratio.SetMarkerColor(kBlue);
  g_electron_dQ_ratio.Write();  
  // -------------- results over all files ---------------------

  // -------------- we need to generate good run list ---------- 
  
  TF1 * f_dN_dQ = new TF1("f_dN_dQ","gaus",8,14);
  h_electron_dQ_ratio[0]->Fit("f_dN_dQ","rq");

  double NSIGMA = 4.00;
  double mu     = f_dN_dQ->GetParameter(1);
  double sigma  = f_dN_dQ->GetParameter(2);
  double left   = mu - NSIGMA*sigma;    
  double right  = mu + NSIGMA*sigma;    
  
  ofstream goodRunList("goodruns.txt",ios::trunc);

  for (int i=0; i<v_e_dQ_ratio.size(); i++)
    {
      if (v_e_dQ_ratio[i] > left && v_e_dQ_ratio[i] < right)
	{
	  goodRunList << Form("/volatile/clas/clas12/dmriser/analysis/e1f_analysis/skim/%d.root",(int)fileNumber[i]) << endl;
	}
    }

  goodRunList.close();
  list.close();
    
  outfile.Write();
  outfile.Close();
  
  return 0;
}
  
  

