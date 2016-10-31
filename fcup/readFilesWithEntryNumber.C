/*

  David Riser 
  University of Connecticut 

  Mar 15, 2016


  readFilesWithEntryNumber.C 

  Program with strange name because I couldn't think of anything better. This
  code loads text files from database that contain the scalar bank FCUP_G2 entries
  and the next head bank event number following the entry. The code below looks at 
  differences in charge dQ between n_i and n_(i+1) as well as dN event difference.


 */

// C++ libs 
#include <fstream>
#include <iostream>
#include <vector>

using namespace std; 

// for class def
#include "readFilesWithEntryNumber.h"

// CERN root libs
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"

int main (int argc, char * argv[])
{

  int nNormal            = 0;
  int nEventsInNormal    = 0;

  int nExceptions        = 0;
  int nEventsInException = 0;

  int lastEvent = 0;

  ifstream runs("allfiles.txt");

  int ifile = 0;
  int nProcessedFiles = 0;

  string ifName;
  int pRunNumber = 0;

  TFile outfile("readFilesWithEntryNumber.root","recreate");

  // ----------------  histograms and graphs --------------------

  TH1D * h1_cdiff        = new TH1D("h1_cdiff","h1_cdiff",100,0,1000000); 
  TH1D * h1_ediff        = new TH1D("h1_ediff","h1_ediff",100,0,50000); 
  TH1D * h1_ediff_cdiff0 = new TH1D("h1_ediff_cdiff0","h1_ediff_cdiff0",100,0,50000); 
  TH1F * h1_ratio  = new TH1F("h1_ratio","h1_ratio",100,0,0.1);

  TH2D * h2_diff = new TH2D("h2_diff","h2_diff",200,0,1000000,200,0,50000);
  h2_diff->GetXaxis()->SetTitle("dQ");
  h2_diff->GetYaxis()->SetTitle("dN");
  // ---------------- end histograms and graphs ---------------

  vector<int> cdiffs;
  vector<int> ediffs;
  vector<int> entries;

  int iEntry = 0;

  // this is a flag, to signify if the events are within a good (dQ > 0)
  // or bad region (dQ = 0)
  bool BAD = false;

  // ------------- we need to have a few things to keep track ---------------
  // ------------ of the fcup charge for each run and make pdf --------------

  int iRunEntry = 1;
  int iCanvas   = 1;

  int xCan = 4;
  int yCan = 4;
  int canMax = xCan*yCan +1;

  vector<int> runStubs;
  vector<int> runEntries;
  vector<int> runCharges;

  vector<int> badEventStart;
  vector<int> badEventEnd;

  TCanvas * c1 = new TCanvas("c1","",800,800);

  // x*y runs per page
  c1->Divide(xCan,yCan);
  c1->Print("accumulation.pdf[");

  // --------------------- end pdf stuff ------------------------------------

  // loop over list of files from command line
  while (getline(runs, ifName) && ifile < atoi(argv[1]))
    {
      string ifTitle = ifName.substr(57,15);
      int cRunNumber = atoi(ifName.substr(62,6).c_str());

      // initialize the fcup reader class 
      FaradayReader fReader(ifTitle);

      // conditional that triggers at the start of 
      // every new run ex: 037801 -> 037802
      // but doesn't trigger on files 
      // 037801.A03 -> 037801.A04
      if (cRunNumber != pRunNumber && pRunNumber != 0)
	{
	  // make a plot for this current run, charge vs. entry 
	  TGraph * g = new TGraph(runEntries.size(),&runEntries[0],&runCharges[0]);
 	  
	  // find out where we are
	  if (iCanvas == canMax)
	    {
	      c1->Print("accumulation.pdf");
	      c1->Clear();
	      c1->Divide(xCan,yCan);
	      iCanvas = 1;
	    }
	  
	  c1->cd(iCanvas);
	  g->SetTitle(Form("run %d",pRunNumber));
	  g->SetMarkerStyle(7);
	  g->Draw();
	  
	  iCanvas++;

	  // this need to be done more carefully so to cut out gaps 
	  // such as .A02 -> .A04 skipping .A03 
	  // I am calling this a 'jump error'
	  int runTotalAccumulation = 0;
	  if (runCharges.size() > 0) 
	    {
	      for (int i=0; i<runCharges.size()-1; i++) 
		{
		  if (runStubs[i+1]-runStubs[i] <= 1)
		    runTotalAccumulation += runCharges[i+1] - runCharges[i];
		}
	    }


	  cout.width(18);
	  cout << pRunNumber;
	  cout.width(18);
	  cout << runTotalAccumulation;
	  cout.width(18);
	  cout << lastEvent << endl;


	  ofstream runAccumulationOut;
	  runAccumulationOut.open(Form("/volatile/clas/clas12/dmriser/analysis/e1f_analysis/fca/%d.fca",pRunNumber),ios::trunc);
	  runAccumulationOut << runTotalAccumulation << " " << lastEvent;
	  runAccumulationOut.close();
	  

	  ofstream e1f_bad_out;
	  e1f_bad_out.open(Form("/volatile/clas/clas12/dmriser/analysis/e1f_analysis/badEvents/%d.txt",pRunNumber),ios::trunc);

	  for (int i=0; i< badEventStart.size(); i++) e1f_bad_out << badEventStart[i] << " " << badEventEnd[i] << endl;

	  e1f_bad_out.close();

	  badEventStart.clear();
	  badEventEnd.clear();

	  runStubs.clear();
	  runEntries.clear();
	  runCharges.clear();
	  iRunEntry = 1;
	} // end conditional on new runfile
      

      // make sure file exists 
      if (fReader.HasFile())
	{

	  int runStub    = atoi(ifTitle.substr(13,2).c_str());
	  nProcessedFiles++;

	  //	  fReader.print();

	  for (int ien = 0; ien < fReader.numberOfEntries()-1; ien++)
	    {
	      // because they are used so much, more efficient to 
	      // calculate once
	      int c = fReader.cdiff(ien, ien+1);
	      int e = fReader.ediff(ien, ien+1);

	      // throw exception for 0 charge 
	      if (c == 0 && fReader.event(ien) > 0 && fReader.event(ien+1) > 0)
		{
		  //		  fReader.print(ien);

		  if (!BAD)
		    {
		      badEventStart.push_back(fReader.event(ien));
		    }

		  BAD = true;

		  nEventsInException += e; 
		  nExceptions++;

		  h1_ediff_cdiff0->Fill(e);
		}
	      
	      else if (c != 0 && fReader.event(ien) > 0 && fReader.event(ien+1) > 0)
		{
		  if (BAD)
		    {
		      badEventEnd.push_back(fReader.event(ien));
		    }

		  BAD = false;
		  nEventsInNormal += e;
		  nNormal++;
		}
	      
	      h1_ratio->Fill((Float_t) e/c);
	      h1_cdiff->Fill(c);
	      h1_ediff->Fill(e);
	      h2_diff->Fill(c,e);

	      runStubs.push_back(runStub);
	      runEntries.push_back(iRunEntry);
	      runCharges.push_back(fReader.charge(ien));
	      
	      lastEvent = fReader.event(ien+1);

	      //  for TGraphs
	      entries.push_back(iEntry);
	      cdiffs.push_back(c);
	      ediffs.push_back(e);
	      iEntry++;
	      iRunEntry++;
	    }  // end loop over entries for a given file
	  
	  if (fReader.numberOfEntries() > 0)
	    {
	      // here we write the file accumulation into the database 
	      int fileAccumulation = 0;
	      for (int i=0; i<fReader.numberOfEntries()-1; i++) fileAccumulation += fReader.cdiff(i,i+1);
	      ofstream accumulationOut;
	      accumulationOut.open(Form("/volatile/clas/clas12/dmriser/analysis/e1f_analysis/fca/%s.fca",ifTitle.c_str()),ios::trunc);
	      accumulationOut << fileAccumulation;
	      accumulationOut.close();
	    }
	}// end conditional on file open 
      
	  ifile++;
	  pRunNumber = cRunNumber;
    }// end loop over files 

  TGraph * g_cdiff = new TGraph(entries.size(),&entries[0],&cdiffs[0]);
  TGraph * g_ediff = new TGraph(entries.size(),&entries[0],&ediffs[0]);

  g_cdiff->Write();
  g_ediff->Write();

  double rException  = (double) nEventsInException/nExceptions;
  double rNormal     = (double) nEventsInNormal/nNormal;

  cout.width(12);
  cout << "exceptions";
  cout.width(16);
  cout << "exeption events";
  cout.width(12);
  cout << "normal";
  cout.width(16);
  cout << "normal events";
  cout.width(12);
  cout << "percentage" << endl;

  cout.width(12);
  cout << nExceptions;

  cout.width(16);
  cout << nEventsInException;

  cout.width(12);
  cout << nNormal;

  cout.width(16);
  cout << nEventsInNormal;

  cout.width(12);
  cout << rException/rNormal *100 << endl;

  cout << " looped on " << nProcessedFiles << " of " << argv[1] << " files" << endl;

  c1->Print("accumulation.pdf");
  c1->Print("accumulation.pdf]");

  runs.close();

  outfile.Write();
  outfile.Close();

  return 0;
}
