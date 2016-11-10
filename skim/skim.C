/*

           ntuple skimmer
   David Riser, University of Connecticut

   April 2, 2016

 */

// C++ Libraries
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>

using namespace std;

// My Libraries
#include "../analysisClasses/h22Event.h"
#include "../analysisClasses/h22Reader.h"
#include "../analysisClasses/Bins.h"
#include "../analysisClasses/FaradayReader.h"
#include "../analysisClasses/ParticleFilter.h"
#include "../analysisClasses/ElasticEvent.h"

// CERN Root Libraries
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TString.h"
#include "TRegexp.h"

int main (int argc, char * argv[])
{
  int GSIM = 0;                          // data (0 - false), gsim (1 - true)

  // get number of files from command line
  if (argc < 2)
    {
      cout << " expected file list as option " << endl;
      exit(0);
    }

  // setup file reader and add files
  h22Reader reader(GSIM);
  ParticleFilter filter;

  for (int ifile = 1; ifile < argc; ifile++) reader.AddFile(argv[ifile]);
  reader.Init(); // set branch addresses

  int nEvents = reader.GetEntries();
  //  int runno = atoi(reader.GetFilenameChunk(68,6).c_str());

  TString filename = reader.fchain->GetFile()->GetName();
  TRegexp runno_regex("[1-9][0-9][0-9][0-9][0-9]");
  TString srunno = filename( runno_regex );
  int runno = srunno.Atoi();

  // -------------- event quality control --------------------
  bool eventStatus = true;
  int iBadWindow = 0;
  ifstream eventQCFile(Form("/volatile/clas/clas12/dmriser/analysis/e1f_analysis/badEvents/%d.txt",runno));

  vector<int> badEventStart;
  vector<int> badEventEnd;

  string line, buffer;

  while(getline(eventQCFile,line))
    {
      istringstream iss(line);
      int istring = 0;

      while(iss)
	{
	  iss >> buffer;
	  if(istring == 0)
	    {
	      int s = atoi(buffer.c_str());
	      badEventStart.push_back(s);
	    }
	  else if (istring == 1)
	    {
	      int e = atoi(buffer.c_str());
	      badEventEnd.push_back(e);
	    }

	  istring++;
	}
    }
  
  eventQCFile.close();
  int NWINDOWS = badEventStart.size();
  
  // print for debuggin 
  cout << " > found " << NWINDOWS << " bad windows " << endl;
  

  // ---------------------------------------------------------


  // -------------- NEW TREE FOR ENTRIES --------------------
  TFile * newNtuple = TFile::Open(Form("/volatile/clas/clas12/dmriser/analysis/e1f_analysis/debug_jobs2/skim_e2n.%d.root",runno),"recreate");
  TChain * newChain = (TChain*) reader.fchain->CloneTree(0);
  TTree * newTree   = newChain->GetTree();
  // -------------------------------------------------------
  
  // loop over events
  for (int iEvent = 0; iEvent < nEvents; iEvent++)
    {
      reader.GetEntry(iEvent);

      // store each event locally for the duration of the loop
      h22Event event = reader.GetEvent();

      // ---------- set event status for QC ------------
      
      if (iBadWindow < NWINDOWS)
	{
	  // catch times when we passed the current window
	  if (event.evntid > badEventEnd[iBadWindow]) iBadWindow++;

	  // we're in a bad event region, dQ = 0
	  if ((event.evntid > badEventStart[iBadWindow]) && (event.evntid < badEventEnd[iBadWindow]))
	    {
	      eventStatus = false;
	    }
	  else 
	    {
	      eventStatus = true;
	    }
	}
      
      // ----------------------------------------------

      filter.loadEvent(event,GSIM,runno);
      int e_index = filter.getIndexByPID(11);
      if (e_index > -123)
	{
	  //	  TLorentzVector recElectron(event.cx[e_index]*event.p[e_index], 
	  //				     event.cy[e_index]*event.p[e_index],
	  //				     event.cz[e_index]*event.p[e_index],
	  //				     event.p[e_index]);

	  //	  ElasticEvent recEvent(recElectron);

      int numberOfNeutrals = 0;
      for (int ipart=1; ipart<event.gpart; ipart++){
	if (event.q[ipart] == 0) numberOfNeutrals++;
      }
      
	  // define the skim options here
	  //	  if (recEvent.getW() > 1.80 && recEvent.getQQ() > 0.95  && eventStatus) newTree->Fill();
	  //	  if (eventStatus) newTree->Fill();
	  if (eventStatus && numberOfNeutrals > 1) newTree->Fill();
	}

    }  // end loop over events

  // save and exit new file
  newTree->AutoSave();
  newNtuple->Close();

  return 0;
}
