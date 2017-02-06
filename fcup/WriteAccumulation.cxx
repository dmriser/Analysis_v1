/*

  David Riser 
  University of Connecticut 

  Mar 15, 2016

  WriteAccumulation.cxx

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
#include "WriteAccumulation.h"

// CERN root libs
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"


class Writer{

public:
  Writer();
  ~Writer();

protected:

  int nNormal, nEventsInNormal, nExceptions, nEventsInException, 
    lastEvent, ifile, nProcessedFiles, pRunNumber, cRunNumber, iEntry, iRunEntry; 

  vector<int> runStubs;
  vector<int> runEntries;
  vector<int> runCharges;
  vector<int> badEventStart;
  vector<int> badEventEnd;
  vector<int> cdiffs;
  vector<int> ediffs;
  vector<int> entries;

  vector<TH1D*> runHistos;

  bool BAD;

  string ifName, ifTitle;

  TFile *outfile;

public:
  void Initialize(); 
  void Execute(int numberOfFiles); 
  void Save(string outputFile);
  void PerformEndOfRunAction();
  void ProcessFile(FaradayReader fReader);
  void EndJob();
  bool FileIsNewRun();
};

Writer::Writer(){

}

Writer::~Writer(){

}

void Writer::Initialize(){
  nNormal            = 0;
  nEventsInNormal    = 0;
  nExceptions        = 0;
  nEventsInException = 0;
  lastEvent          = 0;
  ifile              = 0;
  nProcessedFiles    = 0;
  pRunNumber         = 0;
  cRunNumber         = 0;
  iEntry             = 0;
  iRunEntry          = 1;


  BAD = false;
  outfile = new TFile("WriteAccumulation_Refactoring.root","recreate");
}

bool Writer::FileIsNewRun(){
  return (cRunNumber != pRunNumber && pRunNumber != 0);
}

void Writer::Execute(int numberOfFiles){

  ifstream runs("allfiles.txt"); 
  
  while (getline(runs, ifName) && ifile < numberOfFiles){
    // Very bad practice to use substring by position, 
    // should be done with TRegex or other regex. 
    ifTitle        = ifName.substr(52,15);
    int cRunNumber = atoi(ifName.substr(57,6).c_str());
    FaradayReader fReader(ifTitle);

    if ( FileIsNewRun() ){      
      PerformEndOfRunAction(); 
    } 
    
    if (fReader.HasFile()){ 
      ProcessFile(fReader); 
    }
  
    ifile++;
    pRunNumber = cRunNumber;
  }// end conditional on file open 

  
  runs.close();
}

void Writer::PerformEndOfRunAction(){
  cout << "[Writer::PerformEndOfRunAction] Performing end of run action for file " << ifTitle << " with Run Number " << pRunNumber << endl; 

      // this need to be done more carefully so to cut out gaps 
      // such as .A02 -> .A04 skipping .A03 
      // I am calling this a 'jump error'
      int runTotalAccumulation = 0;
      if (runCharges.size() > 0){
	for (int i=0; i<runCharges.size()-1; i++){
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
      ofstream e1f_bad_out;
      
      runAccumulationOut.open(Form("testOutput/accumulation/%d.fca",pRunNumber),ios::trunc);
      runAccumulationOut << runTotalAccumulation << " " << lastEvent;
      runAccumulationOut.close();

      e1f_bad_out.open(Form("testOutput/badEvents/%d.txt",pRunNumber),ios::trunc);
      
      for (int i=0; i< badEventStart.size(); i++) e1f_bad_out << badEventStart[i] << " " << badEventEnd[i] << endl;

      int numberReadings = runCharges.size();

      TH1D *runHisto = new TH1D(Form("fcup_%d",pRunNumber),"",numberReadings,1,numberReadings); 

      for(int i=0; i<numberReadings-1; i++){
	double chargeDiff = runCharges[i+1]-runCharges[i]; 
	runHisto->SetBinContent(i+1,chargeDiff); 
      }
      runHistos.push_back(runHisto); 

      e1f_bad_out.close();
      badEventStart.clear();
      badEventEnd.clear();
      runStubs.clear();
      runEntries.clear();
      runCharges.clear();
      iRunEntry = 1;
}

void Writer::ProcessFile(FaradayReader fReader){

  int runStub = atoi(ifTitle.substr(13,2).c_str());
  nProcessedFiles++;
  
  cout << "Processing file: " << ifTitle << endl; 

  for (int ien = 0; ien < fReader.numberOfEntries()-1; ien++){
    int c = fReader.cdiff(ien, ien+1);
    int e = fReader.ediff(ien, ien+1);
    
    // throw exception for 0 charge 
    if (c == 0 && fReader.event(ien) > 0 && fReader.event(ien+1) > 0){
      if (!BAD){
	badEventStart.push_back(fReader.event(ien));
      }
      
      BAD = true;
      
      nEventsInException += e; 
      nExceptions++;
    }
    
    else if (c != 0 && fReader.event(ien) > 0 && fReader.event(ien+1) > 0){
      if (BAD){
	badEventEnd.push_back(fReader.event(ien));
      }
      
      BAD = false;
      nEventsInNormal += e;
      nNormal++;
    }
    
    runStubs  .push_back(runStub);
    runEntries.push_back(iRunEntry);
    runCharges.push_back(fReader.charge(ien));
    
    lastEvent = fReader.event(ien+1);
    
    entries.push_back(iEntry);
    cdiffs .push_back(c);
    ediffs .push_back(e);
    iEntry++;
    iRunEntry++;
  }    
}

void Writer::EndJob(){
  
  
  double rException  = (double) nEventsInException/nExceptions;
  double rNormal     = (double) nEventsInNormal/nNormal;

  cout << "------------------------------ EndJob -----------------------------" << endl; 
  cout << "exceptions ";
  cout << nExceptions << endl;
  cout << "exeption events ";
  cout << nEventsInException << endl;
  cout << "normal ";
  cout << nNormal << endl;
  cout << "normal events ";
  cout << nEventsInNormal << endl;
  cout << "percentage ";
  cout << rException/rNormal *100 << endl;
  cout << "[Writer::EndJob] Processed files: " << nProcessedFiles << endl; 
  cout << "------------------------------ EndJob -----------------------------" << endl; 


  for(int h=0; h<runHistos.size(); ++h){
    runHistos[h]->Write(); 
  }
  
  outfile->Write();
  outfile->Close();
  
}

int main (int argc, char * argv[]){

  if (argc > 1){
    Writer writer; 
    writer.Initialize(); 
    writer.Execute(atoi(argv[1])); 
    writer.EndJob(); 
  } else {
    cout << "[WriteAccumulation] Please pass in the number of files to process. " << endl; 
  }

  return 0;
}
