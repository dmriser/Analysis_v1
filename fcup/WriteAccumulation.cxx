/*

  David Riser 
  University of Connecticut 

  Feb 6, 2017

  WriteAccumulation.cxx

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
#include "TH1.h"

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
  vector<int> runEntriesGood;
  vector<int> runChargesGood;
  vector<int> badEventStart;
  vector<int> badEventEnd;
  vector<int> cdiffs;
  vector<int> eventStatus;
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
  long int CalculateTotalCharge();
  bool FileIsNewRun();
  bool EntryIsGood(int entry);
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

  BAD     = false;
  outfile = new TFile("WriteAccumulationTest.root","recreate");
}
  
void Writer::Execute(int numberOfFiles){

  ifstream runs("allfiles.txt"); 
  
  while (getline(runs, ifName) && ifile < numberOfFiles){
    // Very bad practice to use substring by position, 
    // should be done with TRegex or other regex. 
    ifTitle    = ifName.substr(52,15);
    cRunNumber = atoi(ifName.substr(57,6).c_str());
    FaradayReader fReader(ifTitle);

    if ( FileIsNewRun() ){      
      PerformEndOfRunAction(); 
    } 
    
    if (fReader.HasFile()){ 
      ProcessFile(fReader); 
    }
  }
  
  runs.close();
}

bool Writer::FileIsNewRun(){
  return (cRunNumber != pRunNumber && pRunNumber != 0);
}

void Writer::PerformEndOfRunAction(){
  //  cout << "[Writer::PerformEndOfRunAction] Performing end of run action for file " << ifTitle << " with Run Number " << pRunNumber << endl; 

  // this need to be done more carefully so to cut out gaps 
  // such as .A02 -> .A04 skipping .A03 
  // I am calling this a 'jump error'
  long int runTotalAccumulation = 0;
  if (runCharges.size() > 0){
    runTotalAccumulation = CalculateTotalCharge();
  }
  
  ofstream runAccumulationOut;
  ofstream e1f_bad_out;
  
  runAccumulationOut.open(Form("testOutput/accumulation/%d.fca",pRunNumber),ios::trunc);
  runAccumulationOut << runTotalAccumulation << " " << lastEvent;
  runAccumulationOut.close();
  
  e1f_bad_out.open(Form("testOutput/badEvents/%d.txt",pRunNumber),ios::trunc);
  
  for (int i=0; i< badEventStart.size(); i++) e1f_bad_out << badEventStart[i] << " " << badEventEnd[i] << endl;
  
  int numberReadings     = runCharges.size();
  int numberReadingsGood = runChargesGood.size();
  
  if (numberReadings > 1){
    TH1D *runHisto       = new TH1D(Form("chargeDiff_%d",pRunNumber),Form("%d",pRunNumber),numberReadings-1,1,numberReadings-1); 
    TH1D *entryHisto     = new TH1D(Form("entryDiff_%d",pRunNumber),Form("%d",pRunNumber),numberReadings-1,1,numberReadings-1); 
    TH1D *runHistoGood   = new TH1D(Form("chargeDiffGood_%d",pRunNumber),Form("%d",pRunNumber),numberReadingsGood-1,1,numberReadingsGood-1); 
    TH1D *entryHistoGood = new TH1D(Form("entryDiffGood_%d",pRunNumber),Form("%d",pRunNumber),numberReadingsGood-1,1,numberReadingsGood-1); 
    
    for(int i=0; i<numberReadings-1; i++){
      double chargeDiff = (runCharges[i+1]-runCharges[i])/9624000.0; 
      double entryDiff  = runEntries[i+1]-runEntries[i];
      
      runHisto       ->SetBinContent(i+1,chargeDiff); 
      entryHisto     ->SetBinContent(i+1,entryDiff); 

      cout.width(16); cout << runStubs[i]; 
      cout.width(16); cout << runStubs[i+1]; 
      cout.width(16); cout << runCharges[i]/9624000.0; 
      cout.width(16); cout << runCharges[i+1]/9624000.0; 
      cout.width(16); cout << runEntries[i]; 
      cout.width(16); cout << runEntries[i+1]; 
      cout.width(16); cout << (runCharges[i+1]-runCharges[i])/9624000.0 << endl; 
      

    }
    
    for(int i=0; i<numberReadings-1; i++){
      double chargeDiffGood = (runChargesGood[i+1]-runChargesGood[i])/9624000.0; 
      double entryDiffGood  = runEntriesGood[i+1]-runEntriesGood[i];
      
      if(EntryIsGood(i)){
	runHistoGood   ->SetBinContent(i+1,chargeDiffGood); 
	entryHistoGood ->SetBinContent(i+1,entryDiffGood); 
      }
    }
    
    runHistos.push_back(runHisto); 
    runHistos.push_back(entryHisto); 
    runHistos.push_back(runHistoGood); 
    runHistos.push_back(entryHistoGood); 
  }
  
  e1f_bad_out   .close();
  badEventStart .clear();
  badEventEnd   .clear();
  runStubs      .clear();
  runEntries    .clear();
  runCharges    .clear();
  runEntriesGood.clear();
  runChargesGood.clear();
  eventStatus   .clear();
  iRunEntry = 1;
}

long int Writer::CalculateTotalCharge(){
  long int totalCharge = 0;
  for(int entry=0; entry<runCharges.size()-1; ++entry){
    if (EntryIsGood(entry)) { 
      totalCharge += runCharges[entry+1]-runCharges[entry];
    }
  }
  return totalCharge; 
}

bool Writer::EntryIsGood(int entry){
  //  return (eventStatus[entry] && eventStatus[entry+1] && (runStubs[entry]-runStubs[entry-1]) == 0);
  return (eventStatus[entry] && eventStatus[entry+1] && (runStubs[entry]-runStubs[entry-1]) < 2);
}

void Writer::ProcessFile(FaradayReader fReader){
  int runStub = atoi(ifTitle.substr(13,2).c_str());
  nProcessedFiles++;
  
  cout << Form("Processing file (%d): ",ifile) << ifTitle << endl;

  for (int ien = 0; ien < fReader.numberOfEntries()-1; ien++){
    int c = fReader.cdiff(ien, ien+1);
    int e = fReader.ediff(ien, ien+1);
    
    if (c == 0){
      if (!BAD){
	badEventStart.push_back(fReader.event(ien));
      }
      BAD = true;
      nEventsInException += e; 
      nExceptions++;
    }

    else if (c != 0) {
      if (BAD){
	badEventEnd.push_back(fReader.event(ien));
      }
      BAD = false;
      nEventsInNormal += e;
      nNormal++;
    }

    runStubs  .push_back(runStub);
    runEntries.push_back(fReader.event(ien));
    runCharges.push_back(fReader.charge(ien));

    if (!BAD) {
      runEntriesGood.push_back(fReader.event(ien));
      runChargesGood.push_back(fReader.charge(ien));
      eventStatus.push_back(1);

      entries.push_back(iEntry);
      cdiffs .push_back(c);
      ediffs .push_back(e);

    } else {
      eventStatus.push_back(0);
    }

    lastEvent = fReader.event(ien+1);    
    iEntry++;
    iRunEntry++;
  }    
  
  pRunNumber = cRunNumber;
  ifile++;
  
}

void Writer::EndJob(){
  
  
  double rException  = (double) nEventsInException/nExceptions;
  double rNormal     = (double) nEventsInNormal/nNormal;

  cout << endl;
  cout << endl;
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
