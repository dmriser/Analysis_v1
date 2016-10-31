///////////////////////////////////////////////
/*

  David Riser, University of Connecticut 
  August 4, 2016 

  elastic_v2.0.h -> Define Class for 
  elastic second iteration. 

  
*/
///////////////////////////////////////////////

#ifndef ELASTIC_V2_H
#define ELASTIC_V2_H


// c++ Includes 
#include <iostream>
#include <vector>

using namespace std; // using in header because class is also developed here

// My Includes 
#include "../../analysisClasses/h22Event.h"
#include "../../analysisClasses/h22Reader.h"
#include "../../analysisClasses/ParticleFilter.h"

class ElasticAnalysis
{
  // Default Constr. Destr.
 public: 
  ElasticAnalysis();
  ~ElasticAnalysis();

  // Data types 
 private: 
  h22Reader * fReader[2];  //! Two readers, one for GSIM and one for data. 
  
  // User Methods
 public: 
  void init();         //! Initialize the readers
  void add_files(vector<string>,int); //! Pass type, adds list
  void loop(int, int);      //! Pass type, adds list 
  void close();        //! Save everything and close

};

// Definitions
ElasticAnalysis::ElasticAnalysis()
{
  for (int t=0; t<2; t++) fReader[t] = new h22Reader(t);
}

ElasticAnalysis::~ElasticAnalysis()
{
  // Destroy Something 
}

void ElasticAnalysis::init()
{
  for (int t=0; t<2; t++) fReader[t]->Init();
}

void ElasticAnalysis::add_files(std::vector<string> files, int index)
{

  for (int f=0; f<files.size(); f++) { 
    fReader[index]->AddFile(files[f]); 
    cout << " Added: " << files[f] << "to fReader " << index << endl;
  }

}

void ElasticAnalysis::close()
{
  // Do something here. 
}

void ElasticAnalysis::loop(int n_events, int index)
{
  // Determine true number of events 
  if ( fReader[index]->GetEntries() < n_events ) n_events = fReader[index]->GetEntries();

  for (int iev=0; iev<n_events; iev++)
    {
      fReader[index]->GetEntry(iev);
      h22Event event = fReader[index]->GetEvent();

    }
}

