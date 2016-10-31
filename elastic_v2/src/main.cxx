/////////////////////////////////////////
/*
  David Riser, University of Connecticut 
  August 4, 2016 

  main.cxx -> Drives Elastic Analysis 
  

*/
/////////////////////////////////////////

// c++ includes 
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std; 

// My Includes 
#include "elastic_v2.0.h"

int main(int argc, char * argv[])
{
  
  // Getting number of events from command line or doing default number.
  int n_files = 10;
  if (argc > 1) { n_files = atoi(argv[1]); }

  // Getting list of files from database.
  ifstream data("data.dat");   vector<string> data_files;
  ifstream mc("mc.pass11.dat");       vector<string> mc_files;
  string file; 
  int ifile = 0;
  while(!data.eof() && ifile < n_files) {   data >> file; data_files.push_back(file); ifile++; } ifile = 0;
  while(!mc.eof()   && ifile < n_files) {   mc >> file;   mc_files.push_back(file); ifile++;   }
  data.close(); 
  mc.close();

  // Start doing analyis
  ElasticAnalysis elastic;

  elastic.add_files(data_files, 0);
  elastic.add_files(mc_files,   1);
  elastic.set_name("out");
  elastic.init(); // Set branch adresses, initialize histograms 

  elastic.set_mom_corr_status( true );
  elastic.loop(0); // Data 
  elastic.loop(1); // Monte Carlo Events 

  elastic.get_charge(data_files);
  elastic.get_model();
  elastic.calc_xs();
  elastic.fit_w();

  elastic.set_gif_status( false );
  elastic.draw("out");
  elastic.close();

  return 0;
}
