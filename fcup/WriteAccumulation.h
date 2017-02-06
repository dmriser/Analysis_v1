#ifndef READ_FARADAY_WITH_ENTRY_H
#define READ_FARADAY_WITH_ENTRY_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "TString.h"

using namespace std; 

class FaradayReader{

 private: 
  string title;
  bool hasfile;

  vector<int> charges;
  vector<int> events;

 public: 
  FaradayReader(string);
  virtual ~FaradayReader();

  bool HasFile(){return hasfile;};

  int charge(int c){return charges[c];}
  int cdiff(int a, int b){return charges[b]-charges[a];};
  int event(int e){return events[e];}
  int ediff(int a, int b){return events[b]-events[a];};
  int numberOfEntries(){return events.size();}

  string getTitle(){return title;}

  void print();
  void print(int);

};

FaradayReader::FaradayReader(string f){

  ifstream inputFile(Form("/volatile/clas12/dmriser/analysis/e1f_analysis/fcup/%s.fcup",f.c_str()));

  if (!inputFile.is_open()){ cout << " Error opening file: " << f << endl; hasfile = false; return;}
  else hasfile = true;

  string line;

  while(getline(inputFile, line)){
      int ien = 0;
      istringstream iss(line); 
      
      while(iss){
	  string buffer; 
	  iss >> buffer; 

	  if (ien == 0) charges.push_back(atoi(buffer.c_str()));
	  else if (ien == 1) events.push_back(atoi(buffer.c_str()));

	  ien++;
	}
    }

  title = f;

}

FaradayReader::~FaradayReader()
{}

void FaradayReader::print(){
  cout << title << endl;

  for (int ien = 0; ien < events.size(); ien++)
    {
      cout.width(12);
      cout << ien;

      cout.width(12);
      cout << charges[ien];

      cout.width(12);
      cout << events[ien];

      cout.width(12);
      cout << charges[ien+1] - charges[ien];

      cout.width(12);
      cout << events[ien +1] - events[ien] << endl;
      
    }
  cout << endl;
}

void FaradayReader::print(int ien)
{
  cout.width(20);
  cout << title;
  
  cout.width(12);
  cout << ien;
  
  cout.width(12);
  cout << charges[ien];
  
  cout.width(12);
  cout << events[ien];
  
  cout.width(12);
  cout << charges[ien+1] - charges[ien];
  
  cout.width(12);
  cout << events[ien +1] - events[ien] << endl;
  
}


#endif
