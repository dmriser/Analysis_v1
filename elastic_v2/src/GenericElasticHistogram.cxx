#include <iostream>
using namespace std;

#include "TFile.h"
#include "TH2.h"

#include "../../analysisClasses/Bins.h"

#ifndef generic_elastic_histogram_cxx
#define generic_elastic_histogram_cxx

GenericElasticHistogram::GenericElasticHistogram(std::string fname, std::string histTitle, BinStructure thetaBins, BinStructure phiBins, bool recalc = true) : histogramTitle(histTitle), filename(fname) {
  
  if ( recalc ){
    Init();
  }

  else {
    Load();
  }

}

GenericElasticHistogram::~GenericElasticHistogram(){
  
}

void GenericElasticHistogram::Init(){
  
}

void GenericElasticHistogram::Save(){

}

void GenericElasticHistogram::Load(){

}

#endif
