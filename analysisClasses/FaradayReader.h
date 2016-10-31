#ifndef FARADAY_READER_H
#define FARADAY_READER_H


// 
// David Riser 
// Feb 1, 2016
// 
// class to load Faraday cup information from 
// files saved locally.  
//

// cern root includes 
#include "TString.h"

// c++ includes
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
/**
 * FaradayCup is a class that is designed to integrate farday cup information stored 
 * in text files into other programs in an easy way.
 */
class FaradayCup
{
 public: 
  FaradayCup(int); /**< the default constructor acceps the run, 38222 for example */
  FaradayCup(string);/**< opens just one file for ex: clas_038222.A01 */
  ~FaradayCup();

  bool hasFile(){return hasfile;};

  double getFcupCharge(){return fcupAccum;}; /**< returns the accumulation directly from the file after fixing scale from DAQ*/
  double getIntegratedLuminosity(){return integratedLuminosity;}; /**< returns the integreated luminosity, a constant multiple of the accumulation */

  int getNumberOfElectrons(){return nElectrons;}; /**< doesn't work at all. */
  int getLastEvent(){return lastEvent;};

 private:
  double nElectrons;
  int lastEvent;

  double fcupAccum;
  double integratedLuminosity;

  bool hasfile;
};
 
#endif

FaradayCup::FaradayCup(int run)
{

  string value;
  ifstream file( Form("/volatile/clas/clas12/dmriser/analysis/e1f_analysis/fca/%d.fca",run) );

  int fc        = 0;
  int lastEvent = 0;

  if (!file.is_open())
    {
      hasfile = false;
      return;
    }

  hasfile = true;

  file >> fc >> lastEvent;
  //  cout << " raw fc from FaradayReader " << fc << " in uC " << (double) fc/9624000.00 << endl;
  // DAQ Scaling Factor 
  fcupAccum = (double) fc/9624000.00; // uC 
  nElectrons = fcupAccum/(1.6*10e-13);
  integratedLuminosity = nElectrons * 5.00 * 4.22e23; // cm^-2


}

// pass in clas_038222.A14 format
FaradayCup::FaradayCup(string _file)
{

  string value;
 
  ifstream file( Form("/volatile/clas/clas12/dmriser/analysis/e1f_analysis/fca/%s.fca",_file.c_str()) );

  int fc        = 0;
  int lastEvent = 0;

  if (!file.is_open())
    {
      hasfile = false;
      return;
    }

  hasfile = true;

  file >> fc >> lastEvent;

  //  cout << " raw fc from FaradayReader " << fc << " in uC " << (double) fc/9624000.00 << endl;

  fcupAccum = (double) fc/9624000.00;
  integratedLuminosity = fc * 9.128*10e24;
  nElectrons = (int) fcupAccum/(1.6*10e-13);

}

FaradayCup::~FaradayCup()
{}


