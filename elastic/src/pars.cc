#ifndef PARS_CXX
#define PARS_CXX

// My Includes
#include "pars.h"


// C++ Includes 
#include <fstream>
#include <iostream>

using namespace std;

pars::pars()
{
  mu      = {};
  sigma   = {};
  N_SIGMA = {};
}

pars::~pars()
{

}

void pars::read(string _file)
{

}

void pars::write(string _file)
{
  
  ofstream f(file.c_str(),std::ios::trunc);

  // WRITE MU 
  f << "# MU: " << endl;
  for (int t=0; t<3; t++)
    {
      for (int s=0; s<7; s++)
	{
	  f.width(12);
	  f << mu[t][s];
	}
      f << endl;
    }

  // WRITE SIGMA 
  f << "# SIGMA: " << endl;
  for (int t=0; t<3; t++)
    {
      for (int s=0; s<7; s++)
	{
	  f.width(12);
	  f << sigma[t][s];
	}
      f << endl;
    }


  // WRITE N_SIGMA 
  f << "# N_SIGMA: " << endl;
  for (int t=0; t<3; t++)
    {
      for (int s=0; s<7; s++)
	{
	  f.width(12);
	  f << N_SIGMA[t][s];
	}
      f << endl;
    }

}



#endif PARS_CXX
