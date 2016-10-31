#ifndef PARS_H
#define PARS_H

#include <iostream>
using namespace std; 

class pars
{

 public: 
  pars();
  ~pars();

  // Elastic Event W Parameters 
  double mu[3][7];
  double sigma[3][7];
  double N_SIGMA[3][7];

  // Methods 
  void read(string);
  void write(string);

  double lower_bound(int, int);
  double upper_bound(int, int);

};

#endif
