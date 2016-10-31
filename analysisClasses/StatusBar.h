#ifndef STATUSBAR_H
#define STATUSBAR_H

// c++ includes 
#include <iostream>
using namespace std;

class StatusBar
{
 public:
  StatusBar(int);
  ~StatusBar();

  int kWidth;

  void print(int, int);  
};

#endif


#ifndef STATUSBAR_CXX
#define STATUSBAR_CXX

StatusBar::StatusBar(int w)
{
  kWidth = w;
}
StatusBar::~StatusBar(){}

void StatusBar::print(int i, int n)
{
  int perc    = (int) 100*i/n+1;
  int step    = (int) (n/kWidth);  
  int current = floor(i/step);

  std::cout << "\r {";

  for (unsigned int c=0; c<kWidth; c++)
    {
      if (c<current) std::cout << "="; 
      else if (c == current ) std::cout << ">";
      else std::cout << " ";
    }

  std::cout << "} " << perc << "%" << flush;

}
#endif
