#ifndef BINS_H
#define BINS_H

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

/**
 * BinStructure is a class which holds your custom bins
 * and returns important properties. h
 */

class BinStructure 
{
 public: 
  /**
   * Default constructor accepts number of bins, min value, and max value
   */
  BinStructure(int,float,float);
  ~BinStructure();
  
  int number(){return number_of_bins;}; /**< returns the number of bins */
  float max(){return bin_max;};/**< returns the bin with lowest value*/
  float min(){return bin_min;};/**< returns the bin with highest value*/
  float width(){return bin_size;};/**< returns the size of the bins */

  string getName(){return bin_name;};/**< returns the name of this custom bin structure */
  void setName(string name){bin_name = name;};/**< sets the name of this custom bin structure */

  float getBinContent(int);/**< returns the value of the ith bin getBinContent(i) */
  int getBin(float); /**< Finds which bin the passed value belongs in, if value > max or value < min, it returns -999 */

  void print();

  vector<double> getBins();

 private:
  int number_of_bins;
  float bin_size;
  float bin_max;
  float bin_min;
  string bin_name;
  vector<double> bins;

};

BinStructure::BinStructure(int nbins, float bmin, float bmax)
{
  number_of_bins = nbins;
  bin_min = bmin;
  bin_max = bmax;
  bin_size = (bin_max - bin_min)/(number_of_bins-1);

  for (int ibin = 0; ibin < number_of_bins; ibin++) bins.push_back(bin_min + ibin*bin_size);
}

BinStructure::~BinStructure()
{}

float BinStructure::getBinContent(int ibin)
{
  return bin_min + ibin*bin_size;
}

int BinStructure::getBin(float bin_content)
{
  int bin = floor((bin_content-bin_min)/bin_size);

  if (bin >= 0 && bin < number_of_bins ) return bin;

  return -999;
}

vector<double> BinStructure::getBins()
{return bins;}

void BinStructure::print()
{
  // print binning summary 
  cout.width(18);
  cout << bin_name;

  cout.width(12);
  cout << bin_min;

  cout.width(12);
  cout << bin_max;

  cout.width(12);
  cout << bin_size;

  cout.width(12);
  cout << number_of_bins << endl;

}

#endif
