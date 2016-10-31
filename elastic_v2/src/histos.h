#ifndef elast_histos_h
#define elast_histos_h

// c++ includes 
#include <vector>

// my includes 
#include "../../analysisClasses/h22Event.h"
#include "../../analysisClasses/Bins.h"
 
// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"

class elast_histos
{
 public:

  // Default Const/Dest
  elast_histos();
  ~elast_histos();

 private:
  BinStructure * thetaBins;
  BinStructure * phiBins;
  BinStructure * momBins; 
  double w_cut; 

  // initialize or load
 public:
  void init(); /** Here we pass the binning scheme for theta and relative phi */
  //  void load(string);
  void fill(int, int, h22Event, ElasticEvent);
  void fill_gen(int, int, h22Event, ElasticEvent);
  void write_and_close(string);
  void draw(std::string, bool);
  void draw_w(std::string);
  void set_bins(BinStructure*, BinStructure*, BinStructure*);
  void set_w_cut(double x) { w_cut = x; }
  // ------------------------------------------------
  // Histograms contain data and monte carlo. 
  // ------------------------------------------------

  // 1D
  TH1F * h1_W[2][7];
  TH1F * h1_QQ[2][7];
  TH1F * h1_fc_mon;
  TH1F * h1_rc;
  TH1F * h1_xs_model[2];

  // Graphs 
  TGraphErrors * g_mp_from_w[7];
  TGraphErrors * g_dmp_from_w[7];

  // Contain binned 
  std::vector<std::vector<TH1D*> > h1_xs_by_phi;
  std::vector<std::vector<TH1D*> > h1_xs_ratio_by_phi;
  std::vector<std::vector<TH1D*> > h1_rxs_by_phi;
  std::vector<std::vector<TH1D*> > h1_rxs_ratio_by_phi;
  std::vector<std::vector<TH1D*> > h1_acc_by_phi;
  std::vector<std::vector<TH1D*> > h1_gen_by_phi;
  std::vector<std::vector<TH1D*> > h1_rec_by_phi;
  std::vector<std::vector<TH1D*> > h1_hits_by_phi;
  std::vector<std::vector<TH1D*> > h1_w_by_mom;


};

#endif
