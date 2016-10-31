#ifndef elast_histos_h
#define elast_histos_h

// my includes 
#include "../../analysisClasses/h22Event.h"
#include "../../analysisClasses/Bins.h"
 
// ROOT includes
#include "TH1.h"
#include "TH2.h"

class elast_histos
{
 public:

  // Default Const/Dest
  elast_histos();
  ~elast_histos();

  // initialize or load
  void init(BinStructure, BinStructure); /** Here we pass the binning scheme for theta and relative phi */
  void load(string);
  void write_and_close(string);

  // the histograms[ELAST MC, ELAST MC + RAD, DATA] [All or Sector 1-6]

  // 1D
  TH1F * h1_W[3][7];
  TH1F * h1_QQ[3][7];
  TH1F * h1_x[3][7];
  TH1F * h1_EP_angle[3][7];
  TH1F * h1_MM2[3][7];
  TH1F * h1_rec_theta[3][7];
  TH1F * h1_gen_theta[3][7];
  TH1F * h1_ME[3][7];
  TH1F * h1_fc_mon;
  TH1F * h1_rc;
  TH1F * h1_acc[7];

  // (norad, rad) (gen, rec, data), (all,s1-s6)
  TH1F * h1_norm_theta[2][3][7];
  TH1D * h1_xs[7];
  TH1D * h1_xs_raw[7];

  // with and without radiation 
  TH1F * h1_xs_model[2];
  TH1F * h1_xs_ratio[2][7];

  // 2D
  TH2F * h2_rec[3][7];
  TH2F * h2_gen[3][7];
  TH2F * h2_acc[3][7];
  TH2F * h2_W_theta[3][7];

  TH2F * h2_xs[7];
  TH2F * h2_scale;

  // Histograms for profiling BH process in MC + RAD
  // [before, after][all, sectors 1-6]
  TH1F * h1_e_E;
  TH1F * h1_gamma_E;
  TH1F * h1_gamma_theta;
  TH1F * h1_gamma_phistar;

  TH1F * h1_proton_E;
  TH1F * h1_proton_theta;
  TH1F * h1_proton_phistar;


};

#endif
