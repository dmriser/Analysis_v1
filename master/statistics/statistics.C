#define statistics_cxx

#include "statistics.h"
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TTree.h>
#include "TF1.h"
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <cmath>
#include "config.h"

using std::ofstream;

// main
int main(int argc, char *argv[]){

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "================================" << std::endl;
  std::cout << "          statistics.C "          << std::endl;
  std::cout << "================================" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  if(debug){std::cout << "debug mode is active" << std::endl << std::endl;}

  statistics * data = new statistics();

  // load in files from command line
  for(int iarg=1; iarg<argc; iarg++){
    data->AddFile(argv[iarg]);
}

  std::cout << "running " << argc -1 << " files though statistics.C" << std::endl << std::endl;

  data->Init();
  data->processThetaCC();
  data->processECSampling();

  return 0;
}

//
// ==============================================================================
//                SUBROUTINES FOLLOW
// ==============================================================================
//

Float_t statistics::getThetaCC(int ipart){
  Float_t cc_pln[3] = {-0.0007840784063, 0.0, -0.001681461571};
  Float_t d = 1.0;

  Float_t dir[3] = {tl3_cx[ipart], tl3_cy[ipart], tl3_cz[ipart]};

  Float_t cm_factor = 1.0; // this is 10 in the original CLAS_Event.cc script
  Float_t P1[3] = {tl3_x[ipart]/cm_factor, tl3_y[ipart]/cm_factor, tl3_z[ipart]/cm_factor};
  Float_t t = (cc_pln[0]*P1[0] + cc_pln[1]*P1[1] + cc_pln[2]*P1[2] + d)/(cc_pln[0]*dir[0] + cc_pln[1]*dir[1] + cc_pln[2]*dir[2]);

  Float_t CCx = (P1[0] + dir[0]*t)*10;
  Float_t CCy = (P1[1] + dir[1]*t)*10;
  Float_t CCz = (P1[2] + dir[2]*t)*10;

  Float_t thetaCC = atan2(sqrt(CCx*CCx + CCy*CCy), CCz);
  return thetaCC * r2d;
}

// ==============================================================================
// range of atan2() is
// -pi to pi here we have
// 0 to 2pi
// ==============================================================================

Float_t atan3(Float_t y, Float_t x){

  if(atan2(y, x) >= 0)
    {
      return atan2(y, x);
    }
  else if(atan2(y, x) < 0)
    {
      return atan2(y, x) + 2*3.14159265359;
    }
  else
    {
  return -123;
}
}

bool statistics::ecGeometricCut(int ipart){

  Float_t x = ech_x[ipart];
  Float_t y = ech_y[ipart];
  Float_t z = ech_z[ipart];

  float uMax = 405;
  float uMin = 70;
  float vMax = 362;
  float wMax = 395;

  Float_t u, v, w, xi, yi, zi;
  Float_t EC_the = 0.4363323;
  Float_t EC_phi;
  Float_t ylow   = -182.974;
  Float_t yhi    = 189.956;
  Float_t tgrho  = 1.95325;
  Float_t sinrho = 0.8901256;
  Float_t cosrho = 0.455715;
  Float_t rot[3][3];

  Float_t phi = 57.2957795*atan3(y, x);

  phi=phi+30.;
  if (phi>=360.) phi=phi-360.;

  EC_phi = (int)(phi/60.) * 1.0471975;

  rot[0][0] = cos(EC_the)*cos(EC_phi);
  rot[0][1] = -sin(EC_phi);
  rot[0][2] = sin(EC_the)*cos(EC_phi);
  rot[1][0] = cos(EC_the)*sin(EC_phi);
  rot[1][1] = cos(EC_phi);
  rot[1][2] = sin(EC_the)*sin(EC_phi);
  rot[2][0] = -sin(EC_the);
  rot[2][1] = 0.;
  rot[2][2] = cos(EC_the);

  yi=x*rot[0][0]+y*rot[1][0]+z*rot[2][0];
  xi=x*rot[0][1]+y*rot[1][1]+z*rot[2][1];
  zi=x*rot[0][2]+y*rot[1][2]+z*rot[2][2];

  zi=zi-510.32;

  u = (yi-ylow)/sinrho;
  v = (yhi-ylow)/tgrho-xi+(yhi-yi)/tgrho;
  w = ((yhi-ylow)/tgrho+xi+(yhi-yi)/tgrho)/2./cosrho;

  if(u >= uMin && u <= uMax && v <= vMax && w <= wMax) return true;

  return false;
}

// ==============================================================================
// get the particle angle phi
// in the correct range
// which is -30 to 330 degrees
// from -180 to 180 degrees
// ==============================================================================

Float_t statistics::getRelPhi(int ipart){

  Float_t phi = r2d*atan3(cy[ipart],cx[ipart]);

  if(phi < -150) return phi + 180;
  if(phi >= -150 && phi < -90) return phi + 120;
  if(phi >= -90 && phi < -30) return phi + 60;
  if(phi >= -30 && phi < 30) return phi;
  if(phi >= 30 && phi < 90) return phi - 60;
  if(phi >= 90 && phi < 150) return phi - 120;
  if(phi >= 150 && phi < 210) return phi - 180;
  if(phi >= 210 && phi < 270) return phi - 240;
  if(phi >= 270 && phi < 330) return phi - 300;
  if(phi >= 330) return phi - 360;

  return 0;
}

// ==============================================================================
// ==============================================================================

void statistics::processThetaCC(){

  // length of tchain in entries
  int nen = fChain->GetEntries();

  std::cout << std::endl;
  std::cout << "inside of processThetaCC()" << std::endl;
  std::cout << "> processing " << (int) nen/1000000 << " million events" << std::endl;
  std::cout << std::endl;

  // initialization of local variables
  Float_t theta_cc = 0.00;
  int cc_segment   = 0;
  int cc_sector    = 0;

  // histograms
  TH2F * hcct[6];

  for(int qsect=0; qsect<6; qsect++){
    hcct[qsect] =  new TH2F(Form("hcct_%d",qsect),Form("theta_cc vs cc_segm -> sector_%d",qsect),19,0,20,100,0,50);
}

  // this histo tries to debug why we have a large signal ~120 degrees
  TH2F * h_large_thetacc = new TH2F("h_large_thetacc","large theta_cc vs sector",6,0,6,18,0,17);
  TH1F * h_debug_p       = new TH1F("h_debug_p","p for theta_cc > 120 events",1000,0,5);
  TH1F * h_debug_vz      = new TH1F("h_debug_vz","vz for theta_cc > 120 events",1000,-100,100);

  TH1F * h_theta_cc[18][6];                                                             // theta_cc dist for each cc_segm, cc_sect
  TF1  * f_theta_cc[18][6];                                                             // corresponding fits for each cegm

  // generate histogram for each cc_segment (18)
  for (int ihist=0; ihist<18; ihist++){
    for(int isect=0; isect<6; isect++){
      h_theta_cc[ihist][isect] = new TH1F(Form("h_theta_cc_%d_%d",ihist,isect),Form("sector %d segment %d",isect,ihist),400,0,50);
      f_theta_cc[ihist][isect] = new TF1(Form("f_theta_cc_%d_%d",ihist,isect),"gaus");
    }
  }

  // main loop over entries
  for(int ien=0; ien<nen; ien++){
    fChain->GetEntry(ien);                                           // retrieve next event

    if (ien%1000000 == 0){std::cout << "> event number: " << (int) ien/1000000 << " million" << std::endl;}

    // main loop over particles in an event
    for(int ipart=0; ipart<gpart; ipart++){
      theta_cc   = getThetaCC(ipart);
      cc_segment = (cc_segm[ipart]%1000)/10 -1;
      cc_sector  = cc_sect[ipart]-1;

	// dont load seg -1
	if(cc_segment != -1)
	  if(cc_sector != -1)
	    if(theta_cc < 50.00 && theta_cc > 0.00)
	      if(q[ipart]<0)
		if(ec_ei[ipart]>0.06){
		    hcct[cc_sector]->Fill(cc_segment,theta_cc);
		    h_theta_cc[cc_segment][cc_sector]->Fill(theta_cc);
		  }

	// just for filling debug histo with large thetacc
	if(cc_segment != -1)
	  if(cc_sector != -1)
	    if(q[ipart]<0)
	      if(ecGeometricCut(ipart))
		if(ec_ei[ipart] > 0.06)
		  if(theta_cc > 120.00){
		    h_large_thetacc->Fill(cc_sector,cc_segment);
		    h_debug_p->Fill(p[ipart]);
		    h_debug_vz->Fill(vz[ipart]);
		  }
    }                                                               // end part loop
  }                                                                 // end event loop

  // ====================================================================================

  TCanvas * c1 = new TCanvas("c1","",1100,800);

  // open .pdf file
  c1->Print("thetaCC.pdf[");

  h_large_thetacc->Draw("colz");
  c1->Print("thetaCC.pdf");

  h_debug_p->Draw();
  c1->Print("thetaCC.pdf");

  h_debug_vz->Draw();
  c1->Print("thetaCC.pdf");

  // ==========================================
  // start fitting things
  // ==========================================

  // fit the segment, sector with gaussian
  Float_t f_theta_cc_mean[18][6];
  Float_t f_theta_cc_stddev[18][6];
  Float_t cc_seg_num[18];

  for (int khist=0; khist<18; khist++){
    for (int jsect=0; jsect<6; jsect++){

      h_theta_cc[khist][jsect]->Fit(Form("f_theta_cc_%d_%d",khist,jsect),"q");
      f_theta_cc_mean[khist][jsect]   = f_theta_cc[khist][jsect]->GetParameter(1); // mean
      f_theta_cc_stddev[khist][jsect] = f_theta_cc[khist][jsect]->GetParameter(2); // stddev

      h_theta_cc[khist][jsect]->Draw();
      c1->Print("thetaCC.pdf");
    }

    // used to pass to TGraph constructor
    cc_seg_num[khist] = (Float_t) khist;
  }

  // one TGraph for each Sector 0-5
  TGraphErrors * g_theta_cc_mean[6];
  TGraphErrors * g_theta_cc_stddev[6];

  // mean and std fit for each sector
  TF1 * f_mean[6];
  TF1 * f_stddev[6];

  // param, sector
  Float_t f_mean_params[3][6];
  Float_t f_stddev_params[3][6];

  // used to create 1-D arrays to pass to TGraph constructor
  Float_t temp_mean[18];
  Float_t temp_stddev[18];

  Float_t temp_mean_error[18];
  Float_t temp_stddev_error[18];

  Float_t empty_bins[18];

  // initialize and load graphs and fits
  for(int ksect=0; ksect<6; ksect++){

    for(int iseg=0; iseg<18; iseg++){
      temp_mean[iseg]   = f_theta_cc_mean[iseg][ksect];
      temp_stddev[iseg] = f_theta_cc_stddev[iseg][ksect];

      temp_mean_error[iseg] = f_theta_cc[iseg][ksect]->GetParError(1);
      temp_stddev_error[iseg] = f_theta_cc[iseg][ksect]->GetParError(2);

      empty_bins[iseg] = 0.00;
}

    g_theta_cc_mean[ksect]   = new TGraphErrors((Int_t) 18, cc_seg_num,temp_mean,empty_bins,temp_mean_error);
    g_theta_cc_stddev[ksect] = new TGraphErrors((Int_t) 18, cc_seg_num,temp_stddev,empty_bins,temp_stddev_error);

    f_mean[ksect]   = new TF1(Form("f_mean_%d",ksect),"pol2",0,16);
    f_stddev[ksect] = new TF1(Form("f_stddev_%d",ksect),"pol2",0,16);

    g_theta_cc_mean[ksect]->SetTitle(Form("g_theta_cc_mean_%d",ksect));
    g_theta_cc_stddev[ksect]->SetTitle(Form("g_theta_cc_stddev_%d",ksect));

    g_theta_cc_mean[ksect]->Fit(Form("f_mean_%d",ksect),"FRq");
    g_theta_cc_stddev[ksect]->Fit(Form("f_stddev_%d",ksect),"FRq");

    // get a,b,c from ax*x + bx + c,  {c=0, b=1, a=2}
    f_mean_params[0][ksect] = f_mean[ksect]->GetParameter(0);
    f_mean_params[1][ksect] = f_mean[ksect]->GetParameter(1);
    f_mean_params[2][ksect] = f_mean[ksect]->GetParameter(2);

    f_stddev_params[0][ksect] = f_stddev[ksect]->GetParameter(0);
    f_stddev_params[1][ksect] = f_stddev[ksect]->GetParameter(1);
    f_stddev_params[2][ksect] = f_stddev[ksect]->GetParameter(2);
  }

  // below lies the code with authors the cut values determined 
  // into a header file to be used by the analysis code for eID
  // open a file up for saving cut values
  ofstream ofile;
  ofile.open("cv_theta_cc.h", std::ios::out | std::ios::trunc);

  // hold the cut values [# of segment][# of sector][n_strict]
  Float_t cut_max[18][6][n_strict];
  Float_t cut_min[18][6][n_strict];

  // change the 5 if you change n_strict
  ofile << Form("Float_t cut_theta_cc_min[18][6][%d] = {",n_strict);

  // load up cut values mean +/- individually maintained
  // values for stddev by sector/segment
  // this loop also prints the arrays for min cut values
  for(int mseg=0; mseg<18; mseg++){
    ofile << "{" ;
    for(int msect=0; msect<6; msect++){
      ofile << "{" ;
      for(int istrict=0; istrict<n_strict; istrict++){
	cut_min[mseg][msect][istrict]= (f_mean_params[2][msect]*cc_seg_num[mseg]*cc_seg_num[mseg] + f_mean_params[1][msect]*cc_seg_num[mseg] + f_mean_params[0][msect]) - (f_stddev_params[2][msect]*cc_seg_num[mseg]*cc_seg_num[mseg] + f_stddev_params[1][msect]*cc_seg_num[mseg] + f_stddev_params[0][msect])*strict[istrict];
	cut_max[mseg][msect][istrict]= (f_mean_params[2][msect]*cc_seg_num[mseg]*cc_seg_num[mseg] + f_mean_params[1][msect]*cc_seg_num[mseg] + f_mean_params[0][msect]) + (f_stddev_params[2][msect]*cc_seg_num[mseg]*cc_seg_num[mseg] + f_stddev_params[1][msect]*cc_seg_num[mseg] + f_stddev_params[0][msect])*strict[istrict];
	
	ofile << cut_min[mseg][msect][istrict];
	if(istrict != (n_strict-1)){ofile << ",";}
      }
      ofile << "}" ;
      if(msect != 5){ofile << ",";}
    }
    ofile << "}" ;
    if(mseg != 17){ofile << ",";}
  }
  ofile << "};\n";

  // now print the max
  ofile << "Float_t cut_theta_cc_max[18][6][5] = {" ;

  for(int fseg=0; fseg<18; fseg++){
    ofile << "{" ;
    for(int fsect=0; fsect<6; fsect++){
      ofile << "{" ;
      for(int jstrict=0; jstrict<n_strict; jstrict++){
	ofile << cut_max[fseg][fsect][jstrict];
	if(jstrict != (n_strict-1)){ofile << ",";}
      }
      ofile << "}" ;
      if(fsect != 5){ofile << ",";}
    }
    ofile << "}" ;
    if(fseg != 17){ofile << ",";}
  }
  ofile << "}; \n";

  ofile.close();

  // using TGraphs, put example cuts onto the 2D histos
  //[min,max][sector] this will use by default strictness 3
  TGraph * cc_cuts[2][6];

  // used to pass to TGraph constructor
  Float_t temp_upper[18];
  Float_t temp_lower[18];

  // loop over sectors, initialize TGraphs
  for(int csect=0; csect<6; csect++){
    for(int cseg=0; cseg<18; cseg++){
      temp_upper[cseg] = cut_max[cseg][csect][3];
      temp_lower[cseg] = cut_min[cseg][csect][3];
    }

    cc_cuts[0][csect] = new TGraph((Int_t) 18,cc_seg_num,temp_lower);
    cc_cuts[1][csect] = new TGraph((Int_t) 18,cc_seg_num,temp_upper);

    // output to pdf
    hcct[csect]->Draw("colz");
    cc_cuts[0][csect]->Draw("same");
    cc_cuts[1][csect]->Draw("same");
    c1->Print("thetaCC.pdf");

    g_theta_cc_mean[csect]->Draw("AC*");
    c1->Print("thetaCC.pdf");

    g_theta_cc_stddev[csect]->Draw("AC*");
    c1->Print("thetaCC.pdf");
}

  // close the .pdf file
  c1->Print("thetaCC.pdf]");

  c1->Clear();
  c1->Divide(2,3);
  
  int ipad = 1;

  for(int iisect=0; iisect<6; iisect++){
    c1->cd(ipad);
    hcct[iisect]->Draw("colz same");

    cc_cuts[0][iisect]->SetLineColor(2);
    cc_cuts[1][iisect]->SetLineColor(2);

    cc_cuts[0][iisect]->SetLineWidth(3);
    cc_cuts[1][iisect]->SetLineWidth(3);

    cc_cuts[0][iisect]->Draw("same");
    cc_cuts[1][iisect]->Draw("same");

    ipad++;
}
  c1->Print("cc_cuts.png");

  return;
}

//==========================================================================

void statistics::processECSampling(){

  // length of tchain in entries
  int nen = fChain->GetEntries();

  std::cout << std::endl;
  std::cout << "inside of processECSampling()" << std::endl;
  std::cout << "> processing " << (int) nen/1000000 << " million events" << std::endl;
  std::cout << std::endl;

  // uses parameters from config.h to set up binning
  Float_t bin_size_p = (p_max-p_min)/(n_p_bins-1);
  Float_t bins[n_p_bins];

  if(debug){
    std::cout << "binning summary : " << std::endl;
    std::cout << "> bin size      : " << bin_size_p << std::endl;
    std::cout << "> number of bins: " << n_p_bins   << std::endl;
    std::cout << std::endl;
  }

  // delcare local variables
  int ec_sector = 0;

  TH2F * h_etot_vs_p[6];

  for(int qsect=0; qsect<6; qsect++){
    h_etot_vs_p[qsect] = new TH2F(Form("h_etot_vs_p_%d",qsect),Form("etot/p vs p sector_%d",qsect),100,0,5,100,0,0.5);
}

  // [bin number][sector number]
  TH1F * h_etot_p_slice[n_p_bins][6];
  TF1  * f_etot_p_slice[n_p_bins][6];

  for(int ihist=0; ihist<n_p_bins; ihist++){
    for(int isect=0; isect<6; isect++){
      h_etot_p_slice[ihist][isect] =  new TH1F(Form("h_etot_p_slice_%d_%d",ihist,isect),Form("bin_%d sector_%d",ihist,isect),100,0,0.5);
      f_etot_p_slice[ihist][isect] =  new TF1(Form("f_etot_p_slice_%d_%d",ihist,isect),"gaus",0.2,0.45);
      bins[ihist] =  ihist * bin_size_p + p_min;
    }
  }

  if(debug){
    std::cout << "histograms and binning structure have been initialized successfully" << std::endl;
    std::cout << std::endl;
    std::cout << "> first bin     : " << bins[0]          << std::endl;;
    std::cout << "> last bin      : " << bins[n_p_bins-1] << std::endl;;
    std::cout << std::endl;
  }

  // main loop over entries
  for(int ien=0; ien<nen; ien++){
    fChain->GetEntry(ien);                                           // retrieve next event

    if (ien%1000000 == 0){std::cout << "> event number: " << (int) ien/1000000 << " million" << std::endl;}

    // main loop over particles in an event
    for(int ipart=0; ipart<gpart; ipart++){
      ec_sector = ec_sect[ipart] -1;

      if(ecGeometricCut(ipart))
	if(q[ipart]<0)
	  if(ec_ei[ipart]>0.06)
	    if(etot[ipart]>0.01){
	      h_etot_vs_p[ec_sector]->Fill(p[ipart],etot[ipart]/p[ipart]);

	      // load up bins in p
	      if(floor((p[ipart]-p_min)/bin_size_p) >= 0 && floor((p[ipart]-p_min)/bin_size_p) <= (n_p_bins-1)){
		h_etot_p_slice[(Int_t) floor((p[ipart]-p_min)/bin_size_p)][ec_sector]->Fill(etot[ipart]/p[ipart]);
	      }
	    }
    }                                    // end part loop
  }                                     // end event loop

  // ====================================================================================

  // gaussian fit [bin number][ec_sector]
  Float_t f_etot_p_slice_mean[n_p_bins][6];
  Float_t f_etot_p_slice_stddev[n_p_bins][6];

  for (int jhist=0; jhist<n_p_bins; jhist++){
    for(int jsect=0; jsect<6; jsect++){
      h_etot_p_slice[jhist][jsect]->Fit(Form("f_etot_p_slice_%d_%d",jhist,jsect),"Rq");
      f_etot_p_slice_mean[jhist][jsect]   = f_etot_p_slice[jhist][jsect]->GetParameter(1);
      f_etot_p_slice_stddev[jhist][jsect] = f_etot_p_slice[jhist][jsect]->GetParameter(2);
    }
  }

  TGraphErrors * g_etot_p_mean[6];
  TGraphErrors * g_etot_p_stddev[6];

  TF1 * f_mean[6];
  TF1 * f_stddev[6];

  Float_t f_mean_params[3][6];
  Float_t f_stddev_params[3][6];

  Float_t temp_mean[n_p_bins];
  Float_t temp_stddev[n_p_bins];

  Float_t temp_mean_error[n_p_bins];
  Float_t temp_stddev_error[n_p_bins];

  Float_t empty_bins[n_p_bins];

  for(int ksect=0; ksect<6; ksect++){
    for(int ibin=0; ibin<n_p_bins; ibin++){
      temp_mean[ibin]   = f_etot_p_slice_mean[ibin][ksect];
      temp_stddev[ibin] = f_etot_p_slice_stddev[ibin][ksect];
      
      temp_mean_error[ibin]   = f_etot_p_slice[ibin][ksect]->GetParError(1);
      temp_stddev_error[ibin] = f_etot_p_slice[ibin][ksect]->GetParError(2);

      empty_bins[ibin] = 0.00;
}
    g_etot_p_mean[ksect]   = new TGraphErrors((Int_t) n_p_bins, bins, temp_mean,empty_bins,temp_mean_error);
    g_etot_p_stddev[ksect] = new TGraphErrors((Int_t) n_p_bins, bins, temp_stddev,empty_bins,temp_stddev_error);

    f_mean[ksect]   = new TF1(Form("f_mean_%d",ksect),"pol2");
    f_stddev[ksect] = new TF1(Form("f_stddev_%d",ksect),"pol2");

    g_etot_p_mean[ksect]->SetTitle(Form("g_etot_p_mean_%d",ksect));
    g_etot_p_stddev[ksect]->SetTitle(Form("g_etot_p_stddev_%d",ksect));

    g_etot_p_mean[ksect]->Fit(Form("f_mean_%d",ksect),"Fq");
    g_etot_p_stddev[ksect]->Fit(Form("f_stddev_%d",ksect),"Fq");

    f_mean_params[0][ksect]= f_mean[ksect]->GetParameter(0);
    f_mean_params[1][ksect]= f_mean[ksect]->GetParameter(1);
    f_mean_params[2][ksect]= f_mean[ksect]->GetParameter(2);

    f_stddev_params[0][ksect]= f_stddev[ksect]->GetParameter(0);
    f_stddev_params[1][ksect]= f_stddev[ksect]->GetParameter(1);
    f_stddev_params[2][ksect]= f_stddev[ksect]->GetParameter(2);
  }

  // c = 0, b = 1, a = 2
  // [param][sect][strict]
  Float_t cuts_min[3][6][n_strict];
  Float_t cuts_max[3][6][n_strict];

  for(int msect=0; msect<6; msect++){
    for(int istrict=0; istrict<n_strict; istrict++){
      cuts_min[0][msect][istrict]= f_mean_params[0][msect] - f_stddev_params[0][msect]*strict[istrict];
      cuts_min[1][msect][istrict]= f_mean_params[1][msect] - f_stddev_params[1][msect]*strict[istrict];
      cuts_min[2][msect][istrict]= f_mean_params[2][msect] - f_stddev_params[2][msect]*strict[istrict];

      cuts_max[0][msect][istrict]= f_mean_params[0][msect] + f_stddev_params[0][msect]*strict[istrict];
      cuts_max[1][msect][istrict]= f_mean_params[1][msect] + f_stddev_params[1][msect]*strict[istrict];
      cuts_max[2][msect][istrict]= f_mean_params[2][msect] + f_stddev_params[2][msect]*strict[istrict];
    }
  }

  ofstream ofile;
  ofile.open("cv_ec_sampling.h", std::ios::out | std::ios::trunc);

  // change the 5 if you change n_strict 
  ofile << "Float_t ec_sampling_min[3][6][5] = {";

  for(int iparam=0; iparam<3; iparam++){
    ofile << "{";
    for(int fsect=0; fsect<6; fsect++){
      ofile << "{";
      for(int jstrict=0; jstrict<n_strict; jstrict++){
	ofile << cuts_min[iparam][fsect][jstrict];
	if(jstrict != (n_strict-1)){ofile << ",";}
      }
      ofile << "}";
      if(fsect != 5){ofile << ",";}
    }
    ofile << "}";
    if(iparam != 2){ofile << ",";}
  }
  ofile << "}; \n";

  ofile << "Float_t ec_sampling_max[3][6][5] = {";

  for(int jparam=0; jparam<3; jparam++){
    ofile << "{";
    for(int usect=0; usect<6; usect++){
      ofile << "{";
      for(int ustrict=0; ustrict<n_strict; ustrict++){
	ofile << cuts_max[jparam][usect][ustrict];
	if(ustrict != (n_strict-1)){ofile << ",";}
      }
      ofile << "}";
      if(usect != 5){ofile << ",";}
    }
    ofile << "}";
    if(jparam != 2){ofile << ",";}
  }
  ofile << "}; \n";

  ofile.close();

  // ===========================================================

  TCanvas * c1 = new TCanvas("c1","",1100,800);

  // open .pdf files
  c1->Print("ecSampling.pdf[");

  // graphs to show upper and lower cut boundary
  // default value will be strictness 2
  TGraph * ec_max[6];
  TGraph * ec_min[6];

  Float_t temp_max[n_p_bins];
  Float_t temp_min[n_p_bins];
  
  for(int asect=0; asect<6; asect++){ 
    for(int abin=0; abin<n_p_bins; abin++){
      temp_max[abin] = cuts_max[2][asect][2]*bins[abin]*bins[abin] + cuts_max[1][asect][2]*bins[abin] + cuts_max[0][asect][2];
      temp_min[abin] = cuts_min[2][asect][2]*bins[abin]*bins[abin] + cuts_min[1][asect][2]*bins[abin] + cuts_min[0][asect][2];
    }

    ec_max[asect] = new TGraph((Int_t) n_p_bins, bins, temp_max);
    ec_min[asect] = new TGraph((Int_t) n_p_bins, bins, temp_min);
  }

  for (int khist=0; khist<n_p_bins; khist++){
    for(int tsect=0; tsect<6; tsect++){
      h_etot_p_slice[khist][tsect]->Draw();
      c1->Print("ecSampling.pdf");
    }
  }

  for(int psect=0; psect<6; psect++){

    h_etot_vs_p[psect]->Draw("colz");

    c1->SetLineWidth(3);
    c1->SetLineColor(2);

    ec_max[psect]->Draw("same");
    ec_min[psect]->Draw("same");
    c1->Print("ecSampling.pdf");

    g_etot_p_mean[psect]->Draw("AC*");
    c1->Print("ecSampling.pdf");

    g_etot_p_stddev[psect]->Draw("AC*");
    c1->Print("ecSampling.pdf");
  }

  c1->Print("ecSampling.pdf]");
  // close the .pdf files

  c1->Clear();
  c1->Divide(2,3);

  int ipad = 1;
  for(int wsect=0; wsect<6; wsect++){
    c1->cd(ipad);
    h_etot_vs_p[wsect]->Draw("colz same");

    ec_max[wsect]->SetLineWidth(3);
    ec_max[wsect]->SetLineColor(2);
    ec_min[wsect]->SetLineWidth(3);
    ec_min[wsect]->SetLineColor(2);

    ec_max[wsect]->Draw("same");
    ec_min[wsect]->Draw("same");

    ipad++;
}
  c1->Print("ec_samping_cuts.png");

  return;
}
