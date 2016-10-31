#define statistics_cxx
#include "statistics.h"
#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TTree.h>
#include "TF1.h"
#include <TGraph.h>
#include <cmath>

const Float_t r2d = 180/3.14159;
const Float_t d2r = 1/r2d;

int main(int argc, char *argv[]){

  statistics * data = new statistics(); 

  // load in files from command line
  for(int iarg=1; iarg<argc; iarg++){
    data->AddFile(argv[iarg]);
}
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

  Float_t phi = atan3(cy[ipart],cx[ipart]);

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
  std::cout << "processing " << (int) nen/1000000 << " million events" << std::endl;
  std::cout << std::endl;

  // initialization of local variables
  Float_t theta_cc = 0.00;
  int cc_segment   = 0;

  // histograms
  TH2F * hcct = new TH2F("hcct","theta_cc vs cc_segm",19,0,20,100,0,90);
 
  // main loop over entries
  for(int ien=0; ien<nen; ien++){
    fChain->GetEntry(ien);                                           // retrieve next event

    // print every 1 million events
    if (ien%1000000 == 0){std::cout << "event number: " << (int) ien/1000000 << " million" << std::endl;}

    // main loop over particles in an event
    for(int ipart=0; ipart<gpart; ipart++){
      theta_cc   = getThetaCC(ipart);
      cc_segment = (cc_segm[ipart]%1000)/10 -1;
      
      // -1 holds all tracks
      if(cc_segment != -1)
	if(theta_cc < 90.00 && theta_cc > 0.00)
	  if(q[ipart]<0)
	    if(ecGeometricCut(ipart))
	  hcct->Fill(cc_segment,theta_cc);
    }                                                               // end part loop
  }                                                                 // end event loop

  // ===========================================================================================

  // cut up the 2D histo into segment bins
  hcct->FitSlicesY();                                       // creates fits for slices
  TH1F * hcct_1 = (TH1F*) gDirectory->Get("hcct_1");        // mean   (mu)
  TH1F * hcct_2 = (TH1F*) gDirectory->Get("hcct_2");        // stddev (sigma)

  // all parameters of gaussians
  Float_t hcct_stddev[18];
  Float_t hcct_mean[18];

  // get fit parameters from histograms into arrays
  for (int iseg=0; iseg<18; iseg++){
   hcct_stddev[iseg] = hcct_2->GetBinContent(iseg);
   hcct_mean[iseg]   = hcct_1->GetBinContent(iseg);
}

  int n_strictnesses = 5;
  Float_t eidParams[n_strictnesses][18][2]; //[strictness][cc_seg][min(0), max(1)]

  // loop over strictnesses, and segments
  for(int istrict=1; istrict<=n_strictnesses; istrict++){
    for(int jseg=0; jseg<18; jseg++){
      eidParams[istrict-1][jseg][0] = hcct_mean[jseg] - (Float_t) istrict*hcct_stddev[jseg]; 
      eidParams[istrict-1][jseg][1] = hcct_mean[jseg] + (Float_t) istrict*hcct_stddev[jseg]; 
    }
}

  std::cout << "float eidParams[n_strictnesses][18][2] = {";

  for(int jstrict=1; jstrict<=n_strictnesses; jstrict++){
    std::cout << "{";
    for(int kseg=0; kseg<18; kseg++){
      std::cout << Form("{%f, %f}" ,eidParams[jstrict-1][kseg][0], eidParams[jstrict-1][kseg][1]); 
      if(kseg<17){std::cout << ", ";}
    }
    std::cout << "}, " << std::endl;  
}
  std::cout << "};" << std::endl;

  // ======================================================================================

  TCanvas * c1 = new TCanvas("c1","",1100,800);

  // open .pdf files  
  c1->Print("thetaCC.pdf[");

  hcct->Draw("colz");
  c1->Print("thetaCC.pdf");

  hcct_1->Draw();
  c1->Print("thetaCC.pdf");
 
  // close the .pdf files 
  c1->Print("thetaCC.pdf]");


  return;
}

// ==========================================================================================

void statistics::processECSampling(){

  // length of tchain in entries
  int nen = fChain->GetEntries();

  std::cout << std::endl;
  std::cout << "inside of processECSampling()" << std::endl;
  std::cout << "processing " << (int) nen/1000000 << " million events" << std::endl;
  std::cout << std::endl;

  TH2F * h_etot_vs_p = new TH2F("h_etot_vs_p","etot/p vs p",100,0,5,100,0,0.5);
 
  for(int ien=0; ien<nen; ien++){                                    // main loop over entries
    fChain->GetEntry(ien);                                           // retrieve next event
    
    if (ien%1000000 == 0){std::cout << "event number: " << (int) ien/1000000 << " million" << std::endl;}
    
    for(int ipart=0; ipart<gpart; ipart++){                          // main loop over particles in an event
      if(ecGeometricCut(ipart) && q[ipart]<0){                       // very basic eid conditions 
 	if(etot[ipart]>0.01){h_etot_vs_p->Fill(p[ipart],etot[ipart]/p[ipart]);}
      }                     
    }                                    // end part loop
  }                                                                 // end event loop

  // ====================================================================================

  h_etot_vs_p->FitSlicesY();
  TH1F * h_etot_vs_p_1 = (TH1F*) gDirectory->Get("h_etot_vs_p_1");
  TH1F * h_etot_vs_p_2 = (TH1F*) gDirectory->Get("h_etot_vs_p_2");

  TCanvas * c1 = new TCanvas("c1","",1100,800);

  // open .pdf files  
  c1->Print("ecSampling.pdf[");

  h_etot_vs_p->Draw("colz");
  c1->Print("ecSampling.pdf");

  h_etot_vs_p_1->Draw();
  c1->Print("ecSampling.pdf");

  c1->Print("ecSampling.pdf]");
  // close the .pdf files 
  
  return;
}
