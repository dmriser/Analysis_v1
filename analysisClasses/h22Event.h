#include <TROOT.h>
#ifndef h22Event_h
#define h22Event_h

#include "programFiles/functions.C"

#include <cmath>
#include <iostream>
using namespace std;

#include "TVector3.h"

/**
 * h22Event class is the basis for the h22Reader class as well as ParticleFilter class.
 * It defines the general structure of an event in the h22 ntuple.
 */

class h22Event 
{
 public :
 
  // Declaration of leaf types
  string filename;
   UInt_t          evntid;
   UChar_t         ihel;
   Float_t         q_l;
   Float_t         tr_time;
   Int_t           gpart;
   Int_t           q[40];   //[gpart]
   Float_t         p[40];   //[gpart]
   Float_t         b[40];   //[gpart]
   Float_t         cx[40];   //[gpart]
   Float_t         cy[40];   //[gpart]
   Float_t         cz[40];   //[gpart]
   Float_t         vz[40];   //[gpart]
   UChar_t         dc_sect[40];   //[gpart]
   Float_t         tl1_cx[40];   //[gpart]
   Float_t         tl1_cy[40];   //[gpart]
   UChar_t         ec_sect[40];   //[gpart]
   Float_t         ec_r[40];   //[gpart]
   Float_t         ec_t[40];   //[gpart]
   Float_t         ec_ei[40];   //[gpart]
   Float_t         ec_eo[40];   //[gpart]
   Float_t         etot[40];   //[gpart]
   UChar_t         cc_sect[40];   //[gpart]
   Float_t         cc_r[40];   //[gpart]
   Float_t         cc_t[40];   //[gpart]
   UShort_t        nphe[40];   //[gpart]
   Float_t         cc_c2[40];   //[gpart]
   UChar_t         sc_sect[40];   //[gpart]
   Float_t         sc_r[40];   //[gpart]
   Float_t         sc_t[40];   //[gpart]
   Float_t         edep[40];   //[gpart]
   UChar_t         sc_pd[40];   //[gpart]
   UShort_t        cc_segm[40];   //[gpart]
   Float_t         ech_x[40];   //[gpart]
   Float_t         ech_y[40];   //[gpart]
   Float_t         ech_z[40];   //[gpart]
   Float_t         tl1_x[40];   //[gpart]
   Float_t         tl1_y[40];   //[gpart]
   Float_t         tl1_z[40];   //[gpart]
   Float_t         tl3_x[40];   //[gpart]
   Float_t         tl3_y[40];   //[gpart]
   Float_t         tl3_z[40];   //[gpart]
   Float_t         tl3_cx[40];   //[gpart]
   Float_t         tl3_cy[40];   //[gpart]
   Float_t         tl3_cz[40];   //[gpart]
   Int_t           id[40];   //[gpart]
   Float_t         vx[40];   //[gpart]
   Float_t         vy[40];   //[gpart]
   Int_t           nprt;
   Int_t           pidpart[20];   //[nprt]
   Float_t         xpart[20];   //[nprt]
   Float_t         ypart[20];   //[nprt]
   Float_t         zpart[20];   //[nprt]
   Float_t         epart[20];   //[nprt]
   Float_t         pxpart[20];   //[nprt]
   Float_t         pypart[20];   //[nprt]
   Float_t         pzpart[20];   //[nprt]
   Float_t         qpart[20];   //[nprt]
   Int_t           Ipart10[20];   //[nprt]
   Float_t         Rpart11[20];   //[nprt]
   Float_t         Rpart12[20];   //[nprt]
   Int_t           Ipart13[20];   //[nprt]
   Int_t           mcnentr;
   UChar_t         mcnpart;
   Int_t           mcst[40];   //[mcnentr]
   Int_t           mcid[40];   //[mcnentr]
   Int_t           mcpid[40];   //[mcnentr]
   Float_t         mctheta[40];   //[mcnentr]
   Float_t         mcphi[40];   //[mcnentr]
   Float_t         mcp[40];   //[mcnentr]
   Float_t         mcm[40];   //[mcnentr]
   Float_t         mcvx[40];   //[mcnentr]
   Float_t         mcvy[40];   //[mcnentr]
   Float_t         mcvz[40];   //[mcnentr]
   Float_t         mctof[40];   //[mcnentr]

   h22Event();
   virtual ~h22Event();
   void printEvent();
   double getRelativePhi(int);
   double rphi(int);
   double getMCRelativePhi(int);
   double getTheta(int);
   double theta(int);
   double getThetaCC(int);
   TVector3 uvw(int);

   double mcpx(int);
   double mcpy(int);
   double mcpz(int);
   double mcrphi(int);
   int mcsect(int);

};

h22Event::h22Event()
{}


h22Event::~h22Event()
{
  /*  delete [] q;
  delete [] p;
  delete [] b;
  delete [] cx;
  delete [] cy;
  delete [] cz;
  delete [] vz;
  delete [] dc_sect;
  delete [] tl1_cx;
  delete [] tl1_cy;
  delete [] ec_sect;
  delete [] ec_r;
  delete [] ec_t;
  delete [] ec_ei;
  delete [] ec_eo;
  delete [] etot;
  delete [] cc_sect;
  delete [] cc_r;
  delete [] cc_t;
  delete [] nphe;
  delete [] cc_c2;
  delete [] sc_sect;
  delete [] sc_r;
  delete [] sc_t;
  delete [] edep;
  delete [] sc_pd;
  delete [] cc_segm;
  delete [] ech_x;
  delete [] ech_y;
  delete [] ech_z;
  delete [] tl1_x;
  delete [] tl1_y;
  delete [] tl1_z;
  delete [] tl3_x;
  delete [] tl3_y;
  delete [] tl3_z;
  delete [] tl3_cx;
  delete [] tl3_cy;
  delete [] tl3_cz;
  delete [] id;
  delete [] vx;
  delete [] vy;
  delete [] pidpart;
  delete [] xpart;
  delete [] ypart;
  delete [] zpart;
  delete [] epart;
  delete [] pxpart;
  delete [] pypart;
  delete [] pzpart;
  delete [] qpart;
  delete [] Ipart10;
  delete [] Rpart11;
  delete [] Rpart12;
  delete [] Ipart13;
  delete [] mcst;
  delete [] mcid;
  delete [] mcpid;
  delete [] mctheta;
  delete [] mcphi;
  delete [] mcp;
  delete [] mcm;
  delete [] mcvx;
  delete [] mcvy;
  delete [] mcvz;
  delete [] mctof;  
  */
}

double h22Event::getThetaCC(int ipart)
{
  return get_thetaCC(tl3_x[ipart], tl3_y[ipart], tl3_z[ipart], tl3_cx[ipart], tl3_cy[ipart], tl3_cz[ipart]) * 180/3.14159;
}

double h22Event::getTheta(int ipart)
{
  return acos(cz[ipart])*(180/3.14159);
}

double h22Event::theta(int ipart)
{
  return acos(cz[ipart])*(180/3.14159);
}

double h22Event::getRelativePhi(int ipart)
{
  double rphi = (180/3.14159)*atan(cy[ipart]/cx[ipart]);
  if (rphi > 330.00) rphi -= 360.00; 
  rphi = rphi - 60*floor((rphi+30)/60);
  return rphi;
}

double h22Event::rphi(int ipart)
{
  double rphi = (180/3.14159)*atan(cy[ipart]/cx[ipart]);
  if (rphi > 330.00) rphi -= 360.00; 
  rphi = rphi - 60*floor((rphi+30)/60);
  return rphi;
}

double h22Event::getMCRelativePhi(int ipart)
{
  double rphi = mcphi[ipart];
  if (rphi > 330.00) rphi -= 360.00; 
  rphi = rphi - 60*floor((rphi+30)/60);
  return rphi;
}

double h22Event::mcrphi(int ipart)
{
  double rphi = mcphi[ipart];
  if (rphi > 330.00) rphi -= 360.00; 
  rphi = rphi - 60*floor((rphi+30)/60);
  return rphi;
}

void h22Event::printEvent()
{
  for (int ipart=0; ipart<gpart; ipart++)
    {
      cout.width(6); cout << ipart;
      cout.width(6); cout << q[ipart];
      cout.width(6); cout << "1";
      cout.width(6); cout << "-123";
      cout.width(4); cout << "0";
      cout.width(4); cout << "0";
      cout.width(12); cout << cx[ipart]*p[ipart];
      cout.width(12); cout << cy[ipart]*p[ipart];
      cout.width(12); cout << cz[ipart]*p[ipart];
      cout.width(12); cout << p[ipart];
      cout.width(4); cout << "0";
      cout.width(12); cout << vx[ipart];
      cout.width(12); cout << vy[ipart];
      cout.width(12); cout << vz[ipart] << endl;

    }
}

TVector3 h22Event::uvw(int ipart)
{

  Float_t u, v, w, xi, yi, zi;
  Float_t EC_the = 0.4363323;
  Float_t EC_phi;
  Float_t ylow   = -182.974;
  Float_t yhi    = 189.956;
  Float_t tgrho  = 1.95325;
  Float_t sinrho = 0.8901256;
  Float_t cosrho = 0.455715;
  Float_t rot[3][3];

  double phi = (180/3.14159)*atan(cy[ipart]/cx[ipart]);
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

  double x = ech_x[ipart];
  double y = ech_y[ipart];
  double z = ech_z[ipart];

  yi=x*rot[0][0]+y*rot[1][0]+z*rot[2][0];
  xi=x*rot[0][1]+y*rot[1][1]+z*rot[2][1];
  zi=x*rot[0][2]+y*rot[1][2]+z*rot[2][2];

  zi=zi-510.32;

  u = (yi-ylow)/sinrho;
  v = (yhi-ylow)/tgrho-xi+(yhi-yi)/tgrho;
  w = ((yhi-ylow)/tgrho+xi+(yhi-yi)/tgrho)/2./cosrho;

  TVector3 uvw(u,v,w);
  return uvw;
  
}

double h22Event::mcpx(int ipart)
{
  return mcp[ipart]*sin(mctheta[ipart]*3.14159/180)*cos(mcphi[ipart]*3.14159/180);
}

double h22Event::mcpy(int ipart)
{
  return mcp[ipart]*sin(mctheta[ipart]*3.14159/180)*sin(mcphi[ipart]*3.14159/180);
}

double h22Event::mcpz(int ipart)
{
  return mcp[ipart]*cos(mctheta[ipart]*3.14159/180);
}

int h22Event::mcsect(int ipart)
{
  return (int)(floor(mcphi[ipart]/60.0)+1);
}

#endif
