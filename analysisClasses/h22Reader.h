#include "h22Event.h"

#include <iostream>
#include <fstream>

using namespace std;

#ifndef h22Reader_h
#define h22Reader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

/**
 * h22Reader is a class which relies on the h22Event class, it constructs a chain of events from file.
 * h22Reader can be initialized for data or gsim, which contains extra mc banks.
 */

class h22Reader{
public :
   TChain         *fchain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; // current number in chain
   Int_t           GSIM; // -1 (uninitialized) 0 (data) 1 (gsim)

   h22Event event;

   // List of branches
   TBranch        *b_evntid;   //!
   TBranch        *b_ihel;   //!
   TBranch        *b_q_l;   //!
   TBranch        *b_tr_time;   //!
   TBranch        *b_gpart;   //!
   TBranch        *b_q;   //!
   TBranch        *b_p;   //!
   TBranch        *b_b;   //!
   TBranch        *b_cx;   //!
   TBranch        *b_cy;   //!
   TBranch        *b_cz;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_dc_sect;   //!
   TBranch        *b_tl1_cx;   //!
   TBranch        *b_tl1_cy;   //!
   TBranch        *b_ec_sect;   //!
   TBranch        *b_ec_r;   //!
   TBranch        *b_ec_t;   //!
   TBranch        *b_ec_ei;   //!
   TBranch        *b_ec_eo;   //!
   TBranch        *b_etot;   //!
   TBranch        *b_cc_sect;   //!
   TBranch        *b_cc_r;   //!
   TBranch        *b_cc_t;   //!
   TBranch        *b_nphe;   //!
   TBranch        *b_cc_c2;   //!
   TBranch        *b_sc_sect;   //!
   TBranch        *b_sc_r;   //!
   TBranch        *b_sc_t;   //!
   TBranch        *b_edep;   //!
   TBranch        *b_sc_pd;   //!
   TBranch        *b_cc_segm;   //!
   TBranch        *b_ech_x;   //!
   TBranch        *b_ech_y;   //!
   TBranch        *b_ech_z;   //!
   TBranch        *b_tl1_x;   //!
   TBranch        *b_tl1_y;   //!
   TBranch        *b_tl1_z;   //!
   TBranch        *b_tl3_x;   //!
   TBranch        *b_tl3_y;   //!
   TBranch        *b_tl3_z;   //!
   TBranch        *b_tl3_cx;   //!
   TBranch        *b_tl3_cy;   //!
   TBranch        *b_tl3_cz;   //!
   TBranch        *b_id;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_nprt;   //!
   TBranch        *b_pidpart;   //!
   TBranch        *b_xpart;   //!
   TBranch        *b_ypart;   //!
   TBranch        *b_zpart;   //!
   TBranch        *b_epart;   //!
   TBranch        *b_pxpart;   //!
   TBranch        *b_pypart;   //!
   TBranch        *b_pzpart;   //!
   TBranch        *b_qpart;   //!
   TBranch        *b_Ipart10;   //!
   TBranch        *b_Rpart11;   //!
   TBranch        *b_Rpart12;   //!
   TBranch        *b_Ipart13;   //!
   TBranch        *b_mcnentr;   //!
   TBranch        *b_mcnpart;   //!
   TBranch        *b_mcst;   //!
   TBranch        *b_mcid;   //!
   TBranch        *b_mcpid;   //!
   TBranch        *b_mctheta;   //!
   TBranch        *b_mcphi;   //!
   TBranch        *b_mcp;   //!
   TBranch        *b_mcm;   //!
   TBranch        *b_mcvx;   //!
   TBranch        *b_mcvy;   //!
   TBranch        *b_mcvz;   //!
   TBranch        *b_mctof;   //!

   h22Reader(int);
   ~h22Reader();
   virtual void     Init();

   int AddFile(string);
   int GetEntries(){return fchain->GetEntries();};

   h22Event GetEvent(){return event;};

   string GetFilename(){return fchain->GetCurrentFile()->GetName();}
   string GetFilenameChunk(int,int);

   void AddList(string, int);
   void AddList(string, int, int);
   void GetEntry(int ien){fchain->GetEntry(ien);};

};
#endif

h22Reader::h22Reader(int mc) 
{
  GSIM = mc;
  fchain = new TChain("h22");
}

h22Reader::~h22Reader()
{
  fchain->Delete();
}
/**< AddList takes the name of a text file and adds an integer nfiles from that list to the current TChain */
void h22Reader::AddList(string _file, int nfiles)
{
  string type[2] = {"data", "GSIM"};

  ifstream file;
  file.open(_file.c_str());
  
  string line = "";
  int ifile   = 0;

  while ( !file.eof() && (ifile < nfiles) )
    {
      file >> line;
      fchain->AddFile( line.c_str() );
      //      cout << " > Added " << line.c_str() << " to the TChain as type " << type[GSIM] << endl;
      ifile++;
    }

  file.close();

}

void h22Reader::AddList(string _file, int nfiles, int startfile)
{
  string type[2] = {"data", "GSIM"};

  ifstream file;
  file.open(_file.c_str());
  
  string line = "";
  int ifile   = 0;
  int jfile   = 0;

  while ( !file.eof() && (jfile < nfiles) )
    {
      file >> line;

      if (ifile >= startfile)
	{
	  fchain->AddFile( line.c_str() );
	  //	  cout << " > Added " << line.c_str() << " to the TChain as type " << type[GSIM] << endl;
	  jfile++;
	}
      ifile++;
    }

  file.close();

}
/**< AddFile can be used to add 1 file to the current TChain by passing the file path in */
int h22Reader::AddFile(string _fname)
{
  return fchain->AddFile(_fname.c_str());
}
/**< Init() must be run once to link the branches of the TChain to the h22Event class members*/
void h22Reader::Init()
{
   // Set branch addresses and branch pointers
   fchain->SetBranchAddress("evntid", &event.evntid, &b_evntid);
   fchain->SetBranchAddress("ihel", &event.ihel, &b_ihel);
   fchain->SetBranchAddress("q_l", &event.q_l, &b_q_l);
   //   fchain->SetBranchAddress("tr_time", &event.tr_time, &b_tr_time);
   fchain->SetBranchAddress("gpart", &event.gpart, &b_gpart);
   fchain->SetBranchAddress("q", event.q, &b_q);
   fchain->SetBranchAddress("p", event.p, &b_p);
   fchain->SetBranchAddress("b", event.b, &b_b);
   fchain->SetBranchAddress("cx", event.cx, &b_cx);
   fchain->SetBranchAddress("cy", event.cy, &b_cy);
   fchain->SetBranchAddress("cz", event.cz, &b_cz);
   fchain->SetBranchAddress("vz", event.vz, &b_vz);
   fchain->SetBranchAddress("dc_sect", event.dc_sect, &b_dc_sect);
   fchain->SetBranchAddress("tl1_cx", event.tl1_cx, &b_tl1_cx);
   fchain->SetBranchAddress("tl1_cy", event.tl1_cy, &b_tl1_cy);
   fchain->SetBranchAddress("ec_sect", event.ec_sect, &b_ec_sect);
   fchain->SetBranchAddress("ec_r", event.ec_r, &b_ec_r);
   fchain->SetBranchAddress("ec_t", event.ec_t, &b_ec_t);
   fchain->SetBranchAddress("ec_ei", event.ec_ei, &b_ec_ei);
   fchain->SetBranchAddress("ec_eo", event.ec_eo, &b_ec_eo);
   fchain->SetBranchAddress("etot", event.etot, &b_etot);
   fchain->SetBranchAddress("cc_sect", event.cc_sect, &b_cc_sect);
   fchain->SetBranchAddress("cc_r", event.cc_r, &b_cc_r);
   fchain->SetBranchAddress("cc_t", event.cc_t, &b_cc_t);
   fchain->SetBranchAddress("nphe", event.nphe, &b_nphe);
   fchain->SetBranchAddress("cc_c2", event.cc_c2, &b_cc_c2);
   fchain->SetBranchAddress("sc_sect", event.sc_sect, &b_sc_sect);
   fchain->SetBranchAddress("sc_r", event.sc_r, &b_sc_r);
   fchain->SetBranchAddress("sc_t", event.sc_t, &b_sc_t);
   fchain->SetBranchAddress("edep", event.edep, &b_edep);
   fchain->SetBranchAddress("sc_pd", event.sc_pd, &b_sc_pd);
   fchain->SetBranchAddress("cc_segm", event.cc_segm, &b_cc_segm);
   fchain->SetBranchAddress("ech_x", event.ech_x, &b_ech_x);
   fchain->SetBranchAddress("ech_y", event.ech_y, &b_ech_y);
   fchain->SetBranchAddress("ech_z", event.ech_z, &b_ech_z);
   fchain->SetBranchAddress("tl1_x", event.tl1_x, &b_tl1_x);
   fchain->SetBranchAddress("tl1_y", event.tl1_y, &b_tl1_y);
   fchain->SetBranchAddress("tl1_z", event.tl1_z, &b_tl1_z);
   fchain->SetBranchAddress("tl3_x", event.tl3_x, &b_tl3_x);
   fchain->SetBranchAddress("tl3_y", event.tl3_y, &b_tl3_y);
   fchain->SetBranchAddress("tl3_z", event.tl3_z, &b_tl3_z);
   fchain->SetBranchAddress("tl3_cx", event.tl3_cx, &b_tl3_cx);
   fchain->SetBranchAddress("tl3_cy", event.tl3_cy, &b_tl3_cy);
   fchain->SetBranchAddress("tl3_cz", event.tl3_cz, &b_tl3_cz);
   fchain->SetBranchAddress("vx", event.vx, &b_vx);
   fchain->SetBranchAddress("vy", event.vy, &b_vy);

   if (GSIM)
     {
       fchain->SetBranchAddress("mcnentr", &event.mcnentr, &b_mcnentr);
       fchain->SetBranchAddress("mcnpart", &event.mcnpart, &b_mcnpart);
       fchain->SetBranchAddress("mcst", event.mcst, &b_mcst);
       fchain->SetBranchAddress("mcid", event.mcid, &b_mcid);
       fchain->SetBranchAddress("mcpid", event.mcpid, &b_mcpid);
       fchain->SetBranchAddress("mctheta", event.mctheta, &b_mctheta);
       fchain->SetBranchAddress("mcphi", event.mcphi, &b_mcphi);
       fchain->SetBranchAddress("mcp", event.mcp, &b_mcp);
       fchain->SetBranchAddress("mcm", event.mcm, &b_mcm);
       fchain->SetBranchAddress("mcvx", event.mcvx, &b_mcvx);
       fchain->SetBranchAddress("mcvy", event.mcvy, &b_mcvy);
       fchain->SetBranchAddress("mcvz", event.mcvz, &b_mcvz);
       fchain->SetBranchAddress("mctof", event.mctof, &b_mctof);
     }

}

string h22Reader::GetFilenameChunk(int stringStart, int stringLen)
{
  if (fchain->GetEntries() == 0){ cerr << "trying to get chunk when fchain empty" << endl; exit(0);}

  string fname = fchain->GetCurrentFile()->GetName();
  return fname.substr(stringStart,stringLen); ;
}
