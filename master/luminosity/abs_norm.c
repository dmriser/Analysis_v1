#include <iostream> 
#include <fstream>
#include "TTree.h"
#include "TChain.h"
#include "TFile.h" 
#include <string.h>
#include "abs_norm.h" 
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"

void abs_norm(){

  TChain *fchain = new TChain("h22"); // put files into tchain
 
  // open the text file which contains
  // runs to be processed
  // and add them to the chain
  ifstream in;                   
  in.open("runs.txt");            
  string line;                    

 // loop over lines of the text file
  while(std::getline(in,line)){      
    fchain->Add(line.c_str());        // add to our tchian
  }                               
  in.close();                     

  // begin initialization stage
  // declaration of leaves types
   UInt_t          evntid;
   UChar_t         ihel;
   Float_t         q_l;
   Int_t           gpart;
   Int_t           q[40];
   Float_t         p[40];
   Float_t         b[40];
   Float_t         cx[40];
   Float_t         cy[40];
   Float_t         cz[40];
   Float_t         vz[40];
   UChar_t         dc_sect[40];
   Float_t         tl1_cx[40];
   Float_t         tl1_cy[40];
   UChar_t         ec_sect[40];
   Float_t         ec_r[40];
   Float_t         ec_t[40];
   Float_t         ec_ei[40];
   Float_t         ec_eo[40];
   Float_t         etot[40];
   UChar_t         cc_sect[40];
   Float_t         cc_r[40];
   Float_t         cc_t[40];
   UShort_t        nphe[40];
   Float_t         cc_c2[40];
   UChar_t         sc_sect[40];
   Float_t         sc_r[40];
   Float_t         sc_t[40];
   Float_t         edep[40];
   UChar_t         sc_pd[40];
   UShort_t        cc_segm[40];
   Float_t         ech_x[40];
   Float_t         ech_y[40];
   Float_t         ech_z[40];
   Float_t         tl1_x[40];
   Float_t         tl1_y[40];
   Float_t         tl1_z[40];
   Float_t         tl3_x[40];
   Float_t         tl3_y[40];
   Float_t         tl3_z[40];
   Float_t         tl3_cx[40];
   Float_t         tl3_cy[40];
   Float_t         tl3_cz[40];
   Int_t           id[40];
   Float_t         vx[40];
   Float_t         vy[40];
 
 
   // Set branch addresses.
   fchain->SetBranchAddress("evntid",&evntid);
   fchain->SetBranchAddress("ihel",&ihel);
   fchain->SetBranchAddress("q_l",&q_l);
   fchain->SetBranchAddress("gpart",&gpart);
   fchain->SetBranchAddress("q",q);
   fchain->SetBranchAddress("p",p);
   fchain->SetBranchAddress("b",b);
   fchain->SetBranchAddress("cx",cx);
   fchain->SetBranchAddress("cy",cy);
   fchain->SetBranchAddress("cz",cz);
   fchain->SetBranchAddress("vz",vz);
   fchain->SetBranchAddress("dc_sect",dc_sect);
   fchain->SetBranchAddress("tl1_cx",tl1_cx);
   fchain->SetBranchAddress("tl1_cy",tl1_cy);
   fchain->SetBranchAddress("ec_sect",ec_sect);
   fchain->SetBranchAddress("ec_r",ec_r);
   fchain->SetBranchAddress("ec_t",ec_t);
   fchain->SetBranchAddress("ec_ei",ec_ei);
   fchain->SetBranchAddress("ec_eo",ec_eo);
   fchain->SetBranchAddress("etot",etot);
   fchain->SetBranchAddress("cc_sect",cc_sect);
   fchain->SetBranchAddress("cc_r",cc_r);
   fchain->SetBranchAddress("cc_t",cc_t);
   fchain->SetBranchAddress("nphe",nphe);
   fchain->SetBranchAddress("cc_c2",cc_c2);
   fchain->SetBranchAddress("sc_sect",sc_sect);
   fchain->SetBranchAddress("sc_r",sc_r);
   fchain->SetBranchAddress("sc_t",sc_t);
   fchain->SetBranchAddress("edep",edep);
   fchain->SetBranchAddress("sc_pd",sc_pd);
   fchain->SetBranchAddress("cc_segm",cc_segm);
   fchain->SetBranchAddress("ech_x",ech_x);
   fchain->SetBranchAddress("ech_y",ech_y);
   fchain->SetBranchAddress("ech_z",ech_z);
   fchain->SetBranchAddress("tl1_x",tl1_x);
   fchain->SetBranchAddress("tl1_y",tl1_y);
   fchain->SetBranchAddress("tl1_z",tl1_z);
   fchain->SetBranchAddress("tl3_x",tl3_x);
   fchain->SetBranchAddress("tl3_y",tl3_y);
   fchain->SetBranchAddress("tl3_z",tl3_z);
   fchain->SetBranchAddress("tl3_cx",tl3_cx);
   fchain->SetBranchAddress("tl3_cy",tl3_cy);
   fchain->SetBranchAddress("tl3_cz",tl3_cz);
   fchain->SetBranchAddress("id",id);
   fchain->SetBranchAddress("vx",vx);
   fchain->SetBranchAddress("vy",vy);
   // end initialization stage

   // setup variables needed to do 
   // fcup ratio testing
   fchain->GetEntry(0);
   string cfile     = fchain->GetCurrentFile()->GetName();
   string pfile     = cfile;
   string piece     = cfile.substr(63,15);
   int nelec, ifile;
   nelec = ifile    = 0;
   Float_t fca      = 0;
   Float_t coef     = 8.457 * 10e22;                                  // target length * proton density 

   // create some storage for my information to be generated
   // of size number of files
   string flist[fchain->GetListOfFiles()->GetEntries()];
   Float_t ratio[fchain->GetListOfFiles()->GetEntries()];
   Float_t ilum[fchain->GetListOfFiles()->GetEntries()];

   TH1F * hfc  = new TH1F("hfc","nelec/fcup accumulation",100,10,50);   // place to put ratio in a histogram
   TH1F * hlum = new TH1F("hlum","luminosity",100,5,15); // luminosity 

  // loop over all entries in fchain
  for(int ien=0;ien<=fchain->GetEntries();ien++){

    fchain->GetEntry(ien);
    cfile = fchain->GetCurrentFile()->GetName();

    // extremely crude electron cuts
    if(eDepCut(etot[0],p[0])) 
      if(qCut(q[0]))
	nelec++; 

    // check to see if we ran out of events in the current run.AXX file
    // and if so open the corresponding fcup file to calculate 
    // the required ratio 
    if (cfile != pfile){
      piece = pfile.substr(63,15);
      ifstream temp;
      temp.open(Form("/volatile/clas/clas12/dmriser/analysis/ef1_analysis/fcAccum_out/%s.fcup.fcAccum",piece.c_str()));
      temp >> fca;
      fca = fca/9264.0000;                 // scaling factor our DAQ adds to fcup_g2_2: 
      temp.close();

      // load up some info
      flist[ifile] = pfile;
      ratio[ifile] = nelec/fca;
      ilum[ifile]  = fca*coef;

      hfc->Fill(nelec/fca);
      hlum->Fill(ilum[ifile]/10e25);
      pfile = cfile;
      nelec = 0;

      std::cout << Form("> file finished with luminosity: %f",ilum[ifile]/10e25) << std::endl;
      ifile++;
    } // end long if statement
 
 } // end loop over all entries in fchain

  // output to a pdf file
  TCanvas *c1 = new TCanvas("c1","c1",1100,800);
  c1->Print("fc.pdf[");
  hfc->Draw();
  c1->Print("fc.pdf");

  // fit the histogram with a guassian
    TF1 * myfit = new TF1("myfit","gaus");
    
    hfc->Fit("myfit");
    c1->Print("fc.pdf");
    hlum->Draw();
    c1->Print("fc.pdf");
    c1->Print("fc.pdf]");
    
    Float_t mean  = myfit->GetParameter(1);
    Float_t sigma = myfit->GetParameter(2);
    Float_t max   = mean + 3.00*sigma;
    Float_t min   = mean - 3.00*sigma;

    std::cout << std::endl;
    std::cout << "> guassian fit completed" << std::endl;
    std::cout << "> mean: " << mean  << std::endl;
    std::cout << "> std : " << sigma << std::endl;
    std::cout << "> cutting between " << min << " and " << max << std::endl;
    std::cout << " " << std::endl;

    ofstream out;
    out.open("goodruns.txt");

    for(int jfile=0; jfile<fchain->GetListOfFiles()->GetEntries();jfile++){
      if(ratio[jfile]>min)
	if(ratio[jfile]<max)
      out << flist[jfile] << std::endl;
    }
    
    out.close();

} 
