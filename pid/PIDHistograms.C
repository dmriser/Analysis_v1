/* 

   Writing out Electron ID 
   histograms 
   
   March 22, 2016

 */

// C++ Libraries
#include <iostream>
#include <cstdlib>
#include <map>
using namespace std;

// My Libraries
#include "../analysisClasses/h22Event.h"
#include "../analysisClasses/h22Reader.h"
#include "../analysisClasses/Bins.h"
#include "../analysisClasses/FaradayReader.h"
#include "../analysisClasses/ParticleFilter.h"

// CERN Root Libraries
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"

class PIDHistograms{
public:
  PIDHistograms();
 ~PIDHistograms(){
  }
  
  // 1-D 
  TH1D * h1_nphe[11][7];
  TH1F * h1_ec_edep_inner[11][7];
  TH1F * h1_ec_edep_outer[11][7];
  TH1F * h1_p[11][7];
  TH1F * h1_z_vertex[11][7];

  // 2-D 
  TH2F * h2_cc_theta[11][7];
  TH2F * h2_etot_p[11][7];
  TH2F * h2_ang_fid[11][7];
  TH2F * h2_ec_edep[11][7];
  TH2F * h2_dcr1_fid[11][7];
  TH2F * h2_dcr3_fid[11][7];
  TH2F * h2_ec_fid[11][7];
  
  void Fill(h22Event event, int ipart, int cutType);
  void Save(string outputFilename);
};

PIDHistograms::PIDHistograms(){

    string type[11] = {"allNegatives", "cuts","Z_VERTEX","CC_FID","CC_PHI","CC_THETA","DC_R1_FID","DC_R3_FID","EC_FID","EC_IN_OUT","EC_SAMPLING"};
    string sect[7]  = {"all", "s1", "s2", "s3", "s4", "s5", "s6"};
    // initialize 
    for (int itype = 0; itype < 11; itype++)
      for(int isect = 0; isect < 7; isect++)
	{
	  // 1d
	  h1_nphe[itype][isect]          = new TH1D(Form("h1_nphe_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h1_nphe_%s_%s",type[itype].c_str(),sect[isect].c_str()),100,0,100);
	  h1_ec_edep_inner[itype][isect] = new TH1F(Form("h1_ec_edep_inner_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h1_ec_edep_inner_%s_%s",type[itype].c_str(),sect[isect].c_str()),100,0,4);
	  h1_ec_edep_outer[itype][isect] = new TH1F(Form("h1_ec_edep_outer_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h1_ec_edep_outer_%s_%s",type[itype].c_str(),sect[isect].c_str()),100,0,4);
	  h1_p[itype][isect]             = new TH1F(Form("h1_p_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h1_p_%s_%s",type[itype].c_str(),sect[isect].c_str()),100,0,5);
	  h1_z_vertex[itype][isect]      = new TH1F(Form("h1_z_vertex_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h1_z_vertex_%s_%s",type[itype].c_str(),sect[isect].c_str()),100,-35,-15);
	  
	  // 2d
	  h2_cc_theta[itype][isect] = new TH2F(Form("h_cc_theta_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h_cc_theta_%s_%s",type[itype].c_str(),sect[isect].c_str()),18,0,17,100,0,60);
	  h2_etot_p[itype][isect]   = new TH2F(Form("h_etot_p_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h_etot_p_%s_%s",type[itype].c_str(),sect[isect].c_str()),100,0,5,100,0.05,0.5);
	  h2_ang_fid[itype][isect]  = new TH2F(Form("h_ang_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h_ang_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),100,-30,30,100,0,60);
	  h2_ec_edep[itype][isect]  = new TH2F(Form("h_ec_edep_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h_ec_edep_%s_%s",type[itype].c_str(),sect[isect].c_str()),100,0.01,1,100,0.01,1);
	  h2_dcr1_fid[itype][isect] = new TH2F(Form("h_dcr1_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h_dcr1_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),100,-100,100,100,-100,100);
	  h2_dcr3_fid[itype][isect] = new TH2F(Form("h_dcr3_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h_dcr3_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),100,-500,500,100,-500,500);
	  h2_ec_fid[itype][isect]   = new TH2F(Form("h_ec_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),Form("h_ec_fid_%s_%s",type[itype].c_str(),sect[isect].c_str()),100,-500,500,100,-500,500);
	}
  }

void PIDHistograms::Fill(h22Event event, int ipart, int cutType){
  
  int sector = event.dc_sect[ipart];
  
  // filling histograms for all negatives hN_abc[0][x]
  // 1-D 
  h1_nphe[cutType][0]          ->Fill(event.nphe[ipart]/10);
  h1_ec_edep_inner[cutType][0] ->Fill(event.ec_ei[ipart]);
  h1_ec_edep_outer[cutType][0] ->Fill(event.ec_eo[ipart]);
  h1_p[cutType][0]             ->Fill(event.p[ipart]);
  h1_z_vertex[cutType][0]      ->Fill(event.vz[ipart]);
  
  // 2-D 
  h2_cc_theta[cutType][0] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
  h2_etot_p[cutType][0]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
  h2_ang_fid[cutType][0]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
  h2_ec_edep[cutType][0]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
  h2_dcr1_fid[cutType][0] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
  h2_dcr3_fid[cutType][0] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
  h2_ec_fid[cutType][0]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
  
  if(sector > 0)
    {
      h1_nphe[cutType][sector]          ->Fill(event.nphe[ipart]/10);
      h1_ec_edep_inner[cutType][sector] ->Fill(event.ec_ei[ipart]);
      h1_ec_edep_outer[cutType][sector] ->Fill(event.ec_eo[ipart]);
      h1_p[cutType][sector]             ->Fill(event.p[ipart]);
      h1_z_vertex[cutType][sector]      ->Fill(event.vz[ipart]);
      
      h2_cc_theta[cutType][sector] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
      h2_etot_p[cutType][sector]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
      h2_ang_fid[cutType][sector]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
      h2_ec_edep[cutType][sector]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
      h2_dcr1_fid[cutType][sector] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
      h2_dcr3_fid[cutType][sector] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
      h2_ec_fid[cutType][sector]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
    }
  
}

void PIDHistograms::Save(string outputFilename){

  TFile outfile(outputFilename.c_str(),"recreate");  

    for (int itype = 0; itype < 11; itype++)
      for(int isect = 0; isect < 7; isect++)
	{
	  // 1d
	  h1_nphe[itype][isect]->Write();          
	  h1_ec_edep_inner[itype][isect]->Write();
	  h1_ec_edep_outer[itype][isect]->Write();
	  h1_p[itype][isect]->Write();             
	  h1_z_vertex[itype][isect]->Write();      
	  
	  // 2d
	  h2_cc_theta[itype][isect]->Write(); 
	  h2_etot_p[itype][isect]->Write();   
	  h2_ang_fid[itype][isect]->Write();  
	  h2_ec_edep[itype][isect]->Write();  
	  h2_dcr1_fid[itype][isect]->Write(); 
	  h2_dcr3_fid[itype][isect]->Write(); 
	  h2_ec_fid[itype][isect]->Write();   
	}

  outfile.Write(); 
  outfile.Close(); 
}

int main (int argc, char * argv[])
{
  int GSIM = 1;                          // data (0 - false), gsim (1 - true)
  
  // get number of files from command line 
  if (argc < 2){ cout << " expected number of files as option " << endl; exit(0);}

  int nFiles = atoi(argv[1]);
 
  // this cant be the skim files, they have too much thrown out to show effect of cuts
  string files = "mcfiles.txt";

  // setup file reader and add files
  h22Reader reader(GSIM);
  reader.AddList(files,nFiles);
  reader.Init(); // set branch addresses

  int nEvents = reader.GetEntries();

  ParticleFilter filter;
  PIDHistograms histos;
  
  // loop over events
  for (int iEvent = 0; iEvent < nEvents; iEvent++){
      reader.GetEntry(iEvent);

      // getting out local event and PID setup
      h22Event event = reader.GetEvent();
      //      int runno = atoi(reader.GetFilenameChunk(68,6).c_str());
      int runno = 0;
      filter.loadEvent(event,GSIM,runno);


      // loop over all negatives in the event 
      for(int ipart = 0; ipart < event.gpart; ipart++){     
	  if (event.q[ipart] < 0)
	    {

	      //  holds the result of all cuts, intensive 
	      map<string, bool> eID_Status = filter.getEIDStatus(ipart);
 
	      histos.Fill(event, ipart, 0);
	      if(eID_Status["Z_VERTEX"]){    histos.Fill(event, ipart, 2); }
	      if(eID_Status["CC_FID"]){      histos.Fill(event, ipart, 3); }
	      if(eID_Status["CC_PHI"]){      histos.Fill(event, ipart, 4); }
	      if(eID_Status["CC_THETA"]){    histos.Fill(event, ipart, 5); }
	      if(eID_Status["DC_R1_FID"]){   histos.Fill(event, ipart, 6); }
	      if(eID_Status["DC_R3_FID"]){   histos.Fill(event, ipart, 7); }
	      if(eID_Status["EC_FID"]){      histos.Fill(event, ipart, 8); }
	      if(eID_Status["EC_IN_OUT"]){   histos.Fill(event, ipart, 9); }
	      if(eID_Status["EC_SAMPLING"]){ histos.Fill(event, ipart, 10); }
	    }
	} // end ipart loop 
      

      // look for electron in event
      int e_index = filter.getIndexByPID(11);
      if (e_index > -123){ histos.Fill(event, e_index, 1); }

      if (iEvent%1000 == 0) { cout << "\r done " << iEvent << " of " << nEvents << flush; }
    }  // end loop over events

  histos.Save("refactorTest.root"); 


  return 0;
}
