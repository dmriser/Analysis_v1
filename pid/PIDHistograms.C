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

int main (int argc, char * argv[])
{
  int GSIM = 1;                          // data (0 - false), gsim (1 - true)
  
  // get number of files from command line 
  if (argc < 2){ cout << " expected number of files as option " << endl; exit(0);}

  int nFiles = atoi(argv[1]);

  TFile outfile("PIDHistograms__MC_RAD.root","recreate");

  // this cant be the skim files, they have too much thrown out to show effect of cuts
  string files = "mcfiles.txt";

  // setup file reader and add files
  h22Reader reader(GSIM);
  reader.AddList(files,nFiles);
  reader.Init(); // set branch addresses

  int nEvents = reader.GetEntries();

  ParticleFilter filter;
  
  // ------------ begin histograms ---------------
  
  // Electron Identification Histograms [11][7]
  // [0] - without cuts,  [1] - with all cuts, 
  // [3-10] diff. cuts 
  // [0] - all sectors,   [1-6] - sectors 1-6

  string type[11] = {"allNegatives", "cuts","Z_VERTEX","CC_FID","CC_PHI","CC_THETA","DC_R1_FID","DC_R3_FID","EC_FID","EC_IN_OUT","EC_SAMPLING"};
  string sect[7] = {"all", "s1", "s2", "s3", "s4", "s5", "s6"};

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


  // ------------ end histograms ---------------

  // loop over events
  for (int iEvent = 0; iEvent < nEvents; iEvent++)
    {
      reader.GetEntry(iEvent);

      // getting out local event and PID setup
      h22Event event = reader.GetEvent();
      //      int runno = atoi(reader.GetFilenameChunk(68,6).c_str());
      int runno = 0;
      filter.loadEvent(event,GSIM,runno);


      // loop over all negatives in the event 
      for(int ipart = 0; ipart < event.gpart; ipart++)
	{     
	  if (event.q[ipart] < 0)
	    {
	      int sector = event.dc_sect[ipart];
	     
	      //  holds the result of all cuts, intensive 
	      map<string, bool> eID_Status = filter.getEIDStatus(ipart);
 
	      // filling histograms for all negatives hN_abc[0][x]
 	      // 1-D 
	      h1_nphe[0][0]          ->Fill(event.nphe[ipart]/10);
	      h1_ec_edep_inner[0][0] ->Fill(event.ec_ei[ipart]);
	      h1_ec_edep_outer[0][0] ->Fill(event.ec_eo[ipart]);
	      h1_p[0][0]             ->Fill(event.p[ipart]);
	      h1_z_vertex[0][0]      ->Fill(event.vz[ipart]);
	      
	      // 2-D 
	      h2_cc_theta[0][0] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
	      h2_etot_p[0][0]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
	      h2_ang_fid[0][0]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
	      h2_ec_edep[0][0]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
	      h2_dcr1_fid[0][0] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
	      h2_dcr3_fid[0][0] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
	      h2_ec_fid[0][0]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);

	      if(sector > 0)
		{
		  h1_nphe[0][sector]          ->Fill(event.nphe[ipart]/10);
		  h1_ec_edep_inner[0][sector] ->Fill(event.ec_ei[ipart]);
		  h1_ec_edep_outer[0][sector] ->Fill(event.ec_eo[ipart]);
		  h1_p[0][sector]             ->Fill(event.p[ipart]);
		  h1_z_vertex[0][sector]      ->Fill(event.vz[ipart]);

		  h2_cc_theta[0][sector] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		  h2_etot_p[0][sector]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		  h2_ang_fid[0][sector]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		  h2_ec_edep[0][sector]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		  h2_dcr1_fid[0][sector] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		  h2_dcr3_fid[0][sector] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		  h2_ec_fid[0][sector]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		}

	      // --------- CC FID ----------

	      if(eID_Status["CC_FID"])
		{
		  // 1-D 
		  h1_nphe[3][0]          ->Fill(event.nphe[ipart]/10);
		  h1_ec_edep_inner[3][0] ->Fill(event.ec_ei[ipart]);
		  h1_ec_edep_outer[3][0] ->Fill(event.ec_eo[ipart]);
		  h1_p[3][0]             ->Fill(event.p[ipart]);
		  h1_z_vertex[3][0]      ->Fill(event.vz[ipart]);
		  
		  // 2-D 
		  h2_cc_theta[3][0] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		  h2_etot_p[3][0]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		  h2_ang_fid[3][0]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		  h2_ec_edep[3][0]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		  h2_dcr1_fid[3][0] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		  h2_dcr3_fid[3][0] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		  h2_ec_fid[3][0]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		  
		  if(sector > 0)
		    {
		      h1_nphe[3][sector]          ->Fill(event.nphe[ipart]/10);
		      h1_ec_edep_inner[3][sector] ->Fill(event.ec_ei[ipart]);
		      h1_ec_edep_outer[3][sector] ->Fill(event.ec_eo[ipart]);
		      h1_p[3][sector]             ->Fill(event.p[ipart]);
		      h1_z_vertex[3][sector]      ->Fill(event.vz[ipart]);
		      
		      h2_cc_theta[3][sector] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		      h2_etot_p[3][sector]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		      h2_ang_fid[3][sector]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		      h2_ec_edep[3][sector]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		      h2_dcr1_fid[3][sector] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		      h2_dcr3_fid[3][sector] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		      h2_ec_fid[3][sector]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		    }
		}
	      // ----- END CC FID ------------

	      // --------- CC PHI ----------
	      if(eID_Status["CC_PHI"])
		{
		  // 1-D 
		  h1_nphe[4][0]          ->Fill(event.nphe[ipart]/10);
		  h1_ec_edep_inner[4][0] ->Fill(event.ec_ei[ipart]);
		  h1_ec_edep_outer[4][0] ->Fill(event.ec_eo[ipart]);
		  h1_p[4][0]             ->Fill(event.p[ipart]);
		  h1_z_vertex[4][0]      ->Fill(event.vz[ipart]);
		  
		  // 2-D 
		  h2_cc_theta[4][0] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		  h2_etot_p[4][0]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		  h2_ang_fid[4][0]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		  h2_ec_edep[4][0]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		  h2_dcr1_fid[4][0] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		  h2_dcr3_fid[4][0] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		  h2_ec_fid[4][0]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		  
		  if(sector > 0)
		    {
		      h1_nphe[4][sector]          ->Fill(event.nphe[ipart]/10);
		      h1_ec_edep_inner[4][sector] ->Fill(event.ec_ei[ipart]);
		      h1_ec_edep_outer[4][sector] ->Fill(event.ec_eo[ipart]);
		      h1_p[4][sector]             ->Fill(event.p[ipart]);
		      h1_z_vertex[4][sector]      ->Fill(event.vz[ipart]);
		      
		      h2_cc_theta[4][sector] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		      h2_etot_p[4][sector]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		      h2_ang_fid[4][sector]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		      h2_ec_edep[4][sector]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		      h2_dcr1_fid[4][sector] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		      h2_dcr3_fid[4][sector] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		      h2_ec_fid[4][sector]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		    }
		}
	      // ----- END CC PHI ------------
	      
	      // --------- CC THETA ----------
	      if(eID_Status["CC_THETA"])
		{
		  // 1-D 
		  h1_nphe[5][0]          ->Fill(event.nphe[ipart]/10);
		  h1_ec_edep_inner[5][0] ->Fill(event.ec_ei[ipart]);
		  h1_ec_edep_outer[5][0] ->Fill(event.ec_eo[ipart]);
		  h1_p[5][0]             ->Fill(event.p[ipart]);
		  h1_z_vertex[5][0]      ->Fill(event.vz[ipart]);
		  
		  // 2-D 
		  h2_cc_theta[5][0] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		  h2_etot_p[5][0]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		  h2_ang_fid[5][0]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		  h2_ec_edep[5][0]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		  h2_dcr1_fid[5][0] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		  h2_dcr3_fid[5][0] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		  h2_ec_fid[5][0]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		  
		  if(sector > 0)
		    {
		      h1_nphe[5][sector]          ->Fill(event.nphe[ipart]/10);
		      h1_ec_edep_inner[5][sector] ->Fill(event.ec_ei[ipart]);
		      h1_ec_edep_outer[5][sector] ->Fill(event.ec_eo[ipart]);
		      h1_p[5][sector]             ->Fill(event.p[ipart]);
		      h1_z_vertex[5][sector]      ->Fill(event.vz[ipart]);
		      
		      h2_cc_theta[5][sector] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		      h2_etot_p[5][sector]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		      h2_ang_fid[5][sector]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		      h2_ec_edep[5][sector]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		      h2_dcr1_fid[5][sector] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		      h2_dcr3_fid[5][sector] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		      h2_ec_fid[5][sector]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		    }
		}
	      // ----- END CC THETA ------------
	      
	      // --------- DC Region 1 FID ----------
	      if(eID_Status["DC_R1_FID"])
		{
		  // 1-D 
		  h1_nphe[6][0]          ->Fill(event.nphe[ipart]/10);
		  h1_ec_edep_inner[6][0] ->Fill(event.ec_ei[ipart]);
		  h1_ec_edep_outer[6][0] ->Fill(event.ec_eo[ipart]);
		  h1_p[6][0]             ->Fill(event.p[ipart]);
		  h1_z_vertex[6][0]      ->Fill(event.vz[ipart]);
		  
		  // 2-D 
		  h2_cc_theta[6][0] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		  h2_etot_p[6][0]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		  h2_ang_fid[6][0]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		  h2_ec_edep[6][0]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		  h2_dcr1_fid[6][0] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		  h2_dcr3_fid[6][0] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		  h2_ec_fid[6][0]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		  
		  if(sector > 0)
		    {
		      h1_nphe[6][sector]          ->Fill(event.nphe[ipart]/10);
		      h1_ec_edep_inner[6][sector] ->Fill(event.ec_ei[ipart]);
		      h1_ec_edep_outer[6][sector] ->Fill(event.ec_eo[ipart]);
		      h1_p[6][sector]             ->Fill(event.p[ipart]);
		      h1_z_vertex[6][sector]      ->Fill(event.vz[ipart]);
		      
		      h2_cc_theta[6][sector] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		      h2_etot_p[6][sector]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		      h2_ang_fid[6][sector]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		      h2_ec_edep[6][sector]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		      h2_dcr1_fid[6][sector] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		      h2_dcr3_fid[6][sector] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		      h2_ec_fid[6][sector]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		    }
		}
	      // ----- END DC Region 1 Fid.  ------------
	      
	      // --------- dc region 3 fid  ----------
	      if(eID_Status["DC_R3_FID"])
		{
		  // 1-D 
		  h1_nphe[7][0]          ->Fill(event.nphe[ipart]/10);
		  h1_ec_edep_inner[7][0] ->Fill(event.ec_ei[ipart]);
		  h1_ec_edep_outer[7][0] ->Fill(event.ec_eo[ipart]);
		  h1_p[7][0]             ->Fill(event.p[ipart]);
		  h1_z_vertex[7][0]      ->Fill(event.vz[ipart]);
		  
		  // 2-D 
		  h2_cc_theta[7][0] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		  h2_etot_p[7][0]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		  h2_ang_fid[7][0]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		  h2_ec_edep[7][0]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		  h2_dcr1_fid[7][0] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		  h2_dcr3_fid[7][0] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		  h2_ec_fid[7][0]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		  
		  if(sector > 0)
		    {
		      h1_nphe[7][sector]          ->Fill(event.nphe[ipart]/10);
		      h1_ec_edep_inner[7][sector] ->Fill(event.ec_ei[ipart]);
		      h1_ec_edep_outer[7][sector] ->Fill(event.ec_eo[ipart]);
		      h1_p[7][sector]             ->Fill(event.p[ipart]);
		      h1_z_vertex[7][sector]      ->Fill(event.vz[ipart]);
		      
		      h2_cc_theta[7][sector] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		      h2_etot_p[7][sector]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		      h2_ang_fid[7][sector]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		      h2_ec_edep[7][sector]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		      h2_dcr1_fid[7][sector] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		      h2_dcr3_fid[7][sector] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		      h2_ec_fid[7][sector]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		    }
		}
	      // ----- END dc r3 fid ------------
	      
	      // --------- ec fid  ----------
	      if(eID_Status["EC_FID"])
		{
		  // 1-D 
		  h1_nphe[8][0]          ->Fill(event.nphe[ipart]/10);
		  h1_ec_edep_inner[8][0] ->Fill(event.ec_ei[ipart]);
		  h1_ec_edep_outer[8][0] ->Fill(event.ec_eo[ipart]);
		  h1_p[8][0]             ->Fill(event.p[ipart]);
		  h1_z_vertex[8][0]      ->Fill(event.vz[ipart]);
		  
		  // 2-D 
		  h2_cc_theta[8][0] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		  h2_etot_p[8][0]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		  h2_ang_fid[8][0]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		  h2_ec_edep[8][0]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		  h2_dcr1_fid[8][0] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		  h2_dcr3_fid[8][0] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		  h2_ec_fid[8][0]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		  
		  if(sector > 0)
		    {
		      h1_nphe[8][sector]          ->Fill(event.nphe[ipart]/10);
		      h1_ec_edep_inner[8][sector] ->Fill(event.ec_ei[ipart]);
		      h1_ec_edep_outer[8][sector] ->Fill(event.ec_eo[ipart]);
		      h1_p[8][sector]             ->Fill(event.p[ipart]);
		      h1_z_vertex[8][sector]      ->Fill(event.vz[ipart]);
		      
		      h2_cc_theta[8][sector] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		      h2_etot_p[8][sector]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		      h2_ang_fid[8][sector]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		      h2_ec_edep[8][sector]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		      h2_dcr1_fid[8][sector] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		      h2_dcr3_fid[8][sector] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		      h2_ec_fid[8][sector]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		    }
		}
	      // ----- EC fid ------------
	      
	      // --------- EC INNER vs. OUTER  ----------
	      if(eID_Status["EC_IN_OUT"])
		{
		  // 1-D 
		  h1_nphe[9][0]          ->Fill(event.nphe[ipart]/10);
		  h1_ec_edep_inner[9][0] ->Fill(event.ec_ei[ipart]);
		  h1_ec_edep_outer[9][0] ->Fill(event.ec_eo[ipart]);
		  h1_p[9][0]             ->Fill(event.p[ipart]);
		  h1_z_vertex[9][0]      ->Fill(event.vz[ipart]);
		  
		  // 2-D 
		  h2_cc_theta[9][0] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		  h2_etot_p[9][0]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		  h2_ang_fid[9][0]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		  h2_ec_edep[9][0]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		  h2_dcr1_fid[9][0] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		  h2_dcr3_fid[9][0] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		  h2_ec_fid[9][0]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		  
		  if(sector > 0)
		    {
		      h1_nphe[9][sector]          ->Fill(event.nphe[ipart]/10);
		      h1_ec_edep_inner[9][sector] ->Fill(event.ec_ei[ipart]);
		      h1_ec_edep_outer[9][sector] ->Fill(event.ec_eo[ipart]);
		      h1_p[9][sector]             ->Fill(event.p[ipart]);
		      h1_z_vertex[9][sector]      ->Fill(event.vz[ipart]);
		      
		      h2_cc_theta[9][sector] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		      h2_etot_p[9][sector]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		      h2_ang_fid[9][sector]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		      h2_ec_edep[9][sector]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		      h2_dcr1_fid[9][sector] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		      h2_dcr3_fid[9][sector] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		      h2_ec_fid[9][sector]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		    }
		}
	      // ----- END CC PHI ------------
	      
	      // --------- ec sampling ----------
	      if(eID_Status["EC_SAMPLING"])
		{
		  // 1-D 
		  h1_nphe[10][0]          ->Fill(event.nphe[ipart]/10);
		  h1_ec_edep_inner[10][0] ->Fill(event.ec_ei[ipart]);
		  h1_ec_edep_outer[10][0] ->Fill(event.ec_eo[ipart]);
		  h1_p[10][0]             ->Fill(event.p[ipart]);
		  h1_z_vertex[10][0]      ->Fill(event.vz[ipart]);
		  
		  // 2-D 
		  h2_cc_theta[10][0] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		  h2_etot_p[10][0]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		  h2_ang_fid[10][0]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		  h2_ec_edep[10][0]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		  h2_dcr1_fid[10][0] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		  h2_dcr3_fid[10][0] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		  h2_ec_fid[10][0]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		  
		  if(sector > 0)
		    {
		      h1_nphe[10][sector]          ->Fill(event.nphe[ipart]/10);
		      h1_ec_edep_inner[10][sector] ->Fill(event.ec_ei[ipart]);
		      h1_ec_edep_outer[10][sector] ->Fill(event.ec_eo[ipart]);
		      h1_p[10][sector]             ->Fill(event.p[ipart]);
		      h1_z_vertex[10][sector]      ->Fill(event.vz[ipart]);
		      
		      h2_cc_theta[10][sector] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		      h2_etot_p[10][sector]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		      h2_ang_fid[10][sector]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		      h2_ec_edep[10][sector]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		      h2_dcr1_fid[10][sector] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		      h2_dcr3_fid[10][sector] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		      h2_ec_fid[10][sector]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		    }
		}
	      // ----- END EC SAMPLING ------------
	      
	      // --------- CC PHI ----------
	      if(eID_Status["Z_VERTEX"])
		{
		  // 1-D 
		  h1_nphe[2][0]          ->Fill(event.nphe[ipart]/10);
		  h1_ec_edep_inner[2][0] ->Fill(event.ec_ei[ipart]);
		  h1_ec_edep_outer[2][0] ->Fill(event.ec_eo[ipart]);
		  h1_p[2][0]             ->Fill(event.p[ipart]);
		  h1_z_vertex[2][0]      ->Fill(event.vz[ipart]);
		  
		  // 2-D 
		  h2_cc_theta[2][0] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		  h2_etot_p[2][0]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		  h2_ang_fid[2][0]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		  h2_ec_edep[2][0]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		  h2_dcr1_fid[2][0] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		  h2_dcr3_fid[2][0] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		  h2_ec_fid[2][0]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		  
		  if(sector > 0)
		    {
		      h1_nphe[2][sector]          ->Fill(event.nphe[ipart]/10);
		      h1_ec_edep_inner[2][sector] ->Fill(event.ec_ei[ipart]);
		      h1_ec_edep_outer[2][sector] ->Fill(event.ec_eo[ipart]);
		      h1_p[2][sector]             ->Fill(event.p[ipart]);
		      h1_z_vertex[2][sector]      ->Fill(event.vz[ipart]);
		      
		      h2_cc_theta[2][sector] ->Fill((event.cc_segm[ipart]%1000)/10, event.getThetaCC(ipart));
		      h2_etot_p[2][sector]   ->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
		      h2_ang_fid[2][sector]  ->Fill(event.getRelativePhi(ipart), event.getTheta(ipart));
		      h2_ec_edep[2][sector]  ->Fill(event.ec_ei[ipart], event.ec_eo[ipart]);
		      h2_dcr1_fid[2][sector] ->Fill(event.tl1_x[ipart], event.tl1_y[ipart]);
		      h2_dcr3_fid[2][sector] ->Fill(event.tl3_x[ipart], event.tl3_y[ipart]);
		      h2_ec_fid[2][sector]   ->Fill(event.ech_x[ipart], event.ech_y[ipart]);
		    }
		}
	      // ----- END ZVERT ------------
	      

	    }
	} // end ipart loop 
      

      // look for electron in event
      int e_index = filter.getIndexByPID(11);
      
      if (e_index > -123)
	{
	  int sector = event.dc_sect[e_index];
	  
	  // filling histograms for all negatives hN_abc[0][x]
	  // 1-D 
	  h1_nphe[1][0]          ->Fill(event.nphe[e_index]/10);
	  h1_ec_edep_inner[1][0] ->Fill(event.ec_ei[e_index]);
	  h1_ec_edep_outer[1][0] ->Fill(event.ec_eo[e_index]);
	  h1_p[1][0]             ->Fill(event.p[e_index]);
	  h1_z_vertex[1][0]      ->Fill(event.vz[e_index]);
	  
	  // 2-D 
	  h2_cc_theta[1][0] ->Fill((event.cc_segm[e_index]%1000)/10, event.getThetaCC(e_index));
	  h2_etot_p[1][0]   ->Fill(event.p[e_index], event.etot[e_index]/event.p[e_index]);
	  h2_ang_fid[1][0]  ->Fill(event.getRelativePhi(e_index), event.getTheta(e_index));
	  h2_ec_edep[1][0]  ->Fill(event.ec_ei[e_index], event.ec_eo[e_index]);
	  h2_dcr1_fid[1][0] ->Fill(event.tl1_x[e_index], event.tl1_y[e_index]);
	  h2_dcr3_fid[1][0] ->Fill(event.tl3_x[e_index], event.tl3_y[e_index]);
	  h2_ec_fid[1][0]   ->Fill(event.ech_x[e_index], event.ech_y[e_index]);
	  
	  if(sector > 0)
	    {
	      h1_nphe[1][sector]          ->Fill(event.nphe[e_index]/10);
	      h1_ec_edep_inner[1][sector] ->Fill(event.ec_ei[e_index]);
	      h1_ec_edep_outer[1][sector] ->Fill(event.ec_eo[e_index]);
	      h1_p[1][sector]             ->Fill(event.p[e_index]);
	      h1_z_vertex[1][sector]      ->Fill(event.vz[e_index]);
	      
	      h2_cc_theta[1][sector] ->Fill((event.cc_segm[e_index]%1000)/10, event.getThetaCC(e_index));
	      h2_etot_p[1][sector]   ->Fill(event.p[e_index], event.etot[e_index]/event.p[e_index]);
	      h2_ang_fid[1][sector]  ->Fill(event.getRelativePhi(e_index), event.getTheta(e_index));
	      h2_ec_edep[1][sector]  ->Fill(event.ec_ei[e_index], event.ec_eo[e_index]);
	      h2_dcr1_fid[1][sector] ->Fill(event.tl1_x[e_index], event.tl1_y[e_index]);
	      h2_dcr3_fid[1][sector] ->Fill(event.tl3_x[e_index], event.tl3_y[e_index]);
	      h2_ec_fid[1][sector]   ->Fill(event.ech_x[e_index], event.ech_y[e_index]);
	    }
	}

      cout << "\r done " << iEvent << " of " << nEvents;
    }  // end loop over events

  outfile.Write();
  outfile.Close();

  cout << "\n > done for " << nFiles << " files, hasta pronto! " << endl;

  return 0;
}
