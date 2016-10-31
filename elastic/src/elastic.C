/* 

               elastic.C 
       Program to construct elastic cross section for 
       e1f data.


              David Riser, 
        University of Connecticut 
              Mar 4, 2016

 */


// C++ Libraries
#include <iostream>
#include <cstdlib>
#include <vector>
using namespace std;

// My Libraries
#include "../../analysisClasses/Bins.h"
#include "../../analysisClasses/ElasticEvent.h"
#include "../../analysisClasses/FaradayReader.h"
#include "../../analysisClasses/h22Event.h"
#include "../../analysisClasses/h22Reader.h"
#include "../../analysisClasses/ParticleFilter.h"
#include "histos.cc"
#include "histos.h"

// CERN Root Libraries
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TStyle.h"

// e1f momentum corrections by Marco 
#include "../../analysisClasses/momCorr/MomCorr.C"

// model cross section from FORTRAN elaslib.o elaslib.f
extern"C"{
  float elas_(float*,float*);
  float elasrad_(float*,float*,float*,float*);
}

int main (int argc, char * argv[])
{

  // get number of files from command line or die 
  if (argc < 2)
    {
      cout << " expected number of files as option " << endl;
      exit(0);
    }

  // ------- setup constants, conversions, and physics -------
  const double toDegrees = 180.0/3.14159;
  const double toRadians = 1/toDegrees;
  const double toBarn    = 1.0e24;
  const double toOuthouse = 1.0e30;
  const double toCoulomb  = 1.0e-6;

  const double qElectron   = 1.60217e-13;  // uC
  const double NAvagadro   = 6.02214e23;   // part mol-1
  const double densityLH2  = 70.8e-3;      // g cm-3
  const double H2MolarMass = 2*1.00794;    // g mol-1
  const double targetLength = 5.00;        // cm
  const double radLength = 5.31;           // cm 

  const double targetProtonNumberDensity = 2*NAvagadro*densityLH2/H2MolarMass;

  //  Used for model calc. have to be float.
  float targetRadiationLengths = (float) targetLength/radLength;
  float W_CUT = 1.05;

  ParticleFilter filter;
  MomCorr_e1f * momCorr = new MomCorr_e1f();
  BinStructure thetaBins(10,18.0,39.0);
  BinStructure relPhiBins(10,-15,15);

  thetaBins.setName("theta bins");
  relPhiBins.setName("rel. phi bins");

  elast_histos histos;
  histos.init(thetaBins, relPhiBins);

  // ----------- setup data management, file lists, readers -------------
  int GSIM[3] = {1,1,0};                          // data (0 - false), gsim (1 - true)

  // There are a lot more events in data files, so number of MC files is multiplied by 20.
  int nFiles[3]   = {atoi(argv[1])*20,atoi(argv[1])*20,atoi(argv[1])};
  string files[3] = {"/u/home/dmriser/mydoc/analysis/root_scripts/elastic/src/mc_norad.txt","/u/home/dmriser/mydoc/analysis/root_scripts/elastic/src/mc_rad.txt","/u/home/dmriser/mydoc/analysis/root_scripts/fcup/goodruns.txt"};

  // Begin loop over type {GSIM_NORAD, GSIM_RAD, DATA} 
  for (int itype = 0; itype<3; itype++)
    {
      h22Reader fReader(GSIM[itype]);
      fReader.AddList(files[itype],nFiles[itype]);
      fReader.Init(); 
      
      // loop over events
      for (int iEvent = 0; iEvent < fReader.GetEntries(); iEvent++)
	{
	  fReader.GetEntry(iEvent);
	  
	  // store each event locally for the duration of the loop 
	  h22Event event = fReader.GetEvent();
	  
	  int runno = atoi(fReader.GetFilenameChunk(57,5).c_str());
	  
	  // we don't need to pass the real run number for MC 
	  filter.loadEvent(event,GSIM[itype],runno);
	  int e_index    = filter.getIndexByPID(11);
	  int prot_index = filter.getIndexByPID(2212);

	  if (GSIM[itype])
	    {	      
		  // Both are MC, now choose with or without radiative effects.
		  switch(itype)
		    {
		     
		      // Elastic ONLY
		    case 0:
		      {
			
			if (1)//(e_index > -123 && prot_index > -123) // N.H. Convention 
			  {
			    TLorentzVector recElectron(event.cx[e_index]*event.p[e_index],
						       event.cy[e_index]*event.p[e_index],
						       event.cz[e_index]*event.p[e_index],
						       event.p[e_index]);						       
			    

			    int sector = event.dc_sect[e_index]%1000;
			    ElasticEvent recEvent(recElectron);

			    histos.h1_W[itype][sector]    ->Fill(recEvent.getW());			    
			    histos.h1_W[itype][0]         ->Fill(recEvent.getW());

			    if (recEvent.getW() < W_CUT)
			      {
				histos.h2_rec[itype][sector]       ->Fill(event.getRelativePhi(e_index),event.getTheta(e_index));					   
				histos.h1_rec_theta[itype][sector] ->Fill(event.getTheta(e_index));

				histos.h1_x[itype][sector]         ->Fill(recEvent.getx());
				histos.h1_MM2[itype][sector]       ->Fill(recEvent.getMM2());
				histos.h1_QQ[itype][sector]        ->Fill(recEvent.getQQ());
				histos.h1_EP_angle[itype][sector]  ->Fill(recEvent.getEPAngle());
				histos.h1_ME[itype][sector]        ->Fill(recEvent.getME());
				// load ALL histograms
				histos.h2_rec[itype][0]       ->Fill(event.getRelativePhi(e_index),event.getTheta(e_index));			
				histos.h1_rec_theta[itype][0] ->Fill(event.getTheta(e_index));

				histos.h1_x[itype][0]         ->Fill(recEvent.getx());
				histos.h1_MM2[itype][0]       ->Fill(recEvent.getMM2());
				histos.h1_QQ[itype][0]        ->Fill(recEvent.getQQ());
				histos.h1_EP_angle[itype][0]  ->Fill(recEvent.getEPAngle());
				histos.h1_ME[itype][0]        ->Fill(recEvent.getME());
			      }
			  }
			
			TLorentzVector genElectron(event.mcpx(0),
						   event.mcpy(0),
						   event.mcpz(0),
						   event.mcp[0]);

			ElasticEvent genEvent(genElectron);
			
			int mcsect = (int) (floor(event.mcphi[0]/60.0)+1);
			if (1)//(genEvent.getW() < W_CUT)
			  {
			    histos.h2_gen[itype][0]       ->Fill(event.getMCRelativePhi(0),event.mctheta[0]);
			    histos.h1_gen_theta[itype][0] ->Fill(event.mctheta[0]);
			    histos.h2_gen[itype][mcsect]       ->Fill(event.getMCRelativePhi(0),event.mctheta[0]);
			    histos.h1_gen_theta[itype][mcsect] ->Fill(event.mctheta[0]);

			  }
		      }
		      // Elastic + Radiative Processes 
		    case 1:
		      {

			if (1)//(e_index > -123) // N.H. Convention 
			  {
 			    TLorentzVector recElectron(event.cx[e_index]*event.p[e_index],
						      event.cy[e_index]*event.p[e_index],
						      event.cz[e_index]*event.p[e_index],
						      event.p[e_index]);
			    			    
			    ElasticEvent recEvent(recElectron);

			    int sector = event.dc_sect[e_index]%1000;

			    if (sector>0)
			      {					
				histos.h2_W_theta[1][0]      ->Fill(recEvent.getW(),recElectron.Theta()*toDegrees);
				histos.h2_W_theta[1][sector] ->Fill(recEvent.getW(),recElectron.Theta()*toDegrees);
				histos.h1_W[itype][sector]         ->Fill(recEvent.getW());
				histos.h1_W[itype][0]         ->Fill(recEvent.getW());
			      }

			    if (recEvent.getW() < W_CUT) 
			      {
				histos.h2_rec[itype][sector]       ->Fill(event.getRelativePhi(e_index),event.getTheta(e_index));
				histos.h1_rec_theta[itype][sector] ->Fill(event.getTheta(e_index));
				histos.h1_x[itype][sector]         ->Fill(recEvent.getx());
				histos.h1_MM2[itype][sector]       ->Fill(recEvent.getMM2());
				histos.h1_QQ[itype][sector]        ->Fill(recEvent.getQQ());
				histos.h1_EP_angle[itype][sector]  ->Fill(recEvent.getEPAngle());
				histos.h1_ME[itype][sector]        ->Fill(recEvent.getME());
				// load ALL histograms
				histos.h2_rec[itype][0]       ->Fill(event.getRelativePhi(e_index),event.getTheta(e_index));			
				histos.h1_rec_theta[itype][0] ->Fill(event.getTheta(e_index));

				histos.h1_x[itype][0]         ->Fill(recEvent.getx());
				histos.h1_MM2[itype][0]       ->Fill(recEvent.getMM2());
				histos.h1_QQ[itype][0]        ->Fill(recEvent.getQQ());
				histos.h1_EP_angle[itype][0]  ->Fill(recEvent.getEPAngle());
				histos.h1_ME[itype][0]        ->Fill(recEvent.getME());
			      } 
			  }

			TLorentzVector genElectron(event.mcpx(0),
						   event.mcpy(0),
						   event.mcpz(0),
						   event.mcp[0]);
			TLorentzVector genProton(event.mcpx(1),
						   event.mcpy(1),
						   event.mcpz(1),
						   event.mcp[1]);
			
			ElasticEvent genEvent(genElectron);
			
			int mcsect = (int) (floor(event.mcphi[0]/60.0)+1);
			    
			    // do photon stuff here 
 			TLorentzVector target(0,0,0,0.938);
			TLorentzVector beam(0,0,5.498,5.498);
			TLorentzVector BHPhoton = (beam + target) - (genElectron + genProton);

			histos.h2_W_theta[0][0]      ->Fill(genEvent.getW(),genElectron.Theta()*toDegrees);
			histos.h2_W_theta[0][mcsect] ->Fill(genEvent.getW(),genElectron.Theta()*toDegrees);


			// Radiative Effects Analysis 
			histos.h1_e_E->Fill(genElectron.E());
			histos.h1_gamma_E       ->Fill(BHPhoton.E());
			histos.h1_gamma_theta   ->Fill(BHPhoton.Theta()*toDegrees);
			histos.h1_proton_E ->Fill(genProton.E());
			histos.h1_proton_theta   ->Fill(event.mctheta[1]);

			// Go to frame of virtual photon 

			BHPhoton.RotateZ(-1*(beam+target-genElectron).Phi() - 3.14159);
			BHPhoton.RotateY((beam+target-genElectron).Theta());
			histos.h1_gamma_phistar   ->Fill(BHPhoton.Phi()*toDegrees);

			genProton.RotateZ(-1*(beam+target-genElectron).Phi() - 3.14159);
			genProton.RotateY((beam+target-genElectron).Theta());
			histos.h1_proton_phistar->Fill(genProton.Phi()*toDegrees);

			if(1)// (genEvent.getW() < W_CUT)
			  {
			    histos.h2_gen[itype][0]       ->Fill(event.getMCRelativePhi(0),event.mctheta[0]);
			    histos.h1_gen_theta[itype][0] ->Fill(event.mctheta[0]);
			    histos.h2_gen[itype][mcsect]       ->Fill(event.getMCRelativePhi(0),event.mctheta[0]);
			    histos.h1_gen_theta[itype][mcsect] ->Fill(event.mctheta[0]);
			  }
		      }
		    }
	    }

	  // data events fall here
	  else
	    {
	      // good time for cookies is 2-3pm
	      if (e_index > -123 && prot_index > -123)
		{
		  int sector = event.dc_sect[e_index]%1000;

		  // make and correct the electron 4-vector 
		      TLorentzVector electron(event.cx[e_index]*event.p[e_index],
					      event.cy[e_index]*event.p[e_index],
					      event.cz[e_index]*event.p[e_index],
					      event.p[e_index]);

		      electron = momCorr->PcorN(electron,1,11);
		      
		      ElasticEvent eevent(electron);

		      if (sector>0)
			{
			  histos.h2_W_theta[2][0]->Fill(eevent.getW(),electron.Theta()*toDegrees);
			  histos.h2_W_theta[2][sector]->Fill(eevent.getW(),electron.Theta()*toDegrees);
			  histos.h1_W[itype][sector]         ->Fill(eevent.getW());
			  histos.h1_W[itype][0]         ->Fill(eevent.getW());
			}

		      if (eevent.getW() < W_CUT && sector > 0) 
			{
			  histos.h2_rec[itype][sector]       ->Fill(event.getRelativePhi(e_index),event.getTheta(e_index));
			  histos.h2_gen[itype][sector]       ->Fill(event.getRelativePhi(e_index),event.getTheta(e_index));
			  histos.h1_rec_theta[itype][sector] ->Fill(event.getTheta(e_index));
			  histos.h1_gen_theta[itype][sector] ->Fill(event.getTheta(e_index));
			  histos.h1_EP_angle[itype][sector]  ->Fill(eevent.getEPAngle());
			  histos.h1_MM2[itype][sector]       ->Fill(eevent.getMM2());
			  histos.h1_EP_angle[itype][sector]  ->Fill(eevent.getEPAngle());
			  histos.h1_x[itype][sector]         ->Fill(eevent.getx());
			  histos.h1_ME[itype][sector]        ->Fill(eevent.getME());
			  histos.h1_QQ[itype][sector]        ->Fill(eevent.getQQ());

			  // load ALL histograms
			  histos.h2_rec[itype][0]       ->Fill(event.getRelativePhi(e_index),event.getTheta(e_index));
			  histos.h2_gen[itype][0]       ->Fill(event.getRelativePhi(e_index),event.getTheta(e_index));
			  histos.h1_rec_theta[itype][0] ->Fill(event.getTheta(e_index));
			  histos.h1_gen_theta[itype][0] ->Fill(event.getTheta(e_index));
			  histos.h1_EP_angle[itype][0]  ->Fill(eevent.getEPAngle());
			  histos.h1_MM2[itype][0]       ->Fill(eevent.getMM2());
			  histos.h1_EP_angle[itype][0]  ->Fill(eevent.getEPAngle());
			  histos.h1_x[itype][0]         ->Fill(eevent.getx());
			  histos.h1_ME[itype][0]        ->Fill(eevent.getME());
			  histos.h1_QQ[itype][0]        ->Fill(eevent.getQQ());
			}
		}
	    }

	  cout << "\r done " << iEvent << " of " << fReader.GetEntries() << " for type " << itype;
	}  // end loop over events
    } // end loop on types 

  
  // ------------- get accumulation over the runs -------------
  ifstream fListForQ(files[2].c_str());
  string cFile;

  int kFile = 0;

  // returned in uC 
  double Q  = 0;

  cout << "\n > Getting electron flux from file... " << endl; 
  
  while ( !fListForQ.eof() && kFile < nFiles[2])
    {
      fListForQ >> cFile;

      cout << " working on " << atoi(cFile.substr(57,5).c_str()) << endl;

      int runno = atoi( cFile.substr(57,5).c_str() );
      FaradayCup fcup(runno);
      Q += fcup.getFcupCharge();
      histos.h1_fc_mon->Fill(fcup.getFcupCharge());
      kFile++;
}

  fListForQ.close();
  
  double totalNumberElectrons = Q/qElectron;
  double normalizationScale   = toOuthouse/(totalNumberElectrons*targetProtonNumberDensity*targetLength);

  ofstream normScaleFile("elast_scale.dat",ios::trunc);
  normScaleFile << normalizationScale;
  normScaleFile.close();

  // load scaling histogram for making calculations with above info from other codes
  // set errors to 0 for scale and data, not true, revisit 
  for (int p=0; p<relPhiBins.number(); p++) 
    for (int t=0; t<thetaBins.number(); t++)
      {
	histos.h2_scale->SetBinContent(p+1,t+1,1/normalizationScale);
	histos.h2_scale->SetBinError(p+1,t+1,0);
	for(int s=0; s<7; s++){	
	  histos.h2_rec[2][s]->SetBinError(p+1,t+1,0);
	  histos.h2_gen[2][s]->SetBinError(p+1,t+1,0);
	}
      }

  // ----- Getting model from FORTRAN code, pass values as address.
  float beamEnergy     = 5.498;
  vector<double> v_theta_bins = thetaBins.getBins();
  for (int ibin=0; ibin< v_theta_bins.size(); ibin++)
    {
      float angle = (float) v_theta_bins[ibin];

      // It looks like this FORTRAN routine returns micro barns,
      // also known as outhouses. 
      histos.h1_xs_model[0]->SetBinContent(ibin+1,elas_(&beamEnergy,&angle));
      histos.h1_xs_model[1]->SetBinContent(ibin+1,elasrad_(&beamEnergy,&angle,&targetRadiationLengths,&W_CUT));
    }

  // write out the histograms 
  histos.write_and_close("/u/home/dmriser/mydoc/analysis/root_scripts/elastic/out/elastic.root");

  cout << endl;
  cout << " --------------- Summary --------------- " << endl;
  thetaBins.print();
  relPhiBins.print();

  cout << " Faraday cup electrons:       " << Q/qElectron << endl;
  cout << " Faraday cup accum:           " << Q << endl;
  cout << " Cross Section Normalization: " << toBarn/(totalNumberElectrons*targetProtonNumberDensity*targetLength) << endl; 
  cout << " Integrated Luminosity:       " << totalNumberElectrons*targetProtonNumberDensity*targetLength << endl;
  cout << endl;
  

  return 0;
  
}
