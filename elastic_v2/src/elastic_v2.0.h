///////////////////////////////////////////////
/*

  David Riser, University of Connecticut 
  August 4, 2016 

  elastic_v2.0.h -> Define Class for 
  elastic second iteration. 

  
*/
///////////////////////////////////////////////

#ifndef elastic_v2_h
#define elastic_v2_h

// c++ Includes 
#include <iostream>
#include <map>
#include <vector>
using namespace std; // using in header because class is also developed here

// Root Includes 
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TGraphErrors.h"

// My Includes 
#include "../../analysisClasses/Bins.h"
#include "../../analysisClasses/ElasticEvent.h"
#include "../../analysisClasses/h22Event.h"
#include "../../analysisClasses/h22Reader.h"
#include "../../analysisClasses/ParticleFilter.h"
#include "../../analysisClasses/FaradayReader.h"
#include "../../analysisClasses/momCorr/MomCorr.C"

#include "histos.h"
#include "histos.cc"


// Model 
extern"C"{
  float elas_(float*,float*);
  float elasrad_(float*,float*,float*,float*);
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
float W_CUT = 1.10;
float beamEnergy = 5.498;

/////////////////////////////////////////////////////////
/*

  ElasticAnalysis Class 

*/
/////////////////////////////////////////////////////////

class ElasticAnalysis
{
  // Default Constr. Destr.
 public: 
  ElasticAnalysis();
  ~ElasticAnalysis();

  // Data types 
 private: 
  int n_files[2];  // Are these really needed? 
  double normalization_scale;
  bool gif_status, mom_corr_status;
  string output_name; 

  BinStructure *thetaBins, *phiBins, *momBins;
  elast_histos histos;
  h22Reader * fReader[2];  //! Two readers, one for GSIM and one for data. 
  MomCorr_e1f momCorr;
  ParticleFilter filter; 



  // User Methods
 public: 
  void add_files(vector<string>,int); //! Pass type, adds list
  void calc_xs();
  void close();        //! Save everything and close
  void draw(string);
  void fit_w();
  void get_charge(vector<string>);    //! Get charge for the runs we processed. 
  void get_model();
  void init();         //! Initialize the readers
  void loop(int);      //! Pass type, adds list 
  void set_gif_status(bool x){ gif_status = x; }
  void set_mom_corr_status(bool x){ mom_corr_status = x; }
  void set_name(string n){ output_name = n; }
};

#endif

// Definitions
ElasticAnalysis::ElasticAnalysis()
{
  for (int t=0; t<2; t++) fReader[t] = new h22Reader(t);
  thetaBins = new BinStructure(26, 24, 50);
  phiBins   = new BinStructure(20, -15, 15);
  momBins   = new BinStructure(20, 2, 4.5);

  n_files[0] = 0; n_files[1] = 0;
  normalization_scale = 0;

  histos.set_bins(thetaBins, phiBins, momBins);
  gif_status = false;
  mom_corr_status = true;  
}

ElasticAnalysis::~ElasticAnalysis()
{
  // Destroy Something 
}

void ElasticAnalysis::init()
{

  for (int t=0; t<2; t++) fReader[t]->Init();
  histos.init();
}
void ElasticAnalysis::add_files(std::vector<string> files, int index)
{

  for (int f=0; f<files.size(); f++) { 
    fReader[index]->AddFile(files[f]); 
    n_files[index]++;
  }
}

/*
void ElasticAnalysis::load(string file)
{
  histos.load(file);
  
  // Recover binning scheme by looking at histograms 
  phiBins = new BinStructure(histos.h2_rec[0][0]->GetNbinsX(),
			     histos.h2_rec[0][0]->GetXaxis()->GetXmin(),
			     histos.h2_rec[0][0]->GetXaxis()->GetXmax());

  thetaBins = new BinStructure(histos.h1_rec_theta[0][0]->GetNbinsX(),
			       histos.h1_rec_theta[0][0]->GetXaxis()->GetXmin(),
			       histos.h1_rec_theta[0][0]->GetXaxis()->GetXmax());

  cout << " File loaded " << file << endl;
}
*/

void ElasticAnalysis::close()
{
  histos.write_and_close( Form("/u/home/dmriser/mydoc/analysis/root_scripts/elastic_v2.0/out/%s.root",output_name.c_str()) );
}

void ElasticAnalysis::get_model()
{
  // Getting Fortran Model 
  cout << " Getting Model.. " << endl;

  for (int ibin=0; ibin<thetaBins->number(); ibin++)
    {
      float angle = ibin*thetaBins->width() + thetaBins->min();

      // It looks like this FORTRAN routine returns micro barns,
      // also known as outhouses.
      histos.h1_xs_model[0]->SetBinContent(ibin+1,elas_(&beamEnergy,&angle));
      histos.h1_xs_model[1]->SetBinContent(ibin+1,elasrad_(&beamEnergy,&angle,&targetRadiationLengths,&W_CUT));

      cout.width(12); cout << angle;
      cout.width(12); cout << elas_(&beamEnergy, &angle);
      cout.width(12); cout << elasrad_(&beamEnergy, &angle, &targetRadiationLengths, &W_CUT) << endl;

    }  
}

void ElasticAnalysis::get_charge(vector<string> files)
{
  cout << " Getting charge from database... " << endl;

  double Q = 0.00;

  for (int ifile=0; ifile<n_files[0]; ifile++)
    {
      int runno = atoi( files[ifile].substr(57,5).c_str() );
      FaradayCup fcup(runno);
      double dQ = fcup.getFcupCharge();
      Q += dQ;
      cout.width(8); cout << runno;
      cout.width(12); cout << dQ;
      cout.width(12); cout << Q << endl;
      histos.h1_fc_mon->Fill(dQ);      
    }

  double totalNumberElectrons = Q/qElectron;
  double normalizationScale   = toOuthouse/(totalNumberElectrons*targetProtonNumberDensity*targetLength);
  normalization_scale = normalizationScale;

  cout << " Scale: " << normalizationScale << endl; 
  cout << " Done! " << endl;

}

void ElasticAnalysis::draw(string outfile)
{
  // passing filename to drawing member function in histogram class
  histos.draw(output_name, gif_status);
}

void ElasticAnalysis::calc_xs()
{

  string SECT[7] = {"ALL","S1","S2","S3","S4","S5","S6"};

  histos.h1_rc = (TH1F*) histos.h1_xs_model[1]->Clone();
  histos.h1_rc->Divide(histos.h1_xs_model[0]);

  // Setting 0 Errors for Model 
  for (int ibin=0; ibin<thetaBins->number(); ibin++)
    {
      histos.h1_rc->SetBinError(ibin+1, 0.00);
      histos.h1_xs_model[0]->SetBinError(ibin+1, 0.00);
      histos.h1_xs_model[1]->SetBinError(ibin+1, 0.00);
    }

   // Do Cross Section here. 
  for (int s=0; s<7; s++) 
    { 

      // Doing acceptance in phi
      for (int b=0; b< histos.h1_hits_by_phi[s].size(); b++) 
	{ 
	  histos.h1_rec_by_phi[s][b]  ->Sumw2();
	  histos.h1_gen_by_phi[s][b]  ->Sumw2();
	  histos.h1_hits_by_phi[s][b] ->Sumw2();

	  string name = Form("h1_acc_by_phi_bin%d_%s",b,SECT[s].c_str());
	  histos.h1_acc_by_phi[s][b] = (TH1D*) histos.h1_rec_by_phi[s][b]->Clone(); 
	  histos.h1_acc_by_phi[s][b]->Divide( histos.h1_gen_by_phi[s][b] ); 
	  histos.h1_acc_by_phi[s][b]->SetName( name.c_str() );
	  histos.h1_acc_by_phi[s][b]->SetTitle( name.c_str() );

	  name = Form("h1_xs_by_phi_bin%d_%s",b,SECT[s].c_str());
	  histos.h1_xs_by_phi[s][b] = (TH1D*) histos.h1_hits_by_phi[s][b]->Clone();
	  histos.h1_xs_by_phi[s][b]->SetName( name.c_str() );
	  histos.h1_xs_by_phi[s][b]->SetTitle( name.c_str() );
	  histos.h1_xs_by_phi[s][b]->Scale( normalization_scale );
	  histos.h1_xs_by_phi[s][b]->Scale( 1/(phiBins->width()*thetaBins->width()) );
	  //histos.h1_xs_by_phi[s][b]->Scale( 1/thetaBins->width() );
					    
	  // Divide by cos(theta)
	  for (int tb=0; tb <= thetaBins->number(); tb++) {
	    double thisTheta = histos.h1_xs_by_phi[s][b]->GetBinCenter(tb+1); 
	    double thisValue = histos.h1_xs_by_phi[s][b]->GetBinContent(tb+1); 
	    double thisFactor = 1/cos(thisTheta*toRadians);
	    histos.h1_xs_by_phi[s][b]->SetBinContent(tb+1,thisValue*thisFactor);
	  }

	  // Convert to the differential cross section (go to events per unit angle). 
	  if (s == 0) { histos.h1_xs_by_phi[s][b]->Scale(1/(6.0)); }
	  //	  if (s == 0) { histos.h1_xs_by_phi[s][b]->Scale(1/(6*phiBins->number()*phiBins->width())); }
	  //	  else if ( b == 0 && s != 0) { histos.h1_xs_by_phi[s][b]->Scale(1/(phiBins->number()*phiBins->width())); }
	  //	  else { histos.h1_xs_by_phi[s][b]->Scale( 1/(phiBins->number()*phiBins->width())); }

	  // Integrate over phi to match the model 
	  //	  histos.h1_xs_by_phi[s][b]->Scale( 360.0 );
	  histos.h1_xs_by_phi[s][b]->Scale( 12.56637 );	  
	  //	  histos.h1_xs_by_phi[s][b]->Scale( 2*3.14159/6 );
	  
	  name = Form("h1_rxs_by_phi_bin%d_%s",b,SECT[s].c_str());
	  histos.h1_rxs_by_phi[s][b] = (TH1D*) histos.h1_xs_by_phi[s][b]->Clone();
	  histos.h1_rxs_by_phi[s][b]->SetName( name.c_str() );
	  histos.h1_rxs_by_phi[s][b]->SetTitle( name.c_str() );

	  // Doing the radiative correction. (R = N_gen(rad)/N_gen(no rad))
	  histos.h1_xs_by_phi[s][b]->Divide( histos.h1_rc );

	  name = Form("h1_xs_ratio_by_phi_bin%d_%s",b,SECT[s].c_str());
	  histos.h1_xs_ratio_by_phi[s][b] = (TH1D*) histos.h1_xs_by_phi[s][b]->Clone();
	  histos.h1_xs_ratio_by_phi[s][b]->SetName( name.c_str() );
	  histos.h1_xs_ratio_by_phi[s][b]->SetTitle( name.c_str() );
	  histos.h1_xs_ratio_by_phi[s][b]->Divide( histos.h1_xs_model[0] );

	  name = Form("h1_rxs_ratio_by_phi_bin%d_%s",b,SECT[s].c_str());
	  histos.h1_rxs_ratio_by_phi[s][b] = (TH1D*) histos.h1_rxs_by_phi[s][b]->Clone();
	  histos.h1_rxs_ratio_by_phi[s][b]->SetName( name.c_str() );
	  histos.h1_rxs_ratio_by_phi[s][b]->SetTitle( name.c_str() );
	  histos.h1_rxs_ratio_by_phi[s][b]->Divide( histos.h1_xs_model[1] );

	}
    }
      
}

void ElasticAnalysis::fit_w()
{

  // Fitting Distributions in W for each bin in momentum 
  vector<double> bins = momBins->getBins();
  vector<double> x_err;
  for (int b=0; b<bins.size(); b++) { x_err.push_back(0.000); }

  for (int s=0; s<7; s++)
    {
      // Holding mean and errors for one sector, over all the momentum bins 
      vector<double> mean; 
      vector<double> dmean; 
      vector<double> mean_err; 
      for (int b=1; b< 1+momBins->number(); b++)
	{
	  TF1 * w_fit = new TF1("w_fit","gaus",0.85,1.0);
	  w_fit->SetParameter(0, histos.h1_w_by_mom[s][b]->GetMaximum());
	  w_fit->SetParameter(1,0.938);
	  histos.h1_w_by_mom[s][b]->Fit("w_fit","rq");
	  mean    .push_back( w_fit->GetParameter(1) );
	  dmean   .push_back(0.938-w_fit->GetParameter(1));
	  mean_err.push_back( w_fit->GetParError(1) );
	} 

      // Shitty trick to use vectors to initialize 
      histos.g_mp_from_w[s] = new TGraphErrors(bins.size(), &bins[0], &mean[0], &x_err[0], &mean_err[0]);
      histos.g_mp_from_w[s]->SetName( Form("g_mp_from_w_s%d",s) );

      histos.g_dmp_from_w[s] = new TGraphErrors(bins.size(), &bins[0], &dmean[0], &x_err[0], &mean_err[0]);
      histos.g_dmp_from_w[s]->SetName( Form("g_dmp_from_w_s%d",s) );
    } 

  histos.draw_w(output_name);

}

void ElasticAnalysis::loop(int index)
{

  cout << " Looping on " << n_files[index] << " files. " << endl;

  int n_electrons = 0;
  int n_elastic   = 0;
  int n_events = fReader[index]->GetEntries();

  for (int iev=0; iev<n_events; iev++)
    {
      // Getting current event from reader. 
      fReader[index]->GetEntry(iev);
      h22Event event = fReader[index]->GetEvent();

      // Getting runno and doing PID with it. 
      int runno = atoi(fReader[index]->GetFilenameChunk(57,5).c_str());
      filter.loadEvent(event, index, runno);

      // Processing event, if we found electron. 
      int e_index = filter.getIndexByPID(11); 
      //      if (index == 1) e_index = 0;           //! Bypassing PID for Monte Carlo sample (needs debugging ~10% pass rate)

      /*
      // Doing debugging of Electron ID 
      map<string, bool> results = filter.getEIDStatus(0);
      
      // Do simplfied identification 
      if (results["EC_SAMPLING"])
	if (results["Z_VERTEX"]) 
	  if (results["EC_IN_OUT"]) 
	    if(results["CC_THETA"])
	      if(results["CC_PHI"])
		e_index = 0;
      */

      if (e_index > -123) { 
	n_electrons++;

	TLorentzVector electron(event.p[e_index]*event.cx[e_index],
				event.p[e_index]*event.cy[e_index],
				event.p[e_index]*event.cz[e_index],
				event.p[e_index]);
	
	if (index == 0 && mom_corr_status) { electron = momCorr.PcorN(electron,1,11); }

	ElasticEvent elastic_candidate(electron);
	if (elastic_candidate.getW() < W_CUT) n_elastic++;
	histos.fill(e_index, index, event, elastic_candidate);

      }      

      // Fill generated event 
      TLorentzVector gen_electron(event.mcpx(0), event.mcpy(0), event.mcpz(0), event.mcp[0]);
      ElasticEvent elastic_gen(gen_electron);
      histos.fill_gen(0, 1, event, elastic_gen);

      // Tell the user something
      if (iev%25000 == 0) 
	{
	  cout.width(8); cout << " event: "; cout.width(8); cout << iev; 
	  cout.width(8); cout << " elastic: "; cout.width(8); cout << n_elastic << endl;

	}
    }
  
  cout << " Found " << n_elastic << " elastic events " << endl;

}

