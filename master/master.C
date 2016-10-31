#define master_cxx

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <unistd.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TStyle.h>
#include <TTree.h>

#include "master.h"
#include "cv_pp_beta.h"
#include "cv_pm_beta.h"
#include "cv_proton_beta.h"
#include "cv_proton_beta_MC.h"
#include "event_elastic.c"

int main(int argc, char *argv[]){
    
    master * data = new master();
    
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "================================" << std::endl;
    std::cout << "            master.C "            << std::endl;
    std::cout << "================================" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    
    int opt       = 0;
    int nfiles    = 1;
    
    // output file stuff
    bool saveFile = false;
    std::string rootFile;
    TFile * myRootOutput;
    
    data->setMCstatus(false);
    
    // parse command line args
    while((opt = getopt(argc,argv,"dmn:ho:")) != -1){
        switch(opt){
            case 'n':
                nfiles = atoi(optarg);
                break;
                
            case 'd':
                data->setMCstatus(false);
                break;
                
            case 'm':
                data->setMCstatus(true);
                break;
                
            case 'o':
                saveFile = true;
                rootFile = optarg;
                break;
                
            case 'h':
                std::cout << "usage: [-h show this message] [-n n_files] [-o filename.root]" << std::endl;
                break;
        }
    }
    
    // print out current executable parameters to user
    std::cout << "> program parameters" << std::endl;
    std::cout << "     > number of files requested: " << nfiles << std::endl;
    
    if(debug)
        std::cout << "     > global debug mode is active" << std::endl;
    if(validate_eid)
        std::cout << "     > validating eID" << std::endl;
    if(!validate_eid)
        std::cout << "     > not validating eID" << std::endl;
    if(debug_phi_cc)
        std::cout << "     > debugging phi cc" << std::endl;
    if(!debug_phi_cc)
        std::cout << "     > not debugging phi cc" << std::endl;
    if(genListOfGoodRuns)
        std::cout << "     > generating good run list" << std::endl;
    if(!genListOfGoodRuns)
        std::cout << "     > using previously existing good run list" << std::endl;
    
    // print out cut values for eID
    std::cout << std::endl;
    std::cout << "> eID parameters" << std::endl;
    std::cout << "     > ec inner edep sctrictness  : " << eid_stricts[0] << std::endl;
    std::cout << "     > ec geometric sctrictness   : " << eid_stricts[1] << std::endl;
    std::cout << "     > ec sampling strictness     : " << eid_stricts[2] << std::endl;
    std::cout << "     > cc theta segment strictness: " << eid_stricts[3] << std::endl;
    std::cout << "     > vertex position strictness : " << eid_stricts[4] << std::endl;
    std::cout << std::endl;
    
    // we have to load all runs for this subroutine
    // so it uses a seperate instance of master class
    if(genListOfGoodRuns){
        
        std::cout << "> now loading all data files for good run list generation" << std::endl;
        
        master * allData = new master();
        
        // contains all runs
        std::ifstream ifile;
        ifile.open("runs.txt");
        
        // check to make sure that things are open
        if(!ifile.is_open()){
            std::cout << "> error opening runs.txt"                                << std::endl;
            std::cout << "> trying creating a text file called runs.txt"           << std::endl;
            exit(0);
        }
        
        std::string line;
        
        while(getline(ifile,line)){
            allData->AddFile(const_cast<char*>(line.c_str()));
        }
        
        allData->Init( data->getMCstatus() ); // true false for isMC
        allData->getGoodRuns();
    }
    
    // now that goodruns.txt is up-to-date, we can load
    // the desired number of files into memory
    int ifile  = 0;
    std::ifstream goodruns;    

    if(!data->getMCstatus())
    	goodruns.open("goodruns.txt");

    if(data->getMCstatus())
	goodruns.open("mcruns.txt");

    if(!goodruns.is_open()){
        std::cout << "> error opening goodruns.txt" << std::endl;
        return 0;
    }
    
    std::string line1;
    
    while(ifile < nfiles && getline(goodruns,line1)){
        data->AddFile(const_cast<char*>(line1.c_str()));
        ifile++;
    }
    
    std::cout << Form("> %d files loaded successfully",ifile) << std::endl;
    
    data->Init( data->getMCstatus() );
    std::cout << "> tree structure has been initialized" << std::endl;
    
    
    if (saveFile)
        myRootOutput = new TFile(rootFile.c_str(), "recreate");
    
    
    data->Process();
    //    data->generateElectronCutValues();

    if(debug_phi_cc)
        data->showPhiCC();
    
    if(validate_eid)
        data->eidValidation();
    
    if(debug_pp_mass)
        data->ppMissingMass();
    
    if(run_pp_time_corr)
        data->ppTimeCorrelation();
    
    if(run_pp_beta)
        data->pionBeta();
    
    if (run_proton_beta)
        data->runProtonBeta();
    
    if (saveFile)
    {
        myRootOutput->Write();
        myRootOutput->Map();
        std::cout << Form(" > histograms saved to %s ",rootFile.c_str()) << std::endl;
    }
    
    
    return 0;
}

// ==================================================
//
// ==================================================

void master::Process(){
    
    std::cout << std::endl;
    std::cout << "inside Process()" << std::endl;
    std::cout << std::endl;
    
    int nen = fChain->GetEntries();
    int nel = 0;
    int npp = 0;
    int npm = 0;
    int np  = 0;

    Float_t theta_cc = 0;
    Int_t cc_segment = 0;
    
    // histograms
    TH1F * h_p_transf_sq = new TH1F("h_p_transf_sq","Photon Virtuality; Q^2 (GeV/c)^2",100,0,5);
    TH1F * h_W           = new TH1F("h_W","Final Hadronic State Momentum W; Momentum W (GeV/c)",100,-5,5);
    TH1F * h_x           = new TH1F("h_x","Bjorken x; x",100,0,1);
    TH2F * h_QQ_w        = new TH2F("h_QQ_w","Q^2 vs. W;  W (GeV/c); Q^2 (GeV/c)^2",500,0,3.5,500,0,5);
    TH2F * h_QQ_x        = new TH2F("h_QQ_x","Q^2 vs. x;  Bjorken x; Q^2 (GeV/c)^2",500,0,1,500,0,5);
    
    TLorentzVector beam(0,0,5.498,5.498);
    TLorentzVector target(0,0,0,0.938);
    
    std::cout << "> looping over all events" << std::endl;
    
    //loop over entries
    for(int ien=0; ien<nen; ien++){
        fChain->GetEntry(ien);
        printStatusBar(ien,nen);
        
        if(isElectron(0)){
            for(int ipart=1;ipart<gpart;ipart++){
                int found_pp = -1;
                int found_pm = -1;
                
                if(isPP(ipart)){
                    found_pp = ipart;
                    npp++;
                }
                
                if(isPM(ipart)){
                    found_pm = ipart;
                    npm++;
                }
		
		if (isProton(ipart))
		  np++;

            }
            
            nel++;
            
            TLorentzVector electron(cx[0]*p[0],cy[0]*p[0],cz[0]*p[0],p[0]);
            TLorentzVector virtual_photon = beam - electron;
            TLorentzVector hadronic_state = target + virtual_photon;
    
	    ElasticEvent *event = new ElasticEvent(electron);
	    std::cout << Form(" > params from class ") << event->getQQ() << std::endl;
        
            h_p_transf_sq->Fill(-1*virtual_photon.Mag2());
            h_W->Fill(hadronic_state.Mag());
            h_QQ_w->Fill(hadronic_state.Mag(),-1*virtual_photon.Mag2());
            h_x->Fill(-1*virtual_photon.Mag2()/(2*target.E()*virtual_photon.E()));
            h_QQ_x->Fill(-1*virtual_photon.Mag2()/(2*target.E()*virtual_photon.E()),-1*virtual_photon.Mag2());
        }
    }
    
    
    std::cout << std::endl;
    std::cout << "> summary " << std::endl;
    std::cout << "     > total events: " << nen << std::endl;
    std::cout << "     > electrons   : " << nel << std::endl;
    std::cout << "     > pi+         : " << npp << std::endl;
    std::cout << "     > pi-         : " << npm << std::endl;
    std::cout << "     > protons     : " << np  << std::endl;
    
    TCanvas * c1 = new TCanvas("c1","",800,600);
    c1->Print("pTrans.pdf[");
    h_p_transf_sq->Draw();
    c1->Print("pTrans.pdf");
    
    h_W->Draw();
    c1->Print("pTrans.pdf");
    
    h_x->Draw();
    c1->Print("pTrans.pdf");
    
    h_QQ_w->Draw("colz");
    c1->Print("pTrans.pdf");
    c1->Print("QQvsW.png");
    
    h_QQ_x->Draw("colz");
    c1->Print("pTrans.pdf");
    c1->Print("QQvsx.png");
    
    c1->Print("pTrans.pdf]");
    
    return;
}

// ================================
//   a huge massive hugely
// massive subroutine list follows
//   all are defined in master.h
// ================================

Float_t master::getStartTime(){
    
    // cm/ns
    Float_t speed_of_light = 29.9792458;
    
    // get the run number
    std::string run   = (std::string) fChain->GetCurrentFile()->GetName();
    std::string piece = run.substr(68,6);
    Int_t runno       = atoi(piece.c_str());
    
    Float_t start_time = electronScTimeCorr(1,0,runno) - sc_r[0]/speed_of_light;
    
    return start_time;
}

// =============================================================
//
// =============================================================
Float_t master::getThetaCC(int ipart){
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

// =============================================================
//
// =============================================================

bool master::isProton(int ipart)
{
  if(q[ipart]>0)
    if(isProtonBetaP(ipart))
      return true;

    return false;
}


// =============================================================
//
// =============================================================

bool master::isProtonBetaP(int ipart)
{
    // some things needed to get beta
  std::string file  = (std::string) fChain->GetCurrentFile()->GetName();
  std::string piece = file.substr(68,6);
  Int_t runno       = atoi(piece.c_str());
    
  Float_t speed_of_light = 29.9792458;
  Float_t start_time     = getStartTime();
  Float_t time           = hadronScTimeCorr(!getMCstatus(),ipart,runno) - start_time;
  Float_t beta           = (sc_r[ipart]/time)/speed_of_light;
  Int_t sector           = sc_sect[ipart] -1;

  Float_t cut_max;
  Float_t cut_min;

  if (!getMCstatus())
    {
      cut_max = p[ipart]/sqrt(cv_p_beta_p_upper[sector]*cv_p_beta_p_upper[sector] + p[ipart]*p[ipart]);
      cut_min = p[ipart]/sqrt(cv_p_beta_p_lower[sector]*cv_p_beta_p_lower[sector] + p[ipart]*p[ipart]);
    }

  if (getMCstatus())
    {
      cut_max = p[ipart]/sqrt(cv_p_beta_p_mc_upper[sector]*cv_p_beta_p_mc_upper[sector] + p[ipart]*p[ipart]);
      cut_min = p[ipart]/sqrt(cv_p_beta_p_mc_lower[sector]*cv_p_beta_p_mc_lower[sector] + p[ipart]*p[ipart]);
    }

  if(sector > -1)
    if(beta <= cut_max && beta >= cut_min)
      return true;

  
    return false;
}

// =============================================================
//
// =============================================================

void master::runProtonBeta()
{
    std::cout << "\n > inside of runProtonBeta() \n" << std::endl;
    
    
    int nen = fChain->GetEntries();
    
    // our canvas
    TCanvas * c1 = new TCanvas("c1","",800,600);
    c1->Print("protonBeta.pdf[");
    
    
    Float_t p_bin_size = (p_max_proton_beta - p_min_proton_beta)/(n_p_bins_proton_beta -1);
    Float_t bins[n_p_bins_proton_beta];
    
    
    // load bins of momentum
    for(int ibin=0; ibin<n_p_bins_proton_beta; ibin++){
        bins[ibin] = ibin*p_bin_size + p_min_proton_beta;
    }
    
    if(debug){
        std::cout << "> binning structure: " << std::endl;
        std::cout << Form("     > number of p bins: %d",n_p_bins_proton_beta) << std::endl;
        std::cout << Form("     > bin size        : %f",p_bin_size) << std::endl;
        std::cout << Form("     > first bin       : %f",bins[0]) << std::endl;
        std::cout << Form("     > last bin        : %f",bins[n_p_bins_proton_beta -1]) << std::endl;
    }
    
    // [sector][bin]
    // note: in situations such as this, always use the scheme which allows you to load a[i][j] -> a[i][j+1]
    // and not a[i][j] -> a[i+1][j] because higher dim arrays are still stored in c++ as 1d arrays and the
    // elements a[i][j] and a[i][j+1] are next to each other but a[i][j] and a[i+1][j] are not.
    TH2F * h_p_beta[6];
    TH2F * h_p_beta_all_sectors = new TH2F("h_p_beta_all_sectors","Beta vs. P for all sectors",500,0,5,500,0.55,1.1);
    
    TH1F * h_p_beta_slice[6][n_p_bins_pion_beta];
    TF1  * f_p_beta_slice[6][n_p_bins_pion_beta];
    
    
    // initializing our histograms as well as our Gaussian fits
    // with unique fit boundaries depending on momentum
    for(int ihist=0; ihist<6; ihist++){
        h_p_beta[ihist] = new TH2F(Form("h_p_beta_%d",ihist),Form("beta vs. momentum sector %d",ihist),500,0,5,500,0.55,1.1);
        
        for(int ibin=0; ibin<n_p_bins_proton_beta; ibin++){
            h_p_beta_slice[ihist][ibin] = new TH1F(Form("h_p_beta_slice_%d_%d",ihist,ibin),Form("sector %d bin %d",ihist,ibin),100,0.4,1.1);
        }
    }
    
    std::cout << "> histograms initialized successfully" << std::endl;
    std::cout << "> now filling histograms" << std::endl;
    
    // cm ns^-1
    Float_t speed_of_light = 29.9792458;
    
    // fill histograms
    for(int ien=0; ien<nen; ien++){
        fChain->GetEntry(ien);
        printStatusBar(ien,nen);
        
        // some things needed to get beta
        std::string file  = (std::string) fChain->GetCurrentFile()->GetName();
        std::string piece = file.substr(68,6);
        Int_t runno       = atoi(piece.c_str());
        
        Float_t start_time = getStartTime();
        
        for(int ipart=1; ipart<gpart; ipart++){
            int sector   = ec_sect[ipart] -1;
            Float_t time = hadronScTimeCorr(!getMCstatus(),ipart,runno) - start_time;
            Float_t beta = (sc_r[ipart]/time)/speed_of_light;
            
            if(isPos(ipart))
                if(sc_sect[ipart]>0)
                    if(sector > -1){
                        h_p_beta[sector]->Fill(p[ipart],beta);
                        h_p_beta_all_sectors->Fill(p[ipart],beta);
                        
                        if(floor((p[ipart]-p_min_proton_beta)/p_bin_size) >= 0 && floor((p[ipart]-p_min_proton_beta)/p_bin_size) <= (n_p_bins_proton_beta-1)){
                            h_p_beta_slice[sector][(Int_t) floor((p[ipart]-p_min_proton_beta)/p_bin_size)]->Fill(beta);
                        }
                    }
            
        }
    }
    
    
    Float_t cv_p_beta_min[6];
    Float_t cv_p_beta_max[6];
    
    std::cout << std::endl;
    std::cout << "> fitting histogram slices in p" << std::endl;
    
    for(int isect=0; isect<6; isect++){
        for(int ibin=0; ibin<n_p_bins_proton_beta; ibin++){
            Float_t our_mom = p_min_proton_beta + ibin*p_bin_size;
            Float_t mu      = our_mom/sqrt(0.938*0.938 + our_mom*our_mom);
            
            f_p_beta_slice[isect][ibin] = new TF1(Form("f_p_beta_slice_%d_%d",isect,ibin),"gaus",mu-0.03,mu+0.03);
            h_p_beta_slice[isect][ibin]->Fit(Form("f_p_beta_slice_%d_%d",isect,ibin),"Rq");
        }
    }
    
    
    TGraphErrors * g_mean_p[6];
    TGraphErrors * g_stddev_p[6];
    TF1 * f_mean_p[6];
    TF1 * f_stddev_p[6];
    
    Float_t mean_p[n_p_bins_proton_beta];
    Float_t stddev_p[n_p_bins_proton_beta];
    Float_t mean_err_p[n_p_bins_proton_beta];
    Float_t stddev_err_p[n_p_bins_proton_beta];
    
    Float_t empty_bins[n_p_bins_proton_beta];
    
    for(int isect=0; isect<6; isect++){
        for(int ibin=0; ibin<n_p_bins_pion_beta; ibin++){
            empty_bins[ibin]   = 0.00;
            mean_p[ibin]       = f_p_beta_slice[isect][ibin]->GetParameter(1);
            stddev_p[ibin]     = f_p_beta_slice[isect][ibin]->GetParameter(2);
            mean_err_p[ibin]   = f_p_beta_slice[isect][ibin]->GetParError(1);
            stddev_err_p[ibin] = f_p_beta_slice[isect][ibin]->GetParError(2);
          
        }
        
        g_mean_p[isect]   = new TGraphErrors(n_p_bins_proton_beta,bins,mean_p,empty_bins,mean_err_p);
        g_stddev_p[isect] = new TGraphErrors(n_p_bins_proton_beta,bins,stddev_p,empty_bins,stddev_err_p);
        f_mean_p[isect]   = new TF1(Form("f_mean_p_%d",isect),"x/sqrt(x*x + [0])",0.5,2.5);
        f_stddev_p[isect] = new TF1(Form("f_stddev_p_%d",isect),"pol2",0.5,2.5);

        // try to help the fitter
        f_mean_p[isect]->SetParameter(0,0.938*0.938);
        
        g_mean_p[isect]->Fit(Form("f_mean_p_%d",isect),"Frq");
        g_stddev_p[isect]->Fit(Form("f_stddev_p_%d",isect),"Frq");
        
    }
    
    // sector, param, bin
    Float_t cv_p_upper[6][3];
    Float_t cv_p_lower[6][3];
    Float_t temp_upper_p[n_p_bins_proton_beta];
    Float_t temp_lower_p[n_p_bins_proton_beta];
    Float_t p_err[n_p_bins_proton_beta];
    TGraphErrors * upper_p[6];
    TGraphErrors * lower_p[6];
    TF1 * f_upper_p[6];
    TF1 * f_lower_p[6];
    
    Float_t n_sigma_p_top = 4.0;
    Float_t n_sigma_p_bot = 4.0;
    
    for(int isect=0; isect<6; isect++){
        for(int ibin=0; ibin<n_p_bins_proton_beta; ibin++){
            temp_upper_p[ibin] =  bins[ibin]/sqrt(f_mean_p[isect]->GetParameter(0) +bins[ibin]*bins[ibin]) + n_sigma_p_top*f_p_beta_slice[isect][ibin]->GetParameter(2);
            temp_lower_p[ibin] =  bins[ibin]/sqrt(f_mean_p[isect]->GetParameter(0) +bins[ibin]*bins[ibin]) - n_sigma_p_bot*f_p_beta_slice[isect][ibin]->GetParameter(2);
            p_err[ibin]        = f_p_beta_slice[isect][ibin]->GetParError(1);
        }
        
        upper_p[isect]   = new TGraphErrors(n_p_bins_proton_beta,bins,temp_upper_p,empty_bins,p_err);
        lower_p[isect]   = new TGraphErrors(n_p_bins_proton_beta,bins,temp_lower_p,empty_bins,p_err);
     
        f_upper_p[isect] = new TF1(Form("f_upper_p_%d",isect),"x/sqrt(x*x + [0])",0.5,3.0);
        f_lower_p[isect] = new TF1(Form("f_lower_p_%d",isect),"x/sqrt(x*x + [0])",0.5,3.0);
     
        upper_p[isect]->Fit(Form("f_upper_p_%d",isect),"Frq");
        lower_p[isect]->Fit(Form("f_lower_p_%d",isect),"Frq");
     
    }
    
    h_p_beta_all_sectors->Draw("colz");
    c1->Print("protonBetaAllSectors.png");
    c1->Clear();
    
    // 2d histograms for beta vs p
    // SHOWN WITH FITS
    c1->Divide(3,2);
    c1->cd(1);
    gPad->SetLogz();
    h_p_beta[0]->Draw("colz");
    upper_p[0]->Draw("same");
    lower_p[0]->Draw("same");
    
    c1->cd(2);
    gPad->SetLogz();
    h_p_beta[1]->Draw("colz");
    upper_p[1]->Draw("same");
    lower_p[1]->Draw("same");
    
    c1->cd(3);
    gPad->SetLogz();
    h_p_beta[2]->Draw("colz");
    upper_p[2]->Draw("same");
    lower_p[2]->Draw("same");
    
    c1->cd(4);
    gPad->SetLogz();
    h_p_beta[3]->Draw("colz");
    upper_p[3]->Draw("same");
    lower_p[3]->Draw("same");
    
    c1->cd(5);
    gPad->SetLogz();
    h_p_beta[4]->Draw("colz");
    upper_p[4]->Draw("same");
    lower_p[4]->Draw("same");
    
    c1->cd(6);
    gPad->SetLogz();
    h_p_beta[5]->Draw("colz");
    upper_p[5]->Draw("same");
    lower_p[5]->Draw("same");
    
    c1->Print("protonBeta.pdf");
    c1->Print("protonBeta.png");
    c1->Clear();
    
    c1->Divide(3,2);
    for(int ibin=0; ibin<n_p_bins_proton_beta; ibin++){
        c1->cd(1);
        h_p_beta_slice[0][ibin]->Draw("same");
        c1->cd(2);
        h_p_beta_slice[1][ibin]->Draw("same");
        c1->cd(3);
        h_p_beta_slice[2][ibin]->Draw("same");
        c1->cd(4);
        h_p_beta_slice[3][ibin]->Draw("same");
        c1->cd(5);
        h_p_beta_slice[4][ibin]->Draw("same");
        c1->cd(6);
        h_p_beta_slice[5][ibin]->Draw("same");
        
        c1->Print("protonBeta.pdf");
        c1->Clear();
        c1->Divide(3,2);
    }
    
    for(int isect=0; isect<6; isect++){
        c1->cd(isect+1);
        g_mean_p[isect]->Draw("same");
    }
    
    c1->Print("protonBeta.pdf");
    c1->Print("protonBeta.pdf]");
    
    // work is required here
    // re-write the method below to output parameters of the fits for pion beta
    // then fix the PID algorithms which depend on said parameters
    
    
    // write it out to the header file
    std::ofstream ofile;
    std::string filename;
    
    filename = "cv_proton_beta.h";
    
    if ( getMCstatus() )
        filename = "cv_proton_beta_MC.h";
    
    ofile.open(filename.c_str(), std::ios::out | std::ios::trunc);
    
    if (!getMCstatus())
      ofile << "Float_t cv_p_beta_p_upper[6] = {";
    
    if (getMCstatus())
      ofile << "Float_t cv_p_beta_p_mc_upper[6] = {";

    for(int isect=0; isect<6; isect++){
        ofile << f_upper_p[isect]->GetParameter(0);
        if(isect<5)
            ofile << ",";
    }
    
    ofile << "}; \n";
    
    if (!getMCstatus())
      ofile << "Float_t cv_p_beta_p_lower[6] = {";
    
    if (getMCstatus())
      ofile << "Float_t cv_p_beta_p_mc_lower[6] = {";
     
    for(int isect=0; isect<6; isect++){
        ofile << f_lower_p[isect]->GetParameter(0);
        if(isect<5)
            ofile << ",";
    }
    
    ofile << "}; \n";
    
    ofile.close();
    
    
}

// =============================================================
//
// =============================================================

bool master::ecGeometricCut(int ipart, int strict){
    
    Float_t x = ech_x[ipart];
    Float_t y = ech_y[ipart];
    Float_t z = ech_z[ipart];
    
    float uMin[5] = {58, 64, 70, 76, 82}; // loosest, loose, nominal, tight, tightest
    float uMax[5] = {412, 406, 400, 394, 388};
    float vMax[5] = {374, 368, 362, 356, 350};
    float wMax[5] = {407, 401, 395, 389, 383};
    
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
    
    if(u >= uMin[strict] && u <= uMax[strict] && v <= vMax[strict] && w <= wMax[strict])
        return true;
    
    return false;
}

// =============================================================
//
// =============================================================

bool master::nPhotoElectrons(int ipart){
    
    if(nphe[ipart]>25)
        return true;
    
    return false;
}

// =============================================================
//
// =============================================================

Float_t  master::getCorrZ(int ipart)
{
    Int_t s = dc_sect[ipart] -1;
    
    Float_t s0, sp, sv;
    Float_t px, py, pz;
    
    px = p[ipart]*cx[ipart];
    py = p[ipart]*cy[ipart];
    pz = p[ipart]*cz[ipart];
    
    Float_t n[3][6];
    
    for(int abc = 0; abc < 3; abc++) // initialize to zero
    {
        for(int def = 0; def < 6; def++)
        {
            n[abc][def] = 0.0;
        }
    }
    
    n[0][0] = 1.0;
    n[1][0] = 0.0;
    
    n[0][1] = 0.5;
    n[1][1] = 0.866025388;
    
    n[0][2] = -0.5;
    n[1][2] = 0.866025388;
    
    n[0][3] = -1.0;
    n[1][3] = 0.0;
    
    n[0][4] = -0.5;
    n[1][4] = -0.866025388;
    
    n[0][5] = 0.5;
    n[1][5] = -0.866025388;
    
    Float_t x0, y0, z0; // beam position (cm)
    x0 = 0.15;
    y0 = -0.25;
    z0 = 0.0;
    
    Float_t A;
    
    s0 = x0*n[0][s] + y0*n[1][s] + z0*n[2][s];
    sp = px*n[0][s] + py*n[1][s] + pz*n[2][s];
    sv = vx[ipart]*n[0][s] + vy[ipart]*n[1][s] + vz[ipart]*n[2][s];
    
    Float_t cvz;
    
    if(fabs(sp) > 0.0000000001)
    {
        A = (s0-sv)/sp;
        cvz = vz[ipart] + A*pz;
    }
    else
    {
        cvz = vz[ipart];
    }
    
    return cvz;
}

// =============================================================
//
// =============================================================

Float_t master::getRelPhi(int ipart){
    
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

// =============================================================
//
// =============================================================

Int_t master::ccPhiMatching(int ipart){
    Int_t ccpmt = cc_segm[ipart]/1000 - 1; // -1 left pmt; +1 right pmt; 0 both
    Float_t relphi = getRelPhi(ipart);
    
    if(relphi > 0 && ccpmt > 0) return 1;
    if(relphi > 0 && ccpmt < 0) return 2;
    if(relphi < 0 && ccpmt < 0) return -1;
    if(relphi < 0 && ccpmt > 0) return -2;
    if(ccpmt == 0 || relphi == 0) return 0;
    
    return 0;
}

// =============================================================
//
// =============================================================

bool master::phiCCcut(int ipart){
    if(ccPhiMatching(ipart)>-2 && ccPhiMatching(ipart)<2)
        return true;
    
    return false;
}

// =============================================================
//
// =============================================================

bool master::isVertex(int ipart, int strict){
    
    Float_t leftcut[5]  = {-27.7302 - 0.6, -27.7302 - 0.3, -27.7302, -27.7302 + 0.3, -27.7302 + 0.6};
    Float_t rightcut[5] = {-22.6864 + 0.6, -22.6864 + 0.3, -22.6864, -22.6864 - 0.3, -22.6864 - 0.6};
    
    Float_t vertex = getCorrZ(ipart);
    
    if(vertex > leftcut[strict] && vertex < rightcut[strict] )
        return true;
    
    return false;
}

// =============================================================
//
// =============================================================

bool master::thetaCCcut(int ipart, int strict){
    Float_t theta_cc = getThetaCC(ipart);
    Int_t cc_segment = (cc_segm[ipart]%1000)/10 -1;
    Int_t cc_sector  = cc_sect[ipart] -1;
    
    if(cc_segment == -1 || cc_sector == -1)
        return false;
    
    if(!getMCstatus())
    if(theta_cc <= cut_theta_cc_max[cc_segment][cc_sector][strict])
        if(theta_cc >= cut_theta_cc_min[cc_segment][cc_sector][strict])
            return true;
    
    // temporary until writing of mc into eID code
    if(getMCstatus())
        return true;
    
    return false;
}

// =============================================================
//
// =============================================================

bool master::eDepEC(int ipart, int strict){
    int ec_sector = ec_sect[ipart] -1;
    
    if(etot[ipart]/p[ipart] <= (ec_sampling_max[2][ec_sector][strict]*p[ipart]*p[ipart] + ec_sampling_max[1][ec_sector][strict]*p[ipart] + ec_sampling_max[0][ec_sector][strict]))
        if(etot[ipart]/p[ipart] >= (ec_sampling_min[2][ec_sector][strict]*p[ipart]*p[ipart] + ec_sampling_min[1][ec_sector][strict]*p[ipart] + ec_sampling_min[0][ec_sector][strict]))
            return true;
    
    return false;
}

// =============================================================
//
// =============================================================

bool master::isNeg(int ipart){
    if(q[ipart]<0)
        return true;
    
    return false;
}

// =============================================================
//
// =============================================================

bool master::isElectron(int ipart){
   
    if(isNeg(ipart))
        if(nPhotoElectrons(ipart))
            if(eDepEC(ipart,eid_stricts[2]))
                if(isVertex(ipart,eid_stricts[4]))
                    if(thetaCCcut(ipart,eid_stricts[3]))
                        if(ecGeometricCut(ipart,eid_stricts[1]))
                            if(ecInnerCut(ipart,eid_stricts[0]))
                                if(phiCCcut(ipart))
                                    if(ec_sect[ipart]>0)
                                        if(sc_sect[ipart]>0)
                                            return true;
    
    return false;
}

// =============================================================
//
// =============================================================

bool master::ecInnerCut(int ipart, int strict){
    
    Float_t bound[5] = {0.08, 0.07, 0.06, 0.05, 0.04};
    
    if(ec_ei[ipart]>bound[strict])
        return true;
    
    return false;
}

// =============================================================
//
// =============================================================

bool master::isPos(int ipart){
    if(q[ipart]>0)
        return true;
    
    return false;
}

// =============================================================
//
// =============================================================

void master::eidValidation(){
    
    std::cout << std::endl;
    std::cout << "inside eidValidation()" << std::endl;
    
    int nen = fChain->GetEntries();
    
    // [with cut number][histogram number][strictness]
    TH2F * histos[6][4][5];
    std::string histo_names[4] = {"ec_ei_ec_eo","ec_samp","ec_geo","cc_theta"};
    std::string cut_names[6]   = {"ecal_inner_cut","ec_samp_cut","ec_geo_cut","cc_theta_cut","cc_phi_cut","nphe_cut"};
    
    TH1F * vertex_before[6];
    TH1F * vertex_after[6];
    
    for(int isect=0; isect<6; isect++){
        vertex_before[isect] = new TH1F(Form("vertex_before_%d",isect),"vertex before corrections; distance (cm)",100,-35,-15);
        vertex_after[isect]  = new TH1F(Form("vertex_after_%d",isect),"vertex after corrections; distance (cm)",100,-35,-15);
    }
    
    
    Int_t ec_sector  = 0;
    Int_t cc_segment = 0;
    Float_t theta_cc = 0.0;
    
    for(int istrict=0; istrict<5; istrict++){
        for(int icut=0; icut<6; icut++){
            histos[icut][0][istrict] = new TH2F(Form("%s_%s_strict_%d",histo_names[0].c_str(),cut_names[icut].c_str(),istrict),Form("histo %s cut %s strict %d",histo_names[0].c_str(),cut_names[icut].c_str(),istrict),100,0,0.4,100,0,0.4);
            histos[icut][1][istrict] = new TH2F(Form("%s_%s_strict_%d",histo_names[1].c_str(),cut_names[icut].c_str(),istrict),Form("histo %s cut %s strict %d",histo_names[1].c_str(),cut_names[icut].c_str(),istrict),100,0,5,100,0,0.5);
            histos[icut][2][istrict] = new TH2F(Form("%s_%s_strict_%d",histo_names[2].c_str(),cut_names[icut].c_str(),istrict),Form("histo %s cut %s strict %d",histo_names[2].c_str(),cut_names[icut].c_str(),istrict),100,-500,500,100,-500,500);
            histos[icut][3][istrict] = new TH2F(Form("%s_%s_strict_%d",histo_names[3].c_str(),cut_names[icut].c_str(),istrict),Form("histo %s cut %s strict %d",histo_names[3].c_str(),cut_names[icut].c_str(),istrict),18,0,17,100,0,50);
        }
    }
    
    if(debug)
        std::cout << "eID validation histograms initialized" << std::endl << std::endl;
    
    std::cout << std::endl;
    std::cout << "> loading all cut/strictness histograms" << std::endl;
    
    // major loop
    for(int jstrict=0; jstrict<n_strict; jstrict++){
        for(int jcut=0; jcut<6; jcut++){
            std::cout << Form("     > filling histo for strict %d cut %d",jstrict,jcut) << std::endl;
            
            for(int ien=0; ien<nen; ien++){
                fChain->GetEntry(ien);
                
                for(int ipart=0; ipart<gpart; ipart++){
                    if(q[ipart]<0){
                        cc_segment = (cc_segm[ipart]%1000)/10 -1;
                        ec_sector  = ec_sect[ipart] -1;
                        theta_cc   = getThetaCC(ipart);
                        
                        
                        if(ec_sector != -1){
                            vertex_before[ec_sector]->Fill(vz[ipart]);
                            vertex_after[ec_sector]->Fill(getCorrZ(ipart));
                        }
                        
                        
                        switch(jcut){
                            case 0:
                                if(ecInnerCut(ipart,jstrict)){
                                    if(ec_sector != -1){
                                        if(ec_eo[ipart]>0.01 && ec_ei[ipart]>0.01)
                                            histos[0][0][jstrict]->Fill(ec_ei[ipart],ec_eo[ipart]);
                                        histos[0][1][jstrict]->Fill(p[ipart],etot[ipart]/p[ipart]);
                                        histos[0][2][jstrict]->Fill(ech_x[ipart],ech_y[ipart]);
                                    }
                                    if(cc_segment != -1)
                                        histos[0][3][jstrict]->Fill(cc_segment,theta_cc);
                                }
                                break;
                                
                            case 1:
                                if(eDepEC(ipart,jstrict)){
                                    if(ec_sector != -1){
                                        if(ec_eo[ipart]>0.01 && ec_ei[ipart]>0.01)
                                            histos[1][0][jstrict]->Fill(ec_ei[ipart],ec_eo[ipart]);
                                        histos[1][1][jstrict]->Fill(p[ipart],etot[ipart]/p[ipart]);
                                        histos[1][2][jstrict]->Fill(ech_x[ipart],ech_y[ipart]);
                                    }
                                    if(cc_segment != -1)
                                        histos[1][3][jstrict]->Fill(cc_segment,theta_cc);
                                }
                                break;
                                
                            case 2:
                                if(ecGeometricCut(ipart,jstrict)){
                                    if(ec_sector != -1){
                                        if(ec_eo[ipart]>0.01 && ec_ei[ipart]>0.01)
                                            histos[2][0][jstrict]->Fill(ec_ei[ipart],ec_eo[ipart]);
                                        histos[2][1][jstrict]->Fill(p[ipart],etot[ipart]/p[ipart]);
                                        histos[2][2][jstrict]->Fill(ech_x[ipart],ech_y[ipart]);
                                    }
                                    if(cc_segment != -1)
                                        histos[2][3][jstrict]->Fill(cc_segment,theta_cc);
                                }
                                break;
                                
                            case 3:
                                if(thetaCCcut(ipart,jstrict)){
                                    if(ec_sector != -1){
                                        if(ec_eo[ipart]>0.01 && ec_ei[ipart]>0.01)
                                            histos[3][0][jstrict]->Fill(ec_ei[ipart],ec_eo[ipart]);
                                        histos[3][1][jstrict]->Fill(p[ipart],etot[ipart]/p[ipart]);
                                        histos[3][2][jstrict]->Fill(ech_x[ipart],ech_y[ipart]);
                                    }
                                    if(cc_segment != -1)
                                        histos[3][3][jstrict]->Fill(cc_segment,theta_cc);
                                }
                                break;
                                
                            case 4:
                                if(phiCCcut(ipart)){
                                    if(ec_sector != -1){
                                        if(ec_eo[ipart]>0.01 && ec_ei[ipart]>0.01)
                                            histos[4][0][jstrict]->Fill(ec_ei[ipart],ec_eo[ipart]);
                                        histos[4][1][jstrict]->Fill(p[ipart],etot[ipart]/p[ipart]);
                                        histos[4][2][jstrict]->Fill(ech_x[ipart],ech_y[ipart]);
                                    }
                                    if(cc_segment != -1)
                                        histos[4][3][jstrict]->Fill(cc_segment,theta_cc);
                                }
                                break;
                                
                            case 5:
                                if(nPhotoElectrons(ipart)){
                                    if(ec_sector != -1){
                                        if(ec_eo[ipart]>0.01 && ec_ei[ipart]>0.01)
                                            histos[5][0][jstrict]->Fill(ec_ei[ipart],ec_eo[ipart]);
                                        histos[5][1][jstrict]->Fill(p[ipart],etot[ipart]/p[ipart]);
                                        histos[5][2][jstrict]->Fill(ech_x[ipart],ech_y[ipart]);
                                    }
                                    if(cc_segment != -1)
                                        histos[5][3][jstrict]->Fill(cc_segment,theta_cc);
                                }
                                break;
                        } // end switch statement
                    }
                }
            }
        }
    }
    
    std::cout << std::endl;
    
    // print out raw distributions
    TH2F * raw[4];
    
    raw[0] = new TH2F("raw_ec_io","raw_ec_io;Energy Dep Inner EC (GeV); Energy Dep Outer EC (GeV)",500,0,0.4,500,0,0.4);
    raw[1] = new TH2F("raw_ec_samp","raw_ec_samp;Momentum (GeV/c); etot/p",500,0,5,500,0,0.5);
    raw[2] = new TH2F("raw_ec_geo","raw_ec_geo;x;y",1000,-500,500,1000,-500,500);
    raw[3] = new TH2F("raw_theta_cc","raw_theta_cc;cc segment; theta cc",18,0,17,100,0,50);
    
    std::cout << "> loading histograms for raw distrobutions" << std::endl;
    
    for(int jen=0; jen<nen; jen++){
        fChain->GetEntry(jen);
        printStatusBar(jen,nen);
        
        for(int jpart=0; jpart<gpart; jpart++){
            theta_cc   = getThetaCC(jpart);
            cc_segment = (cc_segm[jpart]%1000)/10 -1;
            ec_sector  = ec_sect[jpart] -1;
            
            if(q[jpart]<0){
                if(ec_sector != -1){
                    if(ec_eo[jpart]>0.005 && ec_ei[jpart]>0.005)
                        raw[0]->Fill(ec_ei[jpart],ec_eo[jpart]);
                    raw[1]->Fill(p[jpart],etot[jpart]/p[jpart]);
                    raw[2]->Fill(ech_x[jpart],ech_y[jpart]);
                }
                if(cc_segment != -1)
                    raw[3]->Fill(cc_segment,theta_cc);
            }
        }
    }
    
    std::cout << std::endl;
    
    // final distributions after eID cuts
    TH2F * final[4];
    
    final[0] = new TH2F("final_ec_io","final_ec_io;Energy Dep Inner EC (GeV); Energy Dep Outer EC (GeV)",500,0,0.4,500,0,0.4);
    final[1] = new TH2F("final_ec_samp","final_ec_samp;momentum (GeV/c);etot/p",500,0,5,500,0,0.5);
    final[2] = new TH2F("final_ec_geo","final_ec_geo;x;y",1000,-500,500,1000,-500,500);
    final[3] = new TH2F("final_theta_cc","final_theta_cc;cc segment;theta cc",18,0,17,100,0,50);
    
    std::cout << "> loading histograms for final distributions" << std::endl;
    
    // load up final distributions
    for(int ken=0; ken<nen; ken++){
        fChain->GetEntry(ken);
        printStatusBar(ken,nen);
        
        for(int kpart=0; kpart<gpart; kpart++){
            theta_cc   = getThetaCC(kpart);
            cc_segment = (cc_segm[kpart]%1000)/10 -1;
            ec_sector  = ec_sect[kpart] -1;
            
            if(isElectron(kpart)){
                if(ec_eo[kpart]>0.005 && ec_ei[kpart]>0.005)
                    final[0]->Fill(ec_ei[kpart],ec_eo[kpart]);
                final[1]->Fill(p[kpart],etot[kpart]/p[kpart]);
                final[2]->Fill(ech_x[kpart],ech_y[kpart]);
                final[3]->Fill(cc_segment,theta_cc);
            }
        }
    }
    
    std::cout << std::endl;
    std::cout << "> authoring eID.pdf file" << std::endl;
    
    TCanvas * c1 = new TCanvas("c1"," ",1100,800);
    c1->Divide(2,2);
    
    c1->Print("eID.pdf[");
    
    c1->cd(1);
    raw[0]->Draw("same colz");
    c1->cd(2);
    raw[1]->Draw("same colz");
    c1->cd(3);
    raw[2]->Draw("same colz");
    c1->cd(4);
    raw[3]->Draw("same colz");
    
    c1->Print("eID.pdf");
    c1->Clear();
    
    c1->Divide(2,2);
    
    for(int kstrict=0; kstrict<5; kstrict++){
        for(int kcut=0; kcut<6; kcut++){
            c1->cd(1);
            histos[kcut][0][kstrict]->Draw("same colz");
            c1->cd(2);
            histos[kcut][1][kstrict]->Draw("same colz");
            c1->cd(3);
            histos[kcut][2][kstrict]->Draw("same colz");
            c1->cd(4);
            histos[kcut][3][kstrict]->Draw("same colz");
            
            c1->Print(Form("eIDvalidation_%d_%d.png",kcut,kstrict));
            c1->Print("eID.pdf");
            c1->Clear();
            c1->Divide(2,2);
        }
    }
    
    
    c1->Clear();
    
    raw[0]->Draw("colz");
    c1->Print("eID.pdf");
    c1->Print("raw_ecio.png");
    
    final[0]->Draw("colz");
    c1->Print("eID.pdf");
    c1->Print("final_ecio.png");
    
    raw[1]->Draw("colz");
    c1->Print("eID.pdf");
    c1->Print("raw_ec_samp.png");
    
    final[1]->Draw("colz");
    c1->Print("eID.pdf");
    c1->Print("final_ec_samp.png");
    
    raw[2]->Draw("colz");
    c1->Print("eID.pdf");
    c1->Print("raw_ec_xy.png");
    
    final[2]->Draw("colz");
    c1->Print("eID.pdf");
    c1->Print("final_ec_xy.png");
    
    raw[3]->Draw("colz");
    c1->Print("eID.pdf");
    c1->Print("raw_cc_theta.png");
    
    final[3]->Draw("colz");
    c1->Print("eID.pdf");
    c1->Print("final_cc_theta.png");
    
    
    c1->Clear();
    c1->Divide(2,1);
    
    vertex_before[0]->SetLineColor(1);
    vertex_before[1]->SetLineColor(2);
    vertex_before[2]->SetLineColor(3);
    vertex_before[3]->SetLineColor(4);
    vertex_before[4]->SetLineColor(5);
    vertex_before[5]->SetLineColor(6);
    
    vertex_after[0]->SetLineColor(1);
    vertex_after[1]->SetLineColor(2);
    vertex_after[2]->SetLineColor(3);
    vertex_after[3]->SetLineColor(4);
    vertex_after[4]->SetLineColor(5);
    vertex_after[5]->SetLineColor(6);
    
    c1->cd(1);
    vertex_before[0]->Draw("same");
    vertex_before[1]->Draw("same");
    vertex_before[2]->Draw("same");
    vertex_before[3]->Draw("same");
    vertex_before[4]->Draw("same");
    vertex_before[5]->Draw("same");
    
    c1->cd(2);
    vertex_after[0]->Draw("same");
    vertex_after[1]->Draw("same");
    vertex_after[2]->Draw("same");
    vertex_after[3]->Draw("same");
    vertex_after[4]->Draw("same");
    vertex_after[5]->Draw("same");
    
    c1->Print("vertex.png");
    c1->Print("eID.pdf");
    
    
    c1->Print("eID.pdf]");
    
    
    return;
}

// =============================================================
//
// =============================================================

void master::showPhiCC(){
    
    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    // debug section for phi cc vs pmt
    
    Int_t nen = fChain->GetEntries();
    
    Int_t debug_pmt;
    Float_t debug_phi;
    
    TH1F * h_phi_cc[3];
    
    // init debug histos
    for(int idebug=0; idebug<3; idebug++){
        h_phi_cc[idebug] = new TH1F(Form("h_phi_cc_%d",idebug),Form("h_phi_cc_%d",idebug),1000,-30,330);
    }
    
    // fill debug histos
    for(int jdebug=0; jdebug<nen; jdebug++){
        fChain->GetEntry(jdebug);
        
        for(int kdebug=0; kdebug<gpart; kdebug++){
            debug_pmt = cc_segm[kdebug]/1000;
            debug_phi = getRelPhi(kdebug);
            h_phi_cc[debug_pmt]->Fill(debug_phi);
        }
    }
    
    // define canvas and print debug histograms to .pdf
    TCanvas * debug_canvas = new TCanvas("debug_canvas","",1100,800);
    
    debug_canvas->Print("debug_phi_cc.pdf[");
    
    for(int ldebug=0; ldebug<3; ldebug++){
        h_phi_cc[ldebug]->Draw();
        debug_canvas->Print("debug_phi_cc.pdf");
    }
    
    debug_canvas->Print("debug_phi_cc.pdf]");
    
    // end debug section
    // -=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
}

// =============================================================
//
// =============================================================

void master::ppMissingMass(){
    
    std::cout << std::endl;
    std::cout << "inside of ppMissingMass()" << std::endl;
    std::cout << std::endl;
    
    
    // get number of entries
    int nen = fChain->GetEntries();
    
    TLorentzVector beam(0,0,5.498,5.498);
    TLorentzVector target(0,0,0,0.938);
    
    TH1F * missing_mass = new TH1F("missing_mass","missing mass {e p -> e pi+ X};Mass (GeV/c^2);Number of Occurances",100,0,5);
    
    // some variables needed for calculating missing mass
    int pos_part_num[40];
    int e_part_num[40];
    
    std::cout << "> looping over all entries to fill histogram" << std::endl;
    
    // loop over all entries and particles
    for(int ien=0; ien<nen; ien++){
        fChain->GetEntry(ien);
        printStatusBar(ien,nen);
        
        int nel = 0, electron_number = 0, npos = 0;
        
        for(int ipart=0; ipart<gpart; ipart++){
            if(isElectron(ipart)){
                e_part_num[nel] = ipart;
                nel++;
            }
            
            if(isPos(ipart)){
                pos_part_num[npos] = ipart;
                npos++;
            }
            
        } // end ipart loop
        
        // find the fastest electron candidate
        if(nel == 1){
            electron_number = e_part_num[0];
        }
        else if(nel > 1){
            Float_t b_max = 0.00;
            
            for(int iel=0; iel<nel; iel++){
                if (b[e_part_num[iel]]>b_max){
                    electron_number = iel;
                    b_max = b[e_part_num[iel]];
                }
            }
        }
        
        if(nel != 0 && npos > 0){
            TLorentzVector electron(p[electron_number]*cx[electron_number],p[electron_number]*cy[electron_number],p[electron_number]*cz[electron_number],p[electron_number]);
            for(int ipos=0; ipos<npos; ipos++){
                TLorentzVector pp_candidate(p[pos_part_num[ipos]]*cx[pos_part_num[ipos]],p[pos_part_num[ipos]]*cy[pos_part_num[ipos]],p[pos_part_num[ipos]]*cz[pos_part_num[ipos]],sqrt(p[pos_part_num[ipos]]*p[pos_part_num[ipos]] + 0.135*0.135));
                TLorentzVector missing = beam + target - electron - pp_candidate;
                
                // fill our histo
                if(missing.M2() >= 0)
                    missing_mass->Fill(sqrt(missing.M2()));
            }
        }
    } // end ien loop
    
    
    std::cout << std::endl;
    std::cout << "> histograms filled, authoring ppMissingMass.pdf" << std::endl;
    
    if(missing_mass->GetEntries() != 0){
        TCanvas * c1 = new TCanvas("c1","",1100,800);
        
        c1->Print("ppMissingMass.pdf[");
        missing_mass->Draw();
        c1->Print("ppMissingMass.pdf");
        c1->Print("ppMissingMass.pdf]");
    }
    else
        std::cout << "> histogram not printed, it was empty" << std::endl;
} // end ppMissingMass()

// =============================================================
//
// =============================================================

void master::getGoodRuns(){
    
    std::cout << std::endl;
    std::cout << "inside of getGoodRuns()" << std::endl;
    std::cout << Form("> %d files loaded",fChain->GetListOfFiles()->GetEntries()) << std::endl;
    std::cout << std::endl;
    
    int nen = fChain->GetEntries();
    
    // setup variables needed to do
    // fcup ratio testing
    fChain->GetEntry(0);
    std::string cfile     = fChain->GetCurrentFile()->GetName();
    std::string pfile     = cfile;
    std::string piece     = cfile.substr(63,15);
    int nelec, ifile;
    nelec            = 0;
    ifile            = 0;
    Float_t fca      = 0;
    Float_t coef     = 8.457 * 10e22;                                  // target length * proton number density
    
    // create some storage for my information to be generated
    // of size number of files
    std::string flist[fChain->GetListOfFiles()->GetEntries()];
    Float_t ratio[fChain->GetListOfFiles()->GetEntries()];
    Float_t ilum[fChain->GetListOfFiles()->GetEntries()];
    
    TH1F * hfc  = new TH1F("hfc","nelec/fcup accumulation",100,10,20);   // place to put ratio in a histogram
    
    // loop over all entries in fchain
    for(int ien=0;ien<nen;ien++){
        printStatusBar(ien,nen);
        
        fChain->GetEntry(ien);
        cfile = fChain->GetCurrentFile()->GetName();
        
        // now using full eID algorithm
        if(isElectron(0))
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
            pfile = cfile;
            nelec = 0;
            
            ifile++;
        } // end long if statement
        
    } // end loop over all entries in fchain
    
    std::cout << std::endl;
    
    // output to a pdf file
    TCanvas * c1 = new TCanvas("c1","c1",1100,800);
    c1->Print("fc.pdf[");
    hfc->Draw();
    c1->Print("fc.pdf");
    
    // fit the histogram with a guassian
    TF1 * myfit = new TF1("myfit","gaus");
    
    hfc->Fit("myfit","q");
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
    
    std::ofstream out;
    out.open("goodruns.txt");
    
    for(int jfile=0; jfile<fChain->GetListOfFiles()->GetEntries();jfile++){
        if(ratio[jfile]>min)
            if(ratio[jfile]<max)
                out << flist[jfile] << std::endl;
    }
    
    out.close();
    
    return;
}

// =============================================================
//
// =============================================================

Float_t master::electronScTimeCorr(int ExpOrSim, int ipart, int runno)
{
    
    int sector     = sc_sect[ipart];
    int paddle     = sc_pd[ipart];
    Float_t sctime = sc_t[ipart];
    
    if(ExpOrSim == 0) return sctime;
    
    if(sector == 2 && paddle == 16 && runno >= 37776 && runno <= 38548) return sctime + 0.5;
    if(sector == 2 && paddle == 16 && runno >= 38549) return sctime + 0.9;
    if(sector == 3 && paddle == 11 && runno >= 37777) return sctime - 2.3;
    if(sector == 4 && paddle == 5 && runno >= 38548) return sctime + 2.0;
    if(sector == 5 && paddle == 3 && runno >= 37673 && runno <= 37854) return sctime + 31.0;
    if(sector == 5 && paddle == 3 && runno >= 37855) return sctime - 0.25;
    if(sector == 5 && paddle == 18 && runno >= 38050 && runno <= 38548) return sctime + 1.1;
    if(sector == 5 && paddle == 20 && runno >= 37777) return sctime - 0.5;
    if(sector == 6 && paddle == 18 && runno >= 38050 && runno <= 38548) return sctime - 1.6;
    
    // note: electrons very rarely hit paddle 2, so the below values were copied from the hadron correction function further down (having these vs not having these makes very little difference)
    if(sector == 3 && paddle == 2 && runno >= 0) return sctime - 15.45;
    if(sector == 4 && paddle == 2 && runno <= 37853) return sctime - 1.9;
    if(sector == 5 && paddle == 2 && runno <= 38240) return sctime - 17.9;
    if(sector == 5 && paddle == 2 && runno >= 38241) return sctime - 19.15;
    
    return sctime;
}

// =============================================================
//
// =============================================================

Float_t master::hadronScTimeCorr(int ExpOrSim, int ipart, int runno) // h for hadron
{
    int sector     = sc_sect[ipart];
    int paddle     = sc_pd[ipart];
    Float_t sctime = sc_t[ipart];
    
    if(ExpOrSim == 0) return sctime;
    
    // calibrated using negative tracks: (low paddle numbers)
    if(sector == 6 && paddle == 1 && runno >= 0) return sctime + 18.25;
    if(sector == 3 && paddle == 2 && runno >= 0) return sctime - 15.45;
    if(sector == 4 && paddle == 2 && runno <= 37853) return sctime - 1.9;
    if(sector == 5 && paddle == 2 && runno <= 38240) return sctime - 17.9;
    if(sector == 5 && paddle == 2 && runno >= 38241) return sctime - 19.15;
    if(sector == 5 && paddle == 3 && runno >= 37673 && runno <= 37854) return sctime + 31.0;
    if(sector == 5 && paddle == 3 && runno >= 37855) return sctime - 0.25;
    
    // calibrated using positive tracks:
    if(sector == 1 && paddle == 24 && runno >= 37749) return sctime + 1.13334;
    if(sector == 2 && paddle == 16 && runno >= 37776 && runno <= 38548) return sctime + 0.565033;
    if(sector == 2 && paddle == 16 && runno >= 38549) return sctime + 1.04168;
    if(sector == 2 && paddle == 38 && runno >= 38535) return sctime - 1.89592;
    if(sector == 3 && paddle == 11 && runno >= 37777) return sctime - 2.26126;
    if(sector == 3 && paddle == 24 && runno >= 37855 && runno <= 38546) return sctime - 1.78266;
    if(sector == 3 && paddle == 25 && runno >= 37743 && runno <= 38266) return sctime + 2.44804;
    if(sector == 3 && paddle == 27 && runno >= 37854 && runno <= 38546) return sctime - 1.85815;
    if(sector == 3 && paddle == 28 && runno >= 37854) return sctime + 1.21704;
    if(sector == 4 && paddle == 5 && runno >= 38549) return sctime + 1.91688;
    if(sector == 4 && paddle == 19 && runno >= 37854) return sctime - 0.365798;
    if(sector == 4 && paddle == 34 && runno >= 37854) return sctime - 2.33721;
    if(sector == 4 && paddle == 42 && runno >= 37750) return sctime - 1.4118;
    if(sector == 4 && paddle == 45 && runno >= 38551) return sctime - 3.36406;
    if(sector == 5 && paddle == 18 && runno >= 37854 && runno <= 38545) return sctime + 1.24884;
    if(sector == 5 && paddle == 20 && runno >= 37809) return sctime - 0.468722;
    if(sector == 5 && paddle == 34 && runno <= 37853) return sctime - 1.0;
    if(sector == 5 && paddle == 34 && runno >= 37854) return sctime + 6.0;
    if(sector == 5 && paddle == 36 && runno >= 37748) return sctime + 1.07962;
    if(sector == 6 && paddle == 18 && runno >= 37854 && runno <= 38545) return sctime - 1.69106;
    if(sector == 6 && paddle == 42 && runno >= 37854 && runno <= 38545) return sctime - 6.0;
    
    return sctime;
}

// =============================================================
//
// =============================================================

bool master::goodORbadSCpaddle(int ipart)
{
    int sector = sc_sect[ipart];
    int paddle = sc_pd[ipart];
    
    if(sector == 3 && paddle == 2) return 0;
    if(sector == 5 && paddle == 3) return 0;
    if(sector == 2 && paddle == 16) return 0;
    if(sector == 2 && paddle == 40) return 0;
    if(sector == 3 && paddle == 40) return 0;
    if(sector == 5 && paddle == 40) return 0;
    if(sector == 6 && paddle == 40) return 0;
    if(sector == 1 && paddle == 41) return 0;
    if(sector == 2 && paddle == 41) return 0;
    if(sector == 3 && paddle == 41) return 0;
    if(sector == 5 && paddle == 41) return 0;
    if(sector == 1 && paddle == 42) return 0;
    if(sector == 2 && paddle == 42) return 0;
    if(sector == 3 && paddle == 42) return 0;
    if(sector == 5 && paddle == 42) return 0;
    if(sector == 6 && paddle == 42) return 0;
    if(sector == 2 && paddle == 43) return 0;
    if(sector == 3 && paddle == 43) return 0;
    if(sector == 4 && paddle == 43) return 0;
    if(sector == 5 && paddle == 43) return 0;
    if(sector == 1 && paddle == 44) return 0;
    if(sector == 3 && paddle == 44) return 0;
    if(sector == 5 && paddle == 44) return 0;
    if(sector == 6 && paddle == 44) return 0;
    if(sector == 1 && paddle == 45) return 0;
    if(sector == 2 && paddle == 45) return 0;
    if(sector == 3 && paddle == 45) return 0;
    if(sector == 6 && paddle == 45) return 0;
    if(sector == 1 && paddle == 46) return 0;
    if(sector == 2 && paddle == 46) return 0;
    if(sector == 3 && paddle == 46) return 0;
    if(sector == 4 && paddle == 46) return 0;
    if(sector == 5 && paddle == 46) return 0;
    if(sector == 1 && paddle == 47) return 0;
    if(sector == 5 && paddle == 47) return 0;
    
    return 1;
}

// not being used for now

// =============================================================
//
// =============================================================

void master::ppTimeCorrelation(){
    
    int nen = fChain->GetEntries();
    
    std::cout << std::endl;
    std::cout << "inside of ppTimeCorrelation()" << std::endl;
    std::cout << Form("> running time correlation for %d entries",nen) << std::endl;
    
    TH2F * debug = new TH2F("debug","debug",100,0,5,100,-15,5);
    
    for(int ien=0; ien<nen; ien++){
        fChain->GetEntry(ien);
        printStatusBar(ien,nen);
        
        // run number is important for correcting the SC paddles
        std::string file = (std::string) fChain->GetCurrentFile()->GetName();
        std::string run  = file.substr(68,6);
        int runno        = atoi(run.c_str()); // creates int, possibly outdated way to do it, but it works
        
        // cm ns^-1
        Float_t speed_of_light = 29.979;
        
        for(int ipart=0; ipart<gpart; ipart++){
            if(isElectron(0)){
                int sector = sc_sect[ipart] -1;
                
                // calculate the time difference between experimental and expected for pions
                Float_t time_exp  = hadronScTimeCorr(1,ipart,runno) - getStartTime();
                Float_t beta      = p[ipart]/sqrt(p[ipart]*p[ipart]+0.135*0.135);
                Float_t time_calc = sc_r[ipart]/(beta*speed_of_light);
                Float_t dt        = time_calc - time_exp;
                
                debug->Fill(p[ipart],dt);
                
            }
        }
    }
    
    std::cout << std::endl;
    
    TCanvas * c1 = new TCanvas("c1","",1100,800);
    
    c1->Print("ppTimeCorrelation.pdf[");
    debug->Draw("colz");
    c1->Print("ppTimeCorrelation.pdf");
    c1->Print("ppTimeCorrelation.pdf]");
    
    return;
}

// =============================================================
//
// =============================================================

void master::pionBeta(){
    
    std::cout << std::endl;
    std::cout << "inside of pionBeta()" << std::endl;
    std::cout << std::endl;
    
    int nen = fChain->GetEntries();
    
    // our canvas
    TCanvas * c1 = new TCanvas("c1","",800,600);
    c1->Print("ppBeta.pdf[");
    
    
    Float_t p_bin_size = (p_max_pion_beta - p_min_pion_beta)/(n_p_bins_pion_beta -1);
    Float_t bins[n_p_bins_pion_beta];
    
    
    // load bins of momentum
    for(int ibin=0; ibin<n_p_bins_pion_beta; ibin++){
        bins[ibin] = ibin*p_bin_size + p_min_pion_beta;
    }
    
    if(debug){
        std::cout << "> binning structure: " << std::endl;
        std::cout << Form("     > number of p bins: %d",n_p_bins_pion_beta) << std::endl;
        std::cout << Form("     > bin size        : %f",p_bin_size) << std::endl;
        std::cout << Form("     > first bin       : %f",bins[0]) << std::endl;
        std::cout << Form("     > last bin        : %f",bins[n_p_bins_pion_beta -1]) << std::endl;
    }
    
    // [sector][bin]
    // note: in situations such as this, always use the scheme which allows you to load a[i][j] -> a[i][j+1]
    // and not a[i][j] -> a[i+1][j] because higher dim arrays are still stored in c++ as 1d arrays and the
    // elements a[i][j] and a[i][j+1] are next to each other but a[i][j] and a[i+1][j] are not.
    TH2F * h_pp_beta[6];
    TH2F * h_pp_beta_all_sectors = new TH2F("h_pp_beta_all_sectors","Beta vs. P for all sectors",500,0,5,500,0.55,1.1);
    
    TH1F * h_pp_beta_slice[6][n_p_bins_pion_beta];
    TF1  * f_pp_beta_slice[6][n_p_bins_pion_beta];
    
    TH2F * h_pm_beta[6];
    TH2F * h_pm_beta_all_sectors = new TH2F("h_pm_beta_all_sectors","Beta vs. P for all sectors",500,0,5,500,0.55,1.1);
    TH1F * h_pm_beta_slice[6][n_p_bins_pion_beta];
    TF1  * f_pm_beta_slice[6][n_p_bins_pion_beta];
    
    // initializing our histograms as well as our Gaussian fits
    // with unique fit boundaries depending on momentum
    for(int ihist=0; ihist<6; ihist++){
        h_pp_beta[ihist] = new TH2F(Form("h_pp_beta_%d",ihist),Form("beta vs. momentum sector %d",ihist),500,0,5,500,0.55,1.1);
        h_pm_beta[ihist] = new TH2F(Form("h_pm_beta_%d",ihist),Form("beta vs. momentum sector %d",ihist),500,0,5,500,0.55,1.1);
        
        for(int ibin=0; ibin<n_p_bins_pion_beta; ibin++){
            h_pp_beta_slice[ihist][ibin] = new TH1F(Form("h_pp_beta_slice_%d_%d",ihist,ibin),Form("sector %d bin %d",ihist,ibin),100,0.4,1.1);
            h_pm_beta_slice[ihist][ibin] = new TH1F(Form("h_pm_beta_slice_%d_%d",ihist,ibin),Form("sector %d bin %d",ihist,ibin),100,0.4,1.1);
        }
    }
    
    std::cout << "> histograms initialized successfully" << std::endl;
    std::cout << "> now filling histograms" << std::endl;
    
    // cm ns^-1
    Float_t speed_of_light = 29.9792458;
    
    // fill histograms
    for(int ien=0; ien<nen; ien++){
        fChain->GetEntry(ien);
        printStatusBar(ien,nen);
        
        // some things needed to get beta
        std::string file  = (std::string) fChain->GetCurrentFile()->GetName();
        std::string piece = file.substr(68,6);
        Int_t runno       = atoi(piece.c_str());
        
        Float_t start_time = getStartTime();
        
        for(int ipart=1; ipart<gpart; ipart++){
            int sector   = ec_sect[ipart] -1;
            Float_t time = hadronScTimeCorr(!getMCstatus(),ipart,runno) - start_time;
            Float_t beta = (sc_r[ipart]/time)/speed_of_light;
            
            if(isPos(ipart))
                if(isPionMissingMass(ipart))
                    if(sc_sect[ipart]>0)
                        if(sector > -1){
                            h_pp_beta[sector]->Fill(p[ipart],beta);
                            h_pp_beta_all_sectors->Fill(p[ipart],beta);
                            
                            if(floor((p[ipart]-p_min_pion_beta)/p_bin_size) >= 0 && floor((p[ipart]-p_min_pion_beta)/p_bin_size) <= (n_p_bins_pion_beta-1)){
                                h_pp_beta_slice[sector][(Int_t) floor((p[ipart]-p_min_pion_beta)/p_bin_size)]->Fill(beta);
                            }
                        }
            
            if(isNeg(ipart))
                if(sc_sect[ipart]>0)
                    if(sector > -1){
                        h_pm_beta[sector]->Fill(p[ipart],beta);
                        h_pm_beta_all_sectors->Fill(p[ipart],beta);
                        
                        if(floor((p[ipart]-p_min_pion_beta)/p_bin_size) >= 0 && floor((p[ipart]-p_min_pion_beta)/p_bin_size) <= (n_p_bins_pion_beta-1)){
                            h_pm_beta_slice[sector][(Int_t) floor((p[ipart]-p_min_pion_beta)/p_bin_size)]->Fill(beta);
                        }
                    }
            
            
        }
    }
    
    
    Float_t cv_pp_beta_min[6];
    Float_t cv_pp_beta_max[6];
    Float_t cv_pm_beta_min[6];
    Float_t cv_pm_beta_max[6];
    
    std::cout << std::endl;
    std::cout << "> fitting histogram slices in p" << std::endl;
    
    for(int isect=0; isect<6; isect++){
        for(int ibin=0; ibin<n_p_bins_pion_beta; ibin++){
            Float_t our_mom = p_min_pion_beta + ibin*p_bin_size;
            Float_t mu      = our_mom/sqrt(0.135*0.135 + our_mom*our_mom);
            
            f_pp_beta_slice[isect][ibin] = new TF1(Form("f_pp_beta_slice_%d_%d",isect,ibin),"gaus",mu-0.03,mu+0.03);
            h_pp_beta_slice[isect][ibin]->Fit(Form("f_pp_beta_slice_%d_%d",isect,ibin),"Rq");
            
            f_pm_beta_slice[isect][ibin] = new TF1(Form("f_pm_beta_slice_%d_%d",isect,ibin),"gaus",mu-0.1,mu+0.1);
            h_pm_beta_slice[isect][ibin]->Fit(Form("f_pm_beta_slice_%d_%d",isect,ibin),"Rq");
        }
    }
    
    
    TGraphErrors * g_mean_pp[6];
    TGraphErrors * g_stddev_pp[6];
    TF1 * f_mean_pp[6];
    TF1 * f_stddev_pp[6];
    
    TGraphErrors * g_mean_pm[6];
    TGraphErrors * g_stddev_pm[6];
    TF1 * f_mean_pm[6];
    TF1 * f_stddev_pm[6];
    
    Float_t mean_pp[n_p_bins_pion_beta];
    Float_t stddev_pp[n_p_bins_pion_beta];
    Float_t mean_err_pp[n_p_bins_pion_beta];
    Float_t stddev_err_pp[n_p_bins_pion_beta];
    
    Float_t mean_pm[n_p_bins_pion_beta];
    Float_t stddev_pm[n_p_bins_pion_beta];
    Float_t mean_err_pm[n_p_bins_pion_beta];
    Float_t stddev_err_pm[n_p_bins_pion_beta];
    
    Float_t empty_bins[n_p_bins_pion_beta];
    
    for(int isect=0; isect<6; isect++){
        for(int ibin=0; ibin<n_p_bins_pion_beta; ibin++){
            empty_bins[ibin]    = 0.00;
            mean_pp[ibin]       = f_pp_beta_slice[isect][ibin]->GetParameter(1);
            stddev_pp[ibin]     = f_pp_beta_slice[isect][ibin]->GetParameter(2);
            mean_pm[ibin]       = f_pm_beta_slice[isect][ibin]->GetParameter(1);
            stddev_pm[ibin]     = f_pm_beta_slice[isect][ibin]->GetParameter(2);
            mean_err_pp[ibin]   = f_pp_beta_slice[isect][ibin]->GetParError(1);
            stddev_err_pp[ibin] = f_pp_beta_slice[isect][ibin]->GetParError(2);
            mean_err_pm[ibin]   = f_pm_beta_slice[isect][ibin]->GetParError(1);
            stddev_err_pm[ibin] = f_pm_beta_slice[isect][ibin]->GetParError(2);
        }
        
        g_mean_pp[isect]   = new TGraphErrors(n_p_bins_pion_beta,bins,mean_pp,empty_bins,mean_err_pp);
        g_stddev_pp[isect] = new TGraphErrors(n_p_bins_pion_beta,bins,stddev_pp,empty_bins,stddev_err_pp);
        f_mean_pp[isect]   = new TF1(Form("f_mean_pp_%d",isect),"x/sqrt(x*x + [0])",0.5,2.5);
        f_stddev_pp[isect] = new TF1(Form("f_stddev_pp_%d",isect),"pol2",0.5,2.5);
        g_mean_pp[isect]->Fit(Form("f_mean_pp_%d",isect),"Frq");
        g_stddev_pp[isect]->Fit(Form("f_stddev_pp_%d",isect),"Frq");
        g_mean_pm[isect]   = new TGraphErrors(n_p_bins_pion_beta,bins,mean_pm,empty_bins,mean_err_pm);
        g_stddev_pm[isect] = new TGraphErrors(n_p_bins_pion_beta,bins,stddev_pm,empty_bins,stddev_err_pm);
        f_mean_pm[isect]   = new TF1(Form("f_mean_pm_%d",isect),"x/sqrt(x*x + [0])",0.2,2.5);
        f_stddev_pm[isect] = new TF1(Form("f_stddev_pm_%d",isect),"pol2",0.2,2.5);
        
        // try to help the fitter
        f_mean_pm[isect]->SetParameter(0,sqrt(0.135));
        g_mean_pm[isect]->Fit(Form("f_mean_pm_%d",isect),"Fr");
        g_stddev_pm[isect]->Fit(Form("f_stddev_pm_%d",isect),"Fr");
    }
    
    // sector, param, bin
    Float_t cv_pp_upper[6][3];
    Float_t cv_pp_lower[6][3];
    Float_t temp_upper_pp[n_p_bins_pion_beta];
    Float_t temp_lower_pp[n_p_bins_pion_beta];
    Float_t temp_upper_pm[n_p_bins_pion_beta];
    Float_t temp_lower_pm[n_p_bins_pion_beta];
    Float_t pp_err[n_p_bins_pion_beta];
    Float_t pm_err[n_p_bins_pion_beta];
    TGraphErrors * upper_pp[6];
    TGraphErrors * lower_pp[6];
    TGraphErrors * upper_pm[6];
    TGraphErrors * lower_pm[6];
    TF1 * f_upper_pp[6];
    TF1 * f_lower_pp[6];
    TF1 * f_upper_pm[6];
    TF1 * f_lower_pm[6];
    
    Float_t n_sigma_pm     = 3.5;
    Float_t n_sigma_pp_top = 4.0;
    Float_t n_sigma_pp_bot = 4.0;
    
    for(int isect=0; isect<6; isect++){
        for(int ibin=0; ibin<n_p_bins_pion_beta; ibin++){
            temp_upper_pm[ibin] =  bins[ibin]/sqrt(f_mean_pm[isect]->GetParameter(0) +bins[ibin]*bins[ibin]) + n_sigma_pm*f_pm_beta_slice[isect][ibin]->GetParameter(2);
            temp_lower_pm[ibin] =  bins[ibin]/sqrt(f_mean_pm[isect]->GetParameter(0) +bins[ibin]*bins[ibin]) - n_sigma_pm*f_pm_beta_slice[isect][ibin]->GetParameter(2);
            temp_upper_pp[ibin] =  bins[ibin]/sqrt(f_mean_pp[isect]->GetParameter(0) +bins[ibin]*bins[ibin]) + n_sigma_pp_top*f_pp_beta_slice[isect][ibin]->GetParameter(2);
            temp_lower_pp[ibin] =  bins[ibin]/sqrt(f_mean_pp[isect]->GetParameter(0) +bins[ibin]*bins[ibin]) - n_sigma_pp_bot*f_pp_beta_slice[isect][ibin]->GetParameter(2);
            pp_err[ibin]        = f_pp_beta_slice[isect][ibin]->GetParError(1);
            pm_err[ibin]        = f_pm_beta_slice[isect][ibin]->GetParError(1);
        }
        upper_pp[isect]   = new TGraphErrors(n_p_bins_pion_beta,bins,temp_upper_pp,empty_bins,pp_err);
        lower_pp[isect]   = new TGraphErrors(n_p_bins_pion_beta,bins,temp_lower_pp,empty_bins,pp_err);
        upper_pm[isect]   = new TGraphErrors(n_p_bins_pion_beta,bins,temp_upper_pm,empty_bins,pm_err);
        lower_pm[isect]   = new TGraphErrors(n_p_bins_pion_beta,bins,temp_lower_pm,empty_bins,pm_err);
        f_upper_pp[isect] = new TF1(Form("f_upper_pp_%d",isect),"x/sqrt(x*x + [0])",0.5,3.0);
        f_lower_pp[isect] = new TF1(Form("f_lower_pp_%d",isect),"x/sqrt(x*x + [0])",0.5,3.0);
        f_upper_pm[isect] = new TF1(Form("f_upper_pm_%d",isect),"x/sqrt(x*x + [0])",0.5,3.0);
        f_lower_pm[isect] = new TF1(Form("f_lower_pm_%d",isect),"x/sqrt(x*x + [0])",0.5,3.0);
        upper_pp[isect]->Fit(Form("f_upper_pp_%d",isect),"Frq");
        lower_pp[isect]->Fit(Form("f_lower_pp_%d",isect),"Frq");
        upper_pm[isect]->Fit(Form("f_upper_pm_%d",isect),"Frq");
        lower_pm[isect]->Fit(Form("f_lower_pm_%d",isect),"Frq");
    }
    
    h_pp_beta_all_sectors->Draw("colz");
    c1->Print("ppBetaAllSectors.png");
    c1->Clear();
    
    h_pm_beta_all_sectors->Draw("colz");
    c1->Print("pmBetaAllSectors.png");
    c1->Clear();
    
    // 2d histograms for beta vs p
    // SHOWN WITH FITS
    c1->Divide(3,2);
    c1->cd(1);
    gPad->SetLogz();
    h_pp_beta[0]->Draw("colz");
    upper_pp[0]->Draw("same");
    lower_pp[0]->Draw("same");
    
    c1->cd(2);
    gPad->SetLogz();
    h_pp_beta[1]->Draw("colz");
    upper_pp[1]->Draw("same");
    lower_pp[1]->Draw("same");
    
    c1->cd(3);
    gPad->SetLogz();
    h_pp_beta[2]->Draw("colz");
    upper_pp[2]->Draw("same");
    lower_pp[2]->Draw("same");
    
    c1->cd(4);
    gPad->SetLogz();
    h_pp_beta[3]->Draw("colz");
    upper_pp[3]->Draw("same");
    lower_pp[3]->Draw("same");
    
    c1->cd(5);
    gPad->SetLogz();
    h_pp_beta[4]->Draw("colz");
    upper_pp[4]->Draw("same");
    lower_pp[4]->Draw("same");
    
    c1->cd(6);
    gPad->SetLogz();
    h_pp_beta[5]->Draw("colz");
    upper_pp[5]->Draw("same");
    lower_pp[5]->Draw("same");
    
    c1->Print("ppBeta.pdf");
    c1->Print("ppBeta.png");
    c1->Clear();
    
    c1->Divide(3,2);
    for(int ibin=0; ibin<n_p_bins_pion_beta; ibin++){
        c1->cd(1);
        h_pp_beta_slice[0][ibin]->Draw("same");
        c1->cd(2);
        h_pp_beta_slice[1][ibin]->Draw("same");
        c1->cd(3);
        h_pp_beta_slice[2][ibin]->Draw("same");
        c1->cd(4);
        h_pp_beta_slice[3][ibin]->Draw("same");
        c1->cd(5);
        h_pp_beta_slice[4][ibin]->Draw("same");
        c1->cd(6);
        h_pp_beta_slice[5][ibin]->Draw("same");
        
        c1->Print("ppBeta.pdf");
        c1->Clear();
        c1->Divide(3,2);
    }
    
    for(int isect=0; isect<6; isect++){
        c1->cd(isect+1);
        g_mean_pp[isect]->Draw("same");
    }
    
    c1->Print("ppBeta.pdf");
    c1->Print("ppBeta.pdf]");
    
    // now for negative pions
    c1->Clear();
    c1->Print("pmBeta.pdf[");
    c1->Divide(3,2);
    
    c1->cd(1);
    h_pm_beta[0]->Draw("colz");
    upper_pm[0]->Draw("same");
    lower_pm[0]->Draw("same");
    
    c1->cd(2);
    h_pm_beta[1]->Draw("colz");
    upper_pm[1]->Draw("same");
    lower_pm[1]->Draw("same");
    
    c1->cd(3);
    h_pm_beta[2]->Draw("colz");
    upper_pm[2]->Draw("same");
    lower_pm[2]->Draw("same");
    
    c1->cd(4);
    h_pm_beta[3]->Draw("colz");
    upper_pm[3]->Draw("same");
    lower_pm[3]->Draw("same");
    
    c1->cd(5);
    h_pm_beta[4]->Draw("colz");
    upper_pm[4]->Draw("same");
    lower_pm[4]->Draw("same");
    
    c1->cd(6);
    h_pm_beta[5]->Draw("colz");
    upper_pm[5]->Draw("same");
    lower_pm[5]->Draw("same");
    
    c1->Print("pmBeta.pdf");
    
    c1->Print("pmBeta.png");
    c1->Clear();
    
    c1->Divide(3,2);
    for(int ibin=0; ibin<n_p_bins_pion_beta; ibin++){
        c1->cd(1);
        h_pm_beta_slice[0][ibin]->Draw("same");
        c1->cd(2);
        h_pm_beta_slice[1][ibin]->Draw("same");
        c1->cd(3);
        h_pm_beta_slice[2][ibin]->Draw("same");
        c1->cd(4);
        h_pm_beta_slice[3][ibin]->Draw("same");
        c1->cd(5);
        h_pm_beta_slice[4][ibin]->Draw("same");
        c1->cd(6);
        h_pm_beta_slice[5][ibin]->Draw("same");
        
        c1->Print("pmBeta.pdf");
        c1->Clear();
        c1->Divide(3,2);
    }
    
    c1->Clear();
    c1->Divide(3,2);
    
    c1->Print("pmBeta.pdf]");
    
    // work is required here
    // re-write the method below to output parameters of the fits for pion beta
    // then fix the PID algorithms which depend on said parameters
    
    
    // write it out to the header file
    std::ofstream ofile;
    std::string filename;
    
    filename = "cv_pp_beta.h";
    
    if ( getMCstatus() )
        filename = "cv_pp_beta_MC.h";
    
    ofile.open(filename.c_str(), std::ios::out | std::ios::trunc);
    
    ofile << "Float_t cv_pp_beta_p_upper[6] = {";
    
    for(int isect=0; isect<6; isect++){
        ofile << f_upper_pp[isect]->GetParameter(0);
        if(isect<5)
            ofile << ",";
    }
    
    ofile << "}; \n";
    
    ofile << "Float_t cv_pp_beta_p_lower[6] = {";
    
    for(int isect=0; isect<6; isect++){
        ofile << f_lower_pp[isect]->GetParameter(0);
        if(isect<5)
            ofile << ",";
    }
    
    ofile << "}; \n";
    
    ofile.close();
    
    filename = "cv_pm_beta.h";
    
    if ( getMCstatus() )
        filename = "cv_pm_beta_MC.h";
    
    ofile.open(filename.c_str(), std::ios::out | std::ios::trunc);
    
    ofile << "Float_t cv_pm_beta_p_upper[6] = {";
    
    for(int isect=0; isect<6; isect++){
        ofile << f_upper_pm[isect]->GetParameter(0);
        if(isect<5)
            ofile << ",";
    }
    
    ofile << "}; \n";
    
    ofile << "Float_t cv_pm_beta_p_lower[6] = {";
    
    for(int isect=0; isect<6; isect++){
        ofile << f_lower_pm[isect]->GetParameter(0);
        if(isect<5)
            ofile << ",";
    }
    
    ofile << "}; \n";
    
    ofile.close();
    
    
    return;
}

// =============================================================
//
// =============================================================

// this routine assumes that particle 0 is electron
bool master::isPionMissingMass(int ipart){
    
    TLorentzVector beam(0,0,5.498,5.498);
    TLorentzVector target(0,0,0,0.938);
    TLorentzVector electron(p[0]*cx[0],p[0]*cy[0],p[0]*cz[0],p[0]);
    TLorentzVector pion(p[ipart]*cx[ipart],p[ipart]*cy[ipart],p[ipart]*cz[ipart],sqrt(p[ipart]*p[ipart] + 0.135*0.135));
    TLorentzVector missing = beam + target - electron - pion;
    
    if(missing.M2() >= 0.00)
        if(sqrt(missing.M2() >= 1.50))
            return true;
    
    return false;
}

bool master::isPpBetaP(int ipart){
    
    // some things needed to get beta
    std::string file  = (std::string) fChain->GetCurrentFile()->GetName();
    std::string piece = file.substr(68,6);
    Int_t runno       = atoi(piece.c_str());
    
    Float_t speed_of_light = 29.9792458;
    Float_t start_time     = getStartTime();
    Float_t time           = hadronScTimeCorr(1,ipart,runno) - start_time;
    Float_t beta           = (sc_r[ipart]/time)/speed_of_light;
    Int_t sector           = sc_sect[ipart] -1;
    
    
    Float_t cut_max = p[ipart]/sqrt(cv_pp_beta_p_upper[sector]*cv_pp_beta_p_upper[sector] + p[ipart]*p[ipart]);
    Float_t cut_min = p[ipart]/sqrt(cv_pp_beta_p_lower[sector]*cv_pp_beta_p_lower[sector] + p[ipart]*p[ipart]);
    
    if(sector > -1)
        if(beta <= cut_max && beta >= cut_min)
            return true;
    
    return false;
}

// =============================================================
//
// =============================================================

bool master::isPmBetaP(int ipart){
    
    // some things needed to get beta
    std::string file  = (std::string) fChain->GetCurrentFile()->GetName();
    std::string piece = file.substr(68,6);
    Int_t runno       = atoi(piece.c_str());
    
    Float_t speed_of_light = 29.9792458;
    Float_t start_time = getStartTime();
    Float_t time = hadronScTimeCorr(1,ipart,runno) - start_time;
    Float_t beta = (sc_r[ipart]/time)/speed_of_light;
    
    Int_t sector = sc_sect[ipart] -1;
    
    Float_t cut_max = p[ipart]/sqrt(cv_pm_beta_p_upper[sector]*cv_pm_beta_p_upper[sector] + p[ipart]*p[ipart]);
    Float_t cut_min = p[ipart]/sqrt(cv_pm_beta_p_lower[sector]*cv_pm_beta_p_lower[sector] + p[ipart]*p[ipart]);
    
    if(sector > -1)
        if(beta <= cut_max && beta >= cut_min)
            return true;
    
    return false;
}

// =============================================================
//
// =============================================================

bool master::isPP(int ipart){
    if(isPos(ipart))
        if(ec_sect[ipart] >0 && sc_sect[ipart] >0)
            if(isPionMissingMass(ipart))
                if(isPpBetaP(ipart))
                    return true;
    
    return false;
}

// =============================================================
//
// =============================================================

bool master::isPM(int ipart){
    if(isNeg(ipart))
        if(ec_sect[ipart] >0 && sc_sect[ipart] >0)
            if(isPionMissingMass(ipart))
                if(isPmBetaP(ipart))
                    return true;
    
    return false;
}

void master::generateElectronCutValues()
{
  
  // we always need to consider if we have data or
  // Monte Carlo, this stores true for MC, false for data
  bool isMC = getMCstatus();
  
  int nentries = fChain->GetEntries();
  
  // binning options, constants 
  int n_p_bins_ec_sampling = 50;
  float p_max_ec_sampling = 3.5;
  float p_min_ec_sampling = 0.5;
  int nsect               = 6;
  int nsegm               = 18;

  // electron histograms
  TH2F * h_ec_sampling[nsect];
  TH2F * h_cc_theta[nsect];
  TH1F * h_ec_sampling_slice[nsect][n_p_bins_ec_sampling];
  TH1F * h_cc_theta_slice[nsect][nsegm];

  // electron fits 
  TF1 * f_ec_sampling_slice[nsect][n_p_bins_ec_sampling];
  TF1 * f_cc_theta_slice[nsect][nsegm];
  TF1 * f_ec_sampling_upper[nsect];
  TF1 * f_ec_sampling_lower[nsect];

  // electron graphs 
  TGraphErrors * g_ec_upper[nsect];
  TGraphErrors * g_ec_lower[nsect];

  // initialize the histograms
  for (int isect = 0; isect < nsect; isect++)
    {
      std::string ecSamplingTitle = Form("h_ec_sampling_s%d",isect);
      std::string ccThetaTitle    = Form("h_cc_theta_s%d",isect);

      h_ec_sampling[isect] = new TH2F(ecSamplingTitle.c_str(),ecSamplingTitle.c_str(),500,0,6,500,0,0.5);
      h_cc_theta[isect]    = new TH2F(ccThetaTitle.c_str(),ccThetaTitle.c_str(),18,0,17,500,0,60);

    }

  std::cout << "\n > looping over entries \n" << std::endl;

  for (int ientry = 0; ientry < nentries; ientry++)
    {
      fChain->GetEntry(ientry);

      printStatusBar(ientry,nentries);

      // since electron is put into 0 if it exists
      // no need to loop over particles, just check
      // ipart = 0
      
      int cc_segment = (cc_segm[0]%1000)/10 -1;
      int cc_sector  = (cc_sect[0] -1);

      if (q[0] < 0)
	if ( (cc_segment > -1) && (cc_sector > -1) ) // both have hits
	  {
	    h_ec_sampling[cc_sector]->Fill(p[0],etot[0]/p[0]);
	    h_cc_theta[cc_sector]->Fill(cc_segment,getThetaCC(0));
	  }
    }

  std::cout << std::endl;

} // end master::generateElectronCutValues()
