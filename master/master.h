
#define master_h

#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TH2.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include "config.h"
#include "statistics/cv_theta_cc.h"
#include "statistics/cv_ec_sampling.h"

// constants of importance for fns
const Float_t r2d = 180/3.14159;
const Float_t d2r = 1/r2d;

// Header file for the classes stored in the TTree if any.
// Fixed size dimensions of array or collections stored in the TTree if any.
Float_t atan3(Float_t y, float_t x);
void printStatusBar(int x, int max);

class master {
private:
    bool isMC;
    
    public :
    TChain *fChain;   //!pointer to the analyzed TTree or TChain
    
    void getGoodRuns();
    void eidValidation();
    void showPhiCC();
    void ppMissingMass();
    void ppTimeCorrelation();
    void pionBeta();
    bool goodORbadSCpaddle(int);
    Float_t electronScTimeCorr(int,int,int);
    Float_t hadronScTimeCorr(int,int,int);
    Float_t getStartTime();
    
    // electron ID stuff
    bool isElectron(int);
    bool isNeg(int);
    bool eDepEC(int, int);
    bool isVertex(int,int);
    bool thetaCCcut(int, int);
    bool phiCCcut(int);
    bool ecGeometricCut(int, int);
    bool ecInnerCut(int, int);
    bool nPhotoElectrons(int);
    Float_t getThetaCC(int);
    Float_t getRelPhi(int);
    int ccPhiMatching(int);
    Float_t getCorrZ(int);
    void generateElectronCutValues();
    
    // pion ID stuff
    bool isPP(int);
    bool isPM(int);
    bool isPos(int);
    bool isPionMissingMass(int);
    bool isPpBetaP(int);
    bool isPmBetaP(int);
    
    // proton ID
    bool isProton(int);
    void runProtonBeta();
    bool isProtonBetaP(int);
    
    // Declaration of leaf types
    UInt_t          evntid;
    UChar_t         ihel;
    Float_t         q_l;
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
    
    // List of branches
    TBranch        *b_evntid;   //!
    TBranch        *b_ihel;   //!
    TBranch        *b_q_l;   //!
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
    
    master();
    virtual ~master();
    virtual void     Init(bool);
    virtual void     Process();
    int AddFile(char *);
    int GetEntries(){return fChain->GetEntries();};
    bool getMCstatus(){return isMC;};
    void setMCstatus(bool stat){isMC = stat;};
};


#ifdef master_cxx
master::master()
{
    fChain = new TChain("h22");
}

master::~master()
{
    if (!fChain) return;
    delete fChain;
}

int master::AddFile(char *_fname){
    return fChain->AddFile(_fname);
}

void master::Init(bool isMC)
{
    fChain->SetBranchAddress("evntid", &evntid, &b_evntid);
    fChain->SetBranchAddress("ihel", &ihel, &b_ihel);
    fChain->SetBranchAddress("q_l", &q_l, &b_q_l);
    fChain->SetBranchAddress("gpart", &gpart, &b_gpart);
    fChain->SetBranchAddress("q", q, &b_q);
    fChain->SetBranchAddress("p", p, &b_p);
    fChain->SetBranchAddress("b", b, &b_b);
    fChain->SetBranchAddress("cx", cx, &b_cx);
    fChain->SetBranchAddress("cy", cy, &b_cy);
    fChain->SetBranchAddress("cz", cz, &b_cz);
    fChain->SetBranchAddress("vz", vz, &b_vz);
    fChain->SetBranchAddress("dc_sect", dc_sect, &b_dc_sect);
    fChain->SetBranchAddress("tl1_cx", tl1_cx, &b_tl1_cx);
    fChain->SetBranchAddress("tl1_cy", tl1_cy, &b_tl1_cy);
    fChain->SetBranchAddress("ec_sect", ec_sect, &b_ec_sect);
    fChain->SetBranchAddress("ec_r", ec_r, &b_ec_r);
    fChain->SetBranchAddress("ec_t", ec_t, &b_ec_t);
    fChain->SetBranchAddress("ec_ei", ec_ei, &b_ec_ei);
    fChain->SetBranchAddress("ec_eo", ec_eo, &b_ec_eo);
    fChain->SetBranchAddress("etot", etot, &b_etot);
    fChain->SetBranchAddress("cc_sect", cc_sect, &b_cc_sect);
    fChain->SetBranchAddress("cc_r", cc_r, &b_cc_r);
    fChain->SetBranchAddress("cc_t", cc_t, &b_cc_t);
    fChain->SetBranchAddress("nphe", nphe, &b_nphe);
    fChain->SetBranchAddress("cc_c2", cc_c2, &b_cc_c2);
    fChain->SetBranchAddress("sc_sect", sc_sect, &b_sc_sect);
    fChain->SetBranchAddress("sc_r", sc_r, &b_sc_r);
    fChain->SetBranchAddress("sc_t", sc_t, &b_sc_t);
    fChain->SetBranchAddress("edep", edep, &b_edep);
    fChain->SetBranchAddress("sc_pd", sc_pd, &b_sc_pd);
    fChain->SetBranchAddress("cc_segm", cc_segm, &b_cc_segm);
    fChain->SetBranchAddress("ech_x", ech_x, &b_ech_x);
    fChain->SetBranchAddress("ech_y", ech_y, &b_ech_y);
    fChain->SetBranchAddress("ech_z", ech_z, &b_ech_z);
    fChain->SetBranchAddress("tl1_x", tl1_x, &b_tl1_x);
    fChain->SetBranchAddress("tl1_y", tl1_y, &b_tl1_y);
    fChain->SetBranchAddress("tl1_z", tl1_z, &b_tl1_z);
    fChain->SetBranchAddress("tl3_x", tl3_x, &b_tl3_x);
    fChain->SetBranchAddress("tl3_y", tl3_y, &b_tl3_y);
    fChain->SetBranchAddress("tl3_z", tl3_z, &b_tl3_z);
    fChain->SetBranchAddress("tl3_cx", tl3_cx, &b_tl3_cx);
    fChain->SetBranchAddress("tl3_cy", tl3_cy, &b_tl3_cy);
    fChain->SetBranchAddress("tl3_cz", tl3_cz, &b_tl3_cz);
    fChain->SetBranchAddress("id", id, &b_id);
    fChain->SetBranchAddress("vx", vx, &b_vx);
    fChain->SetBranchAddress("vy", vy, &b_vy);
    
    if (isMC)
    {
        fChain->SetBranchAddress("mcnentr", &mcnentr, &b_mcnentr);
        fChain->SetBranchAddress("mcnpart", &mcnpart, &b_mcnpart);
        fChain->SetBranchAddress("mcst", mcst, &b_mcst);
        fChain->SetBranchAddress("mcid", mcid, &b_mcid);
        fChain->SetBranchAddress("mcpid", mcpid, &b_mcpid);
        fChain->SetBranchAddress("mctheta", mctheta, &b_mctheta);
        fChain->SetBranchAddress("mcphi", mcphi, &b_mcphi);
        fChain->SetBranchAddress("mcp", mcp, &b_mcp);
        fChain->SetBranchAddress("mcm", mcm, &b_mcm);
        fChain->SetBranchAddress("mcvx", mcvx, &b_mcvx);
        fChain->SetBranchAddress("mcvy", mcvy, &b_mcvy);
        fChain->SetBranchAddress("mcvz", mcvz, &b_mcvz);
        fChain->SetBranchAddress("mctof", mctof, &b_mctof);
    }
    
    
}
#endif // #ifdef master_cxx

void printStatusBar(int x, int max){
    
    int width = 50;
    
    if(width > max)
        width=max;
    
    double ratio = (double) x/max;
    int percentage = (int) ratio*100;
    int step = max/width;
    int current = floor(x/step);
    
    std::cout << "\r {";
    
    for(int i=0; i<width; i++){
        if(i <= current)
            std::cout << "=";
        if(i > current)
            std::cout << " ";
    }
    
    std::cout << "} " << floor(ratio*100) << "%";
    
    return;
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
