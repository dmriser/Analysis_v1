#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#include <vector>
#include <map>

// root files 
#include "TLorentzVector.h"
#include "TVector3.h"

// my files
#include "h22Event.h"

// Nathan H. Files for PID
#include "programFiles/functions.C"
#include "programFiles/eID.C"
#include "programFiles/hadronID.C"

/**
 * ParticleFilter is the class which manages particle identification 
 * for the E1F dataset.  All of the core particle ID algorithms are 
 * from Nathan Harrison, and the ParticleFilter class acts as a wrapper.
 * The ParticleFilter class returns the identified particles as TLorentzVectors,
 * or by event particle index number.  The ParticleFilter class also contains 
 * counting algorithms.  
 */
class ParticleFilter
{
 public:
  int runno; // needed for PID
  int GSIM;
  h22Event event;

  // constructor (event, gsim (0,1), run number )
  ParticleFilter();
  virtual ~ParticleFilter();

  int getIndexByPID(int);
  vector<int> getAllHadronIndices();
  map<string,bool> getEIDStatus(int);
  void loadEvent(h22Event, int, int);  
  //  int newGetIndexByPID(int);

};
#endif

ParticleFilter::ParticleFilter()
{
  // Load parameters from the database to do PID.
}

ParticleFilter::~ParticleFilter()
{
  // Close all. In this case nothing to do. 
}

void ParticleFilter::loadEvent(h22Event input, int mc, int run_num)
{
  event = input;
  GSIM  = mc;
  runno = run_num;
}

vector<int> ParticleFilter::getAllHadronIndices()
{
  int e_index = eID(event.gpart, event.q, event.p, event.cc_sect, event.sc_sect, event.ec_sect, event.dc_sect, event.cx, event.cy, 
		  event.cz, event.tl1_x, event.tl1_y, event.tl3_x, event.tl3_y, event.tl3_z, event.tl3_cx, event.tl3_cy, event.tl3_cz, 
		  0, event.vz, event.vy, event.vx, 0, !GSIM, event.etot, 0, event.ec_ei, event.ech_x, event.ech_y, event.ech_z, 0,
		  event.cc_segm, 0, 0, 0, 0, event.sc_pd, 0);   

  // construct the unknown hardonic state as a 4-Vector, v4_H
  TLorentzVector v4_k(0, 0, 5.498, 5.498);
  TLorentzVector v4_kprime(event.p[e_index]*event.cx[e_index], event.p[e_index]*event.cy[e_index], event.p[e_index]*event.cz[e_index], event.p[e_index]);
  TLorentzVector v4_target(0,0,0,0.938);
  TLorentzVector v4_H = (v4_k - v4_kprime) + v4_target;

  return hadronID(event.gpart, e_index, event.q, event.p, event.sc_sect, event.dc_sect, event.sc_t, event.sc_r, event.sc_pd, 0, 0, 0, !GSIM, event.ec_ei, event.ec_sect, event.cc_sect, event.nphe, event.ec_eo, event.cx, event.cy, event.cz, event.b, event.tl1_x, event.tl1_y, v4_H, runno, 0, 0, 0);
}

map <string, bool> ParticleFilter::getEIDStatus(int ipart)
{

    return eID_map(ipart, event.gpart, event.q, event.p, event.cc_sect, event.sc_sect, event.ec_sect, event.dc_sect, event.cx, event.cy, 
		  event.cz, event.tl1_x, event.tl1_y, event.tl3_x, event.tl3_y, event.tl3_z, event.tl3_cx, event.tl3_cy, event.tl3_cz, 
		  0, event.vz, event.vy, event.vx, 0, !GSIM, event.etot, 0, event.ec_ei, event.ech_x, event.ech_y, event.ech_z, 0,
		  event.cc_segm, 0, 0, 0, 0, event.sc_pd, 0);   
}

int ParticleFilter::getIndexByPID(int pid)
{

  int index = -123; // NH convention

  // strictnesses can be included later by altering the class def and adding a SetEIDParameters(), SetHardonIDParameters() method.
  // all set to 0 for now.
  int e_index = eID(event.gpart, event.q, event.p, event.cc_sect, event.sc_sect, event.ec_sect, event.dc_sect, event.cx, event.cy, 
		  event.cz, event.tl1_x, event.tl1_y, event.tl3_x, event.tl3_y, event.tl3_z, event.tl3_cx, event.tl3_cy, event.tl3_cz, 
		  0, event.vz, event.vy, event.vx, 0, !GSIM, event.etot, 0, event.ec_ei, event.ech_x, event.ech_y, event.ech_z, 0,
		  event.cc_segm, 0, 0, 0, 0, event.sc_pd, 0);   

  if (e_index < 0) return index; // hadron ID depends on finding an electron, so don't continue unless we have electrons
 
 if (pid == 11)
    {
      index = e_index;
    }


  if ((pid == 2212) || (pid == 211) || (pid == -211))
    {
      // construct the unknown hardonic state as a 4-Vector, v4_H
      TLorentzVector v4_k(0, 0, 5.498, 5.498);
      TLorentzVector v4_kprime(event.p[e_index]*event.cx[e_index], event.p[e_index]*event.cy[e_index], event.p[e_index]*event.cz[e_index], event.p[e_index]);
      TLorentzVector v4_target(0,0,0,0.938);
      TLorentzVector v4_H = (v4_k - v4_kprime) + v4_target;

      vector<int> indicies = hadronID(event.gpart, e_index, event.q, event.p, event.sc_sect, event.dc_sect, event.sc_t, event.sc_r, event.sc_pd, 0, 0, 0, !GSIM, event.ec_ei, event.ec_sect, event.cc_sect, event.nphe, event.ec_eo, event.cx, event.cy, event.cz, event.b, event.tl1_x, event.tl1_y, v4_H, runno, 0, 0, 0);
	
      if (pid ==  211) index = indicies[0];
      if (pid == -211) index = indicies[1];
      if (pid == 2212) index = indicies[2];
    }

  return index;
}
