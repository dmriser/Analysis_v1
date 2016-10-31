#ifndef ELASTIC_EVENT_H
#define ELASTIC_EVENT_H

#include <iostream>
#include "TLorentzVector.h"

using namespace std;

/**
 * ElasticEvent class accepts a TLorentzVector which is the final state 
 * electron from PID.  It then calculates the 4-vectors for the virtual 
 * photon, and hadronic state (just a proton here).  It gives back 
 * details on the kinematics of the event.
 */
class ElasticEvent
{
 private:
  TLorentzVector electron;
  TLorentzVector proton;
  TLorentzVector virtualPhoton;

  Float_t W;
  Float_t QQ;
  Float_t qq;
  Float_t x;
  Float_t mp;
  Float_t theta_e_p;
  Float_t MM2;
  Float_t ME;

 public:
  // stuff
  ElasticEvent(TLorentzVector elec);
  ElasticEvent(TLorentzVector elec, TLorentzVector prot);
 ~ElasticEvent();
  
 Float_t getW(){return W;} /**< returns the value of W for the event  */
 Float_t getqq(){return qq;}/**< returns the value of q^2 for the event  */
 Float_t getQQ(){return QQ;}/**< returns the value of virtuality for the event  */
 Float_t getx(){return x;}/**< returns the value of Bjorken x for the event  */
 Float_t getEPAngle(){return theta_e_p;}
 Float_t getMM2(){return MM2;}  
 Float_t getME(){return ME;}

  TLorentzVector getProton(); //*< returns a TLorentzVector that contains the hadronic state (just proton) */

  bool passes();
  void printKinematics();

};

#endif

ElasticEvent::ElasticEvent(TLorentzVector elec)
{
  TLorentzVector beam(0,0,5.498,5.498);
  TLorentzVector target(0,0,0,0.938);

  mp = 0.938;

  electron      = elec;
  virtualPhoton = beam - electron;
  proton        = target + virtualPhoton;

  theta_e_p = electron.Angle(proton.Vect())*180/3.14159;

  // kinematics 
  W  = proton.Mag();
  qq = virtualPhoton.Mag();
  QQ = -qq;
  x  = QQ/(2*target.E()*virtualPhoton.E()); 
  MM2 = (electron-(beam + target)).Mag();
  ME = (electron - (beam+target)).E();
}

ElasticEvent::ElasticEvent(TLorentzVector elec, TLorentzVector proton)
{
  TLorentzVector beam(0,0,5.498,5.498);
  TLorentzVector target(0,0,0,0.938);

  mp = 0.938;

  electron      = elec;
  virtualPhoton = beam - electron;

  theta_e_p = electron.Angle(proton.Vect())*180/3.14159;

  // kinematics 
  W  = proton.Mag();
  qq = virtualPhoton.Mag();
  QQ = -qq;
  x  = QQ/(2*target.E()*virtualPhoton.E()); 
  MM2 = ((electron+proton)-(beam + target)).Mag();
  ME = ((electron+proton)-(beam+target)).E();
}

ElasticEvent::~ElasticEvent()
{
 
}

bool ElasticEvent::passes()
{
  // Elastic event pass criteria
  if(W < 1.1) 
    return true;

  return false;
}

void ElasticEvent::printKinematics()
{
  Float_t toDegrees = 180/3.14159;

  cout.width(15);
  cout << W;
  cout.width(15);
  cout << qq;
  cout.width(15);
  cout << QQ;
  cout.width(15);
  cout << x << endl;


}

TLorentzVector ElasticEvent::getProton()
{return proton;}


