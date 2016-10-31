#include "TLorentzVector.h"

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

 public:
  // stuff
  ElasticEvent(TLorentzVector elec);
 ~ElasticEvent();
  
  Float_t getW(){return W;}
  Float_t getqq(){return qq;}
  Float_t getQQ(){return QQ;}
  Float_t getx(){return x;}

};

ElasticEvent::ElasticEvent(TLorentzVector elec)
{
  TLorentzVector beam(0,0,5.498,5.498);
  TLorentzVector target(0,0,0,0.938);

  electron      = elec;
  virtualPhoton = beam - electron;
  proton        = target + virtualPhoton;

  W  = proton.Mag();
  qq = virtualPhoton.Mag();
  QQ = -qq;
  x  = 1; 
}

ElasticEvent::~ElasticEvent()
{
 
}
