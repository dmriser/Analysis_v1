#include "TLorentzVector.h"

class elasticEvent 
{
 private:
  TLorentzVector electron;
  TLorentzVector proton;
  TLorentzVector virtualPhoton;
  
  Float_t W;
  Float_t qq;
  Float_t QQ;
  Float_t x;

 public:
  elasticEvent();
  virtual void ~elasticEvent();

  elasticEvent(TLorentzVector, TLorentzVector);
  void calculateKinematics();

  void setElectron(TLorentzVector electronIn){electron = electronIn;};
  void setProton(TLorentzVector protonIn){proton = protonIn;};

  Float_t getW(){return W;};
  Float_t getqq(){return qq;};
  Float_t getQQ(){return QQ;};
  Float_t getx(){return x;};

  // in case user needs it
  TLorentzVector getVirtualPhoton(){return virtualPhoton;};

}

elasticEvent::elasticEvent()
{
  W = 0;
  qq = 0;
  QQ = -qq;
  x = 0;
}

elasticEvent::~elasticEvent(){}

elasticEvent::elasticEvent(TLorentzVector electronIn, TLorentzVector protonIn)
{
  electron = electronIn;
  proton   = protonIn; 
}

// to do 
void elasticEvent::calculateKinematics()
{
  W = 0;
}
