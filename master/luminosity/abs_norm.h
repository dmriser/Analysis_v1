#include <iostream>

bool eDepCut(Float_t edep, Float_t p){
  if(edep/p >= 0.25)
  return 1;

  return 0;
}

bool qCut(Int_t q){
  if(q<0)
    return 1;

  return 0;
}
