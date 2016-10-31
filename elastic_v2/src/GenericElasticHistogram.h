#include <iostream>
#include "TH2.h"
#include "../../analysisClasses/Bins.h"

class GenericElasticHistogram{
 public:
  GenericElasticHistogram(std::string filename, std::string histogramTitle, BinStructure thetaBins, BinStructure phiBins, bool recalc = true);
  ~GenericElasticHistogram();

  TH2D * thetaVsPhi[7]; 
  std::string filename, histogramTitle; 

  void Init();
  void Load();
  void Save();
};
