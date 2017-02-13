{

  // --------------------------------------------------------
  //      User Parameters 
  // --------------------------------------------------------

  TFile *inputFile = TFile::Open("../WriteAccumulation.root");
  string imagePath = "/volatile/clas12/dmriser/plots/faradayCup/";   

  // --------------------------------------------------------
  //      Business  
  // --------------------------------------------------------
  vector<TH1D*> chargeHistos; 
  vector<TH1D*>  entryHistos; 
  vector<TH1D*> chargeHistosGood; 
  vector<TH1D*>  entryHistosGood; 

  // Load from file
  TList *tableOfContents = inputFile->GetListOfKeys();
  TIter next(tableOfContents);

  while(TObject *objectFromFile = next()){
    TString currentObjectName  = objectFromFile->GetName();

    if (currentObjectName.Contains("chargeDiff_"))          { chargeHistos.push_back((TH1D*) inputFile->Get(currentObjectName)); }
    else if (currentObjectName.Contains("entryDiff_"))      { entryHistos.push_back((TH1D*) inputFile->Get(currentObjectName)); }
    else if (currentObjectName.Contains("chargeDiffGood_")) { chargeHistosGood.push_back((TH1D*) inputFile->Get(currentObjectName)); }
    else if (currentObjectName.Contains("entryDiffGood_"))  { entryHistosGood.push_back((TH1D*) inputFile->Get(currentObjectName)); }

    cout << "[LoadingHistograms] Found " << currentObjectName << endl;
  }  

  cout << chargeHistos.size() << endl; 
  cout << entryHistos.size() << endl; 
  cout << chargeHistosGood.size() << endl; 
  cout << entryHistosGood.size() << endl; 

  TCanvas *c1 = new TCanvas("c1","",800,500); 

  TLatex stats, tit, xtit, ytit;
  stats.SetNDC();
  stats.SetTextSize(0.03);
  stats.SetTextFont(42);

  tit.SetNDC();
  tit.SetTextSize(0.05);
  tit.SetTextFont(42);

  xtit.SetNDC();
  xtit.SetTextSize(0.05);
  xtit.SetTextFont(42);

  ytit.SetNDC();
  ytit.SetTextSize(0.05);
  ytit.SetTextFont(42);
  ytit.SetTextAngle(90.0);

  for (int i=0; i<chargeHistos.size(); ++i) {

    cout << " Charge = " << chargeHistos[i]->Integral() << endl;

    //    c1->Divide(2,2); 

    //    c1->cd();

    //    c1->cd(1); 
    //    gPad->SetMargin(0.2, 0.03, 0.01, 0.2); 
    gPad->SetGridx();
    gPad->SetGridy();
    chargeHistos[i]->SetFillColorAlpha(99, 0.4); 
    chargeHistosGood[i]->SetFillColorAlpha(78, 0.4); 
    //    chargeHistos[i]->Draw(); 
    chargeHistosGood[i]->Draw(); 
    stats.DrawLatex(0.75, 0.9, Form("Q = %4.f #muC, N = %4.E",chargeHistosGood[i]->Integral(),entryHistosGood[i]->Integral()));
    tit.DrawLatex(0.45, 0.9, Form("Run %s",chargeHistosGood[i]->GetTitle())); 
    xtit.DrawLatex(0.45, 0.05, "Scalar Entry"); 
    ytit.DrawLatex(0.05, 0.45, "#DeltaQ [#muC]");
    /*
    c1->cd(2); 
    gPad->SetMargin(0.03, 0.2, 0.01, 0.2); 
    chargeHistosGood[i]->Draw(); 
    stats.DrawLatex(0.5, 0.75, Form("Q = %4.f #muC",chargeHistosGood[i]->Integral()));

    c1->cd(3); 
    gPad->SetMargin(0.2, 0.03, 0.2, 0.01); 
    entryHistos[i]->Draw(); 
    stats.DrawLatex(0.5, 0.95, Form("N = %4.E",entryHistos[i]->Integral()));

    c1->cd(4); 
    gPad->SetMargin(0.03, 0.2, 0.2, 0.01); 
    entryHistosGood[i]->Draw(); 
    stats.DrawLatex(0.5, 0.95, Form("N = %4.E",entryHistosGood[i]->Integral()));
    */

    c1->Print(Form("%scharge_%s.png",imagePath.c_str(),chargeHistos[i]->GetTitle())); 
    c1->Clear(); 
    gPad->SetGridx();
    gPad->SetGridy();
    entryHistosGood[i]->SetFillColorAlpha(78, 0.4); 
    //    chargeHistos[i]->Draw(); 
    entryHistosGood[i]->Draw(); 
    stats.DrawLatex(0.75, 0.9, Form("Q = %4.f #muC, N = %4.E",chargeHistosGood[i]->Integral(),entryHistosGood[i]->Integral()));
    tit.DrawLatex(0.45, 0.9, Form("Run %s",chargeHistosGood[i]->GetTitle())); 
    xtit.DrawLatex(0.45, 0.05, "Scalar Entry"); 
    ytit.DrawLatex(0.05, 0.17, "Number of Entries per Scalar");
    /*
    c1->cd(2); 
    gPad->SetMargin(0.03, 0.2, 0.01, 0.2); 
    chargeHistosGood[i]->Draw(); 
    stats.DrawLatex(0.5, 0.75, Form("Q = %4.f #muC",chargeHistosGood[i]->Integral()));

    c1->cd(3); 
    gPad->SetMargin(0.2, 0.03, 0.2, 0.01); 
    entryHistos[i]->Draw(); 
    stats.DrawLatex(0.5, 0.95, Form("N = %4.E",entryHistos[i]->Integral()));

    c1->cd(4); 
    gPad->SetMargin(0.03, 0.2, 0.2, 0.01); 
    entryHistosGood[i]->Draw(); 
    stats.DrawLatex(0.5, 0.95, Form("N = %4.E",entryHistosGood[i]->Integral()));
    */

    c1->Print(Form("%sentry_%s.png",imagePath.c_str(),chargeHistos[i]->GetTitle())); 
    c1->Clear(); 


  }

}
