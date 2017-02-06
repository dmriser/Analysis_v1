{

  TFile *file = TFile::Open("refactorTest.root");

  TH2F *ECAllHits = (TH2F*) file->Get("h_ec_fid_allNegatives_all");
  TH2F *ECElectrons = (TH2F*) file->Get("h_ec_fid_EC_FID_all");
  
  TH2F *EDepAllHits = (TH2F*) file->Get("h_etot_p_allNegatives_all");
  TH2F *EDepElectrons = (TH2F*) file->Get("h_etot_p_cuts_all");

  TCanvas *c = new TCanvas("c","",800,800);
  TLatex xCaption, yCaption, title; 

  c->SetMargin(0.12,0.12,0.12,0.12);

  xCaption.SetNDC();
  xCaption.SetTextSize(0.025);

  title.SetNDC();
  title.SetTextSize(0.025);

  yCaption.SetNDC();
  yCaption.SetTextSize(0.025);
  yCaption.SetTextAngle(90.0);

  gStyle->SetPalette(62);
  gPad  ->SetGridx();
  gPad  ->SetGridy();
  gPad  ->SetLogz();

  ECAllHits->GetXaxis()->SetLabelSize(0.02);
  ECElectrons->GetXaxis()->SetLabelSize(0.02);

  ECAllHits->GetYaxis()->SetLabelSize(0.02);
  ECElectrons->GetYaxis()->SetLabelSize(0.02);
 
  ECAllHits->Draw();
  ECElectrons->Draw("colzsame");
  xCaption.DrawLatex(0.7, 0.05, "EC X Coord. [cm]");
  yCaption.DrawLatex(0.05, 0.7, "EC Y Coord. [cm]");
  title.DrawLatex(0.57, 0.94, Form("All hits %.1e ",            ECAllHits->GetEntries()));
  title.DrawLatex(0.57, 0.91, Form("Identified Electrons %.1e ",ECElectrons->GetEntries()));

  c->SaveAs("ECFiducial.png");
  c->Clear();

  EDepAllHits->GetXaxis()->SetLabelSize(0.02);
  EDepElectrons->GetXaxis()->SetLabelSize(0.02);
  EDepAllHits->GetYaxis()->SetLabelSize(0.02);
  EDepElectrons->GetYaxis()->SetLabelSize(0.02);

  gPad->SetLogz(0);
  gPad->SetGridx(0);
  gPad->SetGridy(0);

  EDepElectrons->Draw("colz");
  xCaption.DrawLatex(0.7, 0.05, "p [GeV/c]");
  yCaption.DrawLatex(0.05, 0.7, "E_{dep}/p");
  title.DrawLatex(0.43, 0.86, Form("Identified Electrons %.1e ",EDepElectrons->GetEntries()));

  EDepAllHits->Draw("colz");
  xCaption.DrawLatex(0.7, 0.05, "p [GeV/c]");
  yCaption.DrawLatex(0.05, 0.7, "E_{dep}/p");
  title.DrawLatex(0.43, 0.86, Form("All hits %.1e ",            EDepAllHits->GetEntries()));

  c->SaveAs("SamplingFraction.png");
}
