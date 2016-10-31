{


  TFile * _file0 = TFile::Open("readFilesWithEntryNumber.root");
  TFile * _file1 = TFile::Open("ratesFBF.root");

  TH2F * h2_dQ_dN = (TH2F*) _file0.Get("h2_diff");
  TH1F * h1_ratio = (TH1F*) _file1.Get("h_electron_dQ_ratio_all");

  h1_ratio->SetFillColor(kYellow);

  TF1 * f_ratio = new TF1("f_ratio","gaus",2,24);
  h1_ratio->Fit("f_ratio","rq");

 
  TGraph * g_dN_dQ = (TGraph*) _file1.Get("Graph");
  g_dN_dQ->GetXaxis()->SetTitle(" Run Number ");
  /*
  TCanvas * c1 = new TCanvas("c1","",400,400);
  h2_dQ_dN->Draw("colz");
  c1->Print("h2_dQ_dN.png");
  */

  double NSIGMA = 4.00;
  double mu    = f_ratio->GetParameter(1);
  double sigma = f_ratio->GetParameter(2);

  double left  = mu - NSIGMA*sigma;
  double right = mu + NSIGMA*sigma;

  TLine cutLeft(left,0,left,80);
  TLine cutRight(right,0,right,80);

  TLine cutTop(37600,right,38800,right);
  TLine cutBottom(37600,left,38800,left);

  cutLeft.SetLineWidth(3);
  cutRight.SetLineWidth(3);

  cutTop.SetLineWidth(1);
  cutBottom.SetLineWidth(1);

  cutTop.SetLineColor(kRed);
  cutBottom.SetLineColor(kRed);

  TLatex tex;


  TCanvas * c2 = new TCanvas("c2","",400,400);
  h1_ratio->Draw();
  cutLeft.Draw();
  cutRight.Draw();
  tex.DrawLatex(5,100," #Delta N / #Delta Q For All Runs ");
  c2->Print("h1_dQ_dN.png");
  c2->Clear();

  // graph 
  g_dN_dQ->Draw("ap");
  cutTop.Draw();
  cutBottom.Draw();
  tex.DrawLatex(37850,195," #Delta N/#Delta Q For All Runs ");
  c2->Print("g_dN_dQ.png");



}
