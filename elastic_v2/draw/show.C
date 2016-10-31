#include "../src/histos.h"
#include "../src/histos.cc"

void show_rec()
{
  if (VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",1200,800);
      c1->Divide(3,2);
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  string title = Form(" Reconstructed for %s %s",T[TYPE].c_str(),S[s].c_str());
	  histos.h2_rec[TYPE][s]->SetTitle(title.c_str());
	  histos.h2_rec[TYPE][s]->GetYaxis()->SetTitle(" #theta [deg] ");
	  histos.h2_rec[TYPE][s]->GetXaxis()->SetTitle(" Relative #phi [deg] ");
	  histos.h2_rec[TYPE][s]->Draw("colz");
	}
      if (PRINT != "") c1->Print(Form("../img/h2_rec_%s.%s",T[TYPE].c_str(),PRINT.c_str()));
    }

  if (!VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",800,800);
      string title = Form(" Reconstructed for %s %s",T[TYPE].c_str(),S[SECT].c_str());
      histos.h2_rec[TYPE][SECT]->SetTitle(title.c_str());
      histos.h2_rec[TYPE][SECT]->GetYaxis()->SetTitle(" #theta [deg] ");
      histos.h2_rec[TYPE][SECT]->GetXaxis()->SetTitle(" Relative #phi [deg] ");
      histos.h2_rec[TYPE][SECT]->Draw("colz");
    }
}

void show_acc()
{
  if (VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",1200,800);
      c1->Divide(3,2);
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  string title = Form(" Acceptance for %s %s",T[TYPE].c_str(),S[s].c_str());
	  histos.h2_acc[TYPE][s]->SetTitle(title.c_str());
	  histos.h2_acc[TYPE][s]->GetYaxis()->SetTitle(" #theta [deg] ");
	  histos.h2_acc[TYPE][s]->GetXaxis()->SetTitle(" Relative #phi [deg] ");
	  histos.h2_acc[TYPE][s]->Draw("colz");
	}
      if (PRINT != "") c1->Print(Form("../img/h2_acc_%s.%s",T[TYPE].c_str(),PRINT.c_str()));
    }

  if (!VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",800,800);
      string title = Form(" Acceptance for %s %s",T[TYPE].c_str(),S[SECT].c_str());
      histos.h2_acc[TYPE][SECT]->SetTitle(title.c_str());
      histos.h2_acc[TYPE][SECT]->GetYaxis()->SetTitle(" #theta [deg] ");
      histos.h2_acc[TYPE][SECT]->GetXaxis()->SetTitle(" Relative #phi [deg] ");
      histos.h2_acc[TYPE][SECT]->Draw("colz");
    }
}

void show_rc()
{

      TCanvas *c1 = new TCanvas("c1","",800,800);
      histos.h1_rc->SetTitle(" Radiative Correction ");
      histos.h1_rc->SetMarkerStyle(2);
      histos.h1_rc->SetMarkerColor(kRed);
      //histos.h1_rc->SetLineColor(kRed);
      histos.h1_rc->GetXaxis()->SetTitle(" #theta [deg] ");
      histos.h1_rc->GetYaxis()->SetTitle(" #sigma_{rad}/#sigma_{born} ");
      histos.h1_rc->Draw("P");
      //      histos.h1_rc->Draw("L");

      if (PRINT != "") c1->Print(Form("../img/h1_rc.%s",PRINT.c_str()));

}

void show_xs()
{

  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.06);
  gStyle->SetPadRightMargin(0.06);

  // Latex Printer 
  TLatex lab;
  lab.SetNDC();
  lab.SetTextFont(32);
  lab.SetTextSize(0.052);
  lab.SetTextColor(kBlack);

  histos.h1_xs_model[0]->SetLineWidth(2);
  histos.h1_xs_model[0]->SetLineStyle(4);
  histos.h1_xs_model[0]->SetLineColorAlpha(kBlue,1.0);

  histos.h1_xs_model[1]->SetLineWidth(2);
  histos.h1_xs_model[1]->SetLineStyle(4);
  histos.h1_xs_model[1]->SetLineColorAlpha(kRed,1.0);

  for (int s=0; s<7; s++) histos.h1_xs[s]->SetMarkerStyle(2);

  if (VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",1200,800);
      c1->SetTopMargin(0.15);
      c1->Divide(3,2);
      
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  //	  string title = Form(" Cross Section for %s",S[s].c_str());
	  //	  histos.h1_xs[s]->SetTitle(title.c_str());
	  histos.h1_xs[s]->SetTitle("");
	  histos.h1_xs[s]->GetXaxis()->SetTitle(" #theta [deg] ");
	  //	  histos.h1_xs[s]->GetYaxis()->SetTitle(" #frac{d#sigma}{d#theta} [#mub] ");
	  histos.h1_xs[s]->Draw("PE");
	  //	  histos.h1_xs_model[0]->Draw("Psame");
	  //	  histos.h1_xs_model[1]->Draw("Psame");
	  histos.h1_xs_model[0]->Draw("Lsame");
	  histos.h1_xs_model[1]->Draw("Lsame");
	  int nev = histos.h2_rec[2][s]->GetEntries();
	  lab.DrawLatex(0.45, 0.01, Form(" Sector %d ",s) );
	  lab.DrawLatex(0.6, 0.78, Form(" %d events ", nev));
	}

      /*      
      TLegend * leg = new TLegend(0.65,0.7,0.95,0.9);
      leg->AddEntry(histos.h1_xs[0]," Cross Section ","pe");
      leg->AddEntry(histos.h1_xs_model[0]," Model, NO RAD ","L");
      leg->AddEntry(histos.h1_xs_model[1]," Model, RAD ","L");
      leg->Draw();
      */

      c1->cd(1);
      lab.DrawLatex(0.2, 0.9, " Elastic Cross Section (#mubarns) ");

      c1->cd(3);
      lab.SetTextColor(kBlue);
      lab.DrawLatex(0.2, 0.9, " #rightarrow Born  ");
      lab.SetTextColor(kRed);
      lab.DrawLatex(0.4, 0.9, " #rightarrow Born + Bethe-Heitler  ");

      if (PRINT != "") c1->Print(Form("../img/h1_xs.%s",PRINT.c_str()));
    }

  if (!VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",800,800);
      string title = Form(" Cross Section for %s",S[SECT].c_str());
      //      histos.h1_xs[SECT]->SetTitle(title.c_str());	  
      histos.h1_xs[SECT]->GetXaxis()->SetTitle(" #theta [deg] ");
      histos.h1_xs[SECT]->GetYaxis()->SetTitle(" #frac{d#sigma}{d#theta} [#mub] ");
      histos.h1_xs[SECT]->Draw("PE");
      histos.h1_xs_model[0]->Draw("Lsame");
      histos.h1_xs_model[1]->Draw("Lsame");
      TLegend * leg = new TLegend(0.65,0.7,0.95,0.9);
      leg->AddEntry(histos.h1_xs[0]," Cross Section ","pe");
      leg->AddEntry(histos.h1_xs_model[0]," Model, NO RAD ","L");
      leg->AddEntry(histos.h1_xs_model[1]," Model, RAD ","L");
      leg->Draw();

    }
}

void show_xs_raw()
{

  gStyle->SetPadTopMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.06);
  gStyle->SetPadRightMargin(0.06);

  // Latex Printer 
  TLatex lab;
  lab.SetNDC();
  lab.SetTextFont(32);
  lab.SetTextSize(0.052);
  lab.SetTextColor(kBlack);

  histos.h1_xs_model[0]->SetMarkerStyle(4);
  histos.h1_xs_model[1]->SetMarkerStyle(4);
  histos.h1_xs_model[0]->SetMarkerColor(kBlue);
  histos.h1_xs_model[1]->SetMarkerColor(kRed);

  for (int s=0; s<7; s++) histos.h1_xs_raw[s]->SetMarkerStyle(2);

  if (VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",1200,800);
      c1->Divide(3,2);
      
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  string title = Form(" Raw Cross Section for %s",S[s].c_str());
	  histos.h1_xs_raw[s]->SetTitle(title.c_str());
	  histos.h1_xs_raw[s]->GetXaxis()->SetTitle(" #theta [deg] ");
	  histos.h1_xs_raw[s]->GetYaxis()->SetTitle(" #frac{d#sigma}{d#theta} [#mub] ");
	  histos.h1_xs_raw[s]->Draw("P");
	  histos.h1_xs_model[0]->Draw("Psame");
	  histos.h1_xs_model[1]->Draw("Psame");
	  lab.DrawLatex(0.3, 0.01, Form(" Sector %d ", s) );
	}
      
      TLegend * leg = new TLegend(0.65,0.7,0.95,0.9);
      leg->AddEntry(histos.h1_xs_raw[0]," Raw Cross Section ","p");
      leg->AddEntry(histos.h1_xs_model[0]," Model, NO RAD ","p");
      leg->AddEntry(histos.h1_xs_model[1]," Model, RAD ","p");
      leg->Draw();

      if (PRINT != "") c1->Print(Form("../img/h1_xs_raw.%s",PRINT.c_str()));
    }

  if (!VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",800,800);
      string title = Form(" Raw Cross Section for %s",S[SECT].c_str());
      histos.h1_xs_raw[SECT]->SetTitle(title.c_str());	  
      histos.h1_xs_raw[SECT]->GetXaxis()->SetTitle(" #theta [deg] ");
      histos.h1_xs_raw[SECT]->GetYaxis()->SetTitle(" #frac{d#sigma}{d#theta} [#mub] ");
      histos.h1_xs_raw[SECT]->Draw("P");
      histos.h1_xs_model[0]->Draw("Psame");
      histos.h1_xs_model[1]->Draw("Psame");
      TLegend * leg = new TLegend(0.65,0.7,0.95,0.9);
      leg->AddEntry(histos.h1_xs_raw[0]," Raw Cross Section ","p");
      leg->AddEntry(histos.h1_xs_model[0]," Model, NO RAD ","p");
      leg->AddEntry(histos.h1_xs_model[1]," Model, RAD ","p");
      leg->Draw();

    }
}

void show_rec_theta()
{
  if (VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",1200,800);
      c1->Divide(3,2);
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  string title = Form(" Reconstructed for %s %s",T[TYPE].c_str(),S[s].c_str());
	  histos.h1_rec_theta[TYPE][s]->SetTitle(title.c_str());
	  histos.h1_rec_theta[TYPE][s]->GetXaxis()->SetTitle(" #theta ");
	  histos.h1_rec_theta[TYPE][s]->Draw();
	}
      if (PRINT != "") c1->Print(Form("../img/h1_rec_theta_%s.%s",T[TYPE].c_str(),PRINT.c_str()));
    }

  if (!VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",800,800);
      string title = Form(" Reconstructed for %s %s",T[TYPE].c_str(),S[SECT].c_str());
      histos.h1_rec_theta[TYPE][SECT]->SetTitle(title.c_str());
      histos.h1_rec_theta[TYPE][SECT]->GetXaxis()->SetTitle(" #theta ");
      histos.h1_rec_theta[TYPE][SECT]->Draw();
    }
}

void show_gen_theta()
{
  if (VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",1200,800);
      c1->Divide(3,2);
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  string title = Form(" Generated for %s %s",T[TYPE].c_str(),S[s].c_str());
	  histos.h1_gen_theta[TYPE][s]->SetTitle(title.c_str());
	  histos.h1_gen_theta[TYPE][s]->GetXaxis()->SetTitle(" #theta ");
	  histos.h1_gen_theta[TYPE][s]->Draw();
	}
      if (PRINT != "") c1->Print(Form("../img/h1_gen_theta_%s.%s",T[TYPE].c_str(),PRINT.c_str()));
    }


  if (!VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",800,800);
      string title = Form(" Generated for %s %s",T[TYPE].c_str(),S[SECT].c_str());
      histos.h1_gen_theta[TYPE][SECT]->SetTitle(title.c_str());
      histos.h1_gen_theta[TYPE][SECT]->GetXaxis()->SetTitle(" #theta ");
      histos.h1_gen_theta[TYPE][SECT]->Draw();
    }
}

void show_gen()
{
  if (VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",1200,800);
      c1->Divide(3,2);
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  string title = Form(" Generated for %s %s",T[TYPE].c_str(),S[s].c_str());
	  histos.h2_gen[TYPE][s]->SetTitle(title.c_str());
	  histos.h2_gen[TYPE][s]->Draw("colz");
	}
      if (PRINT != "") c1->Print(Form("../img/h2_gen_%s.%s",T[TYPE].c_str(),PRINT.c_str()));
    }

  if (!VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",800,800);
      string title = Form(" Generated for %s %s",T[TYPE].c_str(),S[SECT].c_str());
      histos.h2_gen[TYPE][SECT]->SetTitle(title.c_str());
      histos.h2_gen[TYPE][SECT]->Draw("colz");
    }
}

void show_x()
{
  if (VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",1200,800);
      c1->Divide(3,2);
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  string title = Form(" Bjorken x for %s %s",T[TYPE].c_str(),S[s].c_str());
	  histos.h1_x[TYPE][s]->SetTitle(title.c_str());
	  histos.h1_x[TYPE][s]->Draw();
	}
      if (PRINT != "") c1->Print(Form("../img/h1_x_%s.%s",T[TYPE].c_str(),PRINT.c_str()));      
    }

  if (!VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",800,800);
      string title = Form(" x for %s %s",T[TYPE].c_str(),S[SECT].c_str());
      histos.h1_x[TYPE][SECT]->SetTitle(title.c_str());
      histos.h1_x[TYPE][SECT]->Draw();
    }
}

void show_W()
{
  if (VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",1200,800);
      c1->Divide(3,2);
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  string title = Form(" W for %s %s",T[TYPE].c_str(),S[s].c_str());
	  histos.h1_W[TYPE][s]->SetTitle(title.c_str());
	  histos.h1_W[TYPE][s]->Draw();
	}
      if (PRINT != "") c1->Print(Form("../img/h1_W_%s.%s",T[TYPE].c_str(),PRINT.c_str()));      
    }

  if (!VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",800,800);
      string title = Form(" W for %s %s",T[TYPE].c_str(),S[SECT].c_str());
      histos.h1_W[TYPE][SECT]->SetTitle(title.c_str());
      histos.h1_W[TYPE][SECT]->Draw();
    }
}

void show_QQ()
{
  if (VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",1200,800);
      c1->Divide(3,2);
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  string title = Form(" Q^2 for %s %s",T[TYPE].c_str(),S[s].c_str());
	  histos.h1_QQ[TYPE][s]->SetTitle(title.c_str());
	  histos.h1_QQ[TYPE][s]->Draw();
	}
      if (PRINT != "") c1->Print(Form("../img/h1_QQ_%s.%s",T[TYPE].c_str(),PRINT.c_str()));      
    }

  if (!VIEW)
    {
  TCanvas *c1 = new TCanvas("c1","",800,800);
  string title = Form(" #Q^2 for %s %s",T[TYPE].c_str(),S[SECT].c_str());
  histos.h1_QQ[TYPE][SECT]->SetTitle(title.c_str());
  histos.h1_QQ[TYPE][SECT]->Draw();
    }
}

void show_MM2()
{
  if (VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",1200,800);
      c1->Divide(3,2);
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  string title = Form(" Missing Mass ^2 for %s %s",T[TYPE].c_str(),S[s].c_str());
	  histos.h1_MM2[TYPE][s]->SetTitle(title.c_str());
	  histos.h1_MM2[TYPE][s]->Draw();
	}
      if (PRINT != "") c1->Print(Form("../img/h1_MM2_%s.%s",T[TYPE].c_str(),PRINT.c_str()));      
    }

  if (!VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",800,800);
      string title = Form(" #MM^2 for %s %s",T[TYPE].c_str(),S[SECT].c_str());
      histos.h1_MM2[TYPE][SECT]->SetTitle(title.c_str());
      histos.h1_MM2[TYPE][SECT]->Draw();
    }
}

void show_EP()
{  
  if (VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",1200,800);
      c1->Divide(3,2);
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  string title = Form(" Angle between eP for %s %s",T[TYPE].c_str(),S[s].c_str());
	  histos.h1_EP_angle[TYPE][s]->SetTitle(title.c_str());
	  histos.h1_EP_angle[TYPE][s]->Draw();
	}
      if (PRINT != "") c1->Print(Form("../img/h1_EP_angle_%s.%s",T[TYPE].c_str(),PRINT.c_str()));      
    }

  if (!VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",800,800);
      string title = Form(" EP Angle for %s %s",T[TYPE].c_str(),S[SECT].c_str());
      histos.h1_EP_angle[TYPE][SECT]->SetTitle(title.c_str());
      histos.h1_EP_angle[TYPE][SECT]->Draw();
    }
}

void show_norm_theta()
{
  for (int s=0; s<7; s++)
    { 
      histos.h1_norm_theta[0][2][s]->SetLineColor(kRed);  
      histos.h1_norm_theta[0][0][s]->SetLineColor(kBlue);
      histos.h1_norm_theta[0][1][s]->SetLineColor(kGreen);
      histos.h1_norm_theta[1][0][s]->SetLineColor(kYellow);  
      histos.h1_norm_theta[1][1][s]->SetLineColor(kCyan);
    }

  if (VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",1200,800);
      c1->Divide(3,2);
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  string title = Form(" Normalized Theta for %s",S[s].c_str());
	  histos.h1_norm_theta[0][2][s]->SetTitle(title.c_str());
	  histos.h1_norm_theta[0][2][s]->Draw();
	  histos.h1_norm_theta[0][0][s]->Draw("same");
	  histos.h1_norm_theta[0][1][s]->Draw("same");
	  histos.h1_norm_theta[1][0][s]->Draw("same");
	  histos.h1_norm_theta[1][1][s]->Draw("same");
	}

      TLegend * leg = new TLegend(0.65,0.7,0.95,0.9);
      leg->AddEntry(histos.h1_norm_theta[0][2][0]," Data ","l");
      leg->AddEntry(histos.h1_norm_theta[0][1][0]," Rec (NORAD) ","l");
      leg->AddEntry(histos.h1_norm_theta[1][1][0]," Rec (RAD) ","l");
      leg->AddEntry(histos.h1_norm_theta[0][0][0]," Gen (NORAD) ","l");
      leg->AddEntry(histos.h1_norm_theta[1][0][0]," Gen (RAD) ","l");
      leg->Draw();

      if (PRINT != "") c1->Print(Form("../img/h1_norm_theta.%s",PRINT.c_str()));
    }

  else
    {
      TCanvas *c1 = new TCanvas("c1","",800,800);
      string title = Form(" Normalized Theta for %s",S[SECT].c_str());
      histos.h1_norm_theta[0][2][SECT]->SetTitle(title.c_str());
      histos.h1_norm_theta[0][2][SECT]->Draw();
      histos.h1_norm_theta[0][0][SECT]->Draw("same");
      histos.h1_norm_theta[0][1][SECT]->Draw("same");
      histos.h1_norm_theta[1][0][SECT]->Draw("same");
      histos.h1_norm_theta[1][1][SECT]->Draw("same");

      TLegend * leg = new TLegend(0.65,0.7,0.95,0.9);
      leg->AddEntry(histos.h1_norm_theta[0][2][0]," Data ","l");
      leg->AddEntry(histos.h1_norm_theta[0][1][0]," Rec (NORAD) ","l");
      leg->AddEntry(histos.h1_norm_theta[1][1][0]," Rec (RAD) ","l");
      leg->AddEntry(histos.h1_norm_theta[0][0][0]," Gen (NORAD) ","l");
      leg->AddEntry(histos.h1_norm_theta[1][0][0]," Gen (RAD) ","l");
      leg->Draw();

    }


}

void show_xs_ratio()
{
  TLine line;
  line.SetLineWidth(2);
  line.SetLineStyle(2);

  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.06);
  gStyle->SetPadRightMargin(0.06);

  // Latex Printer 
  TLatex lab;
  lab.SetNDC();
  lab.SetTextFont(32);
  lab.SetTextSize(0.052);
  lab.SetTextColor(kBlack);

  
  for (int s=0; s<7; s++){ 
    histos.h1_xs_ratio[0][s]->SetMarkerStyle(3);  
    histos.h1_xs_ratio[0][s]->SetMarkerColor(kBlue);  
    histos.h1_xs_ratio[1][s]->SetMarkerStyle(3);
    histos.h1_xs_ratio[1][s]->SetMarkerColor(kRed);  
  }

  if (VIEW)
    {
      TCanvas *c1 = new TCanvas("c1","",1200,800);
      c1->Divide(3,2);
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  //	  string title = Form(" Cross Section Ratio for %s",S[s].c_str());
	  string title = "";
	  histos.h1_xs_ratio[0][s]->SetTitle(title.c_str());


	  histos.h1_xs_ratio[0][s]->Draw("P");
	  histos.h1_xs_ratio[1][s]->Draw("Psame");

	  line.DrawLine(thetaBins.min(),1,thetaBins.max(),1);
	  lab.DrawLatex(0.35, 0.01, Form(" Sector %d   #theta [deg]", s) );
	}

      c1->cd(2);
      lab.DrawLatex(0.3, 0.9, " Cross Section Ratio ");

      c1->cd(3);
      lab.SetTextColor(kBlue);
      lab.DrawLatex(0.2, 0.9, " #rightarrow Born  ");
      lab.SetTextColor(kRed);
      lab.DrawLatex(0.4, 0.9, " #rightarrow Born + Bethe-Heitler  ");

      if (PRINT != "") c1->Print(Form("../img/h1_xs_ratio_born.%s",PRINT.c_str()));
    }

  else
    {
      TCanvas *c1 = new TCanvas("c1","",800,800);
      string title = Form(" Cross Section Ratio for %s",S[SECT].c_str());
      histos.h1_xs[SECT]          ->SetTitle(title.c_str());
      histos.h1_xs_ratio[0][SECT] ->Draw("P");
      histos.h1_xs_ratio[1][SECT] ->Draw("Psame");
      line.DrawLine(thetaBins.min(),1,thetaBins.max(),1);
    }

  /*  
  if (VIEW)
    {
      TCanvas *c2 = new TCanvas("c2","",1200,800);
      c1->Divide(3,2);
      for (int s=1; s<7; s++) 
	{
	  c1->cd(s);
	  string title = Form(" Cross Section Ratio for %s",S[s].c_str());
	  histos.h1_xs_ratio[1][s]->SetTitle(title.c_str());
	  histos.h1_xs_ratio[1][s]->Draw("P");
	}

      if (PRINT != "") c1->Print(Form("../img/h1_xs_ratio_rad.%s",PRINT.c_str()));
    }

  else
    {
      TCanvas *c2 = new TCanvas("c2","",800,800);
      string title = Form(" Cross Section Ratio for %s",S[SECT].c_str());
      histos.h1_xs[SECT]->SetTitle(title.c_str());
      histos.h1_xs_ratio[1][SECT]->Draw("PE");
    }
  */
}

void show_fc()
{
  TCanvas *c1 = new TCanvas("c1","",800,800);
  histos.h1_fc_mon->SetTitle(" Faraday Cup Monitor ");
  histos.h1_fc_mon->GetXaxis()->SetTitle(" Charge/run [#muC]");
  histos.h1_fc_mon->SetFillColor(kRed);
  histos.h1_fc_mon->Draw();
}

void show_W_comparison()
{
  TCanvas *c1 = new TCanvas("c1","",1200,800);
  c1->Divide(3,2);
  for (int s=1; s<7; s++) 
    {
      c1->cd(s);
      string title = Form(" W for %s",S[s].c_str());
      histos.h1_W[TYPE][s]->SetTitle(title.c_str());
      histos.h1_W[1][s]->SetLineColor(kRed);
      histos.h1_W[2][s]->SetLineColor(kBlue);
      histos.h1_W[1][s]->Scale(1/histos.h1_W[1][s]->Integral());
      histos.h1_W[2][s]->Scale(1/histos.h1_W[2][s]->Integral());
      histos.h1_W[1][s]->Draw();
      histos.h1_W[2][s]->Draw("same");
    }
  if (PRINT != "") c1->Print(Form("../img/h1_W_comparison.%s",PRINT.c_str()));      
  
}
