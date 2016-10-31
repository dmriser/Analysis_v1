#ifndef HTML_PRINTER_H
#define HTML_PRINTER_H

// c++ libs
#include <fstream>
#include <iostream>
#include <vector>

// CERN root libs
#include "TCanvas.h"
#include "TROOT.h"
#include "TObject.h"

class HistogramWebsite
{
 public: 
  HistogramWebsite(string, int, int, int, int);
  HistogramWebsite(string, int, int);
  virtual ~HistogramWebsite();

  void Print();
  void AddHistogram(TObject*);
  void SetDrawOptions(string);

 private:
  vector<TObject*> histograms;

  string title;  
  string drawOption;

  int nRows;
  int nColumns;
  int imgWidth;
  int imgHeight;

};

HistogramWebsite::HistogramWebsite(string inTitle, int nrows, int ncols, int w, int h)
{
  title     = inTitle;
  nRows     = nrows;
  nColumns  = ncols;
  imgWidth  = w;
  imgHeight = h;
  drawOption = "";
}

HistogramWebsite::HistogramWebsite(string inTitle, int nrows, int ncols)
{
  title     = inTitle;
  nRows     = nrows;
  nColumns  = ncols;

  // default image size
  imgWidth  = 400;
  imgHeight = 400;

  drawOption = "";
}

HistogramWebsite::~HistogramWebsite()
{}

void HistogramWebsite::SetDrawOptions(string opts)
{
  drawOption = opts;
}

void HistogramWebsite::AddHistogram(TObject* histoPointer)
{
  histograms.push_back(histoPointer);
}

void HistogramWebsite::Print()
{

  cout << endl << " > Writing website " << title << endl;
  cout << " > we have " << histograms.size() << " histograms " << endl;
  // first we need to create the images
  TCanvas * HTMLCanvas = new TCanvas("HTMLCanvas","",imgWidth,imgHeight);
  
  // check for consistancy between histograms 
  // and table size.
  if (histograms.size() != nRows*nColumns)
    {
      cout << " > error - number of histos is not same as size of table " << endl;
      return;
    }

  // print pictures 
  for (int irow = 0; irow < nRows; irow++)
    for (int icol = 0; icol < nColumns; icol++)
      {
	HTMLCanvas->Clear();
	int index = icol + irow*nColumns;
	
	histograms[index]->Draw(drawOption.c_str());
	HTMLCanvas->Print(Form("%s_%d_%d.png",title.c_str(),irow,icol));
      }

  // these parts are seperate because later i will probably
  // split this into 2 classes 
  // ImagePrinter and HTMLTableMaker 
  // or something like that 
  // write html code 
  ofstream outfile;
  outfile.open(Form("%s.html",title.c_str()));

  outfile << "<!DOCTYPE html>" << endl;

  outfile << "<script language=JavaScript>" << endl;

  outfile << "function changePicture(x, y)" << endl;
  outfile << "{" << endl;
  outfile << "      var currentName = \"\";" << endl;
  outfile << "      currentName = \""<< Form("%s_",title.c_str()) << "\" + x + \"_\" + y + \".png\";" << endl;
  outfile << "      document.images.IMAGE.src = currentName;" << endl;
  outfile << "}" << endl;
  outfile << "</script>" << endl;

  outfile << "<html>" << endl;
  outfile << "<body>" << endl;
  outfile << "<h1 align=center>" << title << "</h1>" << endl;
  outfile << "" << endl;
  outfile << "<img src=" << Form("%s_0_0.png",title.c_str()) << " name=IMAGE " << " width=" << imgWidth << ">" << endl;
  outfile << "</p>" << endl;
  outfile << "<table>" << endl;

  // draw table 
  for (int irow = 0; irow < nRows; irow++)
    {
      outfile << "      <tr>" << endl;
    for (int icol = 0; icol < nColumns; icol++)
      {
	int index = icol + irow*nColumns;

	outfile << "            <td onMouseOver=changePicture('" << irow << "','" << icol <<"') >" << endl;
	outfile << "                  " << Form("%d_%d",irow,icol) << endl;
	outfile << "            </td>"    << endl;
      }
    outfile << "      </tr>" << endl;
    }
  outfile << "</table>" << endl;
  outfile << "</body>" << endl;
  outfile << "</html>" << endl;

  outfile.close();
}


#endif
