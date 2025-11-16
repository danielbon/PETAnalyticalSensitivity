// Author: Daniel Bonifacio
// usage inside output directory: root -b -l -q './GenerateSpectrum.C("./output/Sensitivity/BGO/0")'

#include "TTree.h"
#include "TFile.h"
//#include "Riostream.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TPolyMarker.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TProfile.h>

void GenerateSpectrum(const char * outdir)
{
 gStyle->SetOptStat(0);
 gStyle->SetOptTitle(0);
 gStyle->SetCanvasColor(kWhite);     // background is no longer mouse-dropping white
 gStyle->SetPalette(1,0);            // blue to red false color palette. Use 9 for b/w
 gStyle->SetCanvasBorderMode(0);     // turn off canvas borders
 gStyle->SetPadBorderMode(0);
 gStyle->SetPaintTextFormat("5.2f");  // What precision to put numbers if plotted with "TEXT"
 gStyle->SetHistFillColor(kWhite);
 gStyle->SetLineWidth(2.);
 gStyle->SetTextSize(0.8);
 gStyle->SetMarkerSize(2); 
 gStyle->SetLabelSize(0.05,"xy");
 gStyle->SetTitleSize(0.05,"xy");
 gStyle->SetPadTopMargin(0.09);
 gStyle->SetPadRightMargin(0.08);
 gStyle->SetPadBottomMargin(0.12);
 gStyle->SetPadLeftMargin(0.14);
 gStyle->SetStripDecimals(kFALSE);

 ofstream outtxtfile;
 cout.setf(ios::fixed);
 cout.setf(ios::left);
 cout<<setprecision(2);
 cout<<setw(8);
 cout<<"output directory = "<< outdir << endl;
 const char* outputtextfile = "/GenerateSpectrum.txt";
 char * fullPath = new char[strlen(outdir) + strlen(outputtextfile) + 1];
 strcpy(fullPath, outdir);
 strcat(fullPath, outputtextfile);
 outtxtfile.open(fullPath, ios::out);
 if (outtxtfile.is_open())
 {
 }
 else cout << "Unable to open file";
 outtxtfile <<"full Path = "<< fullPath << endl;
 outtxtfile << "Showing results...\n" ;
 Int_t layerID, eventID;
 Float_t globalPosX,globalPosY,globalPosZ;
//TFile *fin1 = new TFile(Form("%s/PETdata.NotOpticalSingles.root",outdir));
 TFile *fin1 = new TFile(Form("%s/output.root",outdir));
 cout << "Opening Singles file ..." << endl;
 if (fin1->IsZombie()) 
 {
  cout << "... unable to open root file" << endl;
  return;
 }
 //TTree *TreeSinglesData = (TTree*)gDirectory->Get("tree");
 TTree *TreeSinglesData = (TTree*)gDirectory->Get("Singles");
 TreeSinglesData->SetBranchAddress("layerID",&layerID);
 TreeSinglesData->SetBranchAddress("globalPosX",&globalPosX);
 TreeSinglesData->SetBranchAddress("globalPosY",&globalPosY);
 TreeSinglesData->SetBranchAddress("globalPosZ",&globalPosZ);
 TreeSinglesData->SetBranchAddress("eventID",&eventID);
// Int_t nEntries = TreeSinglesData->GetEntries("layerID==1");
 Int_t nEntries = TreeSinglesData->GetEntries("");
 outtxtfile<<"number of entries = "<< nEntries << endl;
 outtxtfile<<"min Z pos = "<<TreeSinglesData->GetMinimum("globalPosZ")<<endl;
 outtxtfile<<"max Z pos = "<<TreeSinglesData->GetMaximum("globalPosZ")<<endl;

 TCanvas* Canvas1 = new TCanvas("Canvas1","Canvas1",0,0,800,600);
 gStyle->SetPalette(1,0);
 gStyle->SetErrorX(0);
 TGaxis::SetMaxDigits(3);
//TH1F *hglobalPosZ = (TH1F*)gDirectory->Get("hglobalPosZ");
//TreeSinglesData->Draw("globalPosZ>>hglobalPosZ(200)", "layerID==2","same");
// TreeSinglesData->Draw("globalPosZ", "layerID==1","", nEntries,0);
//TreeSinglesData->Draw("globalPosZ>>hglobalPosZ(220)", "layerID==0","");
TreeSinglesData->SetLineWidth(1.);
TreeSinglesData->SetTitle("");
TreeSinglesData->SetLineColor(kBlue);
TreeSinglesData->Draw("globalPosZ>>hglobalPosZ(220, -33.0, 33.0)", "layerID==0",""); // && globalPosZ>0
TH1F *hglobalPosZ = (TH1F*)gDirectory->Get("hglobalPosZ");
hglobalPosZ->SetTitle(";Z (mm);counts");
TreeSinglesData->SetLineColor(kBlue);
TreeSinglesData->Draw("globalPosZ", "layerID==1","same");
TreeSinglesData->SetLineColor(kBlue);
TreeSinglesData->Draw("globalPosZ", "layerID==2","same");
TreeSinglesData->SetLineColor(kBlue);
TreeSinglesData->Draw("globalPosZ", "layerID==3","same");
TreeSinglesData->SetLineColor(kBlue);
TreeSinglesData->Draw("globalPosZ", "layerID==4","same");
TreeSinglesData->SetLineColor(kRed);
TreeSinglesData->Draw("globalPosZ", "layerID==5","same");
TreeSinglesData->SetLineColor(kRed);
TreeSinglesData->Draw("globalPosZ", "layerID==6","same");
TreeSinglesData->SetLineColor(kRed);
TreeSinglesData->Draw("globalPosZ", "layerID==7","same");
TreeSinglesData->SetLineColor(kRed);
TreeSinglesData->Draw("globalPosZ", "layerID==8","same");
TreeSinglesData->SetLineColor(kRed);
TreeSinglesData->Draw("globalPosZ", "layerID==9","same");

/*
 htemp->SetMarkerStyle(21);
 htemp->SetMarkerSize(1);
 Int_t bin1, bin2;
 bin2 = htemp->FindLastBinAbove(htemp->GetMaximum()/2);
 Int_t binMax = htemp->GetMaximumBin();
 bin1 = binMax - (bin2-binMax)-1;
cout<<" bin1 = "<< bin1 << endl;
cout<<" bin2 = "<< bin2 << endl;
 TF1 *fcharge = new TF1("fcharge", "gaus", htemp->GetBinCenter(bin1-1), htemp->GetBinCenter(bin2+1));
 htemp->Fit("fcharge", "RNQ");
 Double_t meanene = fcharge->GetParameter(1);
 Double_t errmeanene = fcharge->GetParError(1);
 Double_t sigene = fcharge->GetParameter(2);
 Double_t errsigene = fcharge->GetParError(2);
 cout<<" *****Results from our simulated data**** "<< endl;
 cout<<" Mean charge = "<< meanene << "(" << errmeanene << ")" << endl;
 cout<<" FWHM charge = "<< 2.355*sigene << "(" << 2.355*errsigene << ")" << endl;
 cout<<" Charge resolution= "<< 100*2.355*sigene/meanene << "(" << 100*2.355*errsigene/meanene << ")%" << endl;
 fcharge->Draw("same");
 TLegend *leg = new TLegend(0.75,0.75,0.9,0.8);
 leg->SetFillColor(0);
 leg->SetBorderSize(1);
 leg->SetTextSize(0.04);
 leg->AddEntry(htemp,"Na-22","lpf");
 leg->Draw("same");
*/
Canvas1->Update();
Canvas1->Print(Form("%s/layerID.pdf",outdir));

TreeSinglesData->SetLineColor(kBlue);
TreeSinglesData->Draw("globalPosZ>>hglobalPosZ511(220, -33.0, 33.0)", "layerID==0&&energy>0.51",""); // && globalPosZ>0
TH1F *hglobalPosZ511 = (TH1F*)gDirectory->Get("hglobalPosZ511");
hglobalPosZ511->SetTitle(";Z (mm);counts");
TreeSinglesData->SetLineColor(kBlue);
TreeSinglesData->Draw("globalPosZ", "layerID==1&&energy>0.51","same");
TreeSinglesData->SetLineColor(kBlue);
TreeSinglesData->Draw("globalPosZ", "layerID==2&&energy>0.51","same");
TreeSinglesData->SetLineColor(kBlue);
TreeSinglesData->Draw("globalPosZ", "layerID==3&&energy>0.51","same");
TreeSinglesData->SetLineColor(kBlue);
TreeSinglesData->Draw("globalPosZ", "layerID==4&&energy>0.51","same");
TreeSinglesData->SetLineColor(kRed);
TreeSinglesData->Draw("globalPosZ", "layerID==5&&energy>0.51","same");
TreeSinglesData->SetLineColor(kRed);
TreeSinglesData->Draw("globalPosZ", "layerID==6&&energy>0.51","same");
TreeSinglesData->SetLineColor(kRed);
TreeSinglesData->Draw("globalPosZ", "layerID==7&&energy>0.51","same");
TreeSinglesData->SetLineColor(kRed);
TreeSinglesData->Draw("globalPosZ", "layerID==8&&energy>0.51","same");
TreeSinglesData->SetLineColor(kRed);
TreeSinglesData->Draw("globalPosZ", "layerID==9&&energy>0.51","same");
Canvas1->Update();
Canvas1->Print(Form("%s/layerID511keV.pdf",outdir));


Double_t einf = 0.;
Double_t esup = 1.3;
TreeSinglesData->SetLineWidth(2.);
TreeSinglesData->SetLineColor(kGray+2);
TreeSinglesData->Draw("energy>>henergyTotal(200, 0., esup)", "","");
TH1F *henergyTotal = (TH1F*)gDirectory->Get("henergyTotal");
henergyTotal->SetTitle(";energy(MeV);counts");
Int_t binmax = henergyTotal->GetMaximumBin();
outtxtfile<<"Energy Bin Center = "<<henergyTotal->GetXaxis()->GetBinCenter(binmax)<<endl;
outtxtfile<<"Counts at the 511 keV peak = "<<henergyTotal->GetMaximum()<<endl;
Canvas1->Update();
Canvas1->Print(Form("%s/energylayerIDTotalSingles.pdf",outdir));
TreeSinglesData->SetLineColor(kBlue+3);
TreeSinglesData->Draw("energy>>henergyBGO(200, 0., esup)", "layerID==0||layerID==1||layerID==2||layerID==3||layerID==4","");
TH1F *henergyBGO = (TH1F*)gDirectory->Get("henergyBGO");
henergyBGO->SetTitle(";energy(MeV);counts");
Canvas1->Update();
Canvas1->Print(Form("%s/energylayerIDBGO.pdf",outdir));

TH1F *zoomBGO = (TH1F*) henergyBGO->Clone("zoomBGO");
zoomBGO->SetTitle(";energy(MeV);counts");
zoomBGO->SetAxisRange(0, 0.03E6, "Y");
zoomBGO->Draw();
Canvas1->Update();
Canvas1->Print(Form("%s/energylayerIDBGOZoom.pdf",outdir));

TreeSinglesData->SetLineColor(kRed+3);
TreeSinglesData->Draw("energy>>henergyBaF2(200, 0., esup)", "layerID==5||layerID==6||layerID==7||layerID==8||layerID==9","");
TH1F *henergyBaF2 = (TH1F*)gDirectory->Get("henergyBaF2");
henergyBaF2->SetTitle(";energy(MeV);counts");
Canvas1->Update();
Canvas1->Print(Form("%s/energylayerIDBaF2.pdf",outdir));

TH1F *zoomBaF2 = (TH1F*) henergyBaF2->Clone("zoomBaF2");
zoomBaF2->SetTitle(";energy(MeV);counts");
zoomBaF2->SetAxisRange(0, 0.01E6, "Y");
zoomBaF2->Draw();
Canvas1->Update();
Canvas1->Print(Form("%s/energylayerIDBaF2Zoom.pdf",outdir));


//TFile *fin2 = new TFile(Form("%s/PETdata.Coincidences.root",outdir));
TFile *fin2 = new TFile(Form("%s/output.root",outdir));
 cout << "Opening Coincidences file ..." << endl;
 if (fin2->IsZombie()) 
 {
  cout << "... unable to open root file" << endl;
  return;
 }
// TTree *TreeCoincData = (TTree*)gDirectory->Get("tree");
 TTree *TreeCoincData = (TTree*)gDirectory->Get("Coincidences");
 Int_t nCoincEntries = TreeCoincData->GetEntries("");
 outtxtfile<<"number of coincidences entries = "<< nCoincEntries << endl;
TreeCoincData->SetLineWidth(2.);
TreeCoincData->SetLineColor(kGray+2);
TreeCoincData->Draw("energy1>>henergyTotalCoinc1(200, 0., esup)", "","");
TH1F *henergyTotalCoinc1 = (TH1F*)gDirectory->Get("henergyTotalCoinc1");
henergyTotalCoinc1->SetTitle(";energy(MeV);counts");
binmax = henergyTotalCoinc1->GetMaximumBin();
outtxtfile<<"Energy Bin Center Coinc 1 = "<<henergyTotalCoinc1->GetXaxis()->GetBinCenter(binmax)<<endl;
outtxtfile<<"Coinc Counts at the 511 keV peak - Channel 1 = "<<henergyTotalCoinc1->GetMaximum()<<endl;
Canvas1->Update();
Canvas1->Print(Form("%s/energylayerIDTotalCoinc1.pdf",outdir));
TreeCoincData->Draw("energy2>>henergyTotalCoinc2(200, 0., esup)", "","");
TH1F *henergyTotalCoinc2 = (TH1F*)gDirectory->Get("henergyTotalCoinc2");
henergyTotalCoinc2->SetTitle(";energy(MeV);counts");
binmax = henergyTotalCoinc2->GetMaximumBin();
outtxtfile<<"Energy Bin Center Coinc 2 = "<<henergyTotalCoinc2->GetXaxis()->GetBinCenter(binmax)<<endl;
outtxtfile<<"Coinc Counts at the 511 keV peak  - Channel 2 = "<<henergyTotalCoinc2->GetMaximum()<<endl;
outtxtfile<<"Total Coinc Counts at the 511 keV peak = "<<henergyTotalCoinc1->GetMaximum()+henergyTotalCoinc2->GetMaximum()<<endl;
Canvas1->Update();
Canvas1->Print(Form("%s/energylayerIDTotalCoinc2.pdf",outdir));
 outtxtfile.close();

}
