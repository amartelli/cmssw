#include "TLegend.h"
#include "TLatex.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TTree.h"
#include "TChain.h"
#include <vector>
#include <fstream>
#include <string>
#include "TROOT.h"
#include "TSystem.h"



void doPlots_Fit_all(){
  gROOT->Macro("/afs/cern.ch/user/a/amartell/public/setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gROOT->Reset();
  gROOT->Macro("~/public/setStyle.C");
  //  gROOT->LoadMacro("~/public/myStyle.C");
  //gROOT->LoadMacro("~/public/setStyle.C");
  //gROOT->LoadMacro("~/public/rootLogon.C");
  //gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  /*
   gStyle->SetTitleOffset(1.1,"x");//X-axis title offset from axis
   gStyle->SetTitleOffset(1.1,"y");//X-axis title offset from axis                                                               
   gStyle->SetTitleSize(0.5,"x");//X-axis title size                         
   gStyle->SetTitleSize(0.5,"y");
   gStyle->SetTitleSize(0.05,"z");
   gStyle->SetLabelOffset(0.025);
  */

  //  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  // gPad->SetGridx();
  // gPad->SetGridy();


    std::string partType = "Pho";
  //  std::string partType = "Pion";

  //  std::cout << " ene = " << ene << std::endl; 

  int iColors[3] = {kGreen+1, kBlue, kRed};


  TGraphErrors* tg[3];
  TGraphErrors* tgR[3];
  TGraphErrors* tgRstd[3];
  for(int iT=0; iT<3; ++iT){
    tg[iT] = new TGraphErrors();
    tg[iT]->SetName(Form("mean_Si%d", (iT+1)*100));
    tg[iT]->SetPoint(0, -1, -1);

    tg[iT]->SetLineColor(iColors[iT]);
    tg[iT]->SetMarkerColor(iColors[iT]);
    tg[iT]->SetMarkerStyle(20);

    tgR[iT] = new TGraphErrors();
    tgR[iT]->SetName(Form("resolution_Si%d", (iT+1)*100));
    tgR[iT]->SetPoint(0, -1, -1);

    tgR[iT]->SetLineColor(iColors[iT]);
    tgR[iT]->SetMarkerColor(iColors[iT]);
    tgR[iT]->SetMarkerStyle(20);

    tgRstd[iT] = new TGraphErrors();
    tgRstd[iT]->SetName(Form("resolution_Si%d", (iT+1)*100));
    tgRstd[iT]->SetPoint(0, -1, -1);

    tgRstd[iT]->SetLineColor(iColors[iT]);
    tgRstd[iT]->SetMarkerColor(iColors[iT]);
    tgRstd[iT]->SetMarkerStyle(20);
  }


  //  TFile *inF = TFile::Open("hydra_CalibOnly_Pho35.root");
  TFile* inF[6];
  inF[0] = TFile::Open(("../test/jobs/Calib_"+partType+"5.root").c_str());
  inF[1] = TFile::Open(("../test/jobs/Calib_"+partType+"10.root").c_str());
  inF[2] = TFile::Open(("../test/jobs/Calib_"+partType+"30.root").c_str());
  inF[3] = TFile::Open(("../test/jobs/Calib_"+partType+"60.root").c_str());
  inF[4] = TFile::Open(("../test/jobs/Calib_"+partType+"100.root").c_str());
  inF[5] = TFile::Open(("../test/jobs/Calib_"+partType+"200.root").c_str());

  //vs gen                                                                                                                                                                      

  TH1F* h_EoP_hitTrack_100[6];
  TH1F* h_EoP_hitTrack_200[6];
  TH1F* h_EoP_hitTrack_300[6];

  TProfile* h_absorberEnergy_vsLayer[6];
  TProfile* h_energy_vsLayer[6];
  TH1F* h_energyAll_vsLayer[6];

  float xMin_a = 0.8;
  float xMin_b = 0.8;
  float xMin_c = 0.8;
  float xMax = 1.4;
  if(partType == "Pion") {
    xMin_a = 0.2;
    xMin_b = 0.8;
    xMin_c = 0.9;
  }

  TF1* fitF_100[6];
  fitF_100[0] = new TF1("Si100_E5", "gaus", xMin_a, xMax );
  fitF_100[1] = new TF1("Si100_E10", "gaus", xMin_a, xMax );
  fitF_100[2] = new TF1("Si100_E30", "gaus", xMin_b, xMax );
  fitF_100[3] = new TF1("Si100_E60", "gaus", xMin_b, xMax );
  fitF_100[4] = new TF1("Si100_E100", "gaus", xMin_b, xMax );
  fitF_100[5] = new TF1("Si100_E200", "gaus", xMin_c, xMax );
  //  fitF_60->SetParameters(600, 0.92, 1.5);
  TF1* fitF_200[6];
  fitF_200[0] = new TF1("Si200_E5", "gaus", xMin_a, xMax );
  fitF_200[1] = new TF1("Si200_E10", "gaus", xMin_a, xMax );
  fitF_200[2] = new TF1("Si200_E30", "gaus", xMin_b, xMax );
  fitF_200[3] = new TF1("Si200_E60", "gaus", xMin_b, xMax );
  fitF_200[4] = new TF1("Si200_E100", "gaus", xMin_b, xMax );
  fitF_200[5] = new TF1("Si200_E200", "gaus", xMin_c, xMax );
  //
  TF1* fitF_300[6];
  fitF_300[0] = new TF1("Si300_E5", "gaus", xMin_a, xMax );
  fitF_300[1] = new TF1("Si300_E10", "gaus", xMin_a, xMax );
  fitF_300[2] = new TF1("Si300_E30", "gaus", xMin_b, xMax );
  fitF_300[3] = new TF1("Si300_E60", "gaus", xMin_b, xMax );
  fitF_300[4] = new TF1("Si300_E100", "gaus", xMin_b, xMax );
  fitF_300[5] = new TF1("Si300_E200", "gaus", xMin_c, xMax );


  int colorsExt[6] = {kRed, kViolet, kCyan+2, kBlue, kGreen+1, kYellow+2};
  //  std::vector<int> colors;

  std::cout << " >>> ora prendo histo " << std::endl;

  for(int iC=0; iC<6; ++iC){
    // h_absorberEnergy_vsLayer[iC] = (TProfile*)(inF[iC]->Get("ana/h_EoP_CPene_1%d_calib_fractionh_absorberEnergy_vsLayer"));
    // h_energy_vsLayer[iC] = (TProfile*)(inF[iC]->Get("HydraCaloCalibrator/h_energy_vsLayer"));
    // h_energyAll_vsLayer[iC] = new TH1F("h_energyAll_vsLayer", "", 50, 0., 50.);
    // h_energyAll_vsLayer[iC]->Add((TH1F*)h_energy_vsLayer[iC], (TH1F*)h_absorberEnergy_vsLayer[iC], 1, 1);


    h_EoP_hitTrack_100[iC] = (TH1F*)(inF[iC]->Get("ana/h_EoP_CPene_100_calib_fraction"));
    h_EoP_hitTrack_200[iC] = (TH1F*)(inF[iC]->Get("ana/h_EoP_CPene_200_calib_fraction"));
    h_EoP_hitTrack_300[iC] = (TH1F*)(inF[iC]->Get("ana/h_EoP_CPene_300_calib_fraction"));

    if(partType == "Pho"){
      h_EoP_hitTrack_100[iC]->Rebin(4);
      h_EoP_hitTrack_200[iC]->Rebin(4);
      h_EoP_hitTrack_300[iC]->Rebin(4);
    }
    else{
      h_EoP_hitTrack_100[iC]->Rebin(10);
      h_EoP_hitTrack_200[iC]->Rebin(10);
      h_EoP_hitTrack_300[iC]->Rebin(10);
    }

    //    colors.push_back(colorsExt[iC]);
    h_EoP_hitTrack_100[iC]->SetLineColor(colorsExt[iC]);
    h_EoP_hitTrack_200[iC]->SetLineColor(colorsExt[iC]);
    h_EoP_hitTrack_300[iC]->SetLineColor(colorsExt[iC]);

    /*
    h_absorberEnergy_vsLayer[iC]->SetLineColor(colorsExt[iC]);
    h_energy_vsLayer[iC]->SetLineColor(colorsExt[iC]);
    h_energyAll_vsLayer[iC]->SetLineColor(colorsExt[iC]);

    h_absorberEnergy_vsLayer[iC]->SetMarkerColor(colorsExt[iC]);
    h_energy_vsLayer[iC]->SetMarkerColor(colorsExt[iC]);
    h_energyAll_vsLayer[iC]->SetMarkerColor(colorsExt[iC]);

    h_absorberEnergy_vsLayer[iC]->SetMarkerStyle(20);
    h_energy_vsLayer[iC]->SetMarkerStyle(20);
    h_energyAll_vsLayer[iC]->SetMarkerStyle(20);
    */

    h_EoP_hitTrack_100[iC]->SetLineWidth(2);
    h_EoP_hitTrack_200[iC]->SetLineWidth(2);
    h_EoP_hitTrack_300[iC]->SetLineWidth(2);

    /*
    h_absorberEnergy_vsLayer[iC]->SetLineWidth(2);
    h_energy_vsLayer[iC]->SetLineWidth(2);
    h_energyAll_vsLayer[iC]->SetLineWidth(2);
    */

    fitF_100[iC]->SetParameters(1000, 1.2, 0.005);
    fitF_200[iC]->SetParameters(1000, 1.2, 0.005);
    fitF_300[iC]->SetParameters(1000, 1.2, 0.005);
    fitF_100[iC]->SetLineColor(colorsExt[iC]);
    fitF_200[iC]->SetLineColor(colorsExt[iC]);
    fitF_300[iC]->SetLineColor(colorsExt[iC]);

    fitF_100[iC]->SetNpx(5000);
    fitF_200[iC]->SetNpx(5000);
    fitF_300[iC]->SetNpx(5000);
    fitF_100[iC]->SetNpx(5000);
    fitF_200[iC]->SetNpx(5000);
    fitF_300[iC]->SetNpx(5000);
  }




  std::cout << " ci sono " << std::endl;

  TLegend *legTGM = new TLegend(0.70,0.30,0.90,0.50,NULL,"brNDC");
  legTGM->SetTextFont(42);
  legTGM->SetFillColor(kWhite);
  legTGM->SetLineColor(kWhite);
  legTGM->SetShadowColor(kWhite);
  legTGM->AddEntry(tg[0], "Si 100 #mum", "p");
  legTGM->AddEntry(tg[1], "Si 200 #mum", "p");
  legTGM->AddEntry(tg[2], "Si 300 #mum", "p");


  TLegend *leg = new TLegend(0.20,0.30,0.35,0.65,NULL,"brNDC");
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  leg->AddEntry(h_EoP_hitTrack_100[0], "#gamma E=5GeV", "l");
  leg->AddEntry(h_EoP_hitTrack_100[1], "#gamma E=10GeV", "l");
  leg->AddEntry(h_EoP_hitTrack_100[2], "#gamma E=30GeV", "l");
  leg->AddEntry(h_EoP_hitTrack_100[3], "#gamma E=60GeV", "l");
  leg->AddEntry(h_EoP_hitTrack_100[4], "#gamma E=100GeV", "l");
  leg->AddEntry(h_EoP_hitTrack_100[5], "#gamma E=200GeV", "l");



  TLatex t10;
  t10.SetNDC();
  t10.SetTextSize(0.03);
  t10.SetTextFont(132);
  t10.SetTextColor(kRed);

  TLatex t30;
  t30.SetNDC();
  t30.SetTextSize(0.03);
  t30.SetTextFont(132);
  t30.SetTextColor(kViolet);


  TLatex t3a;
  t3a.SetNDC();
  t3a.SetTextSize(0.03);
  t3a.SetTextFont(132);
  t3a.SetTextColor(kCyan+2);

  h_EoP_hitTrack_100[0]->Fit("Si100_E5", "R");
  h_EoP_hitTrack_100[1]->Fit("Si100_E10", "R");
  h_EoP_hitTrack_100[2]->Fit("Si100_E30", "R");
  h_EoP_hitTrack_100[3]->Fit("Si100_E60", "R");
  h_EoP_hitTrack_100[4]->Fit("Si100_E100", "R");
  h_EoP_hitTrack_100[5]->Fit("Si100_E200", "R");

  std::cout << " passa a 200 " << std::endl;
  TLatex t1a;
  t1a.SetNDC();
  t1a.SetTextSize(0.03);
  t1a.SetTextFont(132);
  t1a.SetTextColor(kBlue);

  h_EoP_hitTrack_200[0]->Fit("Si200_E5", "R");
  h_EoP_hitTrack_200[1]->Fit("Si200_E10", "R");
  h_EoP_hitTrack_200[2]->Fit("Si200_E30", "R");
  h_EoP_hitTrack_200[3]->Fit("Si200_E60", "R");
  h_EoP_hitTrack_200[4]->Fit("Si200_E100", "R");
  h_EoP_hitTrack_200[5]->Fit("Si200_E200", "R");

  std::cout << " passa a 300 " << std::endl;
  TLatex t3p;
  t3p.SetNDC();
  t3p.SetTextSize(0.03);
  t3p.SetTextFont(132);
  t3p.SetTextColor(kGreen+1);

  TLatex t1p;
  t1p.SetNDC();
  t1p.SetTextSize(0.03);
  t1p.SetTextFont(132);
  t1p.SetTextColor(kYellow+2);

  h_EoP_hitTrack_300[0]->Fit("Si300_E5", "R");
  h_EoP_hitTrack_300[1]->Fit("Si300_E10", "R");
  h_EoP_hitTrack_300[2]->Fit("Si300_E30", "R");
  h_EoP_hitTrack_300[3]->Fit("Si300_E60", "R");
  h_EoP_hitTrack_300[4]->Fit("Si300_E100", "R");
  h_EoP_hitTrack_300[5]->Fit("Si300_E200", "R");


  float energyValues[6] = {5., 10., 30., 60., 100., 200.};

  //  std::string folder = "plots"+std::string(Form("_E%d",ene));
  std::string folder = "plots_Efit";
  folder = folder + "_" + partType;

   TCanvas* ch_100 = new TCanvas();
  ch_100->cd();
  h_EoP_hitTrack_100[0]->GetXaxis()->SetTitle("#Sigma E_{i} / E_{gen} in 100#mum");
  h_EoP_hitTrack_100[0]->GetXaxis()->SetRangeUser(0.6, 1.4);
  h_EoP_hitTrack_100[0]->GetYaxis()->SetRangeUser(0., 1500);
  if(partType == "Pion"){
    h_EoP_hitTrack_100[0]->GetXaxis()->SetRangeUser(0., 1.4);
    h_EoP_hitTrack_100[0]->GetYaxis()->SetRangeUser(0., 600);
  }
  h_EoP_hitTrack_100[0]->Draw();
  h_EoP_hitTrack_100[1]->Draw("same");
  h_EoP_hitTrack_100[2]->Draw("same");
  h_EoP_hitTrack_100[3]->Draw("same");
  h_EoP_hitTrack_100[4]->Draw("same");
  h_EoP_hitTrack_100[5]->Draw("same");


  t10.DrawLatex(0.20,0.9,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_100[0]->GetParameter(1), fitF_100[0]->GetParError(1), 
   			     fitF_100[0]->GetParameter(2), fitF_100[0]->GetParError(2)));
  t30.DrawLatex(0.20,0.85,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_100[1]->GetParameter(1), fitF_100[1]->GetParError(1), 
   			     fitF_100[1]->GetParameter(2), fitF_100[1]->GetParError(2)));
  t3a.DrawLatex(0.20,0.8,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_100[2]->GetParameter(1), fitF_100[2]->GetParError(1), 
   			     fitF_100[2]->GetParameter(2), fitF_100[2]->GetParError(2)));
  t1a.DrawLatex(0.20,0.75,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_100[3]->GetParameter(1), fitF_100[3]->GetParError(1),
   			     fitF_100[3]->GetParameter(2), fitF_100[3]->GetParError(2)));
  t3p.DrawLatex(0.20,0.7,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_100[4]->GetParameter(1), fitF_100[4]->GetParError(1),
			     fitF_100[4]->GetParameter(2), fitF_100[4]->GetParError(2)) );
  t1p.DrawLatex(0.20,0.65,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_100[5]->GetParameter(1), fitF_100[5]->GetParError(1),
			     fitF_100[5]->GetParameter(2), fitF_100[5]->GetParError(2)) );
  //  t1p.DrawLatex(0.2, 0.85, "#sigma");

  for(int iT=0; iT<6; ++iT){
    tg[0]->SetPoint(iT+1, energyValues[iT], fitF_100[iT]->GetParameter(1));
    tg[0]->SetPointError(iT+1, 0., fitF_100[iT]->GetParError(1));

    //    tgR[0]->SetPoint(iT+1, energyValues[iT], fitF_100[iT]->GetParameter(2));
    tgR[0]->SetPoint(iT+1, 1./energyValues[iT], 100.*pow(fitF_100[iT]->GetParameter(2),2));
    tgR[0]->SetPointError(iT+1, 0., 100.*pow(fitF_100[iT]->GetParError(2),2.));

    tgRstd[0]->SetPoint(iT+1, energyValues[iT], fitF_100[iT]->GetParameter(2)/fitF_100[iT]->GetParameter(1));
    //tgRstd[0]->SetPointError(iT+1, 0., fitF_100[iT]->GetParError(2));
    tgRstd[0]->SetPointError(iT+1, 0., sqrt( pow(fitF_100[iT]->GetParError(2)*fitF_100[iT]->GetParameter(1), 2) + 
					     pow(fitF_100[iT]->GetParError(1)*fitF_100[iT]->GetParameter(2), 2)) / pow(fitF_100[iT]->GetParameter(1), 2) );
  }

  ////////////
  // TLatex latexLabel;
  // latexLabel.SetNDC();
  // latexLabel.SetTextSize(0.04);
  // //  latexLabel.DrawLatex(0.18, 0.96, "#font[132]{CMS Preliminary 2010}");
  // // latexLabel.DrawLatex(0.77, 0.955, "#font[132]{#sqrt{s} = 7 TeV}");
  // //  latexLabel.DrawLatex(0.2, 0.85, "#font[132]{#sigma}");


  leg->Draw("same");
  ch_100->Print((folder+"/h_EoP_100_calibratedRH.png").c_str(), "png");
  ch_100->Print((folder+"/h_EoP_100_calibratedRH.pdf").c_str(), "pdf");  
  ch_100->Print((folder+"/h_EoP_100_calibratedRH.root").c_str(), "root");




  TCanvas* ch_200 = new TCanvas();
  ch_200->cd();
  h_EoP_hitTrack_200[0]->GetXaxis()->SetTitle("#Sigma E_{i} / E_{gen} in 200#mum");
  h_EoP_hitTrack_200[0]->GetXaxis()->SetRangeUser(0.6, 1.4);
  //  h_EoP_hitTrack_200[0]->GetYaxis()->SetRangeUser(0., 1600);
  h_EoP_hitTrack_200[0]->GetYaxis()->SetRangeUser(0., 1000);
  if(partType == "Pion"){
    h_EoP_hitTrack_200[0]->GetXaxis()->SetRangeUser(0., 1.4);
    h_EoP_hitTrack_200[0]->GetYaxis()->SetRangeUser(0., 1500);
  }
  h_EoP_hitTrack_200[0]->Draw();
  h_EoP_hitTrack_200[1]->Draw("same");
  h_EoP_hitTrack_200[2]->Draw("same");
  h_EoP_hitTrack_200[3]->Draw("same");
  h_EoP_hitTrack_200[4]->Draw("same");
  h_EoP_hitTrack_200[5]->Draw("same");
  leg->Draw("same");

  t10.DrawLatex(0.20,0.9,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_200[0]->GetParameter(1), fitF_200[0]->GetParError(1), 
			     fitF_200[0]->GetParameter(2), fitF_200[0]->GetParError(2)));
  t30.DrawLatex(0.20,0.85,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_200[1]->GetParameter(1), fitF_200[1]->GetParError(1), 
			     fitF_200[1]->GetParameter(2), fitF_200[1]->GetParError(2)));
  t3a.DrawLatex(0.20,0.8,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_200[2]->GetParameter(1), fitF_200[2]->GetParError(1), 
			     fitF_200[2]->GetParameter(2), fitF_200[2]->GetParError(2)));
  t1a.DrawLatex(0.20,0.75,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_200[3]->GetParameter(1), fitF_200[3]->GetParError(1),
			     fitF_200[3]->GetParameter(2), fitF_200[3]->GetParError(2)));
  t3p.DrawLatex(0.20,0.7,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_200[4]->GetParameter(1), fitF_200[4]->GetParError(1),
			     fitF_200[4]->GetParameter(2), fitF_200[4]->GetParError(2)));
  t1p.DrawLatex(0.20,0.65,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_200[5]->GetParameter(1), fitF_200[5]->GetParError(1),
			     fitF_200[5]->GetParameter(2), fitF_200[5]->GetParError(2)));

  for(int iT=0; iT<6; ++iT){
    tg[1]->SetPoint(iT+1, energyValues[iT], fitF_200[iT]->GetParameter(1));
    tg[1]->SetPointError(iT+1, 0., fitF_200[iT]->GetParError(1));

    //    tgR[1]->SetPoint(iT+1, energyValues[iT], fitF_200[iT]->GetParameter(2));
    tgR[1]->SetPoint(iT+1, 1./energyValues[iT], 100.*pow(fitF_200[iT]->GetParameter(2),2));
    tgR[1]->SetPointError(iT+1, 0., 100.*pow(fitF_200[iT]->GetParError(2),2));

    tgRstd[1]->SetPoint(iT+1, energyValues[iT], fitF_200[iT]->GetParameter(2)/fitF_200[iT]->GetParameter(1));
    //tgRstd[1]->SetPointError(iT+1, 0., fitF_200[iT]->GetParError(2));
    tgRstd[1]->SetPointError(iT+1, 0., sqrt( pow(fitF_200[iT]->GetParError(2)*fitF_200[iT]->GetParameter(1), 2) + 
					     pow(fitF_200[iT]->GetParError(1)*fitF_200[iT]->GetParameter(2), 2)) / pow(fitF_200[iT]->GetParameter(1), 2) );
			     
  }

  ch_200->Print((folder+"/h_EoP_200_calibratedRH.png").c_str(), "png");
  ch_200->Print((folder+"/h_EoP_200_calibratedRH.pdf").c_str(), "pdf");
  ch_200->Print((folder+"/h_EoP_200_calibratedRH.root").c_str(), "root");




  TCanvas* ch_300 = new TCanvas();
  ch_300->cd();
  h_EoP_hitTrack_300[0]->GetXaxis()->SetTitle("#Sigma E_{i} / E_{gen} in 300#mum");
  h_EoP_hitTrack_300[0]->GetXaxis()->SetRangeUser(0.6, 1.4);
  //  h_EoP_hitTrack_300[0]->GetYaxis()->SetRangeUser(0., 400);
  h_EoP_hitTrack_300[0]->GetYaxis()->SetRangeUser(0., 600);
  if(partType == "Pion"){
    h_EoP_hitTrack_300[0]->GetXaxis()->SetRangeUser(0., 1.4);
    h_EoP_hitTrack_300[0]->GetYaxis()->SetRangeUser(0., 900);
  }
  h_EoP_hitTrack_300[0]->Draw();
  h_EoP_hitTrack_300[1]->Draw("same");
  h_EoP_hitTrack_300[2]->Draw("same");
  h_EoP_hitTrack_300[3]->Draw("same");
  h_EoP_hitTrack_300[4]->Draw("same");
  h_EoP_hitTrack_300[5]->Draw("same");
  leg->Draw("same");

  t10.DrawLatex(0.20,0.9,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_300[0]->GetParameter(1), fitF_300[0]->GetParError(1), 
			     fitF_300[0]->GetParameter(2), fitF_300[0]->GetParError(2)));
  t30.DrawLatex(0.20,0.85,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_300[1]->GetParameter(1), fitF_300[1]->GetParError(1), 
			     fitF_300[1]->GetParameter(2), fitF_300[1]->GetParError(2)));
  t3a.DrawLatex(0.20,0.8,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_300[2]->GetParameter(1), fitF_300[2]->GetParError(1), 
			     fitF_300[2]->GetParameter(2), fitF_300[2]->GetParError(2)));
  t1a.DrawLatex(0.20,0.75,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_300[3]->GetParameter(1), fitF_300[3]->GetParError(1),
			     fitF_300[3]->GetParameter(2), fitF_300[3]->GetParError(2)));
  t3p.DrawLatex(0.20,0.7,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_300[4]->GetParameter(1), fitF_300[4]->GetParError(1),
			     fitF_300[4]->GetParameter(2), fitF_300[4]->GetParError(2)));
  t1p.DrawLatex(0.20,0.65,Form("m = %.2e +/- %.2e   #sigma = %.2e +/- %.2e", fitF_300[5]->GetParameter(1), fitF_300[5]->GetParError(1),
			     fitF_300[5]->GetParameter(2), fitF_300[5]->GetParError(2)));


  for(int iT=0; iT<6; ++iT){
    tg[2]->SetPoint(iT+1, energyValues[iT], fitF_300[iT]->GetParameter(1));
    tg[2]->SetPointError(iT+1, 0., fitF_300[iT]->GetParError(1));

//    tgR[2]->SetPoint(iT+1, 1./energyValues[iT], fitF_300[iT]->GetParameter(2));
    tgR[2]->SetPoint(iT+1, 1./energyValues[iT], 100.*pow(fitF_300[iT]->GetParameter(2),2));
    tgR[2]->SetPointError(iT+1, 0., 100.*pow(fitF_300[iT]->GetParError(2),2));

    tgRstd[2]->SetPoint(iT+1, energyValues[iT], fitF_300[iT]->GetParameter(2)/fitF_300[iT]->GetParameter(1));
    tgRstd[2]->SetPointError(iT+1, 0., sqrt( pow(fitF_300[iT]->GetParError(2)*fitF_300[iT]->GetParameter(1), 2) + 
					     pow(fitF_300[iT]->GetParError(1)*fitF_300[iT]->GetParameter(2), 2)) / pow(fitF_300[iT]->GetParameter(1), 2) );
  }

  //    ch_300->Print((folder+"/h_EoP_hitTrack_300.png").c_str(), "png");
  ch_300->Print((folder+"/h_EoP_300_calibratedRH.png").c_str(), "png");
  ch_300->Print((folder+"/h_EoP_300_calibratedRH.pdf").c_str(), "pdf");
  ch_300->Print((folder+"/h_EoP_300_calibratedRH.root").c_str(), "root");


    /////////TGraph
    TCanvas* tgM = new TCanvas();
    tgM->cd();
    tg[0]->GetYaxis()->SetTitle("#Sigma E_{i} / E");
    tg[0]->GetXaxis()->SetTitle("E (GeV)");
    tg[0]->GetXaxis()->SetRangeUser(0, 300.);
    tg[0]->GetYaxis()->SetRangeUser(0.9, 1.3);
    if(partType == "Pion") tg[0]->GetYaxis()->SetRangeUser(0.7, 1.1);
    tg[0]->Draw("ap");
    tg[1]->Draw("p, same");
    tg[2]->Draw("p, same");
    legTGM->Draw("same");
    tgM->Print((folder+"/tgraph_MeanVsEnergy_calibratedRH.png").c_str(), "png");
    tgM->Print((folder+"/tgraph_MeanVsEnergy_calibratedRH.pdf").c_str(), "pdf");
    tgM->Print((folder+"/tgraph_MeanVsEnergy_calibratedRH.root").c_str(), "root");


    /////////TGraph
    TCanvas* tgReso = new TCanvas();
    tgReso->cd();
    tgR[0]->GetYaxis()->SetTitle("#sigma(E) / E)^{2} (%)");
    tgR[0]->GetXaxis()->SetTitle("1/E (GeV)");
    tgR[0]->GetXaxis()->SetRangeUser(0, 0.15);
    tgR[0]->GetYaxis()->SetRangeUser(0., 1.5);
    if(partType == "Pion") tgR[0]->GetYaxis()->SetRangeUser(0., 10.);
    tgR[0]->Draw("ap");
    tgR[1]->Draw("p, same");
    tgR[2]->Draw("p, same");
    legTGM->Draw("same");
    tgReso->Print((folder+"/tgraph_ResolutionVs1overEnergy_calibratedRH.png").c_str(), "png");
    tgReso->Print((folder+"/tgraph_ResolutionVs1overEnergy_calibratedRH.pdf").c_str(), "pdf");
    tgReso->Print((folder+"/tgraph_ResolutionVs1overEnergy_calibratedRH.root").c_str(), "root");


    /////////TGraph
    TCanvas* tgResoStd = new TCanvas();
    tgResoStd->cd();
    tgRstd[0]->GetYaxis()->SetTitle("#sigma(E) / E");
    tgRstd[0]->GetXaxis()->SetTitle("E (GeV)");
    tgRstd[0]->GetXaxis()->SetRangeUser(0, 300.);
    tgRstd[0]->GetYaxis()->SetRangeUser(0., 0.3);
    if(partType == "Pion") tgRstd[0]->GetYaxis()->SetRangeUser(0., 0.8);
    tgRstd[0]->Draw("ap");
    tgRstd[1]->Draw("p, same");
    tgRstd[2]->Draw("p, same");
    legTGM->Draw("same");
    tgResoStd->Print((folder+"/tgraph_ResolutionVsEnergy_calibratedRH.png").c_str(), "png");
    tgResoStd->Print((folder+"/tgraph_ResolutionVsEnergy_calibratedRH.pdf").c_str(), "pdf");
    tgResoStd->Print((folder+"/tgraph_ResolutionVsEnergy_calibratedRH.root").c_str(), "root");

    return;
    /*
    TCanvas* ch_Weights_Abs = new TCanvas();
    gPad->SetLogy();
    ch_Weights_Abs->cd();
    h_absorberEnergy_vsLayer[0]->GetYaxis()->SetTitle("energy weight for absorber");
    h_absorberEnergy_vsLayer[0]->GetXaxis()->SetTitle("layer");
    h_absorberEnergy_vsLayer[0]->GetYaxis()->SetRangeUser(0.01, 30.);
    h_absorberEnergy_vsLayer[0]->Draw("");
    h_absorberEnergy_vsLayer[1]->Draw("same");
    h_absorberEnergy_vsLayer[2]->Draw("same");
    h_absorberEnergy_vsLayer[3]->Draw("same");
    h_absorberEnergy_vsLayer[4]->Draw("same");
    leg->Draw("same");
    ch_Weights_Abs->Print((folder+"/ch_Weights_Abs.png").c_str(), "png");
    

    TCanvas* ch_Weights_Si = new TCanvas();
    gPad->SetLogy();
    ch_Weights_Si->cd();
    h_energy_vsLayer[0]->GetYaxis()->SetTitle("energy weight for Si");
    h_energy_vsLayer[0]->GetXaxis()->SetTitle("layer");
    h_energy_vsLayer[0]->GetYaxis()->SetRangeUser(1e-05, 0.2);
    h_energy_vsLayer[0]->Draw();
    h_energy_vsLayer[1]->Draw("same");
    h_energy_vsLayer[2]->Draw("same");
    h_energy_vsLayer[3]->Draw("same");
    h_energy_vsLayer[4]->Draw("same");
    leg->Draw("same");
    ch_Weights_Si->Print((folder+"/ch_Weights_Si.png").c_str(), "png");


    TCanvas* ch_Weights_All = new TCanvas();
    gPad->SetLogy();
    ch_Weights_All->cd();
    h_energyAll_vsLayer[0]->GetYaxis()->SetTitle("energy weight");
    h_energyAll_vsLayer[0]->GetXaxis()->SetTitle("layer");
    h_energyAll_vsLayer[0]->GetYaxis()->SetRangeUser(0.01, 30.);
    h_energyAll_vsLayer[0]->Draw();
    h_energyAll_vsLayer[1]->Draw("same");
    h_energyAll_vsLayer[2]->Draw("same");
    h_energyAll_vsLayer[3]->Draw("same");
    h_energyAll_vsLayer[4]->Draw("same");
    leg->Draw("same");
    ch_Weights_All->Print((folder+"/ch_Weights_All.png").c_str(), "png");
    */

}
