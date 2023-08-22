// This script overlays all the points of a 2d histogram on a plane.
// I am trying to explore 2d phase spaces to determine SR and CR.

#include <iostream>
#include <vector>
#include <string>
#include <TSystem.h>
#include "TF2.h"
#include "settings_bkgonly.h"
using namespace std;

bool fileExists(const char *filePath){
  TSystem *sys = gSystem;
  if (!sys){
    std::cerr << "Error: ROOT system not available." << std::endl;
    return false;
  }
  if (sys->AccessPathName(filePath) == 0) return true;
  else return false;
}

void explore2d(){

  TString plotname = "region_ht_deta";
  TString jobname = "hst_Aug10_WR";
  TString plottitle = "HT-dEta plane";
  TString tag = "";
  float xlim[2] = {0, 1000};
  float ylim[2] = {0, 6};
  int rebin[2] = {10, 10};

  //Switches:
  bool toDebug = true; //Prints out additional statements.
  bool toSave = true;

  //Don't make changes below:
  //----------------------------------------------------------------------------
  TString inputfolder = "inputs/" + jobname;

  struct hst{
    TString group;
    TString subgr;
    float lumi;
    TH2F *hist;
  };

  vector<hst> bkg = {
      {.group = "QCD_MuEnriched", .subgr = "20to30", .lumi = 23.893},
      {.group = "QCD_MuEnriched", .subgr = "30to50", .lumi = 42.906},
      {.group = "QCD_MuEnriched", .subgr = "50to80", .lumi = 105.880},
      {.group = "QCD_MuEnriched", .subgr = "80to120", .lumi = 508.715},
      {.group = "QCD_MuEnriched", .subgr = "120to170", .lumi = 1802.854},
      {.group = "QCD_MuEnriched", .subgr = "170to300", .lumi = 10265.815},
      {.group = "QCD_MuEnriched", .subgr = "300to470", .lumi = 95249.242},
      {.group = "QCD_MuEnriched", .subgr = "470to600", .lumi = 656872.156},
      {.group = "QCD_MuEnriched", .subgr = "600to800", .lumi = 2060827.812},
      {.group = "QCD_MuEnriched", .subgr = "800to1000", .lumi = 11337379.457},
      {.group = "HTbinnedWJets", .subgr = "70to100", .lumi = 52100.910},
      {.group = "HTbinnedWJets", .subgr = "100to200", .lumi = 41127.174},
      {.group = "HTbinnedWJets", .subgr = "200to400", .lumi = 172265.183},
      {.group = "HTbinnedWJets", .subgr = "400to600", .lumi = 163821.083},
      {.group = "HTbinnedWJets", .subgr = "600to800", .lumi = 708793.021},
      {.group = "HTbinnedWJets", .subgr = "800to1200", .lumi = 1481985.193},
      {.group = "HTbinnedWJets", .subgr = "1200to2500", .lumi = 5602003.457},
      {.group = "HTbinnedWJets", .subgr = "2500toInf", .lumi = 79396214.989},
      {.group = "DYJetsToLL", .subgr = "M10to50", .lumi = 5925.522},
      {.group = "DYJetsToLL", .subgr = "M50", .lumi = 30321.155},
      {.group = "SingleTop", .subgr = "s-channelLeptonDecays", .lumi = 5456748.098},
      {.group = "SingleTop", .subgr = "t-channelAntiTop", .lumi = 1407728.544},
      {.group = "SingleTop", .subgr = "t-channelTop", .lumi = 1572627.866},
      {.group = "SingleTop", .subgr = "tW_AntiTop", .lumi = 238357.428},
      {.group = "SingleTop", .subgr = "tW_Top", .lumi = 245177.196},
      {.group = "TTBar", .subgr = "TTTo2L2Nu", .lumi = 1642541.624},
      {.group = "TTBar", .subgr = "TTToSemiLeptonic", .lumi = 1304012.700},
      {.group = "WW", .subgr = "WWTo1L1Nu2Q", .lumi = 805257.731},
      {.group = "WW", .subgr = "WWTo4Q", .lumi = 773049.853},
      {.group = "WZ", .subgr = "WZTo1L1Nu2Q", .lumi = 805257.731},
      {.group = "WZ", .subgr = "WZTo2Q2L", .lumi = 4499605.731},
      {.group = "WZ", .subgr = "WZTo3LNu", .lumi = 1889798.538},
      {.group = "ZZ", .subgr = "ZZTo2L2Nu", .lumi = 58416512.631},
      {.group = "ZZ", .subgr = "ZZTo2Q2L", .lumi = 7928149.608},
      {.group = "ZZ", .subgr = "ZZTo2Q2Nu", .lumi = 4405016.452},
      {.group = "ZZ", .subgr = "ZZTo4L", .lumi = 74330566.038},
      {.group = "TTW", .subgr = "TTWToLNu", .lumi = 48627268.5}};
  
  //-----------------------------------------------------------------------------
  //Reading histograms:
  if(toDebug) cout<<"No of files in the list = "<<(int)bkg.size()<<endl;
  for(int i = 0; i < (int)bkg.size(); i++){
    TString filename = "hst_" + bkg[i].group + "_" + bkg[i].subgr + ".root";
    TString fullpath = inputfolder + "/" + bkg[i].group + "/" + filename;
    if(fileExists(fullpath)){

      //Opening the file and accessing the histogram:
      TFile *tfile = new TFile(fullpath);
      TH2F *hist = (TH2F *)tfile->Get(plotname);
      float scalefactor = 59700/bkg[i].lumi;
      if(bkg[i].group == "QCD_MuEnriched") scalefactor = scalefactor*0.024164;
      hist->Scale(scalefactor);
      hist->Rebin2D(rebin[0], rebin[1]);

      bkg[i].hist = hist;
    }
    else bkg.erase(bkg.begin() + i);
  }
  if(toDebug) cout<<"No of files read = "<<(int)bkg.size()<<endl;

  //Signal: 
  TFile *sigfile = new TFile("inputs/"+jobname+"/signal/hst_VLL100.root");
  TH2F *sighist = (TH2F *)sigfile->Get(plotname);
  float sigscale = 59700/508761.54;
  //sighist->Scale(sigscale);
  sighist->Rebin2D(rebin[0], rebin[1]);

  cout<<"File reading successful!"<<endl;
  //-----------------------------------------------------------------------------

  //Plotting:
  TCanvas *c1 = new TCanvas(plotname,plotname,800,600);
  //gStyle->SetPalette(kRainBow); // Set a color palette
  c1->SetRightMargin(0.15);

  //Add all the backgrounds together:
  TList *list = new TList;
  for(int i=0; i<(int)bkg.size(); i++){
    if(bkg[i].hist) list->Add(bkg[i].hist);
    else cout<<"Warning : null hist found for "+bkg[i].group+"_"+bkg[i].subgr<<endl; 
  }
  TH2F *totalbkg = (TH2F *)bkg[0].hist->Clone(0);
  totalbkg->Reset();
  totalbkg->Merge(list);

  totalbkg->SetStats(0);
  totalbkg->SetTitle(plottitle);
  totalbkg->GetXaxis()->SetRangeUser(xlim[0], xlim[1]);
  totalbkg->GetYaxis()->SetRangeUser(ylim[0], ylim[1]);

  totalbkg->Draw("COLZ");
  sighist->SetMarkerColor(kRed);
  sighist->SetMarkerStyle(20);
  //sighist->SetMarkerSize(2.0);
  sighist->Draw("SAME P");

  c1->Update();
  
  //-------------------------------------------------------------------
  TString date_stamp = todays_date();
  TString dump_folder = "outputs/"+date_stamp+"_"+jobname+tag+"/2d-plots";
  createFolder(dump_folder);

  if(toSave){
    c1->SaveAs(dump_folder+"/"+plotname+".png");
  }

  cout<<"\nDone!"<<endl;
}
