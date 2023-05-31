#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TF1.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
#include "TMultiGraph.h"
#include "TMath.h"
#include <sys/stat.h>
#include <vector>
#include <string>
#include "settings.h"
#include "decorations.h"
using namespace std;

//This is the bulk of the code.
void plot(TString path, TString jobname, TString plotname, TString plottitle, int binw, float xmin, float xmax);

//#######################################################################################
//This is the main function:
void stackmaker(){

  TString path = "/mnt/d/work/GitHub/StackOverlayMaker/inputs/";
  //TString jobname = "Skimmed_BasicSelection_May30";
  //TString jobname = "Skimmed_IsoMu0_May30";
  TString jobname = "Skimmed_IsoMuons_mu030_May31";
 
  struct varlist{ TString name; TString title; int rebin; float xmin; float xmax;};
  vector<varlist> variables = {
    {.name="SS_mu0_Pt",       .title="mu0 pT",       .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_mu0_Eta",      .title="mu0 Eta",      .rebin = 10, .xmin=-4, .xmax=4},
    {.name="SS_mu0_Phi",      .title="mu0 Phi",      .rebin = 10, .xmin=-4, .xmax=4},
    {.name="SS_mu0_mT",       .title="mu0 mT",       .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_mu0_reliso03", .title="mu0 reliso03", .rebin = 1, .xmin= 0, .xmax=1},
    {.name="SS_mu0_reliso04", .title="mu0 reliso04", .rebin = 1, .xmin= 0, .xmax=1},
    {.name="SS_mu0_sip3d",    .title="mu0 sip3d",    .rebin = 5, .xmin= 0, .xmax=5},
    {.name="SS_mu1_Pt",       .title="mu1 pT",       .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_mu1_Eta",      .title="mu1 Eta",      .rebin = 10, .xmin=-4, .xmax=4},
    {.name="SS_mu1_Phi",      .title="mu1 Phi",      .rebin = 10, .xmin=-4, .xmax=4},
    {.name="SS_mu1_mT",       .title="mu1 mT",       .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_mu1_reliso03", .title="mu1 reliso03", .rebin = 1, .xmin= 0, .xmax=1},
    {.name="SS_mu1_reliso04", .title="mu1 reliso04", .rebin = 1, .xmin= 0, .xmax=1},
    {.name="SS_mu1_sip3d",    .title="mu1 sip3d",    .rebin = 5, .xmin= 0, .xmax=5},
    {.name="SS_dimuon_mass",  .title="dimuon mass (SS pair)",  .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_dEta_muons",   .title="dEta (SS pair)",.rebin = 20, .xmin= 0, .xmax=6},
    {.name="SS_dPhi_muons",   .title="dPhi (SS pair)",.rebin = 20, .xmin= 0, .xmax=6},
    {.name="SS_dR_muons",     .title="dR (SS pair)",  .rebin = 20, .xmin= 0, .xmax=6},
    {.name="SS_ptratio",      .title="pT ratio (SS pair)",.rebin = 10, .xmin= 0, .xmax=1},
    {.name="SS_met",          .title="MET",           .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_metphi",       .title="MET phi",       .rebin = 10, .xmin=-4, .xmax=4},
    {.name="SS_LT",           .title="Sum(mu pT)",    .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_HT",           .title="Sum(Jet pT)",   .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_nJet",         .title="nJet",          .rebin =  1, .xmin= 0, .xmax=10},
    {.name="SS_nbJet",        .title="nbJet",         .rebin =  1, .xmin= 0, .xmax=10}/*,
    {.name="SS_cutflow_obj",  .title="CutFlow (object level)",  .rebin = 1, .xmin= 0, .xmax=15},
    {.name="SS_cutflow_dimuon",.title="CutFlow (dimuon system)",  .rebin = 1, .xmin= 0, .xmax=15},
    {.name="SS_cutflow_evt",  .title="CutFlow (event level)",  .rebin = 1, .xmin= 0, .xmax=15},
    {.name="SS_cutflow_combined",  .title="CutFlow (combined)",  .rebin = 1, .xmin= 0, .xmax=15},*/
  };
  
  for(int i=0; i<(int)variables.size(); i++){
    plot(path,
	 jobname,
	 variables[i].name,
	 variables[i].title,
	 variables[i].rebin,
	 variables[i].xmin,
	 variables[i].xmax);
    cout<<"Hist no."<<i+1<<" ("<<variables[i].name<<") plotted succesfully."<<endl;
  }

  cout<<"\nSucess!!\n"<<endl;
}

//#######################################################################################

void plot(TString path, TString jobname, TString plotname, TString plottitle, int binw, float xmin, float xmax){
  
  //Global settings: 
  bool toStack      = true;
  bool toSave       = true;
  bool toLog        = false;
  bool toSetRange   = true;
  bool toOverlaySig = true;

  float sigscale = 100;
  
  //###################################
  //Reading histograms and sorting them
  //###################################
  
  //Reading histograms from each file, lumiscaling them  and storing them in a vector<struct>
  vector<hst> QCD, WJets, DY, ST, TTbar, WW, WZ, ZZ, Data;
  QCD   = read_files(hst_qcd,       path, jobname, plotname);
  WJets = read_files(hst_wjets,     path, jobname, plotname);
  DY    = read_files(hst_dy,        path, jobname, plotname);
  ST    = read_files(hst_singletop, path, jobname, plotname);
  TTbar = read_files(hst_ttbar,     path, jobname, plotname);
  WW    = read_files(hst_ww,        path, jobname, plotname);
  WZ    = read_files(hst_wz,        path, jobname, plotname);
  ZZ    = read_files(hst_zz,        path, jobname, plotname);
  Data  = read_files(hst_data,      path, jobname, plotname);

  //Push back all the backgrounds into a vector:
  vector<merged> bkg;
  bkg.push_back(merge_and_decorate(DY,    "DY",    plotname, kRed-9));
  bkg.push_back(merge_and_decorate(WJets, "WJets", plotname, kGray+1));
  bkg.push_back(merge_and_decorate(ST,    "ST",    plotname, kCyan-9));
  bkg.push_back(merge_and_decorate(TTbar, "TTBar", plotname, kCyan-10));
  bkg.push_back(merge_and_decorate(WW,    "WW",    plotname, kGreen-4));
  bkg.push_back(merge_and_decorate(WZ,    "WZ",    plotname, kBlue-9));
  bkg.push_back(merge_and_decorate(ZZ,    "ZZ",    plotname, kBlue-7));
  bkg.push_back(merge_and_decorate(QCD,   "QCD",   plotname, kYellow));
  
  //#################################################################
  //Decorating the histograms and setting the bins:
  //Decorating the backgrounds:
  for(int i=0; i<(int)bkg.size(); i++){
    SetHistoStyle(bkg[i].hist, plotname, bkg[i].color);
    bkg[i].hist->Rebin(binw);
    SetOverflowBin(bkg[i].hist);
  }

  //Reading and decorating data:
  merged data = merge_and_decorate(Data, "Data", plotname, kBlack);
  SetHistoStyle(data.hist, plotname, kBlack);
  data.hist->Rebin(binw);
  SetOverflowBin(data.hist);

  //Reading and decorating signal:
  //Only one signal file for now.
  TString sigpath="inputs/"+jobname+"/signal/";
  vector<TString> sigfile = {"hst_VLL100.root"};
  vector<float> lumi = {508761.54}; // 595251 events /1.17 pb
  //Make sure that the sizes of sigfile and lumi are the same.
  vector<TH1F*> sighist;
  for(int i=0; i<sigfile.size(); i++){
    TString fullpath = sigpath+sigfile[i];
    TFile *t = new TFile(fullpath);
    TH1F *h = (TH1F *)t->Get(plotname);
    h->Scale(59700/lumi[i]);
    h->Scale(sigscale);
    SetHistoStyle(h, plotname, kRed);
    h->Rebin(binw);
    SetOverflowBin(h);
    sighist.push_back(h);
  }

  

  //#################################################################
  //Preparing the stack:
  //Sorting the backgrounds according to their integrals (acsending order):
  for(int i=0; i<(int)bkg.size()-1; i++){
    for(int j=i+1; j<(int)bkg.size(); j++){
      float inti = bkg[i].hist->Integral();
      float intj = bkg[j].hist->Integral();
      if(inti > intj) swap(bkg[i], bkg[j]);
    }
  }
  //Making a stack of the sorted backgrounds:
  //cout<<"\nPraparing the stack with "<<(int)bkg.size()<<" no. of backgrounds ...";
  //Preparing fill colors:
  for(int i=0; i<(int)bkg.size(); i++){
    if(toStack) bkg[i].hist->SetFillColor(bkg[i].color);
    else bkg[i].hist->Scale(1/bkg[i].hist->Integral());
  }
  //Filling the stack with histograms:
  THStack *stack = new THStack("Stacked",plotname+";"+plotname+";Events");
  for(int i=0; i<(int)bkg.size(); i++){stack->Add(bkg[i].hist);}
  //cout<<"....Done!"<<endl;

  //#######################################################################
  // Making the ratio hist. Comment this out, if not needed.
  TH1F *hst_ratio = (TH1F *)data.hist->Clone("copy");
  TH1F *total_bkg = (TH1F *)data.hist->Clone("copy"); total_bkg->Reset();
  for(int i=0; i<(int)bkg.size(); i++) total_bkg->Add(bkg[i].hist);
  hst_ratio->Divide(total_bkg);
  float obs = data.hist->Integral();
  float exp = total_bkg->Integral();
  float obsbyexp = obs/exp;
  //cout<<"\nObs(data)/exp(MC) = "<<obsbyexp<<endl;
  
  //######################
  //Setting up the canvas:
  //######################
  TCanvas *c1 = new TCanvas(plotname,plotname,700,600);
  
  TPad *mainPad = create_mainpad();
  mainPad->SetMargin(mainPad->GetLeftMargin(), mainPad->GetRightMargin(), 0.15, 0.1);
  if(toLog) mainPad->SetLogy(1);
  mainPad->Draw();
  
  TPad *ratioPad = create_ratiopad();
  ratioPad->SetMargin(ratioPad->GetLeftMargin(), ratioPad->GetRightMargin(), 0.2, 0.1);
  ratioPad->Draw();

  //Legend:
  TLegend *lg1 = create_legend();
  TString legendheader = "obs/exp = "+to_string(obsbyexp);
  lg1->SetHeader(legendheader);
  TString datacount = "2018 data ["+to_string((int)data.hist->Integral())+"]";
  lg1->AddEntry(data.hist, datacount, "f");
  for(int i = (int)bkg.size()-1; i>-1;  i--){
    TString count = to_string((int)bkg[i].hist->Integral());
    TString name = bkg[i].name + " ["+count+"]";
    lg1->AddEntry(bkg[i].hist, name, "f");
  }
  if(toOverlaySig){
    TString sigcount100 = "VLL100 ["+to_string((int)(sighist[0]->Integral()/sigscale))+"]";
    lg1->AddEntry(sighist[0], sigcount100, "f");
  }
  
  //####################
  //Drawing on mainpad
  //###################
  PutText(jobname);

  mainPad->cd();
  mainPad->SetFillStyle(4000);

  //Normalising the histogrms, if not to stack:
  if(!toStack){
    data.hist->Scale(1/data.hist->Integral());
    if(toOverlaySig) for(int i=0; i<(int)sighist.size(); i++) sighist[i]->Scale(1/sighist[i]->Integral());
    for(int i=0; i<(int)bkg.size(); i++) bkg[i].hist->Scale(1/bkg[i].hist->Integral());
  }

  //####################
  //Finding the 'peaks':
  //####################
  int nbins = data.hist->GetNbinsX();
  float data_peak = 0; float stack_peak=0;
  for(int bin=0; bin<nbins; bin++){
    float ydata = 0;
    float ystack = 0;
    
    ydata = data.hist->GetBinContent(bin);
    for(int i=0; i<(int)bkg.size(); i++){ystack = ystack + bkg[i].hist->GetBinContent(bin);}

    if(ydata  > data_peak)  data_peak  = ydata;
    if(ystack > stack_peak) stack_peak = ystack;
  }
  //Finding the tallest background index:
  float bkg_peak =0; int tallest_bkg = -1;
  for(int i=0; i<(int)bkg.size(); i++){
    for(int bin=0; bin<nbins; bin++){
      int yvalue = bkg[i].hist->GetBinContent(bin);
      if(yvalue > bkg_peak){
	bkg_peak = yvalue;
	tallest_bkg = i;
      }	
    }
  }
  //Finding the tallest signal index:
  float sig_peak =0; int tallest_sig = -1;
  for(int i=0; i<(int)sighist.size(); i++){
    for(int bin=0; bin<nbins; bin++){
      int yvalue = sighist[i]->GetBinContent(bin);
      if(yvalue > sig_peak){
	sig_peak = yvalue;
	tallest_sig = i;
      }	
    }
  }

  //###################
  //Drawing with stack:
  //###################
  
  if(toStack){
    //If the stack is tallest:
    if(stack_peak > data_peak){
      stack->Draw("hist");
      stack->GetYaxis()->SetTitle("Events/sildbv");
      stack->GetXaxis()->SetTitle(plottitle);
      stack->SetTitle("");
      if(toSetRange) stack->GetXaxis()->SetRangeUser(xmin, xmax);
    }
    //If data is tallest:
    else{
      data.hist->GetXaxis()->SetTitle(plottitle);
      data.hist->GetYaxis()->SetTitle("Events");
      data.hist->SetStats(0);
      if(toSetRange) data.hist->GetXaxis()->SetRangeUser(xmin, xmax);
      //data.hist->GetYaxis()->SetRangeUser(1, 10E7);
      data.hist->Draw("ep");
    } 
    stack->Draw("hist same");
    data.hist->Draw("ep same");
    if(toOverlaySig) sighist[0]->Draw("hist same");
  }

  /*
  data.hist->GetXaxis()->SetTitle(plottitle);
  data.hist->GetYaxis()->SetTitle("Events");
  data.hist->SetStats(0);
  if(toSetRange) data.hist->GetXaxis()->SetRangeUser(xmin, xmax);
  //data.hist->GetYaxis()->SetRangeUser(1, 10E7);
  data.hist->Draw("ep");*/
  stack->Draw("hist same");
  if(toOverlaySig) sighist[0]->Draw("hist same");
  data.hist->Draw("ep same");
  
  //################
  // Overlay only:
  //################
  if(!toStack){
    //if data is tallest:
    if(data_peak > bkg_peak && data_peak > sig_peak){
      data.hist->GetXaxis()->SetTitle(plottitle);
      data.hist->GetYaxis()->SetTitle("Events");
      data.hist->SetStats(0);
      if(toSetRange) data.hist->GetXaxis()->SetRangeUser(xmin, xmax);
      data.hist->Draw("hist");
    }
    //if one of the backgroud is tallest:
    else if(bkg_peak > data_peak && bkg_peak > sig_peak){
      bkg[tallest_bkg].hist->GetXaxis()->SetTitle(plottitle);
      bkg[tallest_bkg].hist->GetYaxis()->SetTitle("Events");
      bkg[tallest_bkg].hist->SetStats(0);
      if(toSetRange) bkg[tallest_bkg].hist->GetXaxis()->SetRangeUser(xmin, xmax);
      bkg[tallest_bkg].hist->Draw("hist");
    }
    //if signal is tallest:
    else if(sig_peak > data_peak && sig_peak > bkg_peak && toOverlaySig){
      sighist[tallest_sig]->GetXaxis()->SetTitle(plottitle);
      sighist[tallest_sig]->GetYaxis()->SetTitle("Events");
      sighist[tallest_sig]->SetStats(0);
      if(toSetRange) sighist[tallest_sig]->GetXaxis()->SetRangeUser(xmin, xmax);
      sighist[tallest_sig]->Draw("hist");
    }

    //Drawing the rest:
    data.hist->Draw("hist same");
    if(toOverlaySig) sighist[0]->Draw("hist same");
    for(int i=0; i<(int)bkg.size(); i++) bkg[i].hist->Draw("hist_same");

  }

  lg1->Draw("same");

  //#########################
  // Drawing on the ratiopad
  //######################### 
  ratioPad->cd();

  SetHistoStyle(hst_ratio, plotname, kBlack);
  hst_ratio->GetYaxis()->SetTitle("obs/exp");
  hst_ratio->GetYaxis()->SetRangeUser(0, 2);
  hst_ratio->GetYaxis()->SetNdivisions(5, kTRUE);
  if(toSetRange) hst_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
  hst_ratio->SetStats(0);
  hst_ratio->Draw("ep");

  TLine line;
  line.SetLineColor(kRed); line.SetLineWidth(2); line.SetLineStyle(1);
  float min_x = data.hist->GetXaxis()->GetXmin();
  float max_x = data.hist->GetXaxis()->GetXmax();
  if(toSetRange) line.DrawLine(xmin,1,xmax,1);
  else line.DrawLine(min_x,1,max_x,1);
  
  //######################################################################
  TString date_stamp = todays_date();
  TString dump_folder = "outputs/"+date_stamp+"_"+jobname;
  createFolder(dump_folder);

  if(toSave){
    c1->SaveAs(dump_folder+"/"+plotname+".png");
  }
}
