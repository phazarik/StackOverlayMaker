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
#include <algorithm>
#include <TGraphErrors.h>
#include "settings.h"
#include "decorations.h"
using namespace std;

//This is the bulk of the code.
void plot(TString path, TString jobname, TString plotname, TString plottitle, int binw, float xmin, float xmax);

//#######################################################################################
//This is the main function:
void stackmaker(){

  TString path = "/mnt/d/work/GitHub/StackOverlayMaker/inputs/";
  TString jobname = "hst_Aug28_Basic";
 
  struct varlist{ TString name; TString title; int rebin; float xmin; float xmax;};
  vector<varlist> variables = {
    {.name="SS_STfrac", .title="LT/ST",        .rebin = 5,  .xmin= 0, .xmax=1},   
    //{.name="SS_dilep_mass", .title="dilep mass (SS pair)", .rebin = 10, .xmin= 0, .xmax=200},
    /*
    //Object level plots:
    {.name="SS_llep0_Pt",       .title="llep0 pT",       .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_llep0_Eta",      .title="llep0 Eta",      .rebin = 10, .xmin=-4, .xmax=4},
    {.name="SS_llep0_Phi",      .title="llep0 Phi",      .rebin = 10, .xmin=-4, .xmax=4},
    {.name="SS_llep0_mT",       .title="llep0 mT",       .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_llep0_reliso03", .title="llep0 reliso03", .rebin = 1, .xmin= 0, .xmax=1},
    {.name="SS_llep0_sip3d",    .title="llep0 sip3d",    .rebin = 5, .xmin= 0, .xmax=5},
    {.name="SS_llepss_Pt",       .title="llepss pT",       .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_llepss_Eta",      .title="llepss Eta",      .rebin = 10, .xmin=-4, .xmax=4},
    {.name="SS_llepss_Phi",      .title="llepss Phi",      .rebin = 10, .xmin=-4, .xmax=4},
    {.name="SS_llepss_mT",       .title="llepss mT",       .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_llepss_reliso03", .title="llepss reliso03", .rebin = 1, .xmin= 0, .xmax=1},
    {.name="SS_llepss_sip3d",    .title="llepss sip3d",    .rebin = 5, .xmin= 0, .xmax=5},
    //Dilepton system:
    {.name="SS_dilep_mass", .title="dilep mass (SS pair)", .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_dEta_llep",  .title="dEta (SS pair)",       .rebin = 20, .xmin= 0, .xmax=6},
    {.name="SS_dPhi_llep",  .title="dPhi (SS pair)",       .rebin = 20, .xmin= 0, .xmax=6},
    {.name="SS_dR_llep",    .title="dR (SS pair)",         .rebin = 20, .xmin= 0, .xmax=6},
    {.name="SS_ptratio",    .title="pT ratio (SS pair)",   .rebin = 10, .xmin= 0, .xmax=1},*/
    /*
    //Event level plots:
    {.name="SS_nllep", .title="nllep", 1, 0, 10},
    {.name="SS_nJet",  .title="nJet" , 1, 0, 10},
    {.name="SS_nbJet", .title="nbJet", 1, 0, 10},
    {.name="SS_met",    .title="MET",     .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_metphi", .title="MET phi", .rebin = 10, .xmin=-4, .xmax=4},
    {.name="SS_LT",     .title="Sum(llep pT)", .rebin = 50, .xmin= 0, .xmax=1000},
    {.name="SS_HT",     .title="Sum(Jet pT)",  .rebin = 50, .xmin= 0, .xmax=1000},
    {.name="SS_ST",     .title="LT + HT",      .rebin = 50, .xmin= 0, .xmax=1000},
    {.name="SS_STfrac", .title="LT/ST",        .rebin = 5,  .xmin= 0, .xmax=1},
    {.name="SS_dPhi_met0",  .title="dPhi (llep0, MET)",  .rebin = 20, .xmin= 0, .xmax=6},
    {.name="SS_dPhi_metss", .title="dPhi (llepss, MET)", .rebin = 20, .xmin= 0, .xmax=6},
    {.name="SS_dPhi_met_max", .title="max dPhi (llep, MET)", .rebin = 20, .xmin= 0, .xmax=6},
    {.name="SS_dPhi_met_min", .title="min dPhi (llep, MET)", .rebin = 20, .xmin= 0, .xmax=6},*/
  };
  
  for(int i=0; i<(int)variables.size(); i++){
    plot(path,
	 jobname,
	 variables[i].name,
	 variables[i].title,
	 variables[i].rebin,
	 variables[i].xmin,
	 variables[i].xmax);
    cout<<"Hist no."<<i+1<<" ("<<variables[i].name<<") plotted succesfully.\n"<<endl;
    //break;
  }

  cout<<"Sucess!!\n"<<endl;
}

//#######################################################################################

void plot(TString path, TString jobname, TString plotname, TString plottitle, int binw, float xmin, float xmax){
  
  //Global settings: 
  bool toStack        = true;
  bool toSave         = false;
  bool toLog          = true;
  bool toSetRange     = true;
  bool toOverlaySig   = true;
  bool toOverlayData  = true;
  bool toScaleHT      = false;
  bool toPrintBinInfo = false;
  bool toPlotUncertainty = true;
  bool toAddEgamma    = true;
  float sigscale = 1;
  TString tag = "";
  //TString tag = "";

  if(toSave) gROOT->SetBatch(kTRUE);
  
  //###################################
  //Reading histograms and sorting them
  //###################################
  
  //Reading histograms from each file, lumiscaling them  and storing them in a vector<struct>
  vector<hst> QCD, WJets, DY, ST, TTbar, WW, WZ, ZZ, TTW;
  QCD   = read_files(hst_qcd,       path, jobname, plotname);
  WJets = read_files(hst_wjets,     path, jobname, plotname);
  DY    = read_files(hst_dy,        path, jobname, plotname);
  ST    = read_files(hst_singletop, path, jobname, plotname);
  TTbar = read_files(hst_ttbar,     path, jobname, plotname);
  WW    = read_files(hst_ww,        path, jobname, plotname);
  WZ    = read_files(hst_wz,        path, jobname, plotname);
  ZZ    = read_files(hst_zz,        path, jobname, plotname);
  TTW   = read_files(hst_ttw,       path, jobname, plotname);

  //Push back all the backgrounds into a vector:
  vector<merged> bkg;
  bkg.push_back(merge_and_decorate(QCD,   "QCD",   plotname, kYellow));
  bkg.push_back(merge_and_decorate(DY,    "DY",    plotname, kRed-7));
  bkg.push_back(merge_and_decorate(WJets, "WJets", plotname, kGray+1));
  bkg.push_back(merge_and_decorate(ST,    "ST",    plotname, kCyan-7));
  bkg.push_back(merge_and_decorate(TTbar, "TTBar", plotname, kAzure+1));
  bkg.push_back(merge_and_decorate(TTW,   "TTW",   plotname, kAzure+2));
  bkg.push_back(merge_and_decorate(WW,    "WW",    plotname, kGreen-3));
  bkg.push_back(merge_and_decorate(WZ,    "WZ",    plotname, kGreen-9));
  bkg.push_back(merge_and_decorate(ZZ,    "ZZ",    plotname, kGreen-10));
  

  //Decorating the histograms and setting the bins:
  //Decorating the backgrounds:
  for(int i=0; i<(int)bkg.size(); i++){
    if(toStack) SetHistoStyle(bkg[i].hist, plotname, kBlack);
    else SetHistoStyle_bold(bkg[i].hist, plotname, bkg[i].color);
    bkg[i].hist->Rebin(binw);
    SetOverflowBin(bkg[i].hist);
  }

  //-------------------------------------------------------------------
  //SIGNAL:
  TString sigpath="inputs/"+jobname+"/signal/";

  for(int i=0; i<(int)sig.size(); i++){
    TString fullpath = sigpath+sig[i].filename;
    TFile *t = new TFile(fullpath);
    TH1F *h = (TH1F *)t->Get(plotname);
    h->Scale(59700/sig[i].lumi);
    h->Scale(sigscale);//If futher scaling is required
    SetHistoStyle_bold(h, plotname, sig[i].color);
    h->Rebin(binw);
    SetOverflowBin(h);
    sig[i].hist = h;
  }

  //Data:
  TString datapath = "inputs/"+jobname+"/SingleMuon/";
  vector<TString> datafile = {"hst_SingleMuon_A.root","hst_SingleMuon_B.root","hst_SingleMuon_C.root","hst_SingleMuon_D.root"};
  TH1F * datahist = nullptr;
  if(toOverlayData){
    for(int i=0; i<(int)datafile.size(); i++){
      TString fullpath = datapath+datafile[i];
      TFile *t = new TFile(fullpath);
      TH1F *h = (TH1F *)t->Get(plotname);
      SetHistoStyle(h, plotname, kBlack);
      h->Rebin(binw);
      SetOverflowBin(h);

      if(datahist == nullptr) datahist = (TH1F *)h->Clone(0);
      else datahist->Add(h);      
    }
  }
  //Adding EGamma:
  TString datapath2 = "inputs/"+jobname+"/EGamma/";
  vector<TString> datafile2 = {"hst_EGamma_A.root","hst_EGamma_B.root","hst_EGamma_C.root","hst_EGamma_D.root"};
  TH1F * datahist2 = nullptr;
  if(toOverlayData && toAddEgamma){
    for(int i=0; i<(int)datafile2.size(); i++){
      TString fullpath = datapath2+datafile2[i];
      TFile *t = new TFile(fullpath);
      TH1F *h = (TH1F *)t->Get(plotname);
      SetHistoStyle(h, plotname, kBlack);
      h->Rebin(binw);
      SetOverflowBin(h);

      if(datahist2 == nullptr) datahist2 = (TH1F *)h->Clone(0);
      else datahist2->Add(h);      
    }
    datahist->Add(datahist2);
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

  //##################################################################
  //Finding the drawing order:
  float peak_stk=0; float peak_dat=0; float peak_sig=0;
  //int tallest_bkg_index = -1;
  
  int nbins = bkg[0].hist->GetNbinsX();
  for(int bin=0; bin<nbins; bin++){
    //total count in each bin:
    int ystk = 0; int ysig=0; int ydata=0;
    for(int i=0; i<(int)bkg.size(); i++) ystk = ystk + bkg[i].hist->GetBinContent(bin);
    //if(toOverlaySig)
    //if(toOverlayData)

    if(ystk > peak_stk) peak_stk = ystk;
  }

  //S over root B and obs-by-exp
  float SbyB_global = 0; float obsbyexp_global = 0;
  float nsig=0; float nbkg=0; float ndata=0;
  for(int i=0; i<(int)bkg.size(); i++) nbkg = nbkg + bkg[i].hist->Integral();
  nsig = sig[0].hist->Integral();
  if(toOverlayData) ndata = datahist->Integral();
  
  if(nbkg!=0){
    SbyB_global = nsig/sqrt(nbkg);
    obsbyexp_global = ndata/nbkg;
  }
  else cout<<"Error! nbkg = 0"<<endl;
  
  cout<<"S over root B = "<<SbyB_global<<endl;

  if(toPrintBinInfo){
    cout<<"--------------------------------------------------------------------------"<<endl;
    cout<<"Bin-wise significance:"<<endl;
    int nbins = (int)sig[0].hist->GetNbinsX();
    for(int bin=0; bin<nbins; bin++){
      float nsig = sig[0].hist->GetBinContent(bin);
      float nbkg = 0; for(int i=0; i<(int)bkg.size(); i++) nbkg = nbkg+bkg[i].hist->GetBinContent(bin);
      float sbyb = 0;
      if(nbkg !=0) sbyb = nsig/sqrt(nbkg);
      cout<<"Cut"<<bin-1<<" : nsig = "<<(int)nsig<<"\t nbkg = "<<(int)nbkg<<"\t S/sqrt{B} = "<<sbyb<<endl;
    }
    cout<<"--------------------------------------------------------------------------"<<endl;
  }
  
  //######################
  //Setting up the canvas:
  //######################
  TCanvas *c1 = new TCanvas(plotname,plotname,700,600);
  TPad *mainPad, *ratioPad;
  
  mainPad = create_mainpad();
  mainPad->SetMargin(mainPad->GetLeftMargin(), mainPad->GetRightMargin(), 0.15, 0.1);
  if(toLog) mainPad->SetLogy(1);
  mainPad->Draw();
  
  ratioPad = create_ratiopad();
  ratioPad->SetMargin(ratioPad->GetLeftMargin(), ratioPad->GetRightMargin(), 0.2, 0.1);
  ratioPad->Draw();

  //Legend:
  TLegend *lg1 = create_legend();
  TString legendheader = "";
  if(toOverlayData) legendheader = "obs/exp = "+to_string(obsbyexp_global);
  else legendheader = "S/sqrt(B) = "+to_string(SbyB_global);
  lg1->SetHeader(legendheader);
  if(toOverlayData){
    TString datacount = "2018 data ["+to_string((int)datahist->Integral())+"]";
    lg1->AddEntry(datahist, datacount, "f");
  }
  for(int i = (int)bkg.size()-1; i>-1;  i--){
    TString count = to_string((int)bkg[i].hist->Integral());
    TString name = bkg[i].name + " ["+count+"]";
    lg1->AddEntry(bkg[i].hist, name, "f");
  }
  if(toOverlaySig){
    for(int i=0; i<(int)sig.size(); i++){
      TString sigcount = "VLL"+sig[i].mass+" ["+to_string((int)(sig[i].hist->Integral()/sigscale))+"]";
      lg1->AddEntry(sig[i].hist, sigcount, "f");
    }
  }

  //####################
  //Drawing on mainpad
  //###################
  PutText(jobname);

  mainPad->cd();
  mainPad->SetFillStyle(4000);

  //Normalising the histogrms, if not to stack:
  if(!toStack){
    //data.hist->Scale(1/data.hist->Integral());
    if(toOverlaySig) for(int i=0; i<(int)sig.size(); i++) sig[i].hist->Scale(1/sig[i].hist->Integral());
    for(int i=0; i<(int)bkg.size(); i++) bkg[i].hist->Scale(1/bkg[i].hist->Integral());
  }

  //Drawing the first hist:
  sig[0].hist->GetXaxis()->SetTitle(plottitle);
  sig[0].hist->GetYaxis()->SetTitle("Events");
  sig[0].hist->SetStats(0);
  if(toSetRange) sig[0].hist->GetXaxis()->SetRangeUser(xmin, xmax);
  sig[0].hist->GetYaxis()->SetRangeUser(0.1, 1E8);
  sig[0].hist->Draw("hist");

  //Drawing the rest:
  stack->Draw("hist same");
  if(toOverlayData) datahist->Draw("ep same");
  for(int i =0; i<(int)sig.size(); i++) sig[i].hist->Draw("hist same");
  lg1->Draw("same");


  //#########################
  // Drawing on the ratiopad
  //#########################
  TH1F *rootb, *ratio;
  TGraphErrors *err;
  TH1F *exp = (TH1F*)bkg[0].hist->Clone("copy"); exp->Reset();
  for(int i=0; i<(int)bkg.size(); i++) exp->Add(bkg[i].hist);
  
  //Setting up obsbyexp:
  if(toOverlayData){    
    ratio = (TH1F*)datahist->Clone("copy");
    ratio->Divide(exp);
    //Keep the error in data only:
    // Copy the errors from h1 to the result histogram
    for (int bin = 0; bin <= ratio->GetNbinsX(); bin++) {
      int ndata = datahist->GetBinContent(bin);
      double dataerr = 0;
      if(ndata !=0) dataerr = datahist->GetBinError(bin)/ndata;
      ratio->SetBinError(bin, dataerr);
    }
    SetHistoStyle(ratio, plotname, kBlack);
    ratio->GetYaxis()->SetTitle("obs/exp");
    ratio->GetYaxis()->SetRangeUser(0, 2);
    ratio->GetYaxis()->SetNdivisions(5, kTRUE);
    if(toSetRange) ratio->GetXaxis()->SetRangeUser(xmin, xmax);
    ratio->SetStats(0);
    if(toPlotUncertainty){
      err = get_err(exp);
      err->GetYaxis()->SetNdivisions(5, kTRUE);
      err->SetStats(0);
      if(toSetRange) err->GetXaxis()->SetRangeUser(xmin, xmax);
    }
  }

  //Setting up s-by-rootB:
  rootb = (TH1F *)bkg[0].hist->Clone("copy"); rootb->Reset();
  for(int i=0; i<(int)bkg.size(); i++) rootb->Add(bkg[i].hist);
  for(int bin=0; bin<rootb->GetNbinsX(); bin++){
    double val = rootb->GetBinContent(bin);
    double err = rootb->GetBinError(bin);
    rootb->SetBinContent(bin, std::sqrt(val));
    rootb->SetBinError(bin, err/(2.0*std::sqrt(val)));
  }
  TH1F * srb = (TH1F *)sig[0].hist->Clone("copy");
  srb->Divide(rootb);
  SetHistoStyle(srb, plotname, kBlack);
  srb->GetYaxis()->SetTitle("S/sqrt{B}");
  srb->GetYaxis()->SetRangeUser(0, 1);
  srb->GetYaxis()->SetNdivisions(5, kTRUE);
  if(toSetRange) srb->GetXaxis()->SetRangeUser(xmin, xmax);
  srb->SetStats(0);

  ratioPad->cd();
  ratioPad->SetGrid();

  if(toOverlayData){
    ratio->Draw("ep");
    if(toPlotUncertainty){
      err->Draw("SAME P E2");
    }
    TLine line;
    line.SetLineColor(kRed); line.SetLineWidth(2); line.SetLineStyle(1);
    float min_x = datahist->GetXaxis()->GetXmin();
    float max_x = datahist->GetXaxis()->GetXmax();
    if(toSetRange) line.DrawLine(xmin,1,xmax,1);
    else line.DrawLine(min_x,1,max_x,1);
  }
  
  else srb->Draw("ep");

  //######################################################################
  TString date_stamp = todays_date();
  TString dump_folder = "outputs/"+date_stamp+"_"+jobname+tag;
  createFolder(dump_folder);

  if(toSave){
    c1->SaveAs(dump_folder+"/"+plotname+".png");
    //Deleting pointers to free up memory:
    c1->Clear();
    delete c1;
    delete lg1;
    delete rootb;
    delete ratio;
    delete datahist;
    delete datahist2;
    delete stack;
    delete srb;
    delete err;
    for(auto p : bkg){delete p.hist;} bkg.clear();
    //for(auto q : sig){delete q.hist;} sig.clear();
  }

  
  
}
