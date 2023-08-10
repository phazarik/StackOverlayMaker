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
#include "settings_bkgonly.h"
#include "decorations.h"
using namespace std;

//This is the bulk of the code.
void plot(TString path, TString jobname, TString plotname, TString plottitle, int binw, float xmin, float xmax);

//#######################################################################################
//This is the main function:
void stackbkgonly(){

  TString path = "/mnt/d/work/GitHub/StackOverlayMaker/inputs/";
  TString jobname = "hst_Aug10_Basic";
  //TString jobname = "hst_Aug10_NonIso";
 
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
    /*
    {.name="SS_dimuon_mass",  .title="dimuon mass (SS pair)",  .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_dEta_muons",   .title="dEta (SS pair)",.rebin = 20, .xmin= 0, .xmax=6},
    {.name="SS_dPhi_muons",   .title="dPhi (SS pair)",.rebin = 20, .xmin= 0, .xmax=6},
    {.name="SS_dR_muons",     .title="dR (SS pair)",  .rebin = 20, .xmin= 0, .xmax=6},
    {.name="SS_ptratio",      .title="pT ratio (SS pair)",.rebin = 10, .xmin= 0, .xmax=1},
    {.name="SS_met",          .title="MET",           .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_metphi",       .title="MET phi",       .rebin = 10, .xmin=-4, .xmax=4},
    {.name="SS_LT",           .title="Sum(mu pT)",    .rebin = 10, .xmin= 0, .xmax=200},
    {.name="SS_HT",           .title="Sum(Jet pT)",   .rebin = 100, .xmin= 0, .xmax=1000},
    {.name="SS_nJet",         .title="nJet",          .rebin =  1, .xmin= 0, .xmax=10},
    {.name="SS_nbJet",        .title="nbJet",         .rebin =  1, .xmin= 0, .xmax=10},
    {.name="SS_dphi_mu0_met",   .title="dPhi (mu0, met)",.rebin = 20, .xmin= 0, .xmax=6},
    {.name="SS_dphi_muss_met",   .title="dPhi (mu1, met)",.rebin = 20, .xmin= 0, .xmax=6},
    {.name="SS_maxdphimet",   .title="max(dPhi mu, met)",.rebin = 20, .xmin= 0, .xmax=6},*/
    /*
    {.name="SS_cutflow",        .title="cutflow",       .rebin = 1, .xmin= 0, .xmax=10},
    {.name="SS_cutflow_leading",.title="cutflow (mu0)", .rebin = 1, .xmin= 0, .xmax=20},
    {.name="SS_cutflow_ssmuon", .title="cutflow (muss)",.rebin = 1, .xmin= 0, .xmax=20},*/

    //{.name="SS_dimuon_mass",  .title="dimuon mass (SS pair)",  .rebin = 10, .xmin= 0, .xmax=200},
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

  cout<<"\nSucess!!\n"<<endl;
}

//#######################################################################################

void plot(TString path, TString jobname, TString plotname, TString plottitle, int binw, float xmin, float xmax){
  
  //Global settings: 
  bool toStack        = true;
  bool toSave         = true;
  bool toLog          = true;
  bool toSetRange     = true;
  bool toOverlaySig   = true;
  bool toOverlayData  = false;
  bool toScaleHT      = false;
  bool toPrintBinInfo = false;
  float sigscale = 1;
  TString tag = "_afterQCDscaling";
  
  //###################################
  //Reading histograms and sorting them
  //###################################
  
  //Reading histograms from each file, lumiscaling them  and storing them in a vector<struct>
  vector<hst> QCD, WJets, DY, ST, TTbar, WW, WZ, ZZ, TTW, Data;
  QCD   = read_files(hst_qcd,       path, jobname, plotname);
  WJets = read_files(hst_wjets,     path, jobname, plotname);
  DY    = read_files(hst_dy,        path, jobname, plotname);
  ST    = read_files(hst_singletop, path, jobname, plotname);
  TTbar = read_files(hst_ttbar,     path, jobname, plotname);
  WW    = read_files(hst_ww,        path, jobname, plotname);
  WZ    = read_files(hst_wz,        path, jobname, plotname);
  ZZ    = read_files(hst_zz,        path, jobname, plotname);
  TTW   = read_files(hst_ttw,       path, jobname, plotname);

  if(toOverlayData) Data  = read_files(hst_data,      path, jobname, plotname);

  //Push back all the backgrounds into a vector:
  vector<merged> bkg;
  bkg.push_back(merge_and_decorate(QCD,   "QCD",   plotname, kYellow));
  bkg.push_back(merge_and_decorate(DY,    "DY",    plotname, kRed-9));
  bkg.push_back(merge_and_decorate(WJets, "WJets", plotname, kGray+1));
  bkg.push_back(merge_and_decorate(ST,    "ST",    plotname, kCyan-9));
  bkg.push_back(merge_and_decorate(TTbar, "TTBar", plotname, kCyan-10));
  bkg.push_back(merge_and_decorate(WW,    "WW",    plotname, kGreen-4));
  bkg.push_back(merge_and_decorate(WZ,    "WZ",    plotname, kBlue-9));
  bkg.push_back(merge_and_decorate(ZZ,    "ZZ",    plotname, kBlue-7));
  bkg.push_back(merge_and_decorate(TTW,   "TTW",   plotname, kBlue-7));

  //Decorating the histograms and setting the bins:
  //Decorating the backgrounds:
  for(int i=0; i<(int)bkg.size(); i++){
    SetHistoStyle(bkg[i].hist, plotname, bkg[i].color);
    bkg[i].hist->Rebin(binw);
    SetOverflowBin(bkg[i].hist);
  }

  //SIGNAL:
  TString sigpath="inputs/"+jobname+"/signal/";
  vector<TString> sigfile = {"hst_VLL100.root", "hst_VLL150.root"};
  vector<float> lumi = {508761.54, 2092051.72}; // 595251 events /1.17 pb, 606695 events /0.290 pb
  vector<int> col = {kRed, kBlack};
  //Make sure that the sizes of sigfile and lumi, col are the same.

  vector<TH1F*> sighist;
  for(int i=0; i<(int)sigfile.size(); i++){
    TString fullpath = sigpath+sigfile[i];
    TFile *t = new TFile(fullpath);
    TH1F *h = (TH1F *)t->Get(plotname);
    h->Scale(59700/lumi[i]);
    h->Scale(sigscale);
    SetHistoStyle(h, plotname, col[i]);
    h->Rebin(binw);
    SetOverflowBin(h);
    sighist.push_back(h);
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
  nsig = sighist[0]->Integral();
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
    int nbins = (int)sighist[0]->GetNbinsX();
    for(int bin=0; bin<nbins; bin++){
      float nsig = sighist[0]->GetBinContent(bin);
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
  
  TPad *mainPad = create_mainpad();
  mainPad->SetMargin(mainPad->GetLeftMargin(), mainPad->GetRightMargin(), 0.15, 0.1);
  if(toLog) mainPad->SetLogy(1);
  mainPad->Draw();
  
  TPad *ratioPad = create_ratiopad();
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
    TString sigcount100 = "VLL100 ["+to_string((int)(sighist[0]->Integral()/sigscale))+"]";
    lg1->AddEntry(sighist[0], sigcount100, "f");
    TString sigcount150 = "VLL150 ["+to_string((int)(sighist[1]->Integral()/sigscale))+"]";
    lg1->AddEntry(sighist[1], sigcount150, "f");
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
    if(toOverlaySig) for(int i=0; i<(int)sighist.size(); i++) sighist[i]->Scale(1/sighist[i]->Integral());
    for(int i=0; i<(int)bkg.size(); i++) bkg[i].hist->Scale(1/bkg[i].hist->Integral());
  }

  sighist[0]->GetXaxis()->SetTitle(plottitle);
  sighist[0]->GetYaxis()->SetTitle("Events");
  sighist[0]->SetStats(0);
  if(toSetRange) sighist[0]->GetXaxis()->SetRangeUser(xmin, xmax);
  sighist[0]->GetYaxis()->SetRangeUser(0.1, 1E8);

  sighist[0]->Draw("hist");
  stack->Draw("hist same");
  if(toOverlayData) datahist->Draw("ep same");
  for(int i =0; i<(int)sighist.size(); i++) sighist[i]->Draw("hist same");
  lg1->Draw("same");


  //#########################
  // Drawing on the ratiopad
  //#########################
  TH1F *rootb, *ratio;

  //Setting up obsbyexp:
  if(toOverlayData){
    TH1F *exp = (TH1F*)bkg[0].hist->Clone("copy"); exp->Reset();
    for(int i=0; i<(int)bkg.size(); i++) exp->Add(bkg[i].hist);
    ratio = (TH1F*)datahist->Clone("copy");
    ratio->Divide(exp);
    SetHistoStyle(ratio, plotname, kBlack);
    ratio->GetYaxis()->SetTitle("obs/exp");
    ratio->GetYaxis()->SetRangeUser(0, 2);
    ratio->GetYaxis()->SetNdivisions(5, kTRUE);
    if(toSetRange) ratio->GetXaxis()->SetRangeUser(xmin, xmax);
    ratio->SetStats(0);
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
  TH1F * srb = (TH1F *)sighist[0]->Clone("copy");
  srb->Divide(rootb);
  SetHistoStyle(srb, plotname, kBlack);
  srb->GetYaxis()->SetTitle("S/sqrt{B}");
  srb->GetYaxis()->SetRangeUser(0, 1);
  srb->GetYaxis()->SetNdivisions(5, kTRUE);
  if(toSetRange) srb->GetXaxis()->SetRangeUser(xmin, xmax);
  srb->SetStats(0);

  ratioPad->cd();

  if(toOverlayData){
    ratio->Draw("ep");
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
  }
  
}
