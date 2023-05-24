#include <vector>
#include <iostream>
#include <fstream>

bool ifexists(TString file_){
  std::string filepath = (std::string)file_;
  std::ifstream file(filepath);
  return file.good();
}

struct hst{
  TString group;
  TString subgr;
  float lumi;
  TH1F *hist;
};

struct merged{
  TH1F *hist;
  TString name;
  int color;
};

//Sample list:
vector<hst> hst_qcd = {
  {.group = "QCD_MuEnriched", .subgr="20to30",   .lumi= 23.893},
  {.group = "QCD_MuEnriched", .subgr="30to50",   .lumi= 42.905},
  {.group = "QCD_MuEnriched", .subgr="50to80",   .lumi= 105.88},
  {.group = "QCD_MuEnriched", .subgr="80to120",  .lumi= 508.715},
  {.group = "QCD_MuEnriched", .subgr="120to170", .lumi= 1802.854},
  {.group = "QCD_MuEnriched", .subgr="170to300", .lumi= 10265.816},
  {.group = "QCD_MuEnriched", .subgr="300to470", .lumi= 95249.234},
  {.group = "QCD_MuEnriched", .subgr="470to600", .lumi= 656872.125},
  {.group = "QCD_MuEnriched", .subgr="600to800", .lumi= 2060827.75},
  {.group = "QCD_MuEnriched", .subgr="800to1000",.lumi= 24060648}
};

vector<hst> hst_wjets = {
  {.group = "HTbinnedWJets", .subgr="70to100",   .lumi= 52100.898},
  {.group = "HTbinnedWJets", .subgr="100to200",  .lumi= 41127.176},
  {.group = "HTbinnedWJets", .subgr="200to400",  .lumi= 172265.187},
  {.group = "HTbinnedWJets", .subgr="400to600",  .lumi= 163821.094},
  {.group = "HTbinnedWJets", .subgr="600to800",  .lumi= 708793.812},
  {.group = "HTbinnedWJets", .subgr="800to1200", .lumi= 1481985.25},
  {.group = "HTbinnedWJets", .subgr="1200to2500",.lumi= 5602003.5},
  {.group = "HTbinnedWJets", .subgr="2500toInf", .lumi= 79396216}
};

vector<hst> hst_dy = {
  {.group = "DYJetsToLL", .subgr="M10to50", .lumi=5925.521},
  {.group = "DYJetsToLL", .subgr="M50",     .lumi=30321.148}
};

vector<hst> hst_singletop = {
  {.group = "SingleTop", .subgr="s-channel_LeptonDecays",            .lumi=5456748.5},
  {.group = "SingleTop", .subgr="t-channel_AntiTop_InclusiveDecays", .lumi=1407728.5},
  {.group = "SingleTop", .subgr="t-channel_Top_InclusiveDecays",     .lumi=1558659.25},
  {.group = "SingleTop", .subgr="tW_AntiTop_InclusiceDecays",        .lumi=238357.437},
  {.group = "SingleTop", .subgr="tW_Top_InclusiveDecays",            .lumi=245177.187}
};

vector<hst> hst_ttbar = {
  {.group = "TTBar", .subgr="TTTo2L2Nu",        .lumi=1642541.625},
  {.group = "TTBar", .subgr="TTToSemiLeptonic", .lumi=1304012.875}
};

vector<hst> hst_ww = {
  {.group = "WW", .subgr="WWTo1L1Nu2Q", .lumi=787485.562},
  {.group = "WW", .subgr="WWTo4Q",      .lumi=773049.75}
};

vector<hst> hst_wz = {
  {.group = "WZ", .subgr="WZTo1L1Nu2Q", .lumi=805170},
  {.group = "WZ", .subgr="WZTo2Q2L",    .lumi=4499606.5},
  {.group = "WZ", .subgr="WZTo3LNu",    .lumi=1889798.5}
};

vector<hst> hst_zz = {
  {.group = "ZZ", .subgr="ZZTo2L2Nu", .lumi=56787840},
  {.group = "ZZ", .subgr="ZZTo2Q2L",  .lumi=7928149},
  {.group = "ZZ", .subgr="ZZTo2Q2Nu", .lumi=4405016.5},
  {.group = "ZZ", .subgr="ZZTo4L",    .lumi=74330560}
};

vector<hst> hst_data ={
  {.group = "SingleMuon", .subgr="A", .lumi=59700},
  {.group = "SingleMuon", .subgr="B", .lumi=59700},
  {.group = "SingleMuon", .subgr="C", .lumi=59700},
  {.group = "SingleMuon", .subgr="D", .lumi=59700},
};

//#################################################################################################
vector<hst> read_files(vector<hst> vec, TString path, TString jobname, TString plotname){
  cout<<"Reading "<<vec[0].group<<" files ....";
  vector<hst> newvec;
  for(int i=0; i<(int)vec.size(); i++){
    hst temp;
    temp.group = vec[i].group;
    temp.subgr = vec[i].subgr;
    temp.lumi = vec[i].lumi;
    
    TString filename = path+jobname+"/"+vec[i].group+"/"+vec[i].group+"_"+vec[i].subgr+".root";
    if(ifexists(filename)){
      TFile *tfile = new TFile(filename);
      TH1F *hist = (TH1F *)tfile->Get(plotname);
      
      //Lumiscaling:
      float scalefactor = 59700/vec[i].lumi;
      hist->Scale(scalefactor);
      
      temp.hist = hist;
      newvec.push_back(temp);
    }  
  }
  cout<<".... success!"<<endl;
  return newvec;
}

merged merge_and_decorate(vector<hst> vec, TString samplename, TString plotname, int color){
  merged HIST;
  TList *list = new TList;
  for(int i=0; i<(int)vec.size(); i++) list->Add(vec[i].hist);
  
  TH1F *h = (TH1F *)vec[0].hist->Clone("copy");
  h->Reset();
  h->Merge(list);

  HIST.hist = h;
  HIST.color = color;
  HIST.name = samplename;
  return HIST; 
}

void SetOverflowBin(TH1F *h){
  int lastbin = h->GetNbinsX();
  float lastbincontent = h->GetBinContent(lastbin);
  float overflow = h->GetBinContent(lastbin+1);
  h->SetBinContent(lastbin, lastbincontent+overflow);
}
