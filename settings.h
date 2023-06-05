#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib> // Required for running shell script
#include <array>   
#include <memory> // for keeping the output of the shell script

bool ifexists(TString file_){
  std::string filepath = (std::string)file_;
  std::ifstream file(filepath);
  return file.good();
}

void createFolder(TString foldername){
  std::string processline = "mkdir -p "+(std::string)foldername;
  int result = system(processline.c_str());
  if (result!=0) cout<<"Error : Failed to create folder."<<endl;
}

std::string todays_date(){
  std::string processline = "date +%Y-%m-%d";
  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(processline.c_str(), "r"), pclose);
  if(!pipe) throw std::runtime_error("Failed to run Bash script.");
  //Storing the buffer data in 'result'
  while(fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) result += buffer.data();
  // Remove trailing newline characters
  while(!result.empty() && (result.back() == '\n' || result.back() == '\r')) result.pop_back();
  return result;
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
  {.group = "QCD_MuEnriched", .subgr="30to50",   .lumi= 42.906},
  {.group = "QCD_MuEnriched", .subgr="50to80",   .lumi= 105.880},
  {.group = "QCD_MuEnriched", .subgr="80to120",  .lumi= 508.715},
  {.group = "QCD_MuEnriched", .subgr="120to170", .lumi= 1802.854},
  {.group = "QCD_MuEnriched", .subgr="170to300", .lumi= 10265.815},
  {.group = "QCD_MuEnriched", .subgr="300to470", .lumi= 95249.242},
  {.group = "QCD_MuEnriched", .subgr="470to600", .lumi= 656872.156},
  {.group = "QCD_MuEnriched", .subgr="600to800", .lumi= 2060827.812},
  {.group = "QCD_MuEnriched", .subgr="800to1000",.lumi= 11337379.457}
};

vector<hst> hst_wjets = {
  {.group = "HTbinnedWJets", .subgr="70to100",   .lumi= 52100.910},
  {.group = "HTbinnedWJets", .subgr="100to200",  .lumi= 41127.174},
  {.group = "HTbinnedWJets", .subgr="200to400",  .lumi= 172265.183},
  {.group = "HTbinnedWJets", .subgr="400to600",  .lumi= 163821.083},
  {.group = "HTbinnedWJets", .subgr="600to800",  .lumi= 708793.021},
  {.group = "HTbinnedWJets", .subgr="800to1200", .lumi= 1481985.193},
  {.group = "HTbinnedWJets", .subgr="1200to2500",.lumi= 5602003.457},
  {.group = "HTbinnedWJets", .subgr="2500toInf", .lumi= 79396214.989}
};

vector<hst> hst_dy = {
  {.group = "DYJetsToLL", .subgr="M10to50", .lumi=5925.522},
  {.group = "DYJetsToLL", .subgr="M50",     .lumi=30321.155}
};

vector<hst> hst_singletop = {
  {.group = "SingleTop", .subgr="s-channel_LeptonDecays",            .lumi=5456748.098},
  {.group = "SingleTop", .subgr="t-channel_AntiTop_InclusiveDecays", .lumi=1407728.544},
  {.group = "SingleTop", .subgr="t-channel_Top_InclusiveDecays",     .lumi=1572627.866},
  {.group = "SingleTop", .subgr="tW_AntiTop_InclusiceDecays",        .lumi=238357.428},
  {.group = "SingleTop", .subgr="tW_Top_InclusiveDecays",            .lumi=245177.196}
};

vector<hst> hst_ttbar = {
  {.group = "TTBar", .subgr="TTTo2L2Nu",        .lumi=1642541.624},
  {.group = "TTBar", .subgr="TTToSemiLeptonic", .lumi=1304012.700}
};

vector<hst> hst_ww = {
  {.group = "WW", .subgr="WWTo1L1Nu2Q", .lumi=805257.731},
  {.group = "WW", .subgr="WWTo4Q",      .lumi=773049.853}
};

vector<hst> hst_wz = {
  {.group = "WZ", .subgr="WZTo1L1Nu2Q", .lumi=805257.731},
  {.group = "WZ", .subgr="WZTo2Q2L",    .lumi=4499605.731},
  {.group = "WZ", .subgr="WZTo3LNu",    .lumi=1889798.538}
};

vector<hst> hst_zz = {
  {.group = "ZZ", .subgr="ZZTo2L2Nu", .lumi=58416512.631},
  {.group = "ZZ", .subgr="ZZTo2Q2L",  .lumi=7928149.608},
  {.group = "ZZ", .subgr="ZZTo2Q2Nu", .lumi=4405016.452},
  {.group = "ZZ", .subgr="ZZTo4L",    .lumi=74330566.038}
};

vector<hst> hst_data ={
  {.group = "SingleMuon", .subgr="A", .lumi=59700},
  {.group = "SingleMuon", .subgr="B", .lumi=59700},
  {.group = "SingleMuon", .subgr="C", .lumi=59700},
  {.group = "SingleMuon", .subgr="D", .lumi=59700},
};

//#################################################################################################
vector<hst> read_files(vector<hst> vec, TString path, TString jobname, TString plotname){
  //cout<<"Reading "<<vec[0].group<<" files ....";
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
  //cout<<".... success!"<<endl;
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


//############################################
// Finding QCD HT-binned scale factors:
//############################################
void find_sf(TH1F*qcd, TH1F*data){
  cout<<"\nCalculating HT binned sf for QCD ..."<<endl;
  cout<<"bin\trange\tnqcd\tndata\tsf"<<endl;
  
  vector<float> scaleqcd; scaleqcd.clear();
  
  int nbins = data->GetNbinsX();
  for(int bin=1; bin<nbins+2; bin++){
    float nqcd = qcd->GetBinContent(bin);
    float ndata = data->GetBinContent(bin);
    float xlow = data->GetXaxis()->GetBinLowEdge(bin);
    float xhigh = data->GetXaxis()->GetBinUpEdge(bin);

    TString range = to_string((int)xlow)+"-"+to_string((int)xhigh);
    if(bin>nbins) range = "overflow";
    
    float sf = ndata/nqcd;

    cout<<bin<<"\t"<<range<<"\t"<<(int)nqcd<<"\t"<<(int)ndata<<"\t"<<sf<<endl;

    scaleqcd.push_back(sf);
  }
}
